#pragma once

#include "graph.h"
#include "atomics.h"
#include <memory>
#include <utility>
#include <set>

namespace gtool {

typedef std::ptrdiff_t offset_diff_t;

template<typename T, typename DstT = T>
class StreamGraph {
public:
    typedef T vertex_type;
    typedef DstT destination_type;
    struct Neighborhood {
        T n;
        offset_t *l_offset;
        offset_t *r_offset;
        T *neigh;
        typedef DstT const *iterator;
        iterator begin() { return neigh + l_offset[n]; }
        iterator end() { return neigh + r_offset[n]; }
    };
    StreamGraph(): vertex_number_{}, edge_capacity_{},
        edge_number_{}, edge_number_prev_{},
        out_l_offset_{}, out_r_offset_{}, out_r_offset_prev_{},
        in_l_offset_{}, in_r_offset_{}, in_r_offset_prev_{},
        out_neigh_{}, in_neigh_{},
#if defined(BUILD_WITH_RESTORE)
         out_r_offset_backup_{}, in_r_offset_backup_{},
         out_neigh_backup_{}, in_neigh_backup_{}
#endif
    {}
    StreamGraph(std::size_t const vertex_number,
                std::size_t const edge_capacity,
                std::unique_ptr<offset_t[]> out_l_offset,
                std::unique_ptr<offset_t[]> out_r_offset,
                std::unique_ptr<DstT[]> out_neigh,
                std::unique_ptr<offset_t[]> in_l_offset,
                std::unique_ptr<offset_t[]> in_r_offset,
                std::unique_ptr<DstT[]> in_neigh)
            : vertex_number_{vertex_number}, edge_capacity_{edge_capacity}, edge_number_{}, edge_number_prev_{},
            out_l_offset_(std::move(out_l_offset)), out_r_offset_(std::move(out_r_offset)), out_neigh_(std::move(out_neigh)),
            in_l_offset_(std::move(in_l_offset)), in_r_offset_(std::move(in_r_offset)), in_neigh_(std::move(in_neigh)),
            out_r_offset_prev_(new offset_t[vertex_number + 1]),
            in_r_offset_prev_(new offset_t[vertex_number + 1]),
#if defined (BUILD_WITH_RESTORE)
            out_r_offset_backup_{new offset_t[vertex_number + 1]},
            in_r_offset_backup_{new offset_t[vertex_number + 1]},
            out_neigh_backup_{new DstT[edge_capacity]},
            in_neigh_backup_{new DstT[edge_capacity]}
#endif
    {
        int tmp{};
#pragma omp parallel for reduction(+ : tmp)
        for (int i = 0; i < vertex_number; ++i) {
            tmp += out_r_offset_[i] - out_l_offset_[i];
        }
        edge_number_prev_ = edge_number_ = tmp;
        std::memcpy(out_r_offset_prev_.get(), out_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(in_r_offset_prev_.get(), in_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
#if defined (BUILD_WITH_RESTORE)
        std::memcpy(out_r_offset_backup_.get(), out_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(in_r_offset_backup_.get(), in_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(out_neigh_backup_.get(), out_neigh_.get(), sizeof(DstT) * edge_capacity_);
        std::memcpy(in_neigh_backup_.get(), in_neigh_.get(), sizeof(DstT) * edge_capacity_);
#endif
    }
    ~StreamGraph() = default;
    [[nodiscard]] std::size_t get_vertex_number() const { return vertex_number_; }
    [[nodiscard]] std::size_t get_edge_capacity() const { return edge_capacity_; }
    [[nodiscard]] std::size_t get_edge_number() const { return edge_number_; }
    [[nodiscard]] std::size_t get_edge_number_prev() const { return edge_number_prev_; }
    offset_t out_degree(T n) const { return out_r_offset_[n] - out_l_offset_[n]; }
    offset_t in_degree(T n) const { return in_r_offset_[n] - in_l_offset_[n]; }
    offset_t out_degree_prev(T n) const { return out_r_offset_prev_[n] - out_l_offset_[n]; }
    offset_t in_degree_prev(T n) const { return in_r_offset_prev_[n] - in_l_offset_[n]; }
    offset_diff_t out_degree_delta(T n) const { return static_cast<offset_diff_t>(out_r_offset_[n]) - out_r_offset_prev_[n]; }
    offset_diff_t in_degree_delta(T n) const { return static_cast<offset_diff_t>(in_r_offset_[n]) - in_r_offset_prev_[n]; }
    Neighborhood out_neighbors(T n) const { return { n, out_l_offset_.get(), out_r_offset_.get(), out_neigh_.get() }; }
    Neighborhood in_neighbors(T n) const { return { n, in_l_offset_.get(), in_r_offset_.get(), in_neigh_.get() }; }
    Neighborhood out_neighbors_prev(T n) const { return { n, out_l_offset_.get(), out_r_offset_prev_.get(), out_neigh_.get() }; }
    Neighborhood in_neighbors_prev(T n) const { return { n, in_l_offset_.get(), in_r_offset_prev_.get(), in_neigh_.get() }; }
    Neighborhood out_neighbors_delta(T n) const {
        if (out_degree_delta(n) >= 0)
            return { n, out_r_offset_prev_.get(), out_r_offset_.get(), out_neigh_.get() };
        return { n, out_r_offset_.get(), out_r_offset_prev_.get(), out_neigh_.get() };
    }
    Neighborhood in_neighbors_delta(T n) const {
        if (in_degree_delta(n) >= 0)
            return { n, in_r_offset_prev_.get(), in_r_offset_.get(), in_neigh_.get() };
        return { n, in_r_offset_.get(), in_r_offset_prev_.get(), in_neigh_.get() };
    }
    template<typename Iter,
            typename=std::enable_if<std::is_same<
                    typename std::iterator_traits<Iter>::value_type, std::pair<T, DstT>>::value>::type>
    void stream_addition(Iter first, Iter last) {
        int num = std::distance(first, last);
#pragma omp parallel for default(none) shared(first, num, out_r_offset_, out_neigh_)
        for (int i = 0; i < num; ++i) {
            auto &[from, to] = *(first + i);
            offset_t write_ptr = fetch_and_add(out_r_offset_[from], 1);
            out_neigh_[write_ptr] = to;
        }
#pragma omp parallel for default(none) shared(first, num, out_r_offset_, out_neigh_)
        for (int i = 0; i < num; ++i) {
            auto &[from, to] = *(first + i);
            offset_t write_ptr = fetch_and_add(in_r_offset_[get_dst_id(to)], 1);
            in_neigh_[write_ptr] = from;
        }
        edge_number_ += num;
    }
    template<typename Iter,
            typename=std::enable_if<std::is_same<
                    typename std::iterator_traits<Iter>::value_type, std::pair<T, DstT>>::value>::type>
    void stream_deletion(Iter first, Iter last) {
        std::vector<std::reference_wrapper<typename std::iterator_traits<decltype(first)>::value_type>> deletion_refs(first, last);
        std::sort(deletion_refs.begin(), deletion_refs.end(), [](auto lhs, auto rhs) { return lhs.get().first < rhs.get().first; });
        using RefIter = std::vector<std::reference_wrapper<typename std::iterator_traits<decltype(first)>::value_type>>::iterator;
        auto find_next_group = [](RefIter beg, RefIter end) {
            if (beg == end)
                return std::make_tuple(beg, end);
            RefIter tmp = beg;
            for (++beg; beg != end; ++beg) {
                if ((*beg).get().first != (*(beg - 1)).get().first)
                    break;
            }
            return std::make_tuple(tmp, beg);
        };
        using Task = std::tuple<RefIter, RefIter>;
        std::vector<Task> tasks;
        auto ref_beg = deletion_refs.begin();
        auto ref_end = deletion_refs.end();
        for (auto [beg, end] = find_next_group(ref_beg, ref_end); beg != end; std::tie(beg, end) = find_next_group(end, ref_end)) {
            tasks.emplace_back(beg, end);
        }
#pragma omp parallel for
        for (int task_id = 0; task_id < tasks.size(); ++task_id) {
            auto [beg, end] = tasks[task_id];
            T u = (*beg).get().first;
            std::vector<bool> bmp(out_degree(u));
            for ( ; beg != end; ++beg) {
                DstT v = (*beg).get().second;
                for (offset_t base{out_l_offset_[u]}, offset{}; offset < out_degree(u); ++offset) {
                    if (!bmp[offset] && out_neigh_[base + offset] == v) {
                        bmp[offset] = true;
                        break;
                    }
                }
            }
            std::vector<DstT> new_order; new_order.reserve(out_degree(u));
            int rest_cnt{};
            for (offset_t base{out_l_offset_[u]}, i{}; i < bmp.size(); ++i) {
                if (!bmp[i]) {
                    ++rest_cnt;
                    new_order.emplace_back(out_neigh_[base + i]);
                }
            }
            for (offset_t base{out_l_offset_[u]}, i{}; i < bmp.size(); ++i) {
                if (bmp[i]) {
                    new_order.emplace_back(out_neigh_[base + i]);
                }
            }
            std::memcpy(&out_neigh_[out_l_offset_[u]], new_order.data(), sizeof(DstT) * new_order.size());
            out_r_offset_[u] = out_l_offset_[u] + rest_cnt;
        }
        edge_number_ = 0;
#pragma omp parallel for reduction(+ : edge_number_)
        for (int i = 0; i < vertex_number_; ++i) {
            edge_number_ += out_degree(i);
        }
        std::sort(deletion_refs.begin(), deletion_refs.end(), [](auto lhs, auto rhs) { return lhs.get().second < rhs.get().second; });
        auto find_next_group_in = [](RefIter beg, RefIter end) {
            if (beg == end)
                return std::make_tuple(beg, end);
            RefIter tmp = beg;
            for (++beg; beg != end; ++beg) {
                if ((*beg).get().second != (*(beg - 1)).get().second)
                    break;
            }
            return std::make_tuple(tmp, beg);
        };
        tasks.clear();
        ref_beg = deletion_refs.begin();
        ref_end = deletion_refs.end();
        for (auto [beg, end] = find_next_group_in(ref_beg, ref_end); beg != end; std::tie(beg, end) = find_next_group_in(end, ref_end)) {
            tasks.emplace_back(beg, end);
        }
#pragma omp parallel for
        for (int task_id = 0; task_id < tasks.size(); ++task_id) {
            auto [beg, end] = tasks[task_id];
            T v = (*beg).get().second;
            std::vector<bool> bmp(in_degree(v));
            for ( ; beg != end; ++beg) {
                DstT u = (*beg).get().first;
                for (offset_t base{in_l_offset_[v]}, offset{}; offset < in_degree(v); ++offset) {
                    if (!bmp[offset] && in_neigh_[base + offset] == u) {
                        bmp[offset] = true;
                        break;
                    }
                }
            }
            std::vector<DstT> new_order; new_order.reserve(in_degree(v));
            int rest_cnt{};
            for (offset_t base{in_l_offset_[v]}, i{}; i < bmp.size(); ++i) {
                if (!bmp[i]) {
                    ++rest_cnt;
                    new_order.emplace_back(in_neigh_[base + i]);
                }
            }
            for (offset_t base{in_l_offset_[v]}, i{}; i < bmp.size(); ++i) {
                if (bmp[i]) {
                    new_order.emplace_back(in_neigh_[base + i]);
                }
            }
            std::memcpy(&in_neigh_[in_l_offset_[v]], new_order.data(), sizeof(DstT) * new_order.size());
            in_r_offset_[v] = in_l_offset_[v] + rest_cnt;
        }
    }
    void update() {
        edge_number_prev_ = edge_number_;
        std::memcpy(out_r_offset_prev_.get(), out_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(in_r_offset_prev_.get(), in_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
    }
#if defined(BUILD_WITH_RESTORE)
    [[nodiscard]] offset_t out_degree_base(T n) const { return out_r_offset_backup_[n] - out_l_offset_[n]; }
    [[nodiscard]] offset_t in_degree_base(T n) const { return in_r_offset_backup_[n] - in_l_offset_[n]; }
    Neighborhood out_neighbors_base(T n) const { return { n, out_l_offset_.get(), out_r_offset_backup_.get(), out_neigh_backup_.get() }; }
    Neighborhood in_neighbors_base(T n) const { return { n, in_l_offset_.get(), in_r_offset_backup_.get(), in_neigh_backup_.get() }; }
    void restore() {
        std::memcpy(out_r_offset_.get(), out_r_offset_backup_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(out_r_offset_prev_.get(), out_r_offset_backup_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(in_r_offset_.get(), in_r_offset_backup_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(in_r_offset_prev_.get(), in_r_offset_backup_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(out_neigh_.get(), out_neigh_backup_.get(), sizeof(DstT) * edge_capacity_);
        std::memcpy(in_neigh_.get(), in_neigh_backup_.get(), sizeof(DstT) * edge_capacity_);
        int tmp{};
#pragma omp parallel for reduction(+ : tmp)
        for (int i = 0; i < vertex_number_; ++i) {
            tmp += out_r_offset_[i] - out_l_offset_[i];
        }
        edge_number_prev_ = edge_number_ = tmp;
    }
#endif
private:
    std::size_t const vertex_number_;
    std::size_t const edge_capacity_;
    std::size_t edge_number_;
    std::size_t edge_number_prev_;
    std::unique_ptr<offset_t[]> out_l_offset_;
    std::unique_ptr<offset_t[]> out_r_offset_;
    std::unique_ptr<DstT[]> out_neigh_;
    std::unique_ptr<offset_t[]> in_l_offset_;
    std::unique_ptr<offset_t[]> in_r_offset_;
    std::unique_ptr<DstT[]> in_neigh_;
    std::unique_ptr<offset_t[]> out_r_offset_prev_;
    std::unique_ptr<offset_t[]> in_r_offset_prev_;
#if defined(BUILD_WITH_RESTORE)
    std::unique_ptr<offset_t[]> out_r_offset_backup_;
    std::unique_ptr<offset_t[]> in_r_offset_backup_;
    std::unique_ptr<DstT[]> out_neigh_backup_;
    std::unique_ptr<DstT[]> in_neigh_backup_;
#endif
};

template<typename T>
struct IsRange {
private:
    template<typename U, typename = decltype(std::declval<T>().begin(), std::declval<T>().end())>
    static std::true_type test(void*);
    template<typename>
    static std::false_type test(...);
public:
    using Type = decltype(test<T>(nullptr));
};
template<typename T>
struct IsRangeT : IsRange<T>::Type {};

template<typename T, typename DstT, typename Rng,
typename=std::enable_if<IsRangeT<typename std::decay<Rng>::type>::value>::type>
class Streamer {
public:
    Streamer(StreamGraph<T, DstT> &stream_graph, Rng &&delta, int max_batch_size=1'000)
        : stream_graph_(stream_graph), delta_(std::forward<Rng>(delta)),
        offset_{}, max_batch_size_{max_batch_size} {}
    [[nodiscard]] int max_batch_size() const { return max_batch_size_; }
    void max_batch_size(int batch_size) { max_batch_size_ = batch_size; }
    [[nodiscard]] int suggest_max_batch_size(int batch_number) const {
        return (delta_.size() + batch_number) / batch_number;
    }
    [[nodiscard]] bool has_more_batches() const { return offset_ < delta_.size(); }
    [[nodiscard]] int next_batch_size() const { return std::min(static_cast<int>(delta_.size()) - offset_, max_batch_size_); }
    auto view_next_batch() {
        auto beg = delta_.begin() + offset_;
        int size = next_batch_size();
        return std::make_tuple(beg, beg + size);
    }
    void stream_next_batch(bool with_update = false) {
        stream_next_batch_detail();
        if (with_update) {
            stream_graph_.update();
        }
    }
    void update_graph() {
        stream_graph_.update();
    }
#if defined (BUILD_WITH_RESTORE)
    void reset() {
        offset_ = 0;
    }
#endif
protected:
    virtual void stream_next_batch_detail() {
        auto beg = delta_.begin() + offset_;
        int size = next_batch_size();
        stream_graph_.stream_addition(beg, beg + size);
        offset_ += size;
    }
    StreamGraph<T, DstT> &stream_graph_;
    Rng delta_;
    int offset_;
    int max_batch_size_;
};

template<typename T, typename DstT, typename Rng,
        typename=std::enable_if<IsRangeT<typename std::decay<Rng>::type>::value>::type>
class DStreamer : public Streamer<T, DstT, Rng> {
public:
    DStreamer(StreamGraph<T, DstT> &stream_graph, Rng &&delta, int max_batch_size=1'00):
        Streamer<T, DstT, Rng>(stream_graph, std::forward<Rng>(delta), max_batch_size) {}\
protected:
    void stream_next_batch_detail() override {
        auto beg = this->delta_.begin() + this->offset_;
        int size = this->next_batch_size();
        this->stream_graph_.stream_deletion(beg, beg + size);
        this->offset_ += size;
    }
};

}