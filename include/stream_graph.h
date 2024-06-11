#pragma once

#include "graph.h"
#include "atomics.h"
#include <memory>
#include <utility>

namespace gtool {

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
    StreamGraph(): vertex_number_{}, vertex_number_prev_{}, edge_number_{}, edge_number_prev_{},
        out_l_offset_{}, out_r_offset_{}, out_r_offset_prev_{},
        in_l_offset_{}, in_r_offset_{}, in_r_offset_prev_{}
#if defined(BUILD_WITH_RESTORE)
        , vertex_number_backup_{}, edge_number_backup_{}, out_r_offset_backup_{}, in_r_offset_backup_{}
#endif
    {}
    StreamGraph(std::size_t vertex_number,
                std::size_t edge_number,
                std::unique_ptr<offset_t[]> out_l_offset,
                std::unique_ptr<offset_t[]> out_r_offset,
                std::unique_ptr<DstT[]> out_neigh,
                std::unique_ptr<offset_t[]> in_l_offset,
                std::unique_ptr<offset_t[]> in_r_offset,
                std::unique_ptr<DstT[]> in_neigh)
            : vertex_number_{vertex_number}, edge_number_{edge_number},
            vertex_number_prev_{vertex_number}, edge_number_prev_{edge_number},
            out_l_offset_(std::move(out_l_offset)), out_r_offset_(std::move(out_r_offset)), out_neigh_(std::move(out_neigh)),
            in_l_offset_(std::move(in_l_offset)), in_r_offset_(std::move(in_r_offset)), in_neigh_(std::move(in_neigh)),
            out_r_offset_prev_(new offset_t[vertex_number + 1]),
            in_r_offset_prev_(new offset_t[vertex_number + 1])
#if defined (BUILD_WITH_RESTORE)
            , vertex_number_backup_{vertex_number}, edge_number_backup_{edge_number},
            out_r_offset_backup_{new offset_t[vertex_number + 1]},
            in_r_offset_backup_{new offset_t[vertex_number + 1]}
#endif
    {
        std::memcpy(out_r_offset_prev_.get(), out_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(in_r_offset_prev_.get(), in_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
#if defined (BUILD_WITH_RESTORE)
        std::memcpy(out_r_offset_backup_.get(), out_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(in_r_offset_backup_.get(), in_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
#endif
    }
    ~StreamGraph() = default;
    [[nodiscard]] std::size_t get_vertex_number() const { return vertex_number_; }
    [[nodiscard]] std::size_t get_edge_number() const { return edge_number_; }
    [[nodiscard]] std::size_t get_vertex_number_prev() const { return vertex_number_prev_; }
    [[nodiscard]] std::size_t get_edge_number_prev() const { return edge_number_prev_; }
    offset_t out_degree(T n) const { return out_r_offset_[n] - out_l_offset_[n]; }
    offset_t in_degree(T n) const { return in_r_offset_[n] - in_l_offset_[n]; }
    offset_t out_degree_prev(T n) const { return out_r_offset_prev_[n] - out_l_offset_[n]; }
    offset_t in_degree_prev(T n) const { return in_r_offset_prev_[n] - in_l_offset_[n]; }
    offset_t out_degree_delta(T n) const { return out_r_offset_[n] - out_r_offset_prev_[n]; }
    offset_t in_degree_delta(T n) const { return in_r_offset_[n] - in_r_offset_prev_[n]; }
    Neighborhood out_neighbors(T n) const { return { n, out_l_offset_.get(), out_r_offset_.get(), out_neigh_.get() }; }
    Neighborhood in_neighbors(T n) const { return { n, in_l_offset_.get(), in_r_offset_.get(), in_neigh_.get() }; }
    Neighborhood out_neighbors_prev(T n) const { return { n, out_l_offset_.get(), out_r_offset_prev_.get(), out_neigh_.get() }; }
    Neighborhood in_neighbors_prev(T n) const { return { n, in_l_offset_.get(), in_r_offset_prev_.get(), in_neigh_.get() }; }
    Neighborhood out_neighbors_delta(T n) const { return { n, out_r_offset_prev_.get(), out_r_offset_.get(), out_neigh_.get() }; }
    Neighborhood in_neighbors_delta(T n) const { return { n, in_r_offset_prev_.get(), in_r_offset_.get(), in_neigh_.get() }; }
    template<typename Iter,
            typename=std::enable_if_t<std::is_same_v<
                    typename std::iterator_traits<Iter>::value_type, std::pair<T, DstT>>>>
    void stream_edge(Iter first, Iter last) {
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
    void update() {
        vertex_number_prev_ = vertex_number_;
        edge_number_prev_ = edge_number_;
        std::memcpy(out_r_offset_prev_.get(), out_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(in_r_offset_prev_.get(), in_r_offset_.get(), sizeof(offset_t) * (vertex_number_ + 1));
    }
#if defined(BUILD_WITH_RESTORE)
    [[nodiscard]] offset_t out_degree_base(T n) const { return out_r_offset_backup_[n] - out_l_offset_[n]; }
    [[nodiscard]] offset_t in_degree_base(T n) const { return in_r_offset_backup_[n] - in_l_offset_[n]; }
    Neighborhood out_neighbors_base(T n) const { return { n, out_l_offset_.get(), out_r_offset_backup_.get(), out_neigh_.get() }; }
    Neighborhood in_neighbors_base(T n) const { return { n, in_l_offset_.get(), in_r_offset_backup_.get(), in_neigh_.get() }; }
    void restore() {
        vertex_number_ = vertex_number_prev_ = vertex_number_backup_;
        edge_number_ = edge_number_prev_ = edge_number_backup_;
        std::memcpy(out_r_offset_.get(), out_r_offset_backup_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(out_r_offset_prev_.get(), out_r_offset_backup_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(in_r_offset_.get(), in_r_offset_backup_.get(), sizeof(offset_t) * (vertex_number_ + 1));
        std::memcpy(in_r_offset_prev_.get(), in_r_offset_backup_.get(), sizeof(offset_t) * (vertex_number_ + 1));
    }
#endif
private:
    std::size_t vertex_number_;
    std::size_t edge_number_;
    std::size_t vertex_number_prev_;
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
    std::size_t vertex_number_backup_;
    std::size_t edge_number_backup_;
    std::unique_ptr<offset_t[]> out_r_offset_backup_;
    std::unique_ptr<offset_t[]> in_r_offset_backup_;
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
typename=std::enable_if_t<IsRangeT<std::decay_t<Rng>>::value>>
class Streamer {
public:
    Streamer(StreamGraph<T, DstT> &stream_graph, Rng &&delta, int max_batch_size=1'000)
        : stream_graph_(stream_graph), delta_(std::forward<Rng>(delta)),
        offset_{}, max_batch_size_{max_batch_size}, finish_{false} {}
    [[nodiscard]] int max_batch_size() const { return max_batch_size_; }
    void max_batch_size(int batch_size) { max_batch_size_ = batch_size; }
    [[nodiscard]] int suggest_max_batch_size(int batch_number) const {
        return (delta_.size() + batch_number) / batch_number;
    }
    [[nodiscard]] bool has_more_batches() const { return offset_ < delta_.size(); }
    [[nodiscard]] int next_batch_size() const { return std::min(static_cast<int>(delta_.size()) - offset_, max_batch_size_); }
    auto view_next_batch() {
        if (finish_) {
            throw std::runtime_error("Error: Stream when no more batches exist.");
        }
        auto beg = delta_.begin() + offset_;
        int size = next_batch_size();
        return std::make_tuple(beg, beg + size);
    }
    void stream_next_batch(bool with_update = false) {
        if (finish_) {
            throw std::runtime_error("Error: Stream when no more batches exist.");
        }
        auto beg = delta_.begin() + offset_;
        int size = next_batch_size();
        stream_graph_.stream_edge(beg, beg + size);
        if (with_update) {
            stream_graph_.update();
        }
        offset_ += size;
        if (offset_ == delta_.size()) {
            finish_ = true;
        }
    }
    void update_graph() {
        stream_graph_.update();
    }
#if defined (BUILD_WITH_RESTORE)
    void restore_graph() {
        stream_graph_.restore();
        finish_ = false;
        offset_ = 0;
    }
#endif
private:
    StreamGraph<T, DstT> &stream_graph_;
    Rng delta_;
    int offset_;
    int max_batch_size_;
    bool finish_;
};

}