#pragma once

#include "builder.h"
#include "stream_graph.h"
#include <random>
#include <type_traits>
#include <unordered_set>
#include <utility>

namespace gtool {

template<typename T>
class FixedRange {
public:
    class iterator {
    public:
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = T*;
        using reference = T&;
        explicit iterator(T &val, long num = 0) : val_(val), num_(num) {}
        iterator(iterator const &other): val_(other.val_), num_(other.num_) {}
        iterator(iterator &&other) noexcept: val_(other.val_), num_(other.num_) {}
        iterator &operator=(iterator const &other) {
            if (this == &other) {
                return *this;
            }
            val_ = other.val_;
            num_ = other.num_;
            return *this;
        }
        iterator &operator=(iterator &&other) noexcept {
            if (this == &other) {
                return *this;
            }
            val_ = other.val_;
            num_ = other.num_;
            return *this;
        }
        iterator& operator++() { ++num_; return *this; }
        iterator operator++(int) { iterator old = *this; ++num_; return old; }
        bool operator==(iterator other) const { return num_ == other.num_; }
        bool operator!=(iterator other) const { return !(*this == other); }
        reference operator*() const { return val_; }
    private:
        T &val_;
        long num_;
    };
    explicit FixedRange(T val): val_(val) {}
    iterator begin() { return iterator(val_); }
    iterator end() { return iterator(val_, std::numeric_limits<long>::max()); }
private:
    T val_;
};

template<typename T, typename DstT = T>
class StreamBuilder {
public:
    explicit StreamBuilder(std::string file_name, double stream_ratio = 0.01) :
            graph_file_(std::move(file_name)), stream_ratio_{stream_ratio}, delta_{} {}
    [[nodiscard]] EdgeList<T, DstT> const &delta() const { return delta_; }
    [[nodiscard]] double stream_ratio() const { return stream_ratio_; }
    StreamGraph<T, DstT> build_csr() {
        return build_base_csr(FixedRange<int>(0), FixedRange<int>(0));
    }
    template<typename Rng, typename=std::enable_if<IsRangeT<Rng>::value>::type>
    StreamGraph<T, DstT> build_csr(Rng &&precise_capacity) {
        return build_base_csr(std::forward<Rng>(precise_capacity), std::forward<Rng>(precise_capacity));
    }
    template<typename Rng, typename=std::enable_if<IsRangeT<Rng>::value>::type>
    StreamGraph<T, DstT> build_csr(Rng &&precise_capacity_out, Rng &&precise_capacity_in) {
        return build_base_csr(std::forward<Rng>(precise_capacity_out), std::forward<Rng>(precise_capacity_in));
    }
    StreamGraph<T, DstT> build_csr(int fixed_capacity) {
        return build_base_csr(FixedRange<int>(fixed_capacity), FixedRange<int>(fixed_capacity));
    }
    StreamGraph<T, DstT> build_csr(int fixed_capacity_out, int fixed_capacity_in) {
        return build_base_csr(FixedRange<int>(fixed_capacity_out), FixedRange<int>(fixed_capacity_in));
    }
private:
    template<typename Rng, typename=std::enable_if<IsRangeT<Rng>::value>::type>
    StreamGraph<T, DstT> build_base_csr(Rng &&capacity_out, Rng &&capacity_in) {
        // treated Rng as an input range
        EdgeList<T, DstT> el = read_edge_list<T, DstT>(graph_file_);
        // build offsets with full graph
        T max_idx{};
#pragma omp parallel for default(none) shared(el) reduction(max : max_idx)
        for (int i = 0; i < el.size(); ++i) {
            max_idx = std::max({max_idx, el[i].first, get_dst_id(el[i].second)});
        }
        std::size_t vertex_number = max_idx + 1;
        std::vector<offset_t> out_degrees(vertex_number, 0);
#pragma omp parallel for default(none) shared(el, out_degrees)
        for (int i = 0; i < el.size(); ++i) {
            auto const &edge = el[i];
            fetch_and_add(out_degrees[edge.first], 1);
        }
        std::unique_ptr<offset_t[]> out_l_offset(new offset_t[vertex_number + 1]);
        offset_t curr{};
        auto capacity_out_it = std::begin(capacity_out);
        for (int i = 0; i < vertex_number; ++i, ++capacity_out_it) {
            out_l_offset[i] = curr;
            curr += std::max(out_degrees[i], static_cast<offset_t>(*capacity_out_it));
        }
        out_l_offset[vertex_number] = curr;
        std::size_t edge_number = curr;
        std::vector<offset_t> in_degrees(vertex_number, 0);
#pragma omp parallel for default(none) shared(el, in_degrees)
        for (int i = 0; i < el.size(); ++i) {
            auto const &edge = el[i];
            fetch_and_add(in_degrees[get_dst_id(edge.second)], 1);
        }
        std::unique_ptr<offset_t[]> in_l_offset(new offset_t[vertex_number + 1]);
        curr = 0;
        auto capacity_in_it = std::begin(capacity_in);
        for (int i = 0; i < vertex_number; ++i, ++capacity_in_it) {
            in_l_offset[i] = curr;
            curr += std::max(in_degrees[i], static_cast<offset_t>(*capacity_in_it));
        }
        in_l_offset[vertex_number] = curr;

        // build basic neighbors
        extract_delta(el);
        std::size_t remain = std::distance(el.begin(),
            std::remove_if(el.begin(), el.end(), [](auto const &pair) { return pair.first == -1; }));
        std::unique_ptr<offset_t[]> out_r_offset(new offset_t[vertex_number + 1]);
        std::memcpy(out_r_offset.get(), out_l_offset.get(), sizeof(offset_t) * (vertex_number + 1));
        std::unique_ptr<DstT[]> out_neigh(new DstT[edge_number]);
#pragma omp parallel for default(none) shared(el, out_r_offset, out_neigh, remain)
        for (int i = 0; i < remain; ++i) {
            auto const &edge = el[i];
            int write_pos = fetch_and_add(out_r_offset[edge.first], 1);
            out_neigh[write_pos] = edge.second;
        }
        std::unique_ptr<offset_t[]> in_r_offset(new offset_t[vertex_number + 1]);
        std::memcpy(in_r_offset.get(), in_l_offset.get(), sizeof(offset_t) * (vertex_number + 1));
        std::unique_ptr<DstT[]> in_neigh(new DstT[edge_number]);
#pragma omp parallel for default(none) shared(el, in_r_offset, in_neigh, remain)
        for (int i = 0; i < remain; ++i) {
            auto const &edge = el[i];
            int write_pos = fetch_and_add(in_r_offset[get_dst_id(edge.second)], 1);
            in_neigh[write_pos] = edge.first;
        }
        // construct
        return { vertex_number, edge_number, std::move(out_l_offset), std::move(out_r_offset),
                 std::move(out_neigh), std::move(in_l_offset), std::move(in_r_offset), std::move(in_neigh) };
    }
    // Extract delta edges from raw edge list, push them to `delta_`
    // and set both (src, dst) of these edge to (-1, dst) in raw edge list.
    void extract_delta(EdgeList<T, DstT> &el) {
        int delta_size = static_cast<int>(el.size() * stream_ratio_);
        delta_.reserve(delta_size);
        std::unordered_set<int> set;
        std::random_device rd{};
        std::mt19937 rng{rd()};
        std::uniform_int_distribution<int> dist(0, el.size() - 1);
        int cnt{};
        while (cnt < delta_size) {
            int idx = dist(rng);
            if (!set.contains(idx)) {
                delta_.emplace_back(el[idx]);
                el[idx].first = -1;
                set.emplace(idx);
                ++cnt;
            }
        }
    }
    std::string graph_file_;
    double stream_ratio_;
    EdgeList<T, DstT> delta_;
};

}
