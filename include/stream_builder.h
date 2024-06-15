#pragma once

#include "builder.h"
#include "stream_graph.h"
#include <random>
#include <type_traits>
#include <unordered_set>

namespace gtool {

template<typename T, typename DstT = T>
class StreamBuilder {
public:
    template<typename String>
    requires std::is_convertible_v<std::remove_cvref_t<String>, std::string>
    explicit
    StreamBuilder(String &&file_name, double stream_ratio = 0.01) :
            graph_file_(std::forward<String>(file_name)), stream_ratio_{stream_ratio}, delta_{} {}
    [[nodiscard]] EdgeList<T, DstT> const &delta() const { return delta_; }
    [[nodiscard]] double stream_ratio() const { return stream_ratio_; }
    StreamGraph<T, DstT> build_csr() {
        return build_base_csr();
    }
private:
    StreamGraph<T, DstT> build_base_csr() {
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
        for (int i = 0; i < vertex_number; ++i) {
            out_l_offset[i] = curr;
            curr += out_degrees[i];
        }
        out_l_offset[vertex_number] = curr;
        std::size_t edge_number = curr;
        std::vector<T> in_degrees(vertex_number, 0);
#pragma omp parallel for default(none) shared(el, in_degrees)
        for (int i = 0; i < el.size(); ++i) {
            auto const &edge = el[i];
            fetch_and_add(in_degrees[get_dst_id(edge.second)], 1);
        }
        std::unique_ptr<offset_t[]> in_l_offset(new offset_t[vertex_number + 1]);
        curr = 0;
        for (int i = 0; i < vertex_number; ++i) {
            in_l_offset[i] = curr;
            curr += in_degrees[i];
        }
        in_l_offset[vertex_number] = curr;

        // build basic neighbors
        extract_delta(el);
        auto iter = std::remove_if(el.begin(), el.end(), [](auto const &pair) { return pair.first == -1; });
        std::size_t remain = std::distance(el.begin(), iter);
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