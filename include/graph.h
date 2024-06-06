#pragma once

#include <type_traits>
#include <concepts>
#include <utility>
#include <vector>
#include <algorithm>
#include <tuple>

namespace gtool {

auto &get_dst_id(std::integral auto &dst) {
    return dst;
}

template<typename T, typename DstT = T>
class Graph {
public:
    typedef std::make_unsigned<std::ptrdiff_t>::type offset_t;
    typedef T vertex_type;
    typedef DstT destination_type;
private:
    bool directed_;
    std::size_t vertex_number_;
    std::size_t edge_number_;
    offset_t *out_offset;
    DstT *out_neigh;
    offset_t *in_offset;
    DstT *in_neigh;

    struct Neighborhood {
        T n;
        offset_t *offset;
        T *neigh;

        typedef DstT const *iterator;
        iterator begin() { return neigh + offset[n]; }
        iterator end() { return neigh + offset[n + 1]; }
    };
public:
    Graph(): directed_{}, vertex_number_{}, edge_number_{}, out_offset{nullptr},
        out_neigh{nullptr}, in_offset{nullptr}, in_neigh{nullptr} {};
    Graph(std::size_t vertex_number, offset_t *&out_offset, DstT *&out_neigh);
    Graph(std::size_t vertex_number, offset_t *&out_offset, DstT *&out_neigh,
          offset_t *&in_offset, DstT *&in_neigh);
    Graph(Graph<T, DstT> const &graph) = delete;
    Graph(Graph<T, DstT> &&graph) noexcept;
    ~Graph();

    Graph<T, DstT> &operator=(Graph<T, DstT> const &other) = delete;
    Graph<T, DstT> &operator=(Graph<T, DstT> &&other) noexcept;

    [[nodiscard]] bool directed() const { return directed_; }
    [[nodiscard]] std::size_t get_vertex_number() const { return vertex_number_; }
    [[nodiscard]] std::size_t get_edge_number() const { return edge_number_; }
    [[nodiscard]] offset_t const *get_out_offset() const { return out_offset; }
    [[nodiscard]] offset_t const *get_in_offset() const { return in_offset; }
    [[nodiscard]] DstT const *get_out_neigh() const { return out_neigh; }
    [[nodiscard]] DstT const *get_in_neigh() const { return in_neigh; }
    offset_t out_degree(T n) const { return out_offset[n + 1] - out_offset[n]; }
    offset_t in_degree(T n) const { return in_offset[n + 1] - in_offset[n]; }
    Neighborhood out_neighbors(T n) const { return {n, out_offset, out_neigh}; }
    Neighborhood in_neighbors(T n) const { return {n, in_offset, in_neigh}; }
    template<typename Comp>
    void sort_neighborhood(Comp comp);
    template<typename TT, typename DstTT>
    friend Graph<TT, DstTT> simplify_graph(Graph<TT, DstTT> &raw);
};

template<typename T, typename DstT>
std::tuple<Graph<T, DstT>, std::vector<T>, std::vector<T>> reorder_by_degree(Graph<T, DstT> const &g) {
    typedef typename Graph<T, DstT>::offset_t offset_t;
    typedef std::pair<offset_t, T> DegreeNVertex;
    std::vector<DegreeNVertex> degree_vertex_pairs;
    for (int i = 0; i < g.get_vertex_number(); ++i) {
        degree_vertex_pairs.emplace_back(g.out_degree(i), i);
    }
    std::sort(degree_vertex_pairs.begin(), degree_vertex_pairs.end(), std::greater<DegreeNVertex>());
    std::size_t vertex_number = g.get_vertex_number();
    std::vector<T> out_degrees(vertex_number, 0);
    std::vector<T> new_ids(vertex_number, 0);
    std::vector<T> new_ids_remap(vertex_number, 0);
    for (int i = 0; i < vertex_number; ++i) {
        out_degrees[i] = degree_vertex_pairs[i].first;
        new_ids[degree_vertex_pairs[i].second] = i;
        new_ids_remap[i] = degree_vertex_pairs[i].second;
    }
    auto *out_offset = new offset_t[vertex_number + 1];
    offset_t curr{};
    for (int i = 0; i < vertex_number; ++i) {
        out_offset[i] = curr;
        curr += out_degrees[i];
    }
    out_offset[vertex_number] = curr;
    std::size_t edge_number = curr;
    auto *out_neigh = new DstT[edge_number];
    auto *tmp = new offset_t[vertex_number + 1];
    std::copy(out_offset, out_offset + (vertex_number + 1), tmp);
    for (int u = 0; u < vertex_number; ++u) {
        for (DstT const &v : g.out_neighbors(u)) {
            out_neigh[tmp[new_ids[u]]] = new_ids[get_dst_id(v)];
            tmp[new_ids[u]]++;
        }
        std::sort(&out_neigh[out_offset[new_ids[u]]],
                  &out_neigh[out_offset[new_ids[u]+1]],
                  [](DstT const &lhs, DstT const &rhs) { return get_dst_id(lhs) < get_dst_id(rhs); });
    }
    if (!g.directed()) {
        delete[] tmp;
        return {Graph<T, DstT>{vertex_number, out_offset, out_neigh}, std::move(new_ids), std::move(new_ids_remap)};
    }
    std::vector<T> in_degrees(vertex_number, 0);
    for (int i = 0; i < vertex_number; ++i) {
        in_degrees[new_ids[i]] = g.in_degree(i);
    }
    auto *in_offset = new offset_t[vertex_number + 1];
    curr = 0;
    for (int i = 0; i < vertex_number; ++i) {
        in_offset[i] = curr;
        curr += in_degrees[i];
    }
    in_offset[vertex_number] = curr;
    auto *in_neigh = new DstT[edge_number];
    std::copy(in_offset, in_offset + (vertex_number + 1), tmp);
    for (int u = 0; u < vertex_number; ++u) {
        for (DstT const &v : g.in_neighbors(u)) {
            in_neigh[tmp[new_ids[u]]] = new_ids[get_dst_id(v)];
            tmp[new_ids[u]]++;
        }
        std::sort(&in_neigh[in_offset[new_ids[u]]],
                  &in_neigh[in_offset[new_ids[u]+1]],
                  [](DstT const &lhs, DstT const &rhs) { return get_dst_id(lhs) < get_dst_id(rhs); });
    }
    delete[] tmp;
    return {Graph<T, DstT>{vertex_number, out_offset, out_neigh, in_offset, in_neigh},
            std::move(new_ids),
            std::move(new_ids_remap)};
}

template<typename T, typename DstT>
std::tuple<Graph<T, DstT>, std::vector<T>, std::vector<T>> squeeze_graph(Graph<T, DstT> const &g) {
    std::size_t mapped_src{};
    std::vector<T> vertex_map(g.get_vertex_number(), 0);
    std::vector<T> vertex_remap(g.get_vertex_number(), 0);
    for (int u = 0; u < g.get_vertex_number(); ++u) {
        vertex_map[u] = mapped_src;
        if (g.out_degree(u) != 0 || g.in_degree(u) != 0) {
            vertex_remap[mapped_src] = u;
            mapped_src++;
        }
    }
    typedef typename Graph<T, DstT>::offset_t offset_t;
    std::size_t squeezed_vertex_number = mapped_src;
    std::vector<offset_t> out_degrees(squeezed_vertex_number, 0);
    for (T u = 0; u < squeezed_vertex_number; ++u) {
        out_degrees[u] = g.out_degree(vertex_remap[u]);
    }
    auto *out_offset = new offset_t[mapped_src + 1];
    offset_t curr{};
    for (T u = 0; u < squeezed_vertex_number; ++u) {
        out_offset[u] = curr;
        curr += out_degrees[u];
    }
    out_offset[squeezed_vertex_number] = curr;
    std::size_t edge_number = curr;
    auto *out_neigh = new DstT[edge_number];
    auto *tmp = new offset_t[squeezed_vertex_number + 1];
    std::copy(out_offset, out_offset + (squeezed_vertex_number + 1), tmp);
    for (T u = 0; u < squeezed_vertex_number ; ++u) {
        for (DstT const &v : g.out_neighbors(vertex_remap[u])) {
            out_neigh[tmp[u]] = vertex_map[get_dst_id(v)];
            tmp[u]++;
        }
        std::sort(&out_neigh[out_offset[u]],
                  &out_neigh[out_offset[u+1]],
                  [](DstT const &lhs, DstT const &rhs) { return get_dst_id(lhs) < get_dst_id(rhs); });
    }
    if (!g.directed()) {
        delete[] tmp;
        return {
            Graph<T, DstT>{squeezed_vertex_number, out_offset, out_neigh},
            std::move(vertex_map),
            std::move(vertex_remap)
        };
    }
    std::vector<offset_t> in_degrees(std::move(out_degrees));
    for (T u = 0; u < squeezed_vertex_number; ++u) {
        in_degrees[u] = g.in_degree(vertex_remap[u]);
    }
    auto *in_offset = new offset_t[squeezed_vertex_number + 1];
    curr = 0;
    for (T u = 0; u < squeezed_vertex_number; ++u) {
        in_offset[u] = curr;
        curr += in_degrees[u];
    }
    in_offset[squeezed_vertex_number] = curr;
    auto *in_neigh = new DstT[edge_number];
    std::copy(in_offset, in_offset + (squeezed_vertex_number + 1), tmp);
    for (T u = 0; u < squeezed_vertex_number; ++u) {
        for (DstT const &v : g.in_neighbors(vertex_remap[u])) {
            in_neigh[tmp[u]] = vertex_map[get_dst_id(v)];
            tmp[u]++;
        }
        std::sort(&in_neigh[in_offset[u]],
                  &in_neigh[in_offset[u+1]],
                  [](DstT const &lhs, DstT const &rhs) { return get_dst_id(lhs) < get_dst_id(rhs); });
    }
    delete[] tmp;
    return {
            Graph<T, DstT>{squeezed_vertex_number, out_offset, out_neigh, in_offset, in_neigh},
            std::move(vertex_map),
            std::move(vertex_remap)
    };
}

/**
 * WARNING: The function will change the order of edgelist in raw
*/
template<typename T, typename DstT>
Graph<T, DstT> simplify_graph(Graph<T, DstT> &raw) {
    typedef typename Graph<T, DstT>::offset_t offset_t;
    std::size_t vertex_number = raw.vertex_number_;
    std::vector<offset_t> out_degrees(vertex_number, 0);
    for (T u = 0; u < vertex_number; ++u) {
        out_degrees[u] = std::distance(&raw.out_neigh[raw.out_offset[u]], 
            std::unique(&raw.out_neigh[raw.out_offset[u]], &raw.out_neigh[raw.out_offset[u+1]], 
                [](DstT const &lhs, DstT const &rhs) { return get_dst_id(lhs) == get_dst_id(rhs); }));
    }
    auto *out_offset = new offset_t[vertex_number + 1];
    offset_t curr{};
    for (size_t i = 0; i < vertex_number; ++i) {
        out_offset[i] = curr;
        curr += out_degrees[i];
    }
    out_offset[vertex_number] = curr;
    std::size_t edge_number = out_offset[vertex_number];
    auto *out_neigh = new DstT[edge_number];
    for (T u = 0; u < vertex_number; ++u) {
        std::copy(&raw.out_neigh[raw.out_offset[u]], &raw.out_neigh[raw.out_offset[u]+out_degrees[u]], &out_neigh[out_offset[u]]);
    }
    if (!raw.directed_) {
        return {vertex_number, out_offset, out_neigh};
    }
    std::vector<offset_t> in_degrees(vertex_number, 0);
    for (T u = 0; u < vertex_number; ++u) {
        in_degrees[u] = std::distance(&raw.in_neigh[raw.in_offset[u]], 
            std::unique(&raw.in_neigh[raw.in_offset[u]], &raw.in_neigh[raw.in_offset[u+1]], 
                [](DstT const &lhs, DstT const &rhs) { return get_dst_id(lhs) == get_dst_id(rhs); }));
    }
    auto *in_offset = new offset_t[vertex_number + 1];
    curr = 0;
    for (size_t i = 0; i < vertex_number; ++i) {
        in_offset[i] = curr;
        curr += in_degrees[i];
    }
    in_offset[vertex_number] = curr;
    auto *in_neigh = new DstT[edge_number];
    for (T u = 0; u < vertex_number; ++u) {
        std::copy(&raw.in_neigh[raw.in_offset[u]], &raw.in_neigh[raw.in_offset[u]+in_degrees[u]], &in_neigh[in_offset[u]]);
    }
    return {vertex_number, out_offset, out_neigh, in_offset, in_neigh};
}

}

#include "graph.tpp"