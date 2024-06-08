#pragma once

#include <type_traits>
#include <utility>
#include <vector>
#include <algorithm>
#include <tuple>
#include <memory>

namespace gtool {

typedef std::make_unsigned<std::ptrdiff_t>::type offset_t;

template<typename T, typename=std::enable_if<std::is_integral<T>::value>::type>
T &get_dst_id(T &dst) {
    return dst;
}

template<typename T, typename DstT = T>
class Graph {
public:
    typedef T vertex_type;
    typedef DstT destination_type;
    struct Neighborhood {
        T n;
        offset_t *offset;
        T *neigh;

        typedef DstT const *iterator;
        iterator begin() { return neigh + offset[n]; }
        iterator end() { return neigh + offset[n + 1]; }
    };

    Graph(): directed_{}, vertex_number_{}, edge_number_{}, out_offset_{},
             out_neigh_{}, in_offset_{}, in_neigh_{} {};
    Graph(std::size_t vertex_number, std::shared_ptr<offset_t> out_offset, std::shared_ptr<DstT> out_neigh)
        : vertex_number_{vertex_number}, out_offset_(std::move(out_offset)), out_neigh_(std::move(out_neigh)) {
        in_offset_ = out_offset_;
        in_neigh_ = out_neigh_;
        edge_number_ = out_offset_.get()[vertex_number_];
    }
    Graph(std::size_t vertex_number, std::shared_ptr<offset_t> out_offset, std::shared_ptr<DstT> out_neigh,
          std::shared_ptr<offset_t> in_offset, std::shared_ptr<DstT> in_neigh)
          : vertex_number_{vertex_number}, out_offset_(std::move(out_offset)), out_neigh_(std::move(out_neigh)),
          in_offset_(std::move(in_offset)), in_neigh_(std::move(in_neigh)) {
        edge_number_ = out_offset_.get()[vertex_number_];
    }
    ~Graph() = default;
    [[nodiscard]] bool directed() const { return directed_; }
    [[nodiscard]] std::size_t get_vertex_number() const { return vertex_number_; }
    [[nodiscard]] std::size_t get_edge_number() const { return edge_number_; }
    offset_t out_degree(T n) const { return out_offset_.get()[n + 1] - out_offset_.get()[n]; }
    offset_t in_degree(T n) const { return in_offset_.get()[n + 1] - in_offset_.get()[n]; }
    Neighborhood out_neighbors(T n) const { return {n, out_offset_.get(), out_neigh_.get()}; }
    Neighborhood in_neighbors(T n) const { return {n, in_offset_.get(), in_neigh_.get()}; }
    template<typename Comp>
    void sort_neighborhood(Comp comp);
    template<typename TT, typename DstTT>
    friend Graph<TT, DstTT> simplify_graph(Graph<TT, DstTT> &raw);
private:
    bool directed_;
    std::size_t vertex_number_;
    std::size_t edge_number_;
    std::shared_ptr<offset_t> out_offset_;
    std::shared_ptr<DstT> out_neigh_;
    std::shared_ptr<offset_t> in_offset_;
    std::shared_ptr<DstT> in_neigh_;
};

template<typename T, typename DstT>
template<typename Comp>
void Graph<T, DstT>::sort_neighborhood(Comp comp) {
    for (T u = 0; u < vertex_number_; ++u) {
        std::sort(&out_neigh_.get()[out_offset_.get()[u]],
                  &out_neigh_.get()[out_offset_.get()[u + 1]],
                  [&](T const &lhs, T const &rhs) { return comp(out_degree(lhs), out_degree(rhs)); });
    }
    if (directed_) {
        for (T u = 0; u < vertex_number_; ++u) {
            std::sort(&in_neigh_.get()[in_offset_.get()[u]],
                      &in_neigh_.get()[in_offset_.get()[u + 1]],
                      [&](T const &lhs, T const &rhs) { return comp(out_degree(lhs), out_degree(rhs)); });
        }
    }
}

template<typename T, typename DstT>
std::tuple<Graph<T, DstT>, std::vector<T>, std::vector<T>> reorder_by_degree(Graph<T, DstT> const &g) {
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
    std::shared_ptr<offset_t> out_offset(new offset_t[vertex_number + 1], std::default_delete<offset_t[]>());
    offset_t curr{};
    for (int i = 0; i < vertex_number; ++i) {
        out_offset.get()[i] = curr;
        curr += out_degrees[i];
    }
    out_offset.get()[vertex_number] = curr;
    std::size_t edge_number = curr;
    std::shared_ptr<DstT> out_neigh(new DstT[edge_number], std::default_delete<DstT[]>());
    auto *tmp = new offset_t[vertex_number + 1];
    std::copy(out_offset.get(), out_offset.get() + (vertex_number + 1), tmp);
    for (int u = 0; u < vertex_number; ++u) {
        for (DstT const &v : g.out_neighbors(u)) {
            out_neigh.get()[tmp[new_ids[u]]] = new_ids[get_dst_id(v)];
            tmp[new_ids[u]]++;
        }
        std::sort(&out_neigh.get()[out_offset.get()[new_ids[u]]],
                  &out_neigh.get()[out_offset.get()[new_ids[u]+1]],
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
    std::shared_ptr<offset_t> in_offset(new offset_t[vertex_number + 1], std::default_delete<offset_t[]>());
    curr = 0;
    for (int i = 0; i < vertex_number; ++i) {
        in_offset.get()[i] = curr;
        curr += in_degrees[i];
    }
    in_offset.get()[vertex_number] = curr;
    std::shared_ptr<DstT> in_neigh(new DstT[edge_number], std::default_delete<DstT[]>());
    std::copy(in_offset.get(), in_offset.get() + (vertex_number + 1), tmp);
    for (int u = 0; u < vertex_number; ++u) {
        for (DstT const &v : g.in_neighbors(u)) {
            in_neigh.get()[tmp[new_ids[u]]] = new_ids[get_dst_id(v)];
            tmp[new_ids[u]]++;
        }
        std::sort(&in_neigh.get()[in_offset.get()[new_ids[u]]],
                  &in_neigh.get()[in_offset.get()[new_ids[u]+1]],
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
    std::size_t squeezed_vertex_number = mapped_src;
    std::vector<offset_t> out_degrees(squeezed_vertex_number, 0);
    for (T u = 0; u < squeezed_vertex_number; ++u) {
        out_degrees[u] = g.out_degree(vertex_remap[u]);
    }
    std::shared_ptr<offset_t> out_offset(new offset_t[mapped_src + 1], std::default_delete<offset_t[]>());
    offset_t curr{};
    for (T u = 0; u < squeezed_vertex_number; ++u) {
        out_offset.get()[u] = curr;
        curr += out_degrees[u];
    }
    out_offset.get()[squeezed_vertex_number] = curr;
    std::size_t edge_number = curr;
    std::shared_ptr<DstT> out_neigh(new DstT[edge_number], std::default_delete<DstT[]>());
    auto *tmp = new offset_t[squeezed_vertex_number + 1];
    std::copy(out_offset.get(), out_offset.get() + (squeezed_vertex_number + 1), tmp);
    for (T u = 0; u < squeezed_vertex_number ; ++u) {
        for (DstT const &v : g.out_neighbors(vertex_remap[u])) {
            out_neigh.get()[tmp[u]] = vertex_map[get_dst_id(v)];
            tmp[u]++;
        }
        std::sort(&out_neigh.get()[out_offset.get()[u]],
                  &out_neigh.get()[out_offset.get()[u+1]],
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
    std::shared_ptr<offset_t> in_offset(new offset_t[squeezed_vertex_number + 1], std::default_delete<offset_t[]>());
    curr = 0;
    for (T u = 0; u < squeezed_vertex_number; ++u) {
        in_offset.get()[u] = curr;
        curr += in_degrees[u];
    }
    in_offset.get()[squeezed_vertex_number] = curr;
    std::shared_ptr<DstT> in_neigh(new DstT[edge_number], std::default_delete<DstT[]>());
    std::copy(in_offset.get(), in_offset.get() + (squeezed_vertex_number + 1), tmp);
    for (T u = 0; u < squeezed_vertex_number; ++u) {
        for (DstT const &v : g.in_neighbors(vertex_remap[u])) {
            in_neigh.get()[tmp[u]] = vertex_map[get_dst_id(v)];
            tmp[u]++;
        }
        std::sort(&in_neigh.get()[in_offset.get()[u]],
                  &in_neigh.get()[in_offset.get()[u+1]],
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
    std::size_t vertex_number = raw.vertex_number_;
    std::vector<offset_t> out_degrees(vertex_number, 0);
    for (T u = 0; u < vertex_number; ++u) {
        out_degrees[u] = std::distance(&raw.out_neigh_.get()[raw.out_offset_.get()[u]],
                                       std::unique(&raw.out_neigh_.get()[raw.out_offset_.get()[u]], &raw.out_neigh_.get()[raw.out_offset_.get()[u + 1]],
                                                   [](DstT const &lhs, DstT const &rhs) { return get_dst_id(lhs) == get_dst_id(rhs); }));
    }
    std::shared_ptr<offset_t> out_offset(new offset_t[vertex_number + 1], std::default_delete<offset_t[]>());
    offset_t curr{};
    for (size_t i = 0; i < vertex_number; ++i) {
        out_offset.get()[i] = curr;
        curr += out_degrees[i];
    }
    out_offset.get()[vertex_number] = curr;
    std::size_t edge_number = curr;
    std::shared_ptr<DstT> out_neigh(new DstT[edge_number], std::default_delete<DstT[]>());
    for (T u = 0; u < vertex_number; ++u) {
        std::copy(&raw.out_neigh_.get()[raw.out_offset_.get()[u]], &raw.out_neigh_.get()[raw.out_offset_.get()[u] + out_degrees[u]], &out_neigh.get()[out_offset.get()[u]]);
    }
    if (!raw.directed_) {
        return {vertex_number, out_offset, out_neigh};
    }
    std::vector<offset_t> in_degrees(vertex_number, 0);
    for (T u = 0; u < vertex_number; ++u) {
        in_degrees[u] = std::distance(&raw.in_neigh_.get()[raw.in_offset_.get()[u]],
                                      std::unique(&raw.in_neigh_.get()[raw.in_offset_.get()[u]], &raw.in_neigh_.get()[raw.in_offset_.get()[u + 1]],
                                                  [](DstT const &lhs, DstT const &rhs) { return get_dst_id(lhs) == get_dst_id(rhs); }));
    }
    std::shared_ptr<offset_t> in_offset(new offset_t[vertex_number + 1], std::default_delete<offset_t[]>());
    curr = 0;
    for (size_t i = 0; i < vertex_number; ++i) {
        in_offset.get()[i] = curr;
        curr += in_degrees[i];
    }
    in_offset.get()[vertex_number] = curr;
    std::shared_ptr<DstT> in_neigh(new DstT[edge_number], std::default_delete<DstT[]>());
    for (T u = 0; u < vertex_number; ++u) {
        std::copy(&raw.in_neigh_.get()[raw.in_offset_.get()[u]], &raw.in_neigh_.get()[raw.in_offset_.get()[u] + in_degrees[u]], &in_neigh.get()[in_offset.get()[u]]);
    }
    return {vertex_number, out_offset, out_neigh, in_offset, in_neigh};
}

}