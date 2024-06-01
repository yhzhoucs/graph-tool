#pragma once

#include "graph.h"
#include "atomics.h"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <ranges>
#include <stdexcept>

namespace gtool {

template<typename T, typename DstT = T>
class Builder {
public:
    typedef std::vector<std::pair<T, DstT>> EdgeList;
    template<std::convertible_to<std::string> Str>
    explicit Builder(Str &&graph_file, bool directed = false, bool symmetrize = false)
        : graph_file_{ std::forward<Str>(graph_file) }, symmetrize_{ symmetrize && directed },
        directed_{ directed } {}
    virtual Graph<T, DstT> build_csr() {
        EdgeList el = read_edge_list();
        return build_csr_from_edge_list(el);
    }
protected:
    EdgeList read_edge_list() {
        std::ifstream in(this->graph_file_, std::ios::in);
        if (!in.is_open()) {
            throw std::runtime_error("Error! Opening graph file failed.");
        }
        std::vector<std::pair<T, T>> el;
        el.reserve(1000000);
        std::string line;

        while (std::getline(in, line)) {
            if (line[0] == '#')
                continue;   // skip commented lines
            std::istringstream iss(line);
            T u, v;
            iss >> u >> v;
            el.emplace_back(u, DstT{ v });
            if (symmetrize_ && u != v)
                el.emplace_back(v, DstT{ u });
        }
        return el;
    }
    template<std::ranges::random_access_range Rng>
    requires std::is_same_v<typename Rng::value_type, typename EdgeList::value_type>
    Graph<T, DstT> build_csr_from_edge_list(Rng const &el) {
        typedef typename Graph<T, DstT>::offset_t offset_t;
        T max_idx{};
#pragma omp parallel for default(none) shared(el) reduction(max : max_idx)
        for (int i = 0; i < el.size(); ++i) {
            max_idx = std::max({ max_idx, el[i].first, get_dst_id(el[i].second) });
        }
        T vertex_number = max_idx + 1;
        std::vector<offset_t> out_degrees(vertex_number, 0);
#pragma omp parallel for default(none) shared(el, out_degrees)
        for (int i = 0; i < el.size(); ++i) {
            auto const& edge = el[i];
            fetch_and_add(out_degrees[edge.first], 1);
        }
        auto *out_offset = new offset_t[vertex_number + 1];
        offset_t curr{};
        for (int i = 0; i < vertex_number; ++i) {
            out_offset[i] = curr;
            curr += out_degrees[i];
        }
        out_offset[vertex_number] = curr;
        int64_t edge_number = out_offset[vertex_number];
        auto *out_neigh = new DstT[edge_number];
        auto *tmp = new offset_t[vertex_number + 1];
        std::copy(out_offset, out_offset + (vertex_number + 1), tmp);
#pragma omp parallel for default(none) shared(el, tmp, out_neigh)
        for (int i = 0; i < el.size(); ++i) {
            auto const& edge = el[i];
            int write_pos = fetch_and_add(tmp[edge.first], 1);
            out_neigh[write_pos] = edge.second;
        }
#pragma omp parallel for default(none) shared(vertex_number, out_neigh, out_offset)
        for (int i = 0; i < vertex_number; ++i) {
            std::sort(&out_neigh[out_offset[i]], &out_neigh[out_offset[i + 1]], [](DstT const& lhs, DstT const& rhs) {
                return get_dst_id(lhs) < get_dst_id(rhs);
                });
        }
        if (!directed_ || symmetrize_) {
            // if the input is undirected or the graph is symmetrized
            delete[] tmp;
            return { vertex_number, out_offset, out_neigh };
        }
        std::vector<T> in_degrees(vertex_number, 0);
#pragma omp parallel for default(none) shared(el, in_degrees)
        for (int i = 0; i < el.size(); ++i) {
            auto const& edge = el[i];
            fetch_and_add(in_degrees[get_dst_id(edge.second)], 1);
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
#pragma omp parallel for default(none) shared(el, tmp, in_neigh)
        for (int i = 0; i < el.size(); ++i) {
            auto const& edge = el[i];
            int write_pos = fetch_and_add(tmp[get_dst_id(edge.second)], 1);
            get_dst_id(in_neigh[write_pos]) = edge.first;
        }
#pragma omp parallel for default(none) shared(vertex_number, in_neigh, in_offset)
        for (int i = 0; i < vertex_number; ++i) {
            std::sort(&in_neigh[in_offset[i]], &in_neigh[in_offset[i + 1]], [](DstT const& lhs, DstT const& rhs) {
                return get_dst_id(lhs) < get_dst_id(rhs);
                });
        }
        delete[] tmp;
        return { vertex_number, out_offset, out_neigh, in_offset, in_neigh };
    }
    std::string graph_file_;
    bool symmetrize_;
    bool directed_;
};

}