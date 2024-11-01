#pragma once

#include "graph.h"
#include "atomics.h"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <limits>

namespace gtool {

template<typename T, typename DstT>
using EdgeList = std::vector<std::pair<T, DstT>>;

template<typename T, typename DstT>
EdgeList<T, DstT> read_edge_list_common(std::string const &file) {
    std::ifstream in(file, std::ios::in);
    if (!in.is_open()) {
        throw std::runtime_error("Error! Opening graph file failed.");
    }
    EdgeList<T, DstT> el;
    el.reserve(1000000);
    std::string line;

    while (std::getline(in, line)) {
        if (line[0] == '#')
            continue;   // skip commented lines
        std::istringstream iss(line);
        T u, v;
        iss >> u >> v;
        el.emplace_back(u, DstT{v});
    }
    in.close();
    return el;
}

template<typename T, typename DstT>
EdgeList<T, DstT> read_edge_list_mtx(std::string const &file) {
    std::ifstream in(file, std::ios::in);
    if (!in.is_open()) {
        throw std::runtime_error("Error! Opening graph file failed.");
    }
    EdgeList<T, DstT> el;
    el.reserve(1000000);
    std::string line;

    std::getline(in, line);
    std::istringstream iss(line);
    // %%MatrixMarket matrix coordinate pattern general
    for (int i = 0; i < 3; ++i) {
        iss.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
    }
    std::string value_type;
    std::string format_type;
    iss >> value_type >> format_type;
    if (value_type != "pattern") {
        std::cerr << "Current Matrix Market format is not supported." << std::endl;
        return {};
    }
    int symmetrize = false;
    if (format_type == "symmetric") {
        symmetrize = true;
    } else if (format_type == "general") {
        symmetrize = symmetrize;
    } else {
        std::cerr << "Current Matrix Market format is not supported." << std::endl;
        return {};
    }

    for ( ; !in.eof() && line[0] == '%' ; std::getline(in, line));

    while (std::getline(in, line)) {
        iss.clear();
        iss.str(line);
        T u, v;
        iss >> u >> v;
        --u; --v;
        el.emplace_back(u, DstT{v});
        if (symmetrize && u != v)
            el.emplace_back(v, DstT{u});
    }
    in.close();
    return el;
}

template<typename T, typename DstT>
EdgeList<T, DstT> read_edge_list(std::string const &file) {
    if (file.substr(file.length() - 3) == "mtx") {
        return read_edge_list_mtx<T, DstT>(file);
    } else {
        return read_edge_list_common<T, DstT>(file);
    }
}

template<typename T, typename DstT>
EdgeList<T, DstT> symmetrize_edgelist(EdgeList<T, DstT> const &el) {
    EdgeList<T, DstT> new_el;
    new_el.reserve(el.size() * 2);
    for (auto const &pair : el) {
        T u; DstT v;
        std::tie(u, v) = pair;
        new_el.emplace_back(u, v);
        if (u != v) {
            new_el.emplace_back(v, u);
        }
    }
    return new_el;
}

template<typename T, typename DstT = T, typename Rng>
Graph<T, DstT> build_graph_from_edgelist(Rng&& el)  {
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
#pragma omp parallel for default(none) shared(el, tmp, out_neigh)
    for (int i = 0; i < el.size(); ++i) {
        auto const &edge = el[i];
        int write_pos = fetch_and_add(tmp[edge.first], 1);
        get_dst_id(out_neigh.get()[write_pos]) = edge.second;
    }
    std::vector<T> in_degrees(vertex_number, 0);
#pragma omp parallel for default(none) shared(el, in_degrees)
    for (int i = 0; i < el.size(); ++i) {
        auto const &edge = el[i];
        fetch_and_add(in_degrees[get_dst_id(edge.second)], 1);
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
#pragma omp parallel for default(none) shared(el, tmp, in_neigh)
    for (int i = 0; i < el.size(); ++i) {
        auto const &edge = el[i];
        int write_pos = fetch_and_add(tmp[get_dst_id(edge.second)], 1);
        get_dst_id(in_neigh.get()[write_pos]) = edge.first;
    }
    delete[] tmp;
    return {vertex_number, out_offset, out_neigh, in_offset, in_neigh};
}

template<typename T, typename DstT = T>
Graph<T, DstT> build_graph_from_file(std::string const &graph_file, bool symmetrize = false) {
    auto el = read_edge_list<T, DstT>(graph_file);
    if (symmetrize) {
        return build_graph_from_edgelist<T, DstT>(symmetrize_edgelist(el));
    } else {
        return build_graph_from_edgelist<T, DstT>(el);
    }
}

}