#pragma once

namespace gtool {

template<typename T, typename DstT>
Graph<T, DstT>::Graph(std::size_t vertex_number, offset_t *&out_offset, DstT *&out_neigh)
        : directed_{false}, vertex_number_{vertex_number} {
    this->out_offset = std::exchange(out_offset, nullptr);
    this->out_neigh = std::exchange(out_neigh, nullptr);
    this->in_offset = this->out_offset;
    this->in_neigh = this->out_neigh;
    edge_number_ = this->out_offset[vertex_number_] - this->out_offset[0];
}

template<typename T, typename DstT>
Graph<T, DstT>::Graph(std::size_t vertex_number, offset_t *&out_offset, DstT *&out_neigh,
                      offset_t *&in_offset, DstT *&in_neigh)
        : Graph<T, DstT>{vertex_number, out_offset, out_neigh} {
    this->directed_ = true;
    this->in_offset = std::exchange(in_offset, nullptr);
    this->in_neigh = std::exchange(in_neigh, nullptr);
    edge_number_ = this->out_offset[vertex_number_] - this->out_offset[0];
}

template<typename T, typename DstT>
Graph<T, DstT>::Graph(Graph<T, DstT> &&graph) noexcept
        : directed_{ std::move(graph.directed_) },
          vertex_number_{std::move(graph.vertex_number)},
          edge_number_{std::move(graph.edge_number_)} {
    out_offset = std::exchange(graph.out_offset, nullptr);
    out_neigh = std::exchange(graph.out_neigh, nullptr);
    in_offset = std::exchange(graph.in_offset, nullptr);
    in_neigh = std::exchange(graph.in_neigh, nullptr);
}

template<typename T, typename DstT>
Graph<T, DstT>::~Graph() {
    delete[] out_offset;
    delete[] out_neigh;
    if (directed_) {
        delete[] in_offset;
        delete[] in_neigh;
    }
}

template<typename T, typename DstT>
Graph<T, DstT> &Graph<T, DstT>::operator=(Graph<T, DstT> &&other) noexcept {
    if (this == &other)
        return *this;
    delete[] this->out_offset;
    delete[] this->out_neigh;
    if (this->directed_) {
        delete[] this->in_offset;
        delete[] this->in_neigh;
    }
    this->directed_ = std::exchange(other.directed_, false);
    this->vertex_number_ = std::exchange(other.vertex_number_, 0);
    this->edge_number_ = std::exchange(other.edge_number_, 0);
    this->out_offset = std::exchange(other.out_offset, nullptr);
    this->out_neigh = std::exchange(other.out_neigh, nullptr);
    this->in_offset = std::exchange(other.in_offset, nullptr);
    this->in_neigh = std::exchange(other.in_neigh, nullptr);
    return *this;
}

template<typename T, typename DstT>
template<typename Comp>
void Graph<T, DstT>::sort_neighborhood(Comp comp) {
    for (T u = 0; u < this->vertex_number_; ++u) {
        std::sort(&out_neigh[out_offset[u]],
                  &out_neigh[out_offset[u + 1]],
                  [&](T const &lhs, T const &rhs) { return comp(out_degree(lhs), out_degree(rhs)); });
    }
    if (directed_) {
        for (T u = 0; u < this->vertex_number_; ++u) {
            std::sort(&in_neigh[in_offset[u]],
                      &in_neigh[in_offset[u + 1]],
                      [&](T const &lhs, T const &rhs) { return comp(out_degree(lhs), out_degree(rhs)); });
        }
    }
}

}