#pragma once

#include "builder.h"
#include "graph.h"
#include <random>
#include <type_traits>

namespace gtool {

    template<typename T, typename DstT = T>
    class StreamBuilder : public gtool::Builder<T, DstT> {
    public:
        using typename gtool::Builder<T, DstT>::EdgeList;

        template<typename String>
        requires std::is_convertible_v<std::remove_cvref_t<String>, std::string>
        explicit
        StreamBuilder(String &&file_name, bool directed = false, bool symmetrize = false, double stream_ratio = 0.001) :
                gtool::Builder<T, DstT>{file_name, directed, symmetrize}, stream_ratio_{stream_ratio},
                delta_{} {}

        [[nodiscard]] EdgeList const &delta() const { return delta_; }

        [[nodiscard]] double stream_ratio() const { return stream_ratio_; }

        gtool::Graph<T, DstT> build_csr() override {
            return build_base_csr();
        }

    private:
        gtool::Graph<T, DstT> build_base_csr() {
            EdgeList el = this->read_edge_list();
            extract_delta(el);
            auto const [first, last] = std::ranges::remove_if(el, [](auto const &pair) { return pair.first == -1; });
            return this->build_csr_from_edge_list(std::ranges::subrange(el.begin(), first));
        }

        // Extract delta edges from raw edge list, push them to `delta_`
        // and set both (src, dst) of these edge to (-1, dst) in raw edge list.
        void extract_delta(EdgeList &el) {
            int delta_size = static_cast<int>(el.size() * stream_ratio_);
            delta_.reserve(delta_size);
            std::random_device rd{};
            std::mt19937 rng{rd()};
            std::uniform_int_distribution<int> dist(0, el.size() - 1);
            for ([[maybe_unused]] auto i: std::views::iota(0, delta_size)) {
                int idx = dist(rng);
                delta_.emplace_back(el[idx]);
                el[idx].first = -1;
            }
        }

        double stream_ratio_;
        EdgeList delta_;
    };

}