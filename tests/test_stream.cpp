#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <filesystem>
#include <ranges>
#include "graph.h"
#include "builder.h"
#include "stream_graph.h"
#include "stream_builder.h"

TEST_CASE("stream graph is constructed successfully", "[construction]") {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";
    REQUIRE(fs::exists(graph_file_path));
    gtool::StreamBuilder<Node> builder{graph_file_path.string()};
    gtool::StreamGraph<Node> graph = builder.build_csr();
    REQUIRE(graph.get_vertex_number() == 77360);
    REQUIRE(graph.get_edge_number() == 905468 - static_cast<int>(905468 * 0.01));

    std::size_t edge_number = graph.get_edge_number();
    gtool::Streamer<Node, Node, decltype(builder.delta())> streamer(graph, builder.delta());
    streamer.stream_next_batch();
    REQUIRE(graph.get_edge_number() == edge_number + 1'000);
    REQUIRE(graph.get_edge_number_prev() == edge_number);
    streamer.update_graph();
    REQUIRE(graph.get_edge_number() == edge_number + 1'000);
    REQUIRE(graph.get_edge_number_prev() == edge_number + 1'000);
}

TEST_CASE("additions are handled properly", "[addition][time-consuming]") {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";
    REQUIRE(fs::exists(graph_file_path));
    gtool::Builder<Node> builder(graph_file_path.string());
    gtool::Graph<Node> graph = builder.build_csr();
    gtool::StreamBuilder<Node> stream_builder(graph_file_path.string());
    gtool::StreamGraph<Node> stream_graph = stream_builder.build_csr();
    gtool::Streamer<Node, Node, decltype(stream_builder.delta())> streamer(stream_graph, stream_builder.delta());

    unsigned edge_number = stream_graph.get_edge_number();
    int next_batch_size = streamer.next_batch_size();
    streamer.stream_next_batch();
    REQUIRE(stream_graph.get_edge_number() == edge_number + next_batch_size);

    while (streamer.has_more_batches()) {
        streamer.stream_next_batch();
    }
    streamer.update_graph();

    bool all_equal = true;
    for (int i = 0; i < graph.get_vertex_number(); ++i) {
        all_equal = all_equal && (graph.out_degree(i) == stream_graph.out_degree(i));
        all_equal = all_equal && (graph.in_degree(i) == stream_graph.in_degree(i));
    }
    REQUIRE(all_equal);

    for (int i = 0; i < graph.get_vertex_number(); ++i) {
        std::vector<int> a;
        std::vector<int> b;
        for (auto const &dst : graph.out_neighbors(i)) {
            a.emplace_back(dst);
        }
        for (auto const &dst : stream_graph.out_neighbors(i)) {
            b.emplace_back(dst);
        }
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        REQUIRE_THAT(a, Catch::Matchers::RangeEquals(b));
    }
}

TEST_CASE("deletions are handled properly", "[deletion]") {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";
    REQUIRE(fs::exists(graph_file_path));
    gtool::StreamBuilder<Node> stream_builder(graph_file_path.string());
    gtool::StreamGraph<Node> graph = stream_builder.build_csr();
    REQUIRE(graph.get_vertex_number() == 77360);
    REQUIRE(graph.get_edge_number() == 905468 - static_cast<int>(905468 * 0.01));
    int edge_number = graph.get_edge_number();

    gtool::Streamer<Node, Node, decltype(stream_builder.delta())> streamer(graph, stream_builder.delta());
    auto [next_add_first, next_add_last] = streamer.view_next_batch();
    gtool::EdgeList<Node, Node> deletions(next_add_first, next_add_last);

    REQUIRE(deletions.size() == streamer.next_batch_size());
    streamer.stream_next_batch(true);
    REQUIRE(graph.get_edge_number() == edge_number + 1'000);
    edge_number = graph.get_edge_number();

    gtool::DStreamer<Node, Node, std::add_lvalue_reference_t<decltype(deletions)>> d_streamer(graph, deletions);
    d_streamer.stream_next_batch();
    int delta{};
    for (auto u : std::views::iota(0, static_cast<int>(graph.get_vertex_number()))) {
        delta += graph.out_degree_delta(u);
    }
    int delta2{};
    for (auto u : std::views::iota(0, static_cast<int>(graph.get_vertex_number()))) {
        delta2 += graph.in_degree_delta(u);
    }
    REQUIRE(delta == -100);
    REQUIRE(delta2 == -100);
    REQUIRE(graph.get_edge_number() == edge_number - 100);

    while (d_streamer.has_more_batches()) {
        d_streamer.stream_next_batch();
    }
    d_streamer.update_graph();
    REQUIRE(graph.get_edge_number() == 905468 - static_cast<int>(905468 * 0.01));
}

TEST_CASE("Graph is restored properly", "[restore]") {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";
    REQUIRE(fs::exists(graph_file_path));
    gtool::StreamBuilder<Node> stream_builder(graph_file_path.string(), 0.04);
    gtool::StreamGraph<Node> graph = stream_builder.build_csr();
    gtool::Streamer<Node, Node, decltype(stream_builder.delta())> streamer(graph, stream_builder.delta());
    while (streamer.has_more_batches()) {
        streamer.stream_next_batch();
    }

    gtool::EdgeList<Node, Node> el; el.reserve(500);
    auto const &delta = stream_builder.delta();
    std::mt19937 rng{ std::random_device{}() };
    std::uniform_int_distribution<unsigned> dist(0, delta.size() - 1);
    for (int i : std::views::iota(0, 500)) {
        unsigned idx = dist(rng);
        el.emplace_back(delta[idx]);
    }
    gtool::DStreamer<Node, Node, std::add_lvalue_reference_t<decltype(el)>> d_streamer(graph, el);
    while (d_streamer.has_more_batches()) {
        d_streamer.stream_next_batch();
    }
    d_streamer.update_graph();
    graph.restore();
    streamer.reset();
    for (int i = 0; i < graph.get_vertex_number(); ++i) {
        REQUIRE(graph.out_degree(i) == graph.out_degree_base(i));
        REQUIRE(graph.in_degree(i) == graph.in_degree_base(i));
        {
            auto first1 = graph.out_neighbors(i).begin(), end1 = graph.out_neighbors(i).end();
            auto first2 = graph.out_neighbors_base(i).begin();
            for ( ; first1 != end1; ++first1, ++first2) {
                REQUIRE(*first1 == *first2);
            }
        }
        {
            auto first1 = graph.in_neighbors(i).begin(), end1 = graph.in_neighbors(i).end();
            auto first2 = graph.in_neighbors_base(i).begin();
            for ( ; first1 != end1; ++first1, ++first2) {
                REQUIRE(*first1 == *first2);
            }
        }
    }
}

TEST_CASE("FixRange works properly", "[tiny][fix-range]") {
    gtool::FixedRange<int> fr(10);
    auto first = std::begin(fr);
    for (int i{}; first != std::end(fr) && i < 10; ++i, ++first) {
        REQUIRE(*first == 10);
    }
    REQUIRE(first != std::end(fr));
    first = std::begin(fr);
    for (int i{}; i < 100; ++i)
        ++first;
    int cnt{};
    for (auto tmp = std::begin(fr); tmp != first; ++tmp) {
        ++cnt;
    }
    REQUIRE(cnt == 100);
}