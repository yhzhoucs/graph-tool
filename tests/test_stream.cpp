#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <filesystem>
#include "graph.h"
#include "builder.h"
#include "stream_graph.h"
#include "stream_builder.h"

TEST_CASE( "stream graph constructed successfully", "[stream]" ) {
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

TEST_CASE("stream graph is restored", "[stream]") {
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
    while (streamer.has_more_batches()) {
        streamer.stream_next_batch();
    }
    streamer.restore_graph();
}

TEST_CASE("stream graph reconstructed correctly", "[stream][time-consuming]") {
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