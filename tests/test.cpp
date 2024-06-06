#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iostream>
#include "graph.h"
#include "builder.h"

TEST_CASE( "Graph is generated properly", "[graph]" ) {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";

    REQUIRE(fs::exists(graph_file_path));

    gtool::Builder<Node> builder{graph_file_path.string(), false, true};
    gtool::Graph<Node> graph = builder.build_csr();
    REQUIRE(graph.get_vertex_number() == 77360);
    REQUIRE(graph.get_edge_number() == 905468);
}