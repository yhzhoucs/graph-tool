#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include "gtool/graph.h"
#include "gtool/builder.h"


TEST_CASE( "snap format is generated properly", "[construction]" ) {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";
    REQUIRE(fs::exists(graph_file_path));
    gtool::Builder<Node> builder{graph_file_path.string()};
    gtool::Graph<Node> graph = builder.build_csr();
    REQUIRE(graph.get_vertex_number() == 77360);
    REQUIRE(graph.get_edge_number() == 905468);
}

TEST_CASE( "snap format is generated properly for undirected", "[construction]" ) {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";
    REQUIRE(fs::exists(graph_file_path));
    gtool::Builder<Node> builder{graph_file_path.string(), true, true};
    gtool::Graph<Node> graph = builder.build_csr();
    REQUIRE(graph.get_vertex_number() == 77360);
    REQUIRE(graph.get_edge_number() == 1733629);
}

TEST_CASE( "matrix market format is generated properly", "[construction]" ) {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.mtx";
    REQUIRE(fs::exists(graph_file_path));
    gtool::Builder<Node> builder{graph_file_path.string()};
    gtool::Graph<Node> graph = builder.build_csr();
    REQUIRE(graph.get_vertex_number() == 77360);
    REQUIRE(graph.get_edge_number() == 905468);
}

TEST_CASE( "matrix market format is generated properly for undiredted", "[construction]" ) {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.mtx";
    REQUIRE(fs::exists(graph_file_path));
    gtool::Builder<Node> builder{graph_file_path.string(), true, true};
    gtool::Graph<Node> graph = builder.build_csr();
    REQUIRE(graph.get_vertex_number() == 77360);
    REQUIRE(graph.get_edge_number() == 1733629);
}

TEST_CASE( "graph is sorted properly", "[functional]" ) {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";
    REQUIRE(fs::exists(graph_file_path));
    gtool::Builder<Node> builder{graph_file_path.string()};
    gtool::Graph<Node> graph = builder.build_csr();
    REQUIRE(graph.get_vertex_number() == 77360);
    REQUIRE(graph.get_edge_number() == 905468);
    graph.sort_neighborhood(std::greater<>{});
    bool is_sorted = true;
    for (int i = 0; i < graph.get_vertex_number(); ++i) {
        auto it = graph.out_neighbors(i).begin();
        auto end = graph.out_neighbors(i).end();
        if (it != end) {
            --end;
            for ( ; it != end; ++it) {
                is_sorted = is_sorted && (graph.out_degree(*it) >= graph.out_degree(*(it+1)));
            }
        }
    }
    REQUIRE(is_sorted);
}

TEST_CASE( "graph is reordered properly", "[functional]" ) {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";
    REQUIRE(fs::exists(graph_file_path));
    gtool::Builder<Node> builder{graph_file_path.string()};
    gtool::Graph<Node> graph;
    std::vector<Node> new_ids;
    std::vector<Node> ids_remap;
    {
        gtool::Graph<Node> raw = builder.build_csr();
        std::tie(graph, new_ids, ids_remap) = gtool::reorder_by_degree(raw);
    }

    REQUIRE(graph.get_vertex_number() == 77360);
    REQUIRE(graph.get_edge_number() == 905468);

    bool is_global_sorted = true;
    for (int i = 0; i < graph.get_vertex_number() - 1; ++i) {
        is_global_sorted = is_global_sorted && (graph.out_degree(i) >= graph.out_degree(i+1));
    }
    REQUIRE(is_global_sorted);

    bool is_sorted = true;
    for (int i = 0; i < graph.get_vertex_number(); ++i) {
        auto it = graph.out_neighbors(i).begin();
        auto end = graph.out_neighbors(i).end();
        if (it != end) {
            --end;
            for ( ; it != end; ++it) {
                is_sorted = is_sorted && (graph.out_degree(*it) >= graph.out_degree(*(it+1)));
            }
        }
    }
    REQUIRE(is_sorted);
}

TEST_CASE( "graph is squeezed properly", "[functional]" ) {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";
    REQUIRE(fs::exists(graph_file_path));
    gtool::Builder<Node> builder{graph_file_path.string()};
    gtool::Graph<Node> graph;
    std::vector<Node> new_ids;
    std::vector<Node> ids_remap;
    {
        gtool::Graph<Node> raw = builder.build_csr();
        std::tie(graph, new_ids, ids_remap) = gtool::squeeze_graph(raw);
    }

    REQUIRE(graph.get_vertex_number() == 77360);
    REQUIRE(graph.get_edge_number() == 905468);

    bool is_global_sorted = true;
    for (int i = 0; i < graph.get_vertex_number() - 1; ++i) {
        is_global_sorted = is_global_sorted && (graph.out_degree(i) > 0 || graph.in_degree(i) > 0);
    }
    REQUIRE(is_global_sorted);
}

TEST_CASE( "graph is simplified properly", "[functional][time-consuming]" ) {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";
    REQUIRE(fs::exists(graph_file_path));
    gtool::Builder<Node> builder{graph_file_path.string()};
    gtool::Graph<Node> graph;
    std::vector<Node> new_ids;
    std::vector<Node> ids_remap;
    {
        gtool::Graph<Node> raw = builder.build_csr();
        graph = gtool::simplify_graph(raw);
    }

    REQUIRE(graph.get_vertex_number() == 77360);
    REQUIRE(graph.get_edge_number() == 905468);

    std::vector<Node> helper;
    helper.reserve(50000);
    for (int i = 0; i < graph.get_vertex_number(); ++i) {
        for (auto const &dst : graph.out_neighbors(i)) {
            helper.emplace_back(gtool::get_dst_id(dst));
        }
        std::sort(helper.begin(), helper.end());
        auto unique = std::unique(helper.begin(), helper.end());
        REQUIRE(unique == helper.end());
        helper.clear();
    }
}