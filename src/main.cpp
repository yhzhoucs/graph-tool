#include "graph.h"
#include "builder.h"
#include "stream_builder.h"
#include <filesystem>
#include <iostream>

using Node = int;

int demo_1(int argc, char *argv[]) {
    namespace fs = std::filesystem;

    fs::path graph_file_path;
    if (argc < 2) {
        graph_file_path = fs::path(DATASET_PATH) / "Slashdot0811.txt";
    } else {
        graph_file_path = argv[1];
    }

    if (!fs::exists(graph_file_path)) {
        std::cerr << "Error: graph file not found!" << std::endl;
        return -1;
    }

    gtool::Builder<Node> builder{graph_file_path.string(), false, true}; // add reverse edges
    gtool::Graph<Node> graph;
    {
        gtool::Graph<Node> raw = builder.build_csr();
        graph = gtool::simplify_graph(raw); // remove redundant edges
    }

    std::clog << "Graph: " << graph_file_path.string() << std::endl;
    std::clog << "Vertex: " << graph.get_vertex_number() << std::endl;
    std::clog << "Edge: " << graph.get_edge_number() << std::endl;

    auto const *out_offset = graph.get_out_offset();
    auto const *out_neigh = graph.get_out_neigh();

    for (std::size_t i{}; i < 10/* graph.get_vertex_number() */; ++i) {
        std::cout << out_offset[i] << " ";
    }
    std::cout << std::endl;
    for (std::size_t i{}; i < 10/* graph.get_edge_number() */; ++i) {
        std::cout << gtool::get_dst_id(out_neigh[i]) << " ";
    }
    std::cout << std::endl;
    return 0;
}

int demo_2(int argc, char *argv[]) {
    namespace fs = std::filesystem;

    fs::path graph_file_path;
    if (argc < 2) {
        graph_file_path = fs::path(DATASET_PATH) / "Slashdot0811.txt";
    } else {
        graph_file_path = argv[1];
    }

    if (!fs::exists(graph_file_path)) {
        std::cerr << "Error: graph file not found!" << std::endl;
        return -1;
    }

    gtool::StreamBuilder<Node> stream_builder(graph_file_path.string());
    gtool::Graph<Node> graph = stream_builder.build_csr();
    auto const &delta = stream_builder.delta();

    std::cout << "Streaming: " << std::endl;
    for (auto const &edge: delta | std::views::take(10)) {
        std::cout << edge.first << "->" << edge.second << std::endl;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    demo_2(argc, argv);
    return 0;
}
