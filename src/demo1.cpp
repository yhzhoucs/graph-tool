#include "gtool/graph.h"
#include "gtool/builder.h"
#include <filesystem>
#include <iostream>

using Node = int;

int main(int argc, char *argv[]) {
    namespace fs = std::filesystem;

    fs::path graph_file_path;
    if (argc < 2) {
        graph_file_path = fs::path(DATASET_PATH) / "soc-Slashdot0811.txt";
    } else {
        graph_file_path = argv[1];
    }

    if (!fs::exists(graph_file_path)) {
        std::cerr << "Error: graph file not found!" << std::endl;
        return -1;
    }

    gtool::Graph<Node> graph;
    {
        auto raw = gtool::build_graph_from_file<Node>(graph_file_path.string(), true);
        graph = gtool::simplify_graph(raw); // remove redundant edges
    }

    std::clog << "Graph: " << graph_file_path.string() << std::endl;
    std::clog << "Vertex: " << graph.get_vertex_number() << std::endl;
    std::clog << "Edge: " << graph.get_edge_number() << std::endl;

    return 0;
}