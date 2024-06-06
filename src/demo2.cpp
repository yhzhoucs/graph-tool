#include "graph.h"
#include "builder.h"
#include "stream_builder.h"
#include <filesystem>
#include <iostream>

using Node = int;

int main(int argc, char *argv[]) {
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