#include "gtool/stream_builder.h"
#include "gtool/stream_graph.h"
#include <filesystem>
#include <iostream>

using Node = int;

int main(int argc, char *argv[]) {
    using Node = int;
    namespace fs = std::filesystem;
    fs::path graph_file_path(DATASET_PATH);
    graph_file_path /= "soc-Slashdot0811.txt";
    gtool::StreamBuilder<Node> builder{graph_file_path.string()};
    gtool::StreamGraph<Node> graph = builder.build_csr();
    gtool::Streamer<Node, Node, decltype(builder.delta())> streamer(graph, builder.delta());
    std::cout << "Edge number: " << graph.get_edge_number() << std::endl;

    std::cout << "Max batch size: " << streamer.max_batch_size() << std::endl;
    std::cout << "Suggested batch size: " << streamer.suggest_max_batch_size(100) << std::endl;
    while (streamer.has_more_batches()) {
        std::cout << "Now stream a batch of edges into graph...";
        streamer.stream_next_batch();
        std::cout << "Edge number: " << graph.get_edge_number() << std::endl;
        streamer.update_graph();
    }
    // build with `-DBUILD_WITH_RESTORE=ON`
    graph.restore();
    std::cout << "Graph restored. Edge number: " << graph.get_edge_number() << std::endl;

    return 0;
}