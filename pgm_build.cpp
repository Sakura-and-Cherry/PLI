#include <iostream>
#include <string>
#include <fstream>
#include <pgm_index_build.hpp>

bool read_only = false;

void create_collection_pgm(const std::string collection_basename, const std::string index_basename, const std::string partition_type) {
    typedef pgm_sequence::pgm_builder<> PGM_INDEX_BUILDER;
    PGM_INDEX_BUILDER index;
    if (!read_only){
        index.build_model(collection_basename, partition_type);
        index.save_model(index_basename);
        index.stats_model();
    } else {
        index.load_model(index_basename);
        index.query_random(collection_basename, 10); // for test
        index.stats_model();
    }
}

int main(int argc, const char** argv) {
    std::ios::sync_with_stdio(0);
    int mandatory = 7;
    if (argc < mandatory) {
        std::cerr << "Usage: " << argv[0] << ":\n" << "\t <collection_basename> <index_basename> <read_only> <lambda> <partition_type> <log_path>" << std::endl;
        return 1;
    }

    const std::string collection_basename = argv[1];
    const std::string index_basename = argv[2];
    const std::string read_only_str = argv[3];
    read_only = (read_only_str == "t");
    lambda = stod(argv[4]);
    lambda_coefficient = 2 * lambda / (1 - lambda); // 更新系数
    // std::cerr << "lambda: " << lambda << std::endl;
    const std::string partition_type = argv[5];
    const std::string log_path = argv[6];

    logStream.open(log_path, std::ios::out | std::ios::app);
    output_message("————Partition Type:\t" + partition_type + "\tLambda:\t" + std::to_string(lambda) + "\t——————");

    create_collection_pgm(collection_basename, index_basename, partition_type);
    return 0;
}
