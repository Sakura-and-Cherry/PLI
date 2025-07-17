#include <iostream>
#include <chrono>
#include <string>
#include <fstream>
#include "pgm_index_query.hpp"

std::string decode_type = "";
bool read_only = false;

void query_collection_pgm(const std::string collection_basename, const std::string index_basename) {
    typedef pgm_sequence::pgm_querier<> PGM_INDEX_QUERIER;
    PGM_INDEX_QUERIER index;
    index.load_model(index_basename);
    index.run_query_rounds(collection_basename);
}

int main(int argc, const char** argv)
{
    std::ios::sync_with_stdio(0);
    int mandatory = 5;
    if (argc < mandatory) {
        std::cerr << "Usage: " << argv[0] << ":\n" << "\t <collection_basename> <index_basename> <lambda> <partition_type> <log_path>" << std::endl;
        return 1;
    }

    const std::string collection_basename = argv[1];
    const std::string index_basename = argv[2];
    lambda = stod(argv[3]);
    const std::string partition_type = argv[4];
    const std::string log_path = argv[5];
    logStream.open(log_path, std::ios::app);
    output_message("————Partition Type:\t" + partition_type + "\tLambda:\t" + std::to_string(lambda) + "\t——————");


    query_collection_pgm(collection_basename, index_basename);

    return 0;
}
