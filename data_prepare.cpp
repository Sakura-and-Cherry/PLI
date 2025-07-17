#include <iostream>
#include <string>
#include <tools.hpp>


int main(int argc, const char** argv)
{
    std::ios::sync_with_stdio(0);
    int mandatory = 3;
    if (argc < mandatory) {
        std::cerr << "Usage: " << argv[0] << ":\n" << "\t <collection_basename> <whether_query_generate>" << std::endl;
        return 1;
    }

    const std::string input_basename = argv[1];
    const std::string whether_query_generate = argv[2];

    // cerr << "Size of K " << sizeof(K) << endl;
    if (whether_query_generate == "f") {
        clean_data(input_basename);
    } else {
        generate_queries(input_basename);
    }

    return 0;
}
