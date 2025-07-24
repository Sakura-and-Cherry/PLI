#include <iostream>
#include <vector>
#include <fstream>
#include <string>

/**
 * @brief Translate the data from the file to the vector
 *
 * @param filename
 * @return std::vector<K>
 */
template<typename K>
std::vector<K> load_data(std::string FILE_NAME)
{

    /* Open file. */
    std::ifstream in(FILE_NAME, std::ios::binary);
    if (!in.is_open())
        exit(EXIT_FAILURE);

    /* Read number of keys. */
    K n_keys;
    in.read(reinterpret_cast<char *>(&n_keys), sizeof(K));

    /* Initialize vector. */
    std::vector<K> data;
    data.resize(n_keys);

    /* Read keys. */
    in.read(reinterpret_cast<char *>(data.data()), n_keys * sizeof(K));
    in.close();

    return data;
}


int main() {
    std::vector<u_int64_t> data = load_data<u_int64_t>("./data/SOSD/books_1M_uint64.bin");
    printf("Loaded %zu keys from the file.\n", data.size());
    return 0;
}
