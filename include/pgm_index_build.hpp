#pragma once

#include <iostream>
#include <cassert>
#include <map>
#include <vector>
#include <pgm_index.hpp>
#include <variant>
#include <config.hpp>
#include <../external/mm_file/include/mm_file/mm_file.hpp>

#include "data_partition.hpp"
#include "tools.hpp"

namespace pgm_sequence {
    #define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))
    template <size_t EpsilonRecursive=4, typename Floating=long double> // K is uint32_t or uint64_t
    class pgm_builder{

    public:
        pgm::PGMIndex<EpsilonRecursive, Floating> index;

        void build_model(std::string input_basename, std::string partition_type) {
            std::cerr << "Read File [Build]: " << input_basename << std::endl;

            // 划分数据为页
            std::vector<K> dataset = load_data(input_basename);
            // std::vector<K> dataset = generate_sample_data();
            // sort(dataset.begin(), dataset.end());

            // 开始计时

            DataPartition* partitioner = new DataPartition(dataset);
            if (partition_type == "greedy")
                partitioner -> greedy_partition();
            else if (partition_type == "optimal")
                partitioner -> optimal_partition();
            else if (partition_type == "random")
                partitioner -> random_partition();
            else if (partition_type == "uniform")
                partitioner -> uniform_partition();
            else {
                throw std::invalid_argument("Invalid partition type");
            }
            partitioner -> summarize();
            auto start = std::chrono::high_resolution_clock::now();
            index = pgm::PGMIndex<EpsilonRecursive, Floating>(dataset, partitioner->blocks);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            output_message("Build Time:\t" + std::to_string(duration.count()));
        }

        void stats_model() {
            // 输出pgm的高度，线段数和总大小
            // std::cout << "PGM Index Stats:\t" << "Height:\t" << index.height() << "\tSegments Count:\t" << index.segments_count() << "\tSize in Bytes:\t" << index.size_in_bytes() << "\tFirst Key:\t" << index.first_key << "\tTotal Elements:\t" << index.n << std::endl;
            output_message("Segments Number:\t" + std::to_string(index.segments_count()));
            output_message("Segments Size:\t" + std::to_string(index.size_in_bytes()));
            output_message("Height:\t" + std::to_string(index.height()));
        }

        void save_model(const std::string& filename) {
            std::ofstream out(filename, std::ios::binary);
            if (!out) {
                throw std::runtime_error("Could not open file for writing: " + filename);
            }

            // 写入基本字段
            out.write(reinterpret_cast<const char*>(&index.n), sizeof(index.n));
            out.write(reinterpret_cast<const char*>(&index.first_key), sizeof(index.first_key));

            // 写入 segments
            size_t segments_size = index.segments.size();
            out.write(reinterpret_cast<const char*>(&segments_size), sizeof(segments_size));
            out.write(reinterpret_cast<const char*>(index.segments.data()), segments_size * sizeof(typename pgm::PGMIndex<EpsilonRecursive, Floating>::Segment));

            // 写入 levels_offsets
            size_t offsets_size = index.levels_offsets.size();
            out.write(reinterpret_cast<const char*>(&offsets_size), sizeof(offsets_size));
            out.write(reinterpret_cast<const char*>(index.levels_offsets.data()), offsets_size * sizeof(size_t));

            // 写入 epsilons_vec
            size_t epsilons_size = index.epsilons_vec.size();
            out.write(reinterpret_cast<const char*>(&epsilons_size), sizeof(epsilons_size));
            out.write(reinterpret_cast<const char*>(index.epsilons_vec.data()), epsilons_size * sizeof(std::pair<size_t, size_t>));

            out.close();
        }

        void load_model(const std::string& filename) {
            std::ifstream in(filename, std::ios::binary);
            if (!in) {
                throw std::runtime_error("Could not open file for reading: " + filename);
            }
            // index = new pgm::PGMIndex<EpsilonRecursive, Floating>();

            // 读取基本字段
            in.read(reinterpret_cast<char*>(&index.n), sizeof(index.n));
            in.read(reinterpret_cast<char*>(&index.first_key), sizeof(index.first_key));

            // 读取 segments
            size_t segments_size = 0;
            in.read(reinterpret_cast<char*>(&segments_size), sizeof(segments_size));
            index.segments.resize(segments_size);
            in.read(reinterpret_cast<char*>(index.segments.data()), segments_size * sizeof(typename pgm::PGMIndex<EpsilonRecursive, Floating>::Segment));

            // 读取 levels_offsets
            size_t offsets_size = 0;
            in.read(reinterpret_cast<char*>(&offsets_size), sizeof(offsets_size));
            index.levels_offsets.resize(offsets_size);
            in.read(reinterpret_cast<char*>(index.levels_offsets.data()), offsets_size * sizeof(size_t));

            // 读取 epsilons_vec
            size_t epsilons_size = 0;
            in.read(reinterpret_cast<char*>(&epsilons_size), sizeof(epsilons_size));
            index.epsilons_vec.resize(epsilons_size);
            in.read(reinterpret_cast<char*>(index.epsilons_vec.data()), epsilons_size * sizeof(std::pair<size_t, size_t>));

            in.close();
        }

        // only used for debug
        void query_random(const std::string &input_basename, uint64_t query_num) {
            std::vector<K> dataset = load_data(input_basename);
            for (uint64_t i = 0; i < query_num; ++i) {
                size_t true_pos = rand() % dataset.size();
                uint64_t query = dataset[true_pos];
                auto& pgm_index = index;
                // auto it = pgm_index.segment_for_key(query);
                // std::cout << "Query: " << query << ", Segment Key: " << it->key << std::endl;
                auto approx_pos = pgm_index.search(query);
                std::cout << "Query: " << query
                          << ", Approx Position: " << approx_pos.pos
                          << ", Range: [" << approx_pos.lo
                          << ", " << approx_pos.hi << "]"
                          << ", " << "Epsilon: " << approx_pos.epsilon
                          << ", True Position: " << true_pos
                          << std::endl;
            }
        }

    };
}
