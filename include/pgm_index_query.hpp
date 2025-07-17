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
    class pgm_querier{

    public:
        pgm::PGMIndex<EpsilonRecursive, Floating> index;
        K* dataset; // 用于存储数据集
        size_t data_size; // 数据集大小

        // 加载查询文件
        std::vector<K> load_queries(const std::string& filename) {
            std::ifstream in(filename, std::ios::binary);
            if (!in) {
                throw std::runtime_error("Could not open query file: " + filename);
            }

            K num_queries;
            in.read(reinterpret_cast<char*>(&num_queries), sizeof(K));

            std::vector<K> queries(num_queries);
            in.read(reinterpret_cast<char*>(queries.data()), num_queries * sizeof(K));
            in.close();

            return queries;
        }

        size_t find_exact_pos(const K& key, size_t lo, size_t hi) const {
            // 限制查找范围在 [0, data_size)
            lo = std::max<size_t>(0, lo);
            hi = std::min<size_t>(data_size, hi);

            // 二分查找
            while (lo <= hi) {
                size_t mid = (lo + hi) / 2;
                if (dataset[mid] == key) {
                    return mid;
                } else if (dataset[mid] < key) {
                    lo = mid + 1;
                } else {
                    hi = mid - 1;
                }
            }

            // 如果没找到，返回第一个大于等于 key 的位置（或可改为最接近的）
            return lo < data_size ? lo : data_size - 1;
        }

        // 执行多轮查询并统计耗时
        void run_query_rounds(const std::string& input_basename, int rounds = 10) {
            std::vector<K> tmp_dataset = load_data(input_basename);
            dataset = tmp_dataset.data();
            data_size = tmp_dataset.size();

            using namespace std::chrono;

            uint64_t total_queries = 0;
            uint64_t total_time_us = 0;

            for (int round = 0; round < rounds; ++round) {
                // 加载本轮查询数据
                std::vector<K> queries = load_queries(input_basename + "_queries_" + std::to_string(round) + ".query");

                // 开始计时
                auto start = high_resolution_clock::now();

                // 执行查询
                for (const auto& key : queries) {
                    auto pos = index.search(key);
                    auto extra_pos = find_exact_pos(key, pos.lo, pos.hi);
                    volatile size_t dummy = extra_pos; // 防止编译器优化
                    (void) dummy;
                }

                // 结束计时
                auto end = high_resolution_clock::now();
                auto elapsed = duration_cast<microseconds>(end - start);
                double avg_time_us = static_cast<double>(elapsed.count()) / queries.size();

                // 累计到总计中
                total_queries += queries.size();
                total_time_us += elapsed.count();

                // 输出本轮结果
                std::cerr << "Round " << round + 1 << ": "
                          << "Total time = " << elapsed.count() << " μs, "
                          << "Avg = " << avg_time_us << " μs/query"
                          << std::endl;
            }

            // 最后输出总体平均
            if (total_queries > 0) {
                double overall_avg_us = static_cast<double>(total_time_us) / rounds;
                std::cerr << "[Overall Average] " << "Across\t" << rounds << "\trounds\t"
                          << "Total Queries =\t " << total_queries << "\t"
                          << "Avg 100K Query Time =\t" << overall_avg_us << "\tμs/100K query\t"
                          << "Avg per Query =\t" << (static_cast<double>(total_time_us) / total_queries) << "\tμs/query\n"
                          << std::endl;
                output_message("Avg 100K Query Time:\t"+std::to_string(overall_avg_us));
            }
        }

        void load_model(const std::string& filename) {
            std::ifstream in(filename, std::ios::binary);
            if (!in) {
                throw std::runtime_error("Could not open file for reading: " + filename);
            }

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

    };
}
