// #pragma once
// #include <vector>
// #include <config.hpp>
// #include <random>

// #include "tools.hpp"

// using gap_prefix_type = __uint128_t;

// struct Page {
//     size_t start_idx;
//     size_t end_idx;
//     Page(const int start_idx, const int end_idx): start_idx(start_idx), end_idx(end_idx) {}
// };

// struct Block {
//     size_t epsilon;
//     size_t start_idx, end_idx; // 记录的是0到n-1
//     Block(const std::vector<gap_prefix_type> &gap_prefix_sum, const std::vector<gap_prefix_type> &gap_prefix_squares, size_t begin_idx, size_t end_idx, long double &cost): start_idx(begin_idx), end_idx(end_idx) {
//         // 计算gap_avg和gap_var
//         size_t n = end_idx - begin_idx;
//         if (n < 2) {
//             return;
//         }

//         long double  sum = gap_prefix_sum[end_idx - 1] - gap_prefix_sum[begin_idx];
//         long double  var_sum = gap_prefix_squares[end_idx - 1] - gap_prefix_squares[begin_idx];
//         long double gap_avg = sum / (n - 1);
//         // 计算无偏方差
//         long double gap_var = (var_sum - (n - 1) * gap_avg * gap_avg) / (n - 2); // Bessel's correction 方差用无偏估计
//         // 计算epsilon
//         epsilon = static_cast<size_t>(std::ceil(std::sqrt(gap_var) / gap_avg * std::sqrt(lambda_coefficient * n))); // 采用向上取整得到epsilon
//         epsilon = epsilon > 1 ? epsilon : 1; // 确保epsilon至少为1
//         long double cov = n * gap_var / (gap_avg * gap_avg * epsilon * epsilon) / scale_factor;
//         cost += lambda * cov + (1 - lambda) * log(epsilon); // 更新成本
//     }
// };

// class DataPartition {
// public:
//     std::vector<Page> pages;
//     std::vector<Block> blocks;
//     std::string partition_type = "greedy"; // false表示DP分区，true表示Greedy分区
//     long double optimal_cost = 0; // DP分区的总成本
//     long double greedy_cost = 0; // Greedy分区的总成本
//     long double uniform_cost = 0; // 均匀分区的总成本
//     long double random_cost = 0; // 随机分区的总成本
//     // gap 的前缀和数组
//     std::vector<gap_prefix_type> gap_prefix_sum;     // gap 的前缀和
//     std::vector<gap_prefix_type> gap_prefix_squares; // gap 的平方的前缀和


//     DataPartition() = default;

//     // 按page_size划分为page
//     explicit DataPartition(const std::vector<K> &data) {
//         // 计算每一页的数量
//         size_t page_size = PageSize / sizeof(K); // 每页的元素数量
//         size_t n = data.size();
//         for (size_t i = 0; i < n; i += page_size) {
//             size_t end = std::min(i + page_size, n);
//             pages.emplace_back(i, end);
//         }
//         // 预处理 gap 的前缀和数组
//         gap_prefix_sum.resize(n, 0.0);
//         gap_prefix_squares.resize(n, 0.0);

//         if (n > 1) {
//             gap_prefix_sum[0] = 0.0; // gap[0] 是 data[1]-data[0]，对应 index 0 的 gap 是无效的
//             gap_prefix_squares[0] = 0.0;

//             for (size_t i = 1; i < n; ++i) {
//                 gap_prefix_type gap = data[i] - data[i - 1];
//                 if (gap <= 0)
//                     throw std::invalid_argument("gap cannot be negative");
//                 gap_prefix_sum[i] = gap_prefix_sum[i - 1] + gap;
//                 gap_prefix_squares[i] = gap_prefix_squares[i - 1] + gap * gap;
//             }
//         }

//         // 输出划分的page
//         std::cout << "Split into " << pages.size() << " pages." << std::endl;
//         estimate_scale();
//     }

//     // DP求最优分区
//     void optimal_partition() {
//         optimal_cost = 0;
//         partition_type = "optimal";
//         size_t n = pages.size();
//         if (pages.empty()) return;
//         std::vector<double> dp(n + 1, 1e10); // 初始化为一个很大的数
//         std::vector<int> prev(n + 1, -1);
//         dp[0] = 0;
//         for (size_t i = 1; i <= n; ++i) {
//             for (size_t j = 0; j < i; ++j) {
//                 long double cost = cost_function(pages[j].start_idx, pages[i - 1].end_idx);
//                 if (dp[j] + cost < dp[i]) {
//                     dp[i] = dp[j] + cost;
//                     prev[i] = j;
//                 }
//             }
//         }
//         // 回溯分区
//         blocks.clear();
//         int idx = n;
//         while (idx > 0) {
//             int start = prev[idx];
//             if (start == -1) break; // 如果没有前驱，结束回溯
//             // 将当前分区添加到block中
//             // 将page中的值提取出来
//             blocks.emplace_back(gap_prefix_sum, gap_prefix_squares, pages[start].start_idx, pages[idx - 1].end_idx, optimal_cost);
//             idx = start;
//         }
//         std::reverse(blocks.begin(), blocks.end());
//     }

//     // Greedy求次优分区
//     void greedy_partition() {
//         partition_type = "greedy";
//         greedy_cost = 0;
//         if (pages.empty()) return;
//         blocks.clear();
//         size_t start_idx = pages[0].start_idx, end_idx = pages[0].end_idx;
//         long double last_cost = cost_function(start_idx, end_idx);
//         for (size_t i = 1; i < pages.size(); ++i) {
//             long double single_cost = cost_function(pages[i].start_idx, pages[i].end_idx); // 当前页面的成本
//             long double split_cost = last_cost + single_cost;
//             long double merge_cost = cost_function(start_idx, pages[i].end_idx); // 合并成本
//             if (single_cost < 0 || split_cost < 0 || merge_cost < 0) {
//                 throw std::invalid_argument("invalid split cost, negative exits");
//             }
//             if (merge_cost < split_cost) {
//                 last_cost = merge_cost;
//                 end_idx = pages[i].end_idx; // 合并到当前分区
//             } else { // 分割了一个block
//                 blocks.emplace_back(gap_prefix_sum, gap_prefix_squares, start_idx, end_idx, greedy_cost);
//                 start_idx = pages[i].start_idx;
//                 end_idx = pages[i].end_idx; // 开始新的分区
//                 last_cost = single_cost;
//             }
//         }
//         // 最后一个分区
//         blocks.emplace_back(gap_prefix_sum, gap_prefix_squares, start_idx, end_idx, greedy_cost);
//     }

//     // 均匀划分，每次取1024个page当一个block
//     void uniform_partition() {
//         partition_type = "uniform";
//         uniform_cost = 0; // 可以复用 optimal_cost 存储总成本
//         blocks.clear();

//         const size_t block_page_size = 1024; // 每个 Block 包含的 Page 数量
//         size_t n_pages = pages.size();

//         for (size_t i = 0; i < n_pages; i += block_page_size) {
//             size_t block_end = std::min(i + block_page_size, n_pages);
//             size_t block_start_idx = pages[i].start_idx;
//             size_t block_end_idx = pages[block_end - 1].end_idx;

//             blocks.emplace_back(gap_prefix_sum, gap_prefix_squares, block_start_idx, block_end_idx, uniform_cost);
//         }

//         std::cout << "Uniform partition: " << blocks.size() << " blocks created." << std::endl;

//     }

//     // 每次随机取一定数量的page化为一个block
//     void random_partition() {
//         partition_type = "random";
//         random_cost = 0; // 可以复用 optimal_cost 存储总成本
//         blocks.clear();

//         size_t n_pages = pages.size();
//         size_t current_idx = 0;

//         std::random_device rd;
//         std::mt19937 gen(rd());
//         std::uniform_int_distribution<> block_size_dist(1, 1024); // 每个 Block 包含 [1, 1024] 个 Page

//         while (current_idx < n_pages) {
//             size_t remaining = n_pages - current_idx;
//             size_t block_size = block_size_dist(gen);
//             size_t actual_block_size = std::min(block_size, remaining);

//             size_t block_start_idx = pages[current_idx].start_idx;
//             size_t block_end_idx = pages[current_idx + actual_block_size - 1].end_idx;

//             blocks.emplace_back(gap_prefix_sum, gap_prefix_squares, block_start_idx, block_end_idx, random_cost);

//             current_idx += actual_block_size;
//         }

//         std::cout << "Random partition: " << blocks.size() << " blocks created." << std::endl;
//     }

//     // 定义成本函数，使用公式(5)，高效计算，使用预计算数据，并处理页面之间的 gap
//     long double cost_function(size_t start_idx, size_t end_idx) {
//         // 整个 Block 的 gap 起止索引
//         size_t gap_start = start_idx;
//         size_t gap_end = end_idx - 1;
//         size_t total_gap_n = gap_end - gap_start;

//         if (gap_start > gap_end || total_gap_n <= 0) {
//             return 0.0; // 没有 gap，直接返回
//         }

//         // 使用前缀和快速获取区间内的 gap_sum 和 gap_squares
//         long double total_gap_sum = gap_prefix_sum[gap_end] - gap_prefix_sum[gap_start];
//         long double total_gap_squares = gap_prefix_squares[gap_end] - gap_prefix_squares[gap_start];

//         long double mu_j = total_gap_sum / total_gap_n;
//         long double sigma_j_squared = (total_gap_squares - total_gap_n * mu_j * mu_j) / (total_gap_n - 1);

//         long double D_j = end_idx - start_idx;

//         long double epsilon = static_cast<size_t>(std::ceil(std::sqrt(sigma_j_squared) / mu_j * std::sqrt(lambda_coefficient * D_j))); // 采用向上取整得到epsilon
//         epsilon = epsilon > 1 ? epsilon : 1; // 确保epsilon至少为1
//         long double cov = (D_j * sigma_j_squared / (mu_j * mu_j * epsilon * epsilon)) / scale_factor; // 计算协方差

//         if (epsilon <= 0 || cov <= 0) {
//             std::cerr << "Error: start_idx: " << start_idx << " end_idx: " << end_idx << " gap_avg: " << mu_j << " gap_var: " << sigma_j_squared << " epsilon: " << epsilon << " cov: " << cov << "scale factor: " << scale_factor << std::endl;
//             return 0.0; // 防止除以零或负数
//         }

//         return (lambda * cov + (1 - lambda) * std::log(epsilon)); // 应用公式(5)
//     }

//     void estimate_scale() {
//         size_t n_pages = pages.size();
//         if (n_pages == 0) return;

//         const size_t num_samples = std::min<size_t>(10, n_pages); // 只取前10个page
//         std::vector<long double> epsilons;
//         std::vector<long double> covs;

//         for (size_t i = 0; i < num_samples; ++i) {
//             size_t start_idx = pages[i].start_idx;
//             size_t end_idx = pages[i].end_idx;

//             long double total_gap_start = start_idx;
//             long double total_gap_end = end_idx - 1;
//             size_t total_gap_n = total_gap_end - total_gap_start;

//             if (total_gap_n <= 0) continue;

//             long double total_gap_sum = gap_prefix_sum[total_gap_end] - gap_prefix_sum[total_gap_start];
//             long double total_gap_squares = gap_prefix_squares[total_gap_end] - gap_prefix_squares[total_gap_start];

//             long double mu_j = total_gap_sum / total_gap_n;
//             long double sigma_j_squared = (total_gap_squares - total_gap_n * mu_j * mu_j) / (total_gap_n - 1);

//             long double D_j = end_idx - start_idx;
//             long double epsilon = std::ceil(std::sqrt(sigma_j_squared) / mu_j * std::sqrt(lambda_coefficient * D_j));
//             epsilon = epsilon > 1 ? epsilon : 1;

//             long double cov = D_j * sigma_j_squared / (mu_j * mu_j * epsilon * epsilon);

//             epsilons.push_back(epsilon);
//             covs.push_back(cov);
//         }

//         if (epsilons.empty()) {
//             std::cerr << "No valid samples for comparing epsilon and cov." << std::endl;
//             return;
//         }

//         // 计算采样均值
//         long double avg_epsilon = std::accumulate(epsilons.begin(), epsilons.end(), 0.0L) / epsilons.size();
//         avg_epsilon = std::max(avg_epsilon, 1.0L); // 确保平均epsilon至少为1
//         long double avg_cov = std::accumulate(covs.begin(), covs.end(), 0.0L) / covs.size();

//         // int multiplier = 0;
//         int miu_epsilon = static_cast<int>(std::floor(std::log10(avg_epsilon)));
//         int miu_cov = static_cast<int>(std::floor(std::log10(avg_cov)));

//         scale_factor = std::pow(10, miu_cov - miu_epsilon);
//         // std::cout << "Scale Factor:\t" << scale_factor << std::endl;
//         std::cerr << "Avg Epsilon:\t" << avg_epsilon << ", Avg Cov:\t" << avg_cov << std::endl;
//         output_message("Scale Factor:\t" + std::to_string(scale_factor));
//     }

//     void summarize() const {
//         std::cout << "Total Pages:\t" << pages.size() << std::endl;

//         output_message("Total Blocks:\t" + std::to_string(blocks.size()));
//         if (partition_type == "optimal") {
//             output_message("Optimal Cost:\t" + std::to_string(optimal_cost));
//             for (const auto &block : blocks) {
//                 output_message("Block: [" + std::to_string(block.start_idx) + ", " + std::to_string(block.end_idx) + "), Epsilon: " + std::to_string(block.epsilon));
//             }
//         } else if (partition_type == "greedy") {
//             output_message("Greedy Cost:\t" + std::to_string(greedy_cost));
//             // for (const auto &block : blocks) {
//             //     output_message("Block: [" + std::to_string(block.start_idx) + ", " + std::to_string(block.end_idx) + "), Epsilon: " + std::to_string(block.epsilon));
//             // }
//         } else if (partition_type == "random") {
//             output_message("Random Cost:\t" + std::to_string(random_cost));
//         } else if (partition_type == "uniform") {
//             output_message("Uniform Cost:\t" + std::to_string(uniform_cost));
//         }
//     }
// };


#pragma once
#include <vector>
#include <config.hpp>
#include <random>
#include <queue>
#include <tuple>
#include <cmath>
#include "tools.hpp"

using gap_prefix_type = __uint128_t;

struct Page {
    size_t start_idx;
    size_t end_idx;
    Page(const int start_idx, const int end_idx): start_idx(start_idx), end_idx(end_idx) {}
};

struct Block {
    size_t epsilon;
    size_t start_idx, end_idx;
    long double cost;
    Block(const std::vector<gap_prefix_type> &gap_prefix_sum, const std::vector<gap_prefix_type> &gap_prefix_squares, size_t begin_idx, size_t end_idx, long double &total_cost): start_idx(begin_idx), end_idx(end_idx) {
        size_t n = end_idx - begin_idx;
        if (n < 2) {
            epsilon = 1;
            cost = 0;
            total_cost += cost;
            return;
        }
        long double sum = gap_prefix_sum[end_idx - 1] - gap_prefix_sum[begin_idx];
        long double var_sum = gap_prefix_squares[end_idx - 1] - gap_prefix_squares[begin_idx];
        long double gap_avg = sum / (n - 1);
        long double gap_var = (var_sum - (n - 1) * gap_avg * gap_avg) / (n - 2);
        epsilon = static_cast<size_t>(std::ceil(std::sqrt(gap_var / (gap_avg * gap_avg)) * std::sqrt(lambda_coefficient * n)));
        epsilon = std::max<size_t>(1, epsilon);
        long double cov = n * gap_var / (gap_avg * gap_avg * epsilon * epsilon);
        cost = lambda * cov + (1 - lambda) * std::log(epsilon);
        total_cost += cost;
    }
};


class DataPartition {
public:
    std::vector<Page> pages;
    std::vector<Block> blocks;
    std::string partition_type = "greedy";
    long double optimal_cost = 0;
    long double greedy_cost = 0;
    long double uniform_cost = 0;
    long double random_cost = 0;
    size_t partition_num = 10;
    std::vector<gap_prefix_type> gap_prefix_sum;
    std::vector<gap_prefix_type> gap_prefix_squares;

    DataPartition() = default;

    explicit DataPartition(const std::vector<K> &data, size_t m = 10) : partition_num(m) {
        size_t page_size = PageSize / sizeof(K);
        size_t n = data.size();
        for (size_t i = 0; i < n; i += page_size) {
            size_t end = std::min(i + page_size, n);
            pages.emplace_back(i, end);
        }

        gap_prefix_sum.resize(n, 0);
        gap_prefix_squares.resize(n, 0);

        for (size_t i = 1; i < n; ++i) {
            gap_prefix_type gap = data[i] - data[i - 1];
            if (gap <= 0)
                throw std::invalid_argument("gap cannot be negative or zero");
            gap_prefix_sum[i] = gap_prefix_sum[i - 1] + gap;
            gap_prefix_squares[i] = gap_prefix_squares[i - 1] + gap * gap;
        }
        std::cout << "Split into " << pages.size() << " pages." << std::endl;
        estimate_scale();
    }

    void optimal_partition() {
        partition_type = "optimal";
        size_t n = pages.size();
        std::vector<std::vector<long double>> dp(n + 1, std::vector<long double>(partition_num + 1, 1e18));
        std::vector<std::vector<int>> prev(n + 1, std::vector<int>(partition_num + 1, -1));
        dp[0][0] = 0;

        for (size_t j = 1; j <= n; ++j) {
            for (size_t k = 1; k <= partition_num; ++k) {
                for (size_t i = k - 1; i < j; ++i) {
                    long double tmp_cost = cost_function(pages[i].start_idx, pages[j - 1].end_idx);
                    if (dp[i][k - 1] + tmp_cost < dp[j][k]) {
                        dp[j][k] = dp[i][k - 1] + tmp_cost;
                        prev[j][k] = i;
                    }
                }
            }
        }

        blocks.clear();
        int j = n, k = partition_num;
        while (k > 0 && j > 0) {
            int i = prev[j][k];
            if (i == -1) break;
            blocks.emplace_back(gap_prefix_sum, gap_prefix_squares, pages[i].start_idx, pages[j - 1].end_idx, optimal_cost);
            j = i;
            k--;
        }
        std::reverse(blocks.begin(), blocks.end());
    }

    void greedy_partition() {
        partition_type = "greedy";
        greedy_cost = 0;
        blocks.clear();

        using Segment = std::tuple<long double, size_t, size_t>;  // cost, start_idx, end_idx
        std::priority_queue<Segment> pq;

        // Push initial segment (entire dataset)
        size_t n_pages = pages.size();
        if (n_pages == 0) return;
        pq.emplace(cost_function(pages.front().start_idx, pages.back().end_idx), 0, n_pages);

        // Split until we have enough total partitions (in blocks + in queue)
        while (blocks.size() + pq.size() < partition_num && !pq.empty()) {
            auto [cur_cost, s, t] = pq.top(); pq.pop();
            if (t - s <= 1) {
                // Too small to split further
                blocks.emplace_back(gap_prefix_sum, gap_prefix_squares, pages[s].start_idx, pages[t - 1].end_idx, greedy_cost);
                continue;
            }

            size_t best_split = s + 1;
            long double best_combined = std::numeric_limits<long double>::max();

            for (size_t i = s + 1; i < t; ++i) {
                long double left = cost_function(pages[s].start_idx, pages[i - 1].end_idx, 2);
                long double right = cost_function(pages[i].start_idx, pages[t - 1].end_idx, 2);
                if (left + right - 1e-5 < best_combined) {
                    best_combined = left + right;
                    best_split = i;
                }
            }

            // Decide whether to split or keep as is
            if (best_combined < cur_cost) {
                pq.emplace(cost_function(pages[s].start_idx, pages[best_split - 1].end_idx), s, best_split);
                pq.emplace(cost_function(pages[best_split].start_idx, pages[t - 1].end_idx), best_split, t);
            } else {
                blocks.emplace_back(gap_prefix_sum, gap_prefix_squares, pages[s].start_idx, pages[t - 1].end_idx, greedy_cost);
            }
        }

        // Add remaining segments in priority queue as final blocks
        while (!pq.empty()) {
            auto [_, s, t] = pq.top(); pq.pop();
            blocks.emplace_back(gap_prefix_sum, gap_prefix_squares, pages[s].start_idx, pages[t - 1].end_idx, greedy_cost);
        }

        // Ensure blocks are sorted by their start index
        std::sort(blocks.begin(), blocks.end(), [](const Block &a, const Block &b) {
            return a.start_idx < b.start_idx;
        });
    }


    void uniform_partition() {
        partition_type = "uniform";
        uniform_cost = 0;
        blocks.clear();

        size_t n_pages = pages.size();
        size_t pages_per_block = std::max<size_t>(1, (n_pages + partition_num - 1) / partition_num);

        for (size_t i = 0; i < n_pages; i += pages_per_block) {
            size_t block_end = std::min(i + pages_per_block, n_pages);
            size_t block_start_idx = pages[i].start_idx;
            size_t block_end_idx = pages[block_end - 1].end_idx;
            blocks.emplace_back(gap_prefix_sum, gap_prefix_squares, block_start_idx, block_end_idx, uniform_cost);
        }
    }

    void random_partition() {
        partition_type = "random";
        random_cost = 0;
        blocks.clear();

        size_t n_pages = pages.size();
        if (partition_num > n_pages) partition_num = n_pages;

        std::vector<size_t> indices;
        for (size_t i = 1; i < n_pages; ++i) indices.push_back(i);
        std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});

        indices.resize(partition_num - 1);
        indices.push_back(0);
        indices.push_back(n_pages);
        std::sort(indices.begin(), indices.end());

        for (size_t i = 1; i < indices.size(); ++i) {
            size_t s = indices[i - 1], t = indices[i];
            size_t block_start_idx = pages[s].start_idx;
            size_t block_end_idx = pages[t - 1].end_idx;
            blocks.emplace_back(gap_prefix_sum, gap_prefix_squares, block_start_idx, block_end_idx, random_cost);
        }
    }

    long double cost_function(size_t start_idx, size_t end_idx, int par_num = 1) {
        if (end_idx - start_idx < 2) return 0.0;
        size_t gap_start = start_idx, gap_end = end_idx - 1;
        long double first_gap = gap_prefix_sum[gap_start];
        long double last_gap = gap_prefix_sum[gap_end];
        long double sum = gap_prefix_sum[gap_end] - gap_prefix_sum[gap_start];
        long double sum_sq = gap_prefix_squares[gap_end] - gap_prefix_squares[gap_start];
        size_t n = gap_end - gap_start;
        if (n <= 1) return 0.0;
        long double mu = sum / n;
        long double var = (sum_sq - n * mu * mu) / (n - 1);
        if (mu <= 0 || var <= 0) return 0.0;
        long double eps = std::max(1.0L, std::ceil(std::sqrt(var / (mu * mu)) * std::sqrt(lambda_coefficient * (end_idx - start_idx))));
        long double cov = (end_idx - start_idx) * var / (mu * mu * eps * eps);
        return lambda * cov + (1 - lambda) * std::log(eps)/ par_num;
    }

    void estimate_scale() {

        if (pages.empty()) return;
        size_t count = std::min<size_t>(10, pages.size());
        long double total_eps = 0, total_cov = 0;
        for (size_t i = 0; i < count; ++i) {
            size_t start_idx = pages[i].start_idx;
            size_t end_idx = pages[i].end_idx;
            size_t gap_start = start_idx, gap_end = end_idx - 1;
            long double sum = gap_prefix_sum[gap_end] - gap_prefix_sum[gap_start];
            long double sum_sq = gap_prefix_squares[gap_end] - gap_prefix_squares[gap_start];
            size_t n = gap_end - gap_start;
            if (n < 2) continue;
            long double mu = sum / n;
            long double var = (sum_sq - n * mu * mu) / (n - 1);
            long double eps = std::max(1.0L, std::ceil(std::sqrt(var / (mu * mu)) * std::sqrt(lambda_coefficient * (end_idx - start_idx))));
            long double cov = (end_idx - start_idx) * var / (mu * mu * eps * eps);
            total_eps += eps;
            total_cov += cov;
        }
        long double avg_eps = std::max(1.0L, total_eps / count);
        long double avg_cov = total_cov / count;
        int log_eps = static_cast<int>(std::floor(std::log10(avg_eps)));
        int log_cov = static_cast<int>(std::floor(std::log10(avg_cov)));
        scale_factor = std::pow(10, log_cov - log_eps);
        output_message("Scale Factor:\t" + std::to_string(scale_factor));
    }

    void summarize() const {
        output_message("Partition type: " + partition_type);
        output_message("Total Blocks: " + std::to_string(blocks.size()));
        long double total_cost = 0;
        for (const auto &block : blocks) {
            total_cost += block.cost;
            output_message("Block: [" + std::to_string(block.start_idx) + ", " + std::to_string(block.end_idx) + "), Epsilon: " + std::to_string(block.epsilon) + ", Cost: " + std::to_string(block.cost));
        }
        output_message("Total Cost: " + std::to_string(total_cost));
    }
};