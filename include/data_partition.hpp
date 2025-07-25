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
        long double cov = n * gap_var / (gap_avg * gap_avg * epsilon * epsilon) / scale_factor ;
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
        long double cov = (end_idx - start_idx) * var / (mu * mu * eps * eps) / scale_factor;
        return lambda * cov + (1 - lambda) * std::log(eps) / par_num;
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