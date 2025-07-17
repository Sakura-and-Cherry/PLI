#pragma once
#include <cstdint>



static constexpr uint32_t PageSize = 4096; // 4KB
static long double lambda = 0.3; // space and time cost平衡参数
static long double lambda_coefficient = 2 * lambda / (1 - lambda); // 用于计算公式(5)中的系数
static long double scale_factor = 1; // 用于公式3和公式4的数量级缩放
static std::string log_path = "./log/";
// using K = uint64_t; // Key type
#ifdef USE_UINT32_KEY
using K = uint32_t;
#else
// 默认使用 uint64_t
using K = uint64_t;
#endif