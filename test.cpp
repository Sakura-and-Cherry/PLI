#include <iostream>
#include <vector>

int main() {
    int start_idx = 0, end_idx = 10; // 假设这是划分的起始和结束索引
    std::vector<int> data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    std::vector<int> block_data(data.begin() + start_idx, data.begin() + end_idx);
    for (const auto& val : block_data) {
        std::cout << val << " ";
    }
}
