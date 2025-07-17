#pragma once
#include <vector>
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <string>
#include <config.hpp>
using namespace std;
ofstream logStream;

// 加载数据返回data vector
inline vector<K> load_data(const string &filename) {
    /* Open file. */
    ifstream in(filename, ios::binary);
    if (!in.is_open())
        exit(EXIT_FAILURE);

    /* Read number of keys. */
    K n_keys;
    in.read(reinterpret_cast<char*>(&n_keys), sizeof(K));
    /* Initialize vector. */
    vector<K> data(n_keys);

    /* Read keys. */
    in.read(reinterpret_cast<char*>(data.data()), n_keys * sizeof(K));  
    in.close();

    if (data.back() == std::numeric_limits<K>::max()) {
        std::cerr << data.back() << " is the last element, which is reserved for infinity." << std::endl;
        while (data.back() == std::numeric_limits<K>::max()) {
            data.pop_back();
            n_keys--;
        }
    }


    std::cout << "Loaded " << n_keys << " elements from " << filename << std::endl;

    return data;

}

void clean_data(const string &filename) {
    /* Open file. */
    ifstream in(filename, ios::binary);
    if (!in.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    /* Read number of keys. */
    K n_keys;
    in.read(reinterpret_cast<char*>(&n_keys), sizeof(K));

    /* Initialize vector. */
    vector<K> data(n_keys);

    /* Read keys. */
    in.read(reinterpret_cast<char*>(data.data()), n_keys * sizeof(K));
    in.close();

    /* 原始数据大小 */
    size_t original_size = data.size();

    /* 去重：排序 + unique */
    sort(data.begin(), data.end());
    auto it = unique(data.begin(), data.end());
    data.erase(it, data.end());

    if (data.back() == std::numeric_limits<K>::max()) {
        std::cerr << data.back() << " is the last element, which is reserved for infinity." << std::endl;
        while (data.back() == std::numeric_limits<K>::max()) {
            data.pop_back();
            n_keys--;
        }
    }

    /* 输出去重结果 */
    size_t removed_count = original_size - data.size();
    cout << "Loaded " << original_size << " elements, removed " << removed_count << " duplicates." << endl;

    /* 保存去重后的数据到新文件 */
    string output_filename = filename + "_unique";
    ofstream out(output_filename, ios::binary);
    if (!out.is_open()) {
        cerr << "Failed to create output file: " << output_filename << endl;
        exit(EXIT_FAILURE);
    }

    // 写入去重后的元素数量
    K unique_n_keys = static_cast<K>(data.size());
    out.write(reinterpret_cast<const char*>(&unique_n_keys), sizeof(K));

    // 写入去重后的数据
    out.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(K));
    out.close();

    cout << "Unique data saved to: " << output_filename << endl;
}

inline vector<K> generate_sample_data() {
    vector<K> data;
    default_random_engine gen(random_device{}());

    // 第一段：gap 均值 ~2，方差 ~1
    normal_distribution<> gap_dist1(2.0, 1.0);  // mean=2, stddev=1
    int64_t current = 0;

    for (int i = 0; i < 1024; ++i) {
        int gap = static_cast<int>(round(gap_dist1(gen)));
        current += max(1, gap);
        data.push_back(current);
    }

    // 第二段：gap 均值 ~100，方差 ~5000
    normal_distribution<> gap_dist2(1000.0, sqrt(5000.0));  // mean=100, stddev=~70.7
    for (int i = 0; i < 1024; ++i) {
        int gap = static_cast<int>(round(gap_dist2(gen)));
        current += max(1, gap);
        data.push_back(current);
    }

    // 第三段：gap 均值 ~5000，方差 ~1,000,000
    normal_distribution<> gap_dist3(50000.0, sqrt(1000000.0));  // mean=5000, stddev=1000
    for (int i = 0; i < 1024; ++i) {
        int gap = static_cast<int>(round(gap_dist3(gen)));
        current += max(1, gap);
        data.push_back(current);
    }

    cout << "Generated " << data.size() << " elements with EXTREME gap variance differences." << endl;
    return data;
}

        // 生成查询文件，只需要调用一次
inline void generate_queries(std::string input_basename, int total_rounds = 10, int queries_per_round = 100000) {
    std::vector<K> dataset = load_data(input_basename);
    if (dataset.empty()) {
        std::cerr << "No data loaded for generating queries." << std::endl;
        return;
    }

    // 使用系统时间作为随机种子
    unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
    std::default_random_engine rng(seed);

    std::cout << "Generating " << total_rounds << " rounds of queries, " << queries_per_round << " per round." << std::endl;

    for (int r = 0; r < total_rounds; ++r) {
        std::vector<K> queries;
        queries.reserve(queries_per_round);

        // 生成 1000 个查询键
        for (int q = 0; q < queries_per_round; ++q) {
            std::uniform_int_distribution<size_t> dist(0, dataset.size() - 1);
            size_t true_pos = dist(rng);
            K query = dataset[true_pos];
            queries.push_back(query);
        }

        // 构造文件名
        std::string filename = input_basename + "_queries_" + std::to_string(r) + ".query";

        // 写入文件
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open()) {
            std::cerr << "Failed to create file: " << filename << std::endl;
            continue;
        }

        // 写入查询数量
        K num_queries = static_cast<K>(queries.size());
        out.write(reinterpret_cast<const char*>(&num_queries), sizeof(K));

        // 写入查询内容
        out.write(reinterpret_cast<const char*>(queries.data()), queries.size() * sizeof(K));
        out.close();

        std::cout << "Saved round " << r << " to " << filename << " (" << queries.size() << " queries)." << std::endl;
    }

    std::cout << "All query files generated." << std::endl;
}

inline string to_string_int128(__int128_t n) {
    if (n == 0) return "0";  
    string res;
    bool negative = n < 0;  
    if (negative) n = -n;  

    while (n > 0) {  
        res += (n % 10) + '0'; // 把数字转换为字符  
        n /= 10;  
    }  

    if (negative) res += '-';  
    reverse(res.begin(), res.end());
    return res;  
}

inline string formatted_time(int hours_offset = 8){
    // set current time as log file name
    chrono::seconds offset_seconds(hours_offset * 3600);
    auto now = chrono::system_clock::now();
    now = now + offset_seconds;
    auto in_time_t = chrono::system_clock::to_time_t(now);
    // decode time
    tm buf;
    localtime_r(&in_time_t, &buf); // Linux
//     localtime_s(&buf, &in_time_t); // Windows

    // get time
    ostringstream ss;
    ss << put_time(&buf, "%Y%m%d_%H%M%S"); // 格式化为YYYYMMDD_HHMMSS
    return ss.str();
}

inline void output_message(string message){
    cout << message << endl;
    logStream << message << endl;
}

// const string create_time = formatted_time();
inline void delete_files_in_directory(string code_out_path) {
    filesystem::path directoryPath = code_out_path; // 替换为实际的目录路径
    if (!filesystem::exists(directoryPath)) {
        cerr << "Directory does not exist: " << directoryPath << endl;
        return;
    }
    try {
        for (const auto& entry : filesystem::directory_iterator(directoryPath)) {
            if (filesystem::is_regular_file(entry.path())) {
                filesystem::remove(entry.path());
            }
        }
        output_message("Files in " + code_out_path + " have been cleared.");
    } catch (const filesystem::filesystem_error& ex) {
        cerr << "Error: " << ex.what() << endl;
    }
}

inline void create_log_path(const string home_path, string dataset_name, string create_time=formatted_time(), string exp_name=""){
    // delete_files_in_directory(home_path + "code/out/");
    // cout << "create_time " << create_time << endl;
    error_code ec;
    filesystem::path dir_path(home_path + "result/" + exp_name + "/" + create_time);
    filesystem::create_directories(dir_path, ec);
    if(ec){
        cerr << "create log directory failed: " << ec.message() << endl;
        return;
    }
    string log_path = home_path + "result/" +  exp_name + "/" + create_time + "/" + exp_name + "_" + dataset_name + "_" + create_time  + ".txt";
    output_message("——————————————Dataset:" + dataset_name + "——————————————");
    logStream.open(log_path);
    if(!logStream.is_open()){
        cerr << "open log file failed" << endl;
    }
}

std::vector<std::string> split_str(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        // 去除前后空格
        token.erase(0, token.find_first_not_of(" \t"));
        token.erase(token.find_last_not_of(" \t") + 1);
        tokens.push_back(token);
    }
    return tokens;
}

// 转换为 double 的函数
inline std::vector<double> parseDoubles(const std::string &line, char delimiter) {
    std::vector<double> values;
    std::vector<std::string> tokens = split_str(line, delimiter);

    for (const auto &token : tokens) {
        try {
            // 尝试将每个 token 转换为 double
            values.push_back(std::stod(token));
        } catch (const std::exception &e) {
            // 如果无法转换，跳过或处理异常
            std::cerr << "Error parsing token: " << token << " (" << e.what() << ")\n";
        }
    }
    return values;
}