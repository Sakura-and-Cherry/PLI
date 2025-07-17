// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2018 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include "piecewise_linear_model.hpp"
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <config.hpp>

#include "data_partition.hpp"

namespace pgm {

#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))
#define PGM_ADD_EPS(x, epsilon, size) ((x) + (epsilon) + 2 >= (size) ? (size) : (x) + (epsilon) + 2)

/**
 * A struct that stores the result of a query to a @ref PGMIndex, that is, a range [@ref lo, @ref hi)
 * centered around an approximate position @ref pos of the sought key.
 */
struct ApproxPos {
    size_t pos; ///< The approximate position of the key.
    size_t lo;  ///< The lower bound of the range.
    size_t hi;  ///< The upper bound of the range.
    size_t epsilon;
};

/**
 * A space-efficient index that enables fast search operations on a sorted sequence of numbers.
 *
 * A search returns a struct @ref ApproxPos containing an approximate position of the sought key in the sequence and
 * the bounds of a range where the sought key is guaranteed to be found if present.
 * If the key is not present, a @ref std::lower_bound search on the range finds a key that is greater or equal to the
 * sought key, if any.
 * In the case of repeated keys, the index finds the position of the first occurrence of a key.
 *
 * The @p Epsilon template parameter should be set according to the desired space-time trade-off. A smaller value
 * makes the estimation more precise and the range smaller but at the cost of increased space usage.
 *
 * Internally the index uses a succinct piecewise linear mapping from keys to their position in the sorted order.
 * This mapping is represented as a sequence of linear models (segments) which, if @p EpsilonRecursive is not zero, are
 * themselves recursively indexed by other piecewise linear mappings.
 *
 * @tparam K the type of the indexed keys
 * @tparam Floating the floating-point type to use for slopes
 */
template<size_t EpsilonRecursive = 4, typename Floating = long double>
class PGMIndex {
public:

    struct Segment;

    size_t n;                           ///< The number of elements this index was built on.
    K first_key;                        ///< The smallest element.
    std::vector<Segment> segments;      ///< The segments composing the index.
    std::vector<size_t> levels_offsets; ///< The starting position of each level in segments[], in reverse order.

    std::vector<std::pair<size_t, size_t>> epsilons_vec; // the epsilon values for each segment

    /// Sentinel value to avoid bounds checking.
    static constexpr K sentinel = std::numeric_limits<K>::has_infinity ? std::numeric_limits<K>::infinity(): std::numeric_limits<K>::max();

    template<typename RandomIt>
    static void build(RandomIt first, RandomIt last,
                      std::vector<Block> blocks,
                      size_t epsilon_recursive,
                      std::vector<Segment> &segments,
                      std::vector<size_t> &levels_offsets,
                      std::vector<std::pair<size_t, size_t>> &epsilons_vec) {
        auto n = (size_t) std::distance(first, last);
        if (n == 0)
            return;

        levels_offsets.push_back(0);
        // segments.reserve(n / (epsilon * epsilon));

        if (*std::prev(last) == sentinel)
            throw std::invalid_argument("The value " + std::to_string(sentinel) + " is reserved as a sentinel.");

        auto build_level = [&](auto epsilon, auto in_fun, auto out_fun, size_t last_n, size_t start_idx = 0, bool is_last_block = true) {
            auto n_segments = internal::make_segmentation_par(last_n, epsilon, in_fun, out_fun);
            if (is_last_block) {
                if (segments.back() == sentinel)
                    --n_segments;
                else {
                    if (segments.back()(sentinel - 1) < last_n + start_idx)
                        segments.emplace_back(*std::prev(last) + 1, 0, last_n + start_idx); // Ensure keys > last are mapped to last_n
                    segments.emplace_back(sentinel, 0, last_n + start_idx);
                }
            }
            return n_segments;
        };

        // Build first level divid{e by blocks

        for (size_t i = 0; i < blocks.size(); ++i) {
            auto & block = blocks[i];
            auto begin_it = first + block.start_idx;
            auto block_size = block.end_idx - block.start_idx;
            auto in_fun = [&](auto i) { return K(begin_it[i]); };
            auto out_fun = [&](auto cs) { segments.emplace_back(cs, block.start_idx); }; // intercept += block.start_idx
            auto last_n = build_level(block.epsilon, in_fun, out_fun, block_size, block.start_idx, (i == blocks.size() - 1)); // 有多少个segment
            epsilons_vec.push_back({last_n, block.epsilon});
        }
        auto last_n = segments.size() - 1; // Last segment is sentinel, so we do not count it

        levels_offsets.push_back(segments.size());

        // Build upper levels
        while (epsilon_recursive && last_n > 1) {
            auto offset = levels_offsets[levels_offsets.size() - 2];
            auto in_fun_rec = [&](auto i) { return segments[offset + i].key; };
            auto out_fun_rec = [&](auto cs) { segments.emplace_back(cs); }; // intercept += block.start_idx
            last_n = build_level(epsilon_recursive, in_fun_rec, out_fun_rec, last_n);
            levels_offsets.push_back(segments.size());
        }
    }

public:
    /**
 * Returns the segment responsible for a given key, that is, the rightmost segment having key <= the sought key.
 * @param key the value of the element to search for
 * @return an iterator to the segment responsible for the given key
 */
    auto segment_for_key(const K &key) const {
        if constexpr (EpsilonRecursive == 0) {
            return std::prev(std::upper_bound(segments.begin(), segments.begin() + segments_count(), key));
        }

        auto it = segments.begin() + *(levels_offsets.end() - 2);
        for (auto l = int(height()) - 2; l >= 0; --l) {
            auto level_begin = segments.begin() + levels_offsets[l];
            auto pos = std::min<size_t>((*it)(key), std::next(it)->intercept);
            auto lo = level_begin + PGM_SUB_EPS(pos, EpsilonRecursive + 1);

            static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Segment);
            if constexpr (EpsilonRecursive <= linear_search_threshold) {
                for (; std::next(lo)->key <= key; ++lo)
                    continue;
                it = lo;
            } else {
                auto level_size = levels_offsets[l + 1] - levels_offsets[l] - 1;
                auto hi = level_begin + PGM_ADD_EPS(pos, EpsilonRecursive, level_size);
                it = std::prev(std::upper_bound(lo, hi, key));
            }
        }
        return it;
    }

    size_t get_data_size() const {
        return n;
    }

    /**
     * Constructs an empty index.
     */
    PGMIndex() = default;

    /**
     * Constructs the index on the given sorted vector.
     * @param data the vector of keys to be indexed, must be sorted
     */
    explicit PGMIndex(const std::vector<K> &data, const std::vector<Block> &blocks) : PGMIndex(data.begin(), data.end(), blocks) {}

    /**
     * Constructs the index on the sorted keys in the range [first, last).
     * @param first, last the range containing the sorted keys to be indexed
     */
    template<typename RandomIt>
    PGMIndex(RandomIt first, RandomIt last, std::vector<Block> blocks)
        : n(std::distance(first, last)),
          first_key(n ? *first : K(0)),
          segments(),
          levels_offsets(),
          epsilons_vec(){
        build(first, last, blocks, EpsilonRecursive, segments, levels_offsets, epsilons_vec);
    }

    /**
     * Returns the approximate position and the range where @p key can be found.
     * @param key the value of the element to search for
     * @return a struct with the approximate position and bounds of the range
     */
    ApproxPos search(const K &key) const {
        auto k = std::max(first_key, key);
        auto it = segment_for_key(k);
        auto idx = (size_t) std::distance(segments.begin(), it) + 1;
        size_t Epsilon;
        for (auto epsilon : epsilons_vec) {
            if (idx <= epsilon.first) {
                Epsilon = epsilon.second;
                break;
            }
            idx -= epsilon.first;
        }
        auto pos = std::min<size_t>((*it)(k), std::next(it)->intercept);
        auto lo = PGM_SUB_EPS(pos, Epsilon);
        auto hi = PGM_ADD_EPS(pos, Epsilon, n);
        return {pos, lo, hi, Epsilon};
    }

    /**
     * Returns the number of segments in the last level of the index.
     * @return the number of segments
     */
    size_t segments_count() const { return segments.empty() ? 0 : levels_offsets[1] - 1; }

    /**
     * Returns the number of levels of the index.
     * @return the number of levels of the index
     */
    size_t height() const { return levels_offsets.size() - 1; }

    /**
     * Returns the size of the index in bytes.
     * @return the size of the index in bytes
     */
    size_t size_in_bytes() const { return segments.size() * sizeof(Segment) + levels_offsets.size() * sizeof(size_t) + epsilons_vec.size() * sizeof(std::pair<size_t, size_t>); }
};

#pragma pack(push, 1)

template<size_t EpsilonRecursive, typename Floating>
struct PGMIndex<EpsilonRecursive, Floating>::Segment {
    K key;              ///< The first key that the segment indexes.
    Floating slope;     ///< The slope of the segment.
    uint32_t intercept; ///< The intercept of the segment.

    Segment() = default;

    Segment(K key, Floating slope, uint32_t intercept) : key(key), slope(slope), intercept(intercept) {};

    explicit Segment(const typename internal::OptimalPiecewiseLinearModel<K, size_t>::CanonicalSegment &cs, size_t start_idx = 0)
        : key(cs.get_first_x()) {
        auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(key);
        if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
            throw std::overflow_error("Change the type of Segment::intercept to uint64");
        if (cs_intercept < 0)
            throw std::overflow_error("Unexpected intercept < 0");
        slope = cs_slope;
        intercept = cs_intercept + start_idx;
    }

    friend inline bool operator<(const Segment &s, const K &k) { return s.key < k; }
    friend inline bool operator<(const K &k, const Segment &s) { return k < s.key; }
    friend inline bool operator<(const Segment &s, const Segment &t) { return s.key < t.key; }

    operator K() { return key; };

    /**
     * Returns the approximate position of the specified key.
     * @param k the key whose position must be approximated
     * @return the approximate position of the specified key
     */
    inline size_t operator()(const K &k) const {
        size_t pos = size_t(slope * double(k - key));
        return pos + intercept;
    }
};

#pragma pack(pop)

}
