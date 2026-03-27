// #######################################
// # Copyright (C) 2020-2023 Otmar Ertl. #
// # All rights reserved.                #
// #######################################

#ifndef _RABBITSKETCH_TMH_WEIGHTED_MINWISE_HPP_
#define _RABBITSKETCH_TMH_WEIGHTED_MINWISE_HPP_

#include "tmh_bitstream.hpp"

#include <algorithm>
#include <cstring>
#include <cmath>
#include <numeric>
#include <vector>

namespace tmh {

template <typename I, typename W> I calculateMaxBoundIdx() {
  static_assert(sizeof(I) == sizeof(W),
                "index type and weight type do not have same size");
  I i = 0;
  const W w = std::numeric_limits<W>::max();
  memcpy(&i, &w, sizeof(I));
  return i;
}

template <typename W, typename I> class WeightDiscretization {
public:
  typedef I index_type;
  typedef W weight_type;

  static const index_type maxBoundIdx;

  static weight_type getBound(index_type boundIdx) {
    W f;
    static_assert(std::numeric_limits<W>::is_iec559,
                  "weight_type is not iec559");
    static_assert(sizeof(weight_type) == sizeof(index_type),
                  "weight_type and index_type do not have same size");
    memcpy(&f, &boundIdx, sizeof(index_type));
    return f;
  }
};

template <typename W, typename I>
const typename WeightDiscretization<W, I>::index_type
    WeightDiscretization<W, I>::maxBoundIdx =
        calculateMaxBoundIdx<WeightDiscretization<W, I>::index_type,
                             WeightDiscretization<W, I>::weight_type>();

typedef WeightDiscretization<float, uint32_t> FloatWeightDiscretization;

typedef WeightDiscretization<double, uint64_t> DoubleWeightDiscretization;

struct UnaryWeightFunction {
  template <typename X> constexpr double operator()(X) const { return 1; }
};

struct Node {
  double lowerBound;
  double invRate;
  double ratio;

  Node(double lowerBound, double midBound, double upperBound)
      : lowerBound(lowerBound), invRate(1. / (upperBound - lowerBound)),
        ratio(invRate * (midBound - lowerBound)) {}
};

static std::vector<Node> preCalculateTree(double factor, double max) {
  assert(max > 0);
  assert(factor > 0);
  assert(factor < 1);

  struct TmpNode {
    uint32_t minIdx;
    uint32_t midIdx;
    uint32_t maxIdx;
  };

  std::vector<double> boundaries;

  boundaries.push_back(max);
  do {
    boundaries.push_back(std::min(factor * boundaries.back(),
                                  std::nexttoward(boundaries.back(), 0)));
  } while (boundaries.back() > 0);
  assert(boundaries.back() == 0.);
  std::reverse(boundaries.begin(), boundaries.end());
  uint32_t numNodes = boundaries.size() - 1;

  std::vector<uint32_t> counts(2 * numNodes - 1);
  for (uint32_t i = 0; i < numNodes; ++i)
    counts[numNodes - 1 + i] = 1;
  for (uint32_t idx = numNodes - 2; idx != std::numeric_limits<uint32_t>::max();
       --idx) {
    uint32_t leftChildIdx = 2 * idx + 1;
    uint32_t rightChildIdx = 2 * idx + 2;
    counts[idx] = counts[leftChildIdx] + counts[rightChildIdx];
  }
  assert(counts[0] == numNodes);

  std::vector<TmpNode> tmpNodes(2 * numNodes - 1);
  tmpNodes[0] = {0, 0, numNodes};

  for (uint32_t idx = 0; idx < numNodes - 1; ++idx) {
    uint32_t leftChildIdx = 2 * idx + 1;
    uint32_t rightChildIdx = 2 * idx + 2;
    assert(tmpNodes[idx].maxIdx - tmpNodes[idx].minIdx ==
           counts[leftChildIdx] + counts[rightChildIdx]);
    tmpNodes[idx].midIdx = tmpNodes[idx].minIdx + counts[leftChildIdx];
    tmpNodes[leftChildIdx] = {tmpNodes[idx].minIdx, 0, tmpNodes[idx].midIdx};
    tmpNodes[rightChildIdx] = {tmpNodes[idx].midIdx, 0, tmpNodes[idx].maxIdx};
  }

  for (uint32_t idx = 0; idx < numNodes - 1; ++idx) {
    uint32_t leftChildIdx = 2 * idx + 1;
    uint32_t rightChildIdx = 2 * idx + 2;
    assert(tmpNodes[leftChildIdx].minIdx == tmpNodes[idx].minIdx);
    assert(tmpNodes[rightChildIdx].maxIdx == tmpNodes[idx].maxIdx);
    assert(tmpNodes[leftChildIdx].maxIdx == tmpNodes[rightChildIdx].minIdx);
  }

  for (uint32_t idx = numNodes - 1; idx < 2 * numNodes - 1; ++idx) {
    assert(tmpNodes[idx].minIdx + 1 == tmpNodes[idx].maxIdx);
    assert(tmpNodes[idx].midIdx == 0);
  }

  std::vector<Node> tree;
  tree.reserve(2 * numNodes - 1);
  for (uint32_t idx = 0; idx < 2 * numNodes - 1; ++idx) {
    tree.emplace_back(
        boundaries[tmpNodes[idx].minIdx],
        (tmpNodes[idx].midIdx == 0) ? -1. : boundaries[tmpNodes[idx].midIdx],
        boundaries[tmpNodes[idx].maxIdx]);
  }

  return tree;
}

template <typename RF> class TreeMinHash {

  typedef typename RF::RngType R;

  const uint32_t m;
  const RF rngFunction;
  const std::vector<Node> tree;
  const uint32_t numNonLeafNodes;
  const double initialLimitFactor;
  std::vector<std::pair<double, uint32_t>> buffer;
  PermutationStream permutationStream;
  std::vector<double> factors;

public:
  TreeMinHash(const uint32_t m, const RF &rngFunction,
              double max = std::numeric_limits<double>::max(),
              double factor = 0.5, double successProbabilityFirstRun = 0.9)
      : m(m), rngFunction(rngFunction), tree(preCalculateTree(factor, max)),
        numNonLeafNodes(tree.size() - (tree.size() + 1) / 2),
        initialLimitFactor(
            -std::log(-std::expm1(std::log(successProbabilityFirstRun) / m)) *
            m),
        permutationStream(m), factors(m - 1) {
    buffer.reserve(numNonLeafNodes);
    for (uint32_t i = 0; i < m - 1; ++i)
      factors[i] = static_cast<double>(m) / static_cast<double>(m - i - 1);
  }

  // In-place fill: avoids per-call allocation of the result vector.
  // Uses an unfilled counter (O(1)) instead of a none_of scan (O(m)).
  void fill(const std::vector<std::pair<uint64_t, double>> &data,
            std::vector<std::pair<uint64_t, double>> &result,
            uint64_t *iterationCounter = nullptr) {

    const double weightSum =
        std::accumulate(data.begin(), data.end(), 0.,
                        [](double x, const std::pair<uint64_t,double> &d) {
                          return x + d.second;
                        });

    const double limitIncrement = initialLimitFactor / weightSum;
    double limit = limitIncrement;

    result.assign(m, {UINT64_C(0), limit});
    uint32_t unfilled = m;   // slots still at 'limit'

    if (iterationCounter != nullptr)
      *iterationCounter = 1;

    while (true) {

      for (const std::pair<uint64_t,double>& dw : data) {
        const uint64_t d = dw.first;
        const double w = dw.second;

        if (w == 0)
          continue;
        buffer.clear();
        uint32_t nodeIdx = 0;
        auto rng = rngFunction(d, nodeIdx);
        double point = getExponential1(rng) * tree[nodeIdx].invRate;
        if (!(point < limit))
          continue;

        while (true) {
          while (nodeIdx < numNonLeafNodes &&
                 tree[nodeIdx].lowerBound < w) {
            const bool inheritToLeft = getBernoulli(tree[nodeIdx].ratio, rng);
            nodeIdx <<= 1;
            const uint32_t siblingIdx = nodeIdx + 1 + inheritToLeft;
            nodeIdx += 2 - inheritToLeft;
            double siblingPoint =
                point + getExponential1(rng) * tree[siblingIdx].invRate;
            if (siblingPoint < limit && tree[siblingIdx].lowerBound < w) {
              buffer.emplace_back(siblingPoint, siblingIdx);
            }
          }
          if (tree[nodeIdx].lowerBound < w) {
            const double invRate = tree[nodeIdx].invRate;
            const double acceptanceProbability =
                (w - tree[nodeIdx].lowerBound) * invRate;

            permutationStream.reset();
            for (uint32_t k = 0; k < m; ++k) {
              double nextPoint =
                  (k < factors.size())
                      ? point + getExponential1(rng) * invRate * factors[k]
                      : std::numeric_limits<double>::infinity();
              const uint32_t idx = permutationStream.next(rng);
              if (point < result[idx].second) {
                const bool was_unfilled = (result[idx].second >= limit);
                bool accepted = false;
                if (!(acceptanceProbability < 1)) {
                  result[idx] = {d, point};
                  accepted = true;
                } else {
                  auto rng2 = rngFunction(
                      d, (static_cast<uint64_t>(nodeIdx) << 32) | idx);
                  do {
                    if (getUniformDouble(rng2) < acceptanceProbability) {
                      result[idx] = {d, point};
                      accepted = true;
                      break;
                    }
                    point += getExponential1(rng2) * invRate * m;
                  } while (point < result[idx].second);
                }
                if (was_unfilled && accepted) --unfilled;
              }
              if (!(nextPoint < limit))
                break;
              point = nextPoint;
            }
          }
          if (buffer.empty())
            break;
          std::tie(point, nodeIdx) = buffer.back();
          buffer.pop_back();
          rng = rngFunction(d, nodeIdx);
        }
      }

      if (unfilled == 0)
        return;

      if (iterationCounter != nullptr)
        (*iterationCounter) += 1;
      double oldLimit = limit;
      limit += limitIncrement;
      for (std::pair<uint64_t,double>& r : result) {
        if (r.second == oldLimit) r.second = limit;
      }
    }
  }

  // Legacy returning version (calls fill internally)
  std::vector<std::pair<uint64_t, double>>
  operator()(const std::vector<std::pair<uint64_t, double>> &data,
             uint64_t *iterationCounter = nullptr) {
    std::vector<std::pair<uint64_t, double>> result;
    result.resize(m);
    fill(data, result, iterationCounter);
    return result;
  }
};

} // namespace tmh

#endif // _RABBITSKETCH_TMH_WEIGHTED_MINWISE_HPP_
