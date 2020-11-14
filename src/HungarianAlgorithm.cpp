/**
 * @file    HungrianAlgorithm.cpp
 *
 * @author  btran
 *
 */

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <vector>

#include <hungarian_algorithm/HungarianAlgorithm.hpp>

namespace hungarian
{
namespace
{
std::vector<double> padSquare(const double* costs, int numRows, int numCols);

/**
 *  @brief Hungarian Assignment Algorithm using Maximum Biparte Matching
 *
 *  Supposed U, V are two partitions of the bipartite graph G=(V,E)
 *
 */
class HungarianAlgorithm
{
 public:
    HungarianAlgorithm(int dim)
        : m_dim(dim)
        , m_minSlack(dim)
        , m_lu(dim)
        , m_lv(dim, 0)
    {
    }

    const auto& Mu() const
    {
        return m_Mu;
    }

    void run(const std::vector<double>& costs) const
    {
        // initial feasible labeling
#pragma omp parallel for
        for (int i = 0; i < m_dim; ++i) {
            m_lu[i] = *std::max_element(costs.data() + i * m_dim, costs.data() + (i + 1) * m_dim);
        }

        while (static_cast<int>(m_Mu.size()) < m_dim) {
            int u0 = -1;

            // pick a free vertex u0 belonging to U
            for (int u = 0; u < m_dim; ++u) {
                if (m_Mu.find(u) == m_Mu.end()) {
                    u0 = u;
                    break;
                }
            }
            std::vector<int> S{u0};

#pragma omp parallel for
            for (int v = 0; v < m_dim; ++v) {
                m_minSlack[v] = std::make_pair(this->slack(u0, v, costs), u0);
            }

            std::unordered_map<int, int> T;
            this->augment(costs, S, T);
        }
    }

 private:
    double slack(int u, int v, const std::vector<double>& costs) const
    {
        return m_lu[u] + m_lv[v] - costs[u * m_dim + v];
    }

    int argMinSlack(const std::unordered_map<int, int>& T) const
    {
        int minV = -1;
        double minW = std::numeric_limits<double>::max();
        for (int v = 0; v < m_dim; ++v) {
            if (T.find(v) != T.end()) {
                continue;
            }
            if (m_minSlack[v].first < minW) {
                minV = v;
                minW = m_minSlack[v].first;
            }
        }

        return minV;
    }

    void augment(const std::vector<double>& costs, std::vector<int>& S, std::unordered_map<int, int>& T) const
    {
        while (true) {
            int curV = this->argMinSlack(T);
            auto minElem = m_minSlack[curV];

            if (minElem.first > 0) {
                this->improveLabels(minElem.first, S, T);
            }

            T[curV] = minElem.second;

            if (m_Mv.find(curV) != m_Mv.end()) {
                int u1 = m_Mv[curV];
                S.emplace_back(u1);
                for (int v = 0; v < m_dim; ++v) {
                    if (T.find(v) != T.end()) {
                        continue;
                    }
                    double slackU1 = this->slack(u1, v, costs);
                    if (m_minSlack[v].first > slackU1) {
                        m_minSlack[v] = std::make_pair(slackU1, u1);
                    }
                }
            } else {
                this->improveMatching(curV, T);
                return;
            }
        }
    }

    void improveMatching(int v, const std::unordered_map<int, int>& T) const
    {
        int u = T.at(v);
        if (m_Mu.find(u) != m_Mu.end()) {
            this->improveMatching(m_Mu.at(u), T);
        }
        m_Mu[u] = v;
        m_Mv[v] = u;
    }

    void improveLabels(double val, const std::vector<int>& S, const std::unordered_map<int, int>& T) const
    {
        for (int u : S) {
            m_lu[u] -= val;
        }

#pragma omp parallel for
        for (int v = 0; v < m_dim; ++v) {
            if (T.find(v) != T.end()) {
                m_lv[v] += val;
            } else {
                m_minSlack[v].first -= val;
            }
        }
    }

 private:
    int m_dim;

    using SlackPair = std::pair<double, int>;
    mutable std::vector<SlackPair> m_minSlack;

    mutable std::vector<double> m_lu;  // label row vertices
    mutable std::vector<double> m_lv;  // label column vertices

    // mapping
    mutable std::unordered_map<int, int> m_Mu;
    mutable std::unordered_map<int, int> m_Mv;
};
}  // namespace

std::vector<int> solve(const double* costs, int numRows, int numCols)
{
    auto paddedCosts = padSquare(costs, numRows, numCols);
    int maxDim = std::max<int>(numRows, numCols);
    HungarianAlgorithm solver(maxDim);
    solver.run(paddedCosts);

    std::vector<int> assignments(numRows);

    int minDim = std::min<int>(numRows, numCols);

#pragma omp parallel for
    for (int i = 0; i < numRows; ++i) {
        int val = solver.Mu().at(i);
        assignments[i] = val < minDim ? val : -1;
    }

    return assignments;
}

namespace
{
std::vector<double> padSquare(const double* costs, int numRows, int numCols)
{
    double minElem = *std::min_element(costs, costs + numRows * numCols);

    int maxDim = std::max<int>(numRows, numCols);
    std::vector<double> paddedCosts(maxDim * maxDim, minElem - 1);

#pragma omp parallel for
    for (int i = 0; i < numRows; ++i) {
        std::copy(costs + i * numCols, costs + (i + 1) * numCols, paddedCosts.data() + i * maxDim);
    }

    return paddedCosts;
}
}  // namespace
}  // namespace hungarian
