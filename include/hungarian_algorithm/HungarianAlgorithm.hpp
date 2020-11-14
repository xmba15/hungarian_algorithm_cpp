/**
 * @file    HungarianAlgorithm.hpp
 *
 * @author  btran
 *
 */

#pragma once

#include <vector>

namespace hungarian
{
/**
 *  @brief solve matching problems of numRows agents to numCols jobs to maximize total outputs
 *
 *  @return indices of matched jobs to each row agent. index -1 means no match
 *
 */
std::vector<int> solve(const double* costMatrix, int numRows, int numCols);
}  // namespace hungarian
