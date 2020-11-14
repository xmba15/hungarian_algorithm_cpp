/**
 * @file    App.cpp
 *
 * @author  btran
 *
 */

#include <iostream>
#include <vector>

#include <hungarian_algorithm/HungarianAlgorithm.hpp>

int main(int argc, char* argv[])
{
    std::vector<double> cost = {
        7,   53,  183, 439, 863,  //
        497, 383, 563, 79,  973,  //
        287, 63,  343, 169, 583,  //
        627, 343, 773, 959, 943,  //
        767, 473, 103, 699, 303   //
    };

    int numRows = 5;
    int numCols = 5;
    auto assignments = hungarian::solve(cost.data(), numRows, numCols);

    for (int i = 0; i < numRows; ++i) {
        std::cout << "match: " << i << " " << assignments[i] << "\n";
    }

    return EXIT_SUCCESS;
}
