
/*! \file test_colored_range_counting_wm_parameter.cpp
    \brief A more advance test of the CRC WM creation + range search.
    \author Fabrizio Barisione
*/
#include<iostream>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>
#include <unordered_map>
#include "../crc_arrays.hpp"

int main(int argc, char* argv[])
{
    //1.
    std::cout << ">> Starting test_wm" << std::endl << " loading CRC array of " << argv[1] << std::endl;
    crc_arrays crc_arrays;
    crc_arrays.load(string(argv[1]));
    crc_arrays.print_arrays();
    //2. Quering it.
    //! This function counts distinct values on a range. It's based on Muthukrishnan's Colored range counting algorithm.
    //! See algorithm 2 here: https://users.dcc.uchile.cl/~gnavarro/ps/tcs11.2.pdf

    //For instance, compare it with sorted data found in dat/wikidata-filtered-enumerated-reduced.dat.
    std::cout << "SPO tests " << std::endl;
    auto num_dist_values = crc_arrays.spo_BWT_S->get_number_distinct_values(1, 9);
    std::cout << " L_s # of Distinct values in range [1, 9] : " << num_dist_values << std::endl;
    num_dist_values = crc_arrays.spo_BWT_S->get_number_distinct_values(1, 21);
    std::cout << " L_s # of Distinct values in range [1, 21] : " << num_dist_values << std::endl;
    num_dist_values = crc_arrays.spo_BWT_S->get_number_distinct_values(1, 65);
    std::cout << " L_s # of Distinct values in range [1, 65] : " << num_dist_values << std::endl;
    num_dist_values = crc_arrays.spo_BWT_S->get_number_distinct_values(39, 65);
    std::cout << " L_s # of Distinct values in range [39, 65] : " << num_dist_values << std::endl;

    num_dist_values = crc_arrays.spo_BWT_P->get_number_distinct_values(1, 20);
    std::cout << " L_p # of Distinct values in range [1, 20] : " << num_dist_values << std::endl;

    num_dist_values = crc_arrays.spo_BWT_O->get_number_distinct_values(1, 20);
    std::cout << " L_o # of Distinct values in range [1, 20] : " << num_dist_values << std::endl;
}