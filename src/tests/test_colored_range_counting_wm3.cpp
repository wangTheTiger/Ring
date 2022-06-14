
/*! \file test_colored_range_counting_wm2.cpp
    \brief A more advance test of the CRC WM creation + range search.
    \author Fabrizio Barisione
*/
#include <iostream>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>
#include <unordered_map>
int main()
{
    //1.
    std::cout << ">> Starting test_wm" << std::endl;
    sdsl::int_vector<> S = {1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            3,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            5,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            7,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9,
                            9};
    sdsl::wm_int<sdsl::bit_vector> L;
    std::cout << "Buiding original WM based on S = " << S << std::endl;
    construct_im(L, S);

    //2.
    sdsl::int_vector<> C(S.size());
    std::cout << "Buiding CRC WM with size = " << C.size() << std::endl;
    uint64_t x_s = 0;
    uint64_t x_e = C.size();
    //O ( (x_e - x_s) * log sigma)
    {
        std::unordered_map<uint64_t, uint64_t> hash_map;
        for (uint64_t i = x_s; i < x_e; i++)
        {
            auto it = hash_map.find(L[i]);
            if (it == hash_map.end())
            {
                hash_map.insert({L[i], i});
                C[i - x_s] = 0;
                std::cout << C[i] << " ";
                fflush(stdout);
            }
            else
            {
                C[i - x_s] = it->second + 1;
                std::cout << C[i] << " ";
                fflush(stdout);
                it->second = i;
            }
        }
    }
    std::cout << "" << std::endl;
    //3.
    std::cout << "C = " << C << std::endl;
    //vector<uint64_t> C_aux; // El arreglo "C" del WT.
    sdsl::wm_int<sdsl::bit_vector> wm;
    //sdsl::int_vector<> src;
    construct_im(wm, C);
    std::cout << "  size = " << wm.size() << std::endl;
    std::cout << "  sigma = " << wm.sigma << std::endl;
    std::cout << "  max_level = " << wm.max_level << std::endl; //efectivamente este WT tiene 4 niveles, porque son 11 symbolos, 2⁴ = 16 permite contener los 11 carácteres.
    std::cout << "  rank de 0s = " << wm.rank(wm.size(), 0) << std::endl;
    std::cout << "  rank de 8s = " << wm.rank(wm.size(), 8) << std::endl;
    //! This function counts distinct values on a range. It's based on Muthukrishnan's Colored range counting algorithm.
    //! See algorithm 2 here: https://users.dcc.uchile.cl/~gnavarro/ps/tcs11.2.pdf
    auto res = wm.range_search_2d(1, 10, 0, 0, false);
    auto count = get<0>(res);
    std::cout << " # of Distinct values in range [1, 10] : " << count << std::endl;
    res = wm.range_search_2d(1, 100, 0, 0, false);//! Important: you must start from 1, not from 0.
    count = get<0>(res);
    std::cout << " # of Distinct values in range [1, 100] : " << count << std::endl;
}