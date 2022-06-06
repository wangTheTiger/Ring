/*
 * crc_array.hpp
 * Copyright (C) 2022
 *
 * This is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CRC_ARRAY
#define CRC_ARRAY

#include "Config.hpp"
#include <exception> // std::exception
#include <unordered_map>

class crc
{
private:
    bwt_type L;
    bwt_type crc_L;
    bool build_crc_wm(uint64_t x_s, uint64_t x_e);

public:
    crc() = default;
    crc(bwt_type &wm_l);
    void load(string filename);
    void save(string filename);
    uint64_t get_number_distinct_values_on_range(uint64_t x_s, uint64_t x_e, uint64_t rng_s, uint64_t rng_e);
    uint64_t get_number_distinct_values(uint64_t l, uint64_t r);
};

void crc::save(string filename)
{
    sdsl::store_to_file(crc_L, filename + ".crc");
}

void crc::load(string filename)
{
    sdsl::load_from_file(crc_L, filename + ".crc");
}
crc::crc(bwt_type &wm_l)
{
    L = wm_l;
    build_crc_wm(0, L.size() - 1);
}
//! Buiding Colored Range Counting Array (CRC WM) based on the given Wavelet Matrix.
/*!
 * \returns boolean depending whether the CRC WM is created successfully or not. If true then the CRC WM are available as the 'crc_L' member.
 */
bool crc::build_crc_wm(uint64_t x_s, uint64_t x_e)
{
    sdsl::int_vector<> C(x_e - x_s);
    std::cout << "L.sigma : " << L.sigma << ", L.size() : " << L.size() << " x_s : " << x_s << " x_e : " << x_e << std::endl;
    std::cout << "Building int vector to store CRC (size = " << C.size() << ")." << std::endl;
    // O ( (x_e - x_s) * log sigma)
    //  CORE >>
    {
        std::unordered_map<uint64_t, uint64_t> hash_map;
        for (uint64_t i = x_s; i < x_e; i++)
        {
            auto it = hash_map.find(L[i]);
            if (it == hash_map.end())
            {
                hash_map.insert({L[i], i});
                // C positions must start from 0 until x_e - x_s.
                C[i - x_s] = 0;
                // std::cout << C[i - x_s] << " ";fflush(stdout);
            }
            else
            {
                C[i - x_s] = it->second + 1;
                it->second = i;
                // std::cout << C[i - x_s] << " ";fflush(stdout);
            }
        }
    }
    // CORE <<
    // std::cout << "C = " << C << std::endl;
    // std::cout << "Building the CRC WM based on the CRC int vector." << std::endl;
    construct_im(crc_L, C);
    try
    {
        // std::cout << "CRC WM size : " << crc_L.size() << std::endl;
        if (crc_L.size() > 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    catch (std::exception &e)
    {
        return false;
    }
}
//! TODO:
/*!
 * \param x_s TODO:
 * \param x_e TODO:
 * \param rng_s TODO:
 * \param rng_e TODO:
 * \returns TODO:
 */
uint64_t crc::get_number_distinct_values_on_range(uint64_t x_s, uint64_t x_e, uint64_t rng_s, uint64_t rng_e)
{
    // std::cout << "get_number_distinct_values_on_range [" << x_s << ", " << x_e << "]" << std::endl;
    auto res = crc_L.range_search_2d(x_s, x_e, rng_s, rng_e, false);
    return get<0>(res);
}

//! TODO: fix comment -> First it builds a compact Colored Range Counting representation (CRC) of the BWT.L and then, it calculates the number of different values. Both CRC and BWT.L are Wavelet Matrices.
//! In a three dimension ring this function should be called 3 times.
/*!
 * \param uint64_t l : left
 * \param uint64_t r : right
 * \returns number of distinct values in the WMs. It also sets the value in the member 'num_dist_values' for future references.
 */
uint64_t crc::get_number_distinct_values(uint64_t l, uint64_t r)
{
    // std::cout << "Calling get_number_distinct_values with range : [" << l << ", " << r << "]." << std::endl;
    uint64_t num_dist_values = 0;
    // Build the crc wm for the entire original WT TODO: in the future this will be part of an adaptive algorithm.
    //bool result = build_crc_wm(l, r);
    uint64_t rng_s = 0;
    uint64_t rng_e = (r == 0) ? 0 : r - 1;

    num_dist_values = get_number_distinct_values_on_range(r, l, rng_s, rng_e);
    //num_dist_values = get_number_distinct_values_on_range(0, crc_L.size() - 1, 0, 0);

    // std::cout << "Num of distinct values : " << num_dist_values << std::endl;
    return num_dist_values;
}
#endif
