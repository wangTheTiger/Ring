/*
 * bwt.hpp
 * Copyright (C) 2020 Author removed for double-blind evaluation
 * 
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

#ifndef BWT_T
#define BWT_T

#include "Config.hpp"
#include <exception> // std::exception
using namespace std;

class bwt
{
    bwt_type L;
    bwt_type crc_L;
    uint64_t num_dist_values;

    C_type *C_bv;
    C_rank_type C_rank;
    C_select_type C_select;
    C_select0_type C_select0;

    std::vector<uint64_t> v1_aux;
    std::vector<uint64_t> v2_aux;

private:
    //! Buiding Colored Range Counting Wavelet Matrix (CRC WM) representation of each Ring's BWT.
    /*!
        * \author Fabrizio Barisione
        * \returns boolean depending whether the CRC WM is created successfully or not. If true then the CRC WM are available as the 'crc_L' member.
        */
    bool build_crc_wm(uint64_t x_s, uint64_t x_e)
    {
        std::cout << "L.sigma : " << L.sigma << std::endl;
        sdsl::int_vector<> C(x_e - x_s);
        std::cout << "Building int vector to store CRC (size = " << C.size() << ")." << std::endl;
        //O ( (x_e - x_s) * log sigma)
        // CORE >>
        std::map<uint64_t, uint64_t> aux_map;
        C[0] = 0;
        {
            aux_map[L[0]] = 0;
            for (uint64_t i = x_s + 1; i < x_e; i++)
            {
                uint64_t l_i = L[i];
                bool found = false;
                for (uint64_t j = aux_map[l_i]; j < i; j++)
                {
                    uint64_t l_j = L[j];
                    if (l_j == l_i)
                    {
                        C[i] = j + 1;
                        aux_map[l_i] = j + 1;
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    //Handling first time found symbols.
                    C[i] = 0;
                    aux_map[l_i] = 0; //aux_map[i] = 0;
                }
            }
        }
        // CORE <<
        std::cout << "Building the CRC WM based on the CRC int vector." << std::endl;
        construct_im(crc_L, C);
        try
        {
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
        * \author Fabrizio Barisione
        * \param x_s TODO:
        * \param x_e TODO:
        * \param rng_s TODO:
        * \param rng_e TODO:
        * \returns TODO:
        */
    uint64_t calculate_number_distinct_values_on_range(uint64_t x_s, uint64_t x_e, uint64_t rng_s, uint64_t rng_e)
    {
        auto res = crc_L.range_search_2d(x_s, x_e, rng_s, rng_e, false);
        return get<0>(res);
    }

public:
    std::vector<uint64_t> C;

    bwt() { ; }

    bwt(int_vector<> &_L, vector<uint64_t> &_C)
    {
        construct_im(L, _L);

        bit_vector C_aux = bit_vector(_C[_C.size() - 1] + 1 + _C.size(), 0);

        for (uint64_t i = 0; i < _C.size(); i++)
        {
            C_aux[_C[i] + i] = 1;
        }

        C_bv = new C_type(C_aux);
        util::init_support(C_rank, C_bv);
        util::init_support(C_select, C_bv);
        util::init_support(C_select0, C_bv);
    }

    ~bwt() { ; }

    uint64_t size()
    {
        return sdsl::size_in_bytes(L) + sdsl::size_in_bytes(*C_bv) + sdsl::size_in_bytes(C_rank) + sdsl::size_in_bytes(C_select) + sdsl::size_in_bytes(C_select0);
    }

    void save(string filename)
    {
        sdsl::store_to_file(L, filename + ".L");
        sdsl::store_to_file(*C_bv, filename + ".C");
    }

    void load(string filename)
    {
        sdsl::load_from_file(L, filename + ".L");
        C_bv = new C_type; //bit_vector;
        sdsl::load_from_file(*C_bv, filename + ".C");
        util::init_support(C_rank, C_bv);
        util::init_support(C_select, C_bv);
        util::init_support(C_select0, C_bv);
        v1_aux = std::vector<uint64_t>(L.sigma);
        v2_aux = std::vector<uint64_t>(L.sigma);
    }

    uint64_t LF(uint64_t i)
    {
        uint64_t s = L[i];
        return get_C(s) + L.rank(i, s) - 1;
    }

    uint64_t nElems(uint64_t val)
    {
        return get_C(val + 1) - get_C(val);
    }

    pair<uint64_t, uint64_t>
    backward_step(uint64_t left_end, uint64_t right_end, uint64_t value)
    {
        uint64_t s = L.rank(left_end, value);
        uint64_t e = L.rank(right_end + 1, value) - 1;
        return pair<uint64_t, uint64_t>(s, e);
    }

    inline uint64_t bsearch_C(uint64_t value)
    {
        uint64_t r = C_rank(C_select0(value + 1));
        return r;
    }

    inline uint64_t get_C(uint64_t v) const
    {
        return C_select(v + 1) - v;
    }

    inline uint64_t ranky(uint64_t pos, uint64_t val)
    {
        return L.rank(pos, val);
    }

    inline uint64_t rank(uint64_t pos, uint64_t val)
    {
        return L.rank(get_C(pos), val);
    }

    inline uint64_t select(uint64_t _rank, uint64_t val)
    {
        return L.select(_rank, val);
    }

    inline std::pair<uint64_t, uint64_t> select_next(uint64_t pos, uint64_t val, uint64_t n_elems)
    {
        return L.select_next(get_C(pos), val, n_elems);
    }

    inline uint64_t min_in_range(uint64_t l, uint64_t r)
    {
        return L.range_minimum_query(l, r);
    }

    inline uint64_t range_next_value(uint64_t x, uint64_t l, uint64_t r)
    {
        return L.range_next_value(x, l, r);
    }

    std::vector<pair<uint64_t, uint64_t>>
    //inline void
    values_in_range(uint64_t pos_min, uint64_t pos_max, uint64_t sigma /*, std::vector<uint64_t> & values, uint64_t & k*/)
    {
        //interval_symbols(L, pos_min, pos_max+1, k, values, r_i, r_j);
        return L.range_search_2d(pos_min, pos_max, 1, sigma).second;
    }

    // backward search for pattern of length 1
    pair<uint64_t, uint64_t> backward_search_1_interval(uint64_t P) const
    {
        return pair<uint64_t, uint64_t>(get_C(P), get_C(P + 1) - 1);
    }

    // backward search for pattern of length 1
    pair<uint64_t, uint64_t> backward_search_1_rank(uint64_t P, uint64_t S) const
    {
        uint64_t s = L.rank(get_C(P), S);
        uint64_t e = L.rank(get_C(P + 1), S);
        return pair<uint64_t, uint64_t>(s, e);
    }

    // backward search for pattern PQ of length 2
    // returns an empty interval if search is unsuccessful
    pair<uint64_t, uint64_t>
    backward_search_2_interval(uint64_t P, pair<uint64_t, uint64_t> &I) const
    {
        return pair<uint64_t, uint64_t>(get_C(P) + I.first, get_C(P) + I.second - 1);
    }

    pair<uint64_t, uint64_t>
    backward_search_2_rank(uint64_t P, uint64_t S, pair<uint64_t, uint64_t> &I) const
    {
        uint64_t c = get_C(P);
        uint64_t s = L.rank(c + I.first, S);
        uint64_t e = L.rank(c + I.second, S);
        return pair<uint64_t, uint64_t>(s, e);
    }

    //! TODO:
    /*!
        * \author Fabrizio Barisione
        * \returns number of distinct values in the WMs. It also sets the value in the member 'num_dist_values' for future references.
        */
    uint64_t calculate_gao()
    {
        num_dist_values = 0;
        //Build the crc wm for the entire original WT TODO: in the future this will be part of an adaptive algorithm.
        if (!build_crc_wm(0, L.size()))
        {
            num_dist_values = calculate_number_distinct_values_on_range(0, L.size(), 0, 0);
        }
        std::cout << "Num of distinct values : " << num_dist_values << std::endl;
        return num_dist_values;
    }
};
#endif
