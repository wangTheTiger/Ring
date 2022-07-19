#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <set>
#include "Triple.h"
#include "ring_spo.hpp"
#include "ring_sop.hpp"
#include <unordered_map>
#include "crc_arrays.hpp"

uint64_t get_size_interval(Triple * triple_pattern, ring_spo & graph) {
    if (triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph.open_SPO();
        return open_interval.size();
    } else if (triple_pattern->s->isVariable && !triple_pattern->p->isVariable && triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph.open_PSO();
        uint64_t cur_p = graph.min_P(open_interval);
        cur_p = graph.next_P(open_interval, triple_pattern->p->constant);
        if (cur_p == 0 || cur_p != triple_pattern->p->constant) {
            return 0;
        } else{
            bwt_interval i_p = graph.down_P(cur_p);
            return i_p.size();
        }
    } else if (triple_pattern->s->isVariable && triple_pattern->p->isVariable && !triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph.open_OSP();
        uint64_t cur_o = graph.min_O(open_interval);
        cur_o = graph.next_O(open_interval, triple_pattern->o->constant);
        if (cur_o == 0 || cur_o != triple_pattern->o->constant) {
            return 0;
        } else{
            bwt_interval i_s = graph.down_O(cur_o);
            return i_s.size();
        }
    } else if (!triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph.open_SPO();
        uint64_t cur_s = graph.min_S(open_interval);
        cur_s = graph.next_S(open_interval, triple_pattern->s->constant);
        if (cur_s == 0 || cur_s != triple_pattern->s->constant) {
            return 0;
        } else{
            bwt_interval i_s = graph.down_S(cur_s);
            return i_s.size();
        }
    } else if (!triple_pattern->s->isVariable && !triple_pattern->p->isVariable && triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph.open_SPO();
        uint64_t cur_s = graph.min_S(open_interval);
        cur_s = graph.next_S(open_interval, triple_pattern->s->constant);
        if (cur_s == 0 || cur_s != triple_pattern->s->constant) {
            return 0;
        } else{
            bwt_interval i_s = graph.down_S(cur_s);
            uint64_t cur_p = graph.min_P_in_S(i_s, cur_s);
            cur_p = graph.next_P_in_S(i_s, cur_s, triple_pattern->p->constant);
            if (cur_p == 0 || cur_p != triple_pattern->p->constant) {
              return 0;
            }
            bwt_interval i_p = graph.down_S_P(i_s, cur_s, cur_p);
            return i_p.size();
        }
    } else if (!triple_pattern->s->isVariable && triple_pattern->p->isVariable && !triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph.open_SOP();
        uint64_t cur_s = graph.min_S(open_interval);
        cur_s = graph.next_S(open_interval, triple_pattern->s->constant);
        if (cur_s == 0 || cur_s != triple_pattern->s->constant) {
            return 0;
        } else{
            bwt_interval i_s = graph.down_S(cur_s);
            uint64_t cur_o = graph.min_O_in_S(i_s);
            cur_o = graph.next_O_in_S(i_s, triple_pattern->o->constant);
            if (cur_o == 0 || cur_o != triple_pattern->o->constant) {
              return 0;
            }
            bwt_interval i_o = graph.down_S_O(i_s, cur_o);
            return i_o.size();
        }
    } else if (triple_pattern->s->isVariable && !triple_pattern->p->isVariable && !triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph.open_POS();
        uint64_t cur_p = graph.min_P(open_interval);
        cur_p = graph.next_P(open_interval, triple_pattern->p->constant);
        if (cur_p == 0 || cur_p != triple_pattern->p->constant) {
            return 0;
        } else{
            bwt_interval i_p = graph.down_P(cur_p);
            uint64_t cur_o = graph.min_O_in_P(i_p, cur_p);
            cur_o = graph.next_O_in_P(i_p, cur_p, triple_pattern->o->constant);
            if (cur_o == 0 || cur_o != triple_pattern->o->constant) {
              return 0;
            }
            bwt_interval i_o = graph.down_P_O(i_p, cur_p, cur_o);
            return i_o.size();
        }
    }
    return 0;
}

//! TODO:
/*!
 * \returns TODO:
 */
std::unordered_map<string, uint64_t> get_num_diff_values(Triple *triple_pattern, ring_spo &graph_spo, crc_arrays& crc_arrays)
{
    std::unordered_map<string, uint64_t> hash_map;
    const uint64_t nTriples = graph_spo.get_n_triples();
    //First case: ?S ?P ?O
    if (triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable)
    {
        bwt_interval open_interval = graph_spo.open_SPO();
        // TODO: do something different here.
        hash_map.insert({triple_pattern->s->varname, open_interval.size()});
        hash_map.insert({triple_pattern->p->varname, open_interval.size()});
        hash_map.insert({triple_pattern->o->varname, open_interval.size()});
        //std::cout << "num_distinct_values S = " << open_interval.size() << " vs. interval size = " << open_interval.size() << std::endl;
        //std::cout << "num_distinct_values P = " << open_interval.size() << " vs. interval size = " << open_interval.size() << std::endl;
        //std::cout << "num_distinct_values O = " << open_interval.size() << " vs. interval size = " << open_interval.size() << std::endl;
        return hash_map;
    }
    //Second case : ?S P ?O
    else if (triple_pattern->s->isVariable && !triple_pattern->p->isVariable && triple_pattern->o->isVariable)
    {
        //we need to do graph_spo.get_number_distinct_values_BWT_S() and for ?O we need the inverse ring, that is graph_sop.get_number_distinct_values_BWT_O().
        bwt_interval open_interval = graph_spo.open_PSO();
        //First we get the minimum value of the range
        uint64_t cur_p = graph_spo.min_P(open_interval);
        //Then the next value >= triple_pattern->p->constant.
        cur_p = graph_spo.next_P(open_interval, triple_pattern->p->constant);
        if (cur_p == 0 || cur_p != triple_pattern->p->constant)
        {
            hash_map.insert({triple_pattern->s->varname, 0});
            hash_map.insert({triple_pattern->o->varname, 0});
            return hash_map;//return 0
        }
        else
        {
            bwt_interval i_p = graph_spo.down_P(cur_p); // Range in C array pointing to the S value.
            // Ring => Going from P to S: values must be shifted back to the left, by subtracting nTriples. Remember P is between nTriples + 1 and 2*nTriples.
            uint64_t num_distinct_values_s = crc_arrays.get_number_distinct_values_spo_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples);
            // Reverse Ring => Going from P to O: values must be shifted back to the left, by subtracting nTriples. Remember P is between 2 * nTriples + 1 and 3*nTriples.
            // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
            uint64_t num_distinct_values_o = crc_arrays.get_number_distinct_values_sop_BWT_O(i_p.left() - nTriples, i_p.right() - nTriples);
            //std::cout << "num_distinct_values S = " << num_distinct_values_s << " vs. interval size = " << i_p.size() << std::endl;
            //std::cout << "num_distinct_values O = " << num_distinct_values_o << " vs. interval size = " << i_p.size() << std::endl;
            hash_map.insert({triple_pattern->s->varname, num_distinct_values_s});
            hash_map.insert({triple_pattern->o->varname, num_distinct_values_o});
            return hash_map;
        }
    }
    //Third case ?S ?P O
    else if (triple_pattern->s->isVariable && triple_pattern->p->isVariable && !triple_pattern->o->isVariable)
    {
        bwt_interval open_interval = graph_spo.open_OSP();
        uint64_t cur_o = graph_spo.min_O(open_interval);
        cur_o = graph_spo.next_O(open_interval, triple_pattern->o->constant);
        if (cur_o == 0 || cur_o != triple_pattern->o->constant)
        {
            hash_map.insert({triple_pattern->p->varname, 0});
            hash_map.insert({triple_pattern->s->varname, 0});
            return hash_map;//return 0;
        }
        else
        {
            bwt_interval i_o = graph_spo.down_O(cur_o); // Range in C array pointing to the O value.
            // Ring => Going from O to P: values must be shifted back to the left, by subtracting nTriples. Remember O is between 2* nTriples + 1 and 3*nTriples.
            uint64_t num_distinct_values_p = crc_arrays.get_number_distinct_values_spo_BWT_P(i_o.left() - nTriples, i_o.right() - nTriples);
            // Reverse Ring => Going from O to S: values must be shifted back to the left, by subtracting nTriples. Remember O is between nTriples + 1 and 2*nTriples.
            // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
            uint64_t num_distinct_values_s = crc_arrays.get_number_distinct_values_sop_BWT_S(i_o.left() - nTriples, i_o.right() - nTriples);
            //std::cout << "num_distinct_values S = " << num_distinct_values_s << " vs. interval size = " << i_o.size() << std::endl;
            //std::cout << "num_distinct_values P = " << num_distinct_values_p << " vs. interval size = " << i_o.size() << std::endl;
            hash_map.insert({triple_pattern->p->varname, num_distinct_values_p});
            hash_map.insert({triple_pattern->s->varname, num_distinct_values_s});
            return hash_map;
        }
    }
    //Fourth case S ?P ?O
    else if (!triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable)
    {
        bwt_interval open_interval = graph_spo.open_SPO();
        uint64_t cur_s = graph_spo.min_S(open_interval);
        cur_s = graph_spo.next_S(open_interval, triple_pattern->s->constant);
        if (cur_s == 0 || cur_s != triple_pattern->s->constant)
        {
            hash_map.insert({triple_pattern->o->varname, 0});
            hash_map.insert({triple_pattern->p->varname, 0});
            return hash_map;//return 0;
        }
        else
        {
            bwt_interval i_s = graph_spo.down_S(cur_s);// Range in C array pointing to the S value.
            // Ring => Going from S to O: values must be shifted to the right, by adding 2 * nTriples. Remember S is between  1 and nTriples.
            uint64_t num_distinct_values_o = crc_arrays.get_number_distinct_values_spo_BWT_O(i_s.left() + 2 * nTriples, i_s.right() + 2 * nTriples);
            // Reverse Ring => Going from S to P:  values must be shifted to the right, by adding 2 * nTriples. Remember S is between  1 and nTriples.
            // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
            uint64_t num_distinct_values_p = crc_arrays.get_number_distinct_values_sop_BWT_P(i_s.left() + 2 * nTriples, i_s.right() + 2 * nTriples);
            //std::cout << "num_distinct_values P = " << num_distinct_values_p << " vs. interval size = " << i_s.size() << std::endl;
            //std::cout << "num_distinct_values O = " << num_distinct_values_o << " vs. interval size = " << i_s.size() << std::endl;
            hash_map.insert({triple_pattern->o->varname, num_distinct_values_o});
            hash_map.insert({triple_pattern->p->varname, num_distinct_values_p});
            return hash_map;
        }
    }
    //Fifth case S P ?O
    else if (!triple_pattern->s->isVariable && !triple_pattern->p->isVariable && triple_pattern->o->isVariable)
    {
        bwt_interval open_interval = graph_spo.open_SPO();
        uint64_t cur_s = graph_spo.min_S(open_interval);
        cur_s = graph_spo.next_S(open_interval, triple_pattern->s->constant);
        if (cur_s == 0 || cur_s != triple_pattern->s->constant)
        {
            hash_map.insert({triple_pattern->o->varname, 0});
            return hash_map;//return 0;
        }
        else
        {
            bwt_interval i_s = graph_spo.down_S(cur_s);
            uint64_t cur_p = graph_spo.min_P_in_S(i_s, cur_s);
            cur_p = graph_spo.next_P_in_S(i_s, cur_s, triple_pattern->p->constant);
            if (cur_p == 0 || cur_p != triple_pattern->p->constant)
            {
                hash_map.insert({triple_pattern->o->varname, 0});
                return hash_map;//return 0;
            }
            bwt_interval i_p = graph_spo.down_S_P(i_s, cur_s, cur_p);
            // Ring => Going from P to O: values must be shifted to the right, by adding nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
            uint64_t num_distinct_values_o = crc_arrays.get_number_distinct_values_spo_BWT_O(i_p.left() + nTriples, i_p.right() + nTriples);
            //std::cout << "num_distinct_values O = " << num_distinct_values_o << " vs. interval size = " << i_p.size() << std::endl;
            hash_map.insert({triple_pattern->o->varname, num_distinct_values_o});
            return hash_map;
        }
    }
    //Sixth case S ?P O
    else if (!triple_pattern->s->isVariable && triple_pattern->p->isVariable && !triple_pattern->o->isVariable)
    {
        bwt_interval open_interval = graph_spo.open_SOP();
        uint64_t cur_s = graph_spo.min_S(open_interval);
        cur_s = graph_spo.next_S(open_interval, triple_pattern->s->constant);
        if (cur_s == 0 || cur_s != triple_pattern->s->constant)
        {
            hash_map.insert({triple_pattern->p->varname, 0});
            return hash_map;//return 0;
        }
        else
        {
            bwt_interval i_s = graph_spo.down_S(cur_s);
            uint64_t cur_o = graph_spo.min_O_in_S(i_s);
            cur_o = graph_spo.next_O_in_S(i_s, triple_pattern->o->constant);
            if (cur_o == 0 || cur_o != triple_pattern->o->constant)
            {
                hash_map.insert({triple_pattern->p->varname, 0});
                return hash_map;//return 0;
            }
            bwt_interval i_o = graph_spo.down_S_O(i_s, cur_o);
            // Ring => Going from O to P: values must be shifted to the left, by subtracting nTriples. Remember O is between 2 * nTriples + 1 and 3 * nTriples.
            uint64_t num_distinct_values_p = crc_arrays.get_number_distinct_values_spo_BWT_P(i_o.left() - nTriples, i_o.right() - nTriples);
            //std::cout << "num_distinct_values P = " << num_distinct_values_p << " vs. interval size = " << i_o.size() << std::endl;
            hash_map.insert({triple_pattern->p->varname, num_distinct_values_p});
            return hash_map;
        }
    }
    //Seventh case ?S P O
    else if (triple_pattern->s->isVariable && !triple_pattern->p->isVariable && !triple_pattern->o->isVariable)
    {
        bwt_interval open_interval = graph_spo.open_POS();
        uint64_t cur_p = graph_spo.min_P(open_interval);
        cur_p = graph_spo.next_P(open_interval, triple_pattern->p->constant);
        if (cur_p == 0 || cur_p != triple_pattern->p->constant)
        {
            hash_map.insert({triple_pattern->s->varname, 0});
            return hash_map;//return 0;
        }
        else
        {
            bwt_interval i_p = graph_spo.down_P(cur_p);
            uint64_t cur_o = graph_spo.min_O_in_P(i_p, cur_p);
            cur_o = graph_spo.next_O_in_P(i_p, cur_p, triple_pattern->o->constant);
            if (cur_o == 0 || cur_o != triple_pattern->o->constant)
            {
                hash_map.insert({triple_pattern->s->varname, 0});
                return hash_map;//return 0;
            }
            i_p = graph_spo.down_P_O(i_p, cur_p, cur_o);
            // Ring => Going from P to S: values must be shifted to the left, by subtracting nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
            uint64_t num_distinct_values_s = crc_arrays.get_number_distinct_values_spo_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples);
            //std::cout << "num_distinct_values S = " << num_distinct_values_s << " vs. interval size = " << i_p.size() << std::endl;
            hash_map.insert({triple_pattern->s->varname, num_distinct_values_s});
            return hash_map;
        }
    }
    return hash_map;//return 0;
}
bool compare_by_second(pair<string, int> a, pair<string, int> b)
{
    return a.second < b.second;
}

//! TODO:
/*!
 * \param vector<Triple*> : The Basic Graph Pattern (query)
 * \param ring_spo&  : The Ring succinct structure
 * \param ring_sop&  : The reverse Ring.
 * \returns Vector of strings representing the global attribute order.
 */
vector<string> get_gao(vector<Triple *> query, ring_spo &graph_spo, crc_arrays &crc_arrays)
{
    map<string, vector<uint64_t>> triple_values;
    map<string, vector<Triple *>> triples_var;
    for (Triple *triple_pattern : query)
    {
        // OLD: getting the triple_size and storing it into each variable of the triple pattern. Which means every variable has the same weight.
        //uint64_t triple_size = get_size_interval(triple_pattern, graph_spo);
        // OLD: Storing the triple_size of each variable for all triple_patterns.
        auto hash_map = get_num_diff_values(triple_pattern, graph_spo, crc_arrays);
        if (triple_pattern->s->isVariable)
        {
            auto var_aux = hash_map.find(triple_pattern->s->varname);
            if(var_aux != hash_map.end()){
                auto var_weight = var_aux->second;
                triple_values[triple_pattern->s->varname].push_back(var_weight);
                triples_var[triple_pattern->s->varname].push_back(triple_pattern);
            }

        }
        if (triple_pattern->p->isVariable)
        {
            auto var_aux = hash_map.find(triple_pattern->p->varname);
            if(var_aux != hash_map.end()){
                auto var_weight = var_aux->second;
                triple_values[triple_pattern->p->varname].push_back(var_weight);
                triples_var[triple_pattern->p->varname].push_back(triple_pattern);
            }
        }
        if (triple_pattern->o->isVariable)
        {
            auto var_aux = hash_map.find(triple_pattern->o->varname);
            if(var_aux != hash_map.end()){
                auto var_weight = var_aux->second;
                triple_values[triple_pattern->o->varname].push_back(var_weight);
                triples_var[triple_pattern->o->varname].push_back(triple_pattern);
            }
        }
    }

    vector<string> gao;
    vector<pair<string, uint64_t>> varmin_pairs;
    vector<string> single_vars;
    map<string, bool> selectable_vars;
    map<string, bool> selected_vars;
    map<string, set<string>> related_vars;

    for (auto it = (triples_var).cbegin(); it != (triples_var).cend(); ++it)
    {
        for (Triple *triple_pattern : triples_var[it->first])
        {
            if (triple_pattern->s->isVariable && it->first.compare(triple_pattern->s->varname) != 0)
            {
                related_vars[it->first].insert(triple_pattern->s->varname);
            }
            if (triple_pattern->p->isVariable && it->first.compare(triple_pattern->p->varname) != 0)
            {
                related_vars[it->first].insert(triple_pattern->p->varname);
            }
            if (triple_pattern->o->isVariable && it->first.compare(triple_pattern->o->varname) != 0)
            {
                related_vars[it->first].insert(triple_pattern->o->varname);
            }
        }
    }

    for (auto it = (triple_values).cbegin(); it != (triple_values).cend(); ++it)
    {
        if (triple_values[it->first].size() == 1)
        {
            single_vars.push_back(it->first);
        }
        else
        {
            varmin_pairs.push_back(pair<string, uint64_t>(it->first, *min_element(triple_values[it->first].begin(), triple_values[it->first].end())));
            selectable_vars[it->first] = false;
        }
    }

    sort((varmin_pairs).begin(), (varmin_pairs).end(), compare_by_second);
    if (varmin_pairs.size() > 0)
    {
        selectable_vars[varmin_pairs[0].first] = true;
    }
    for (int i = 0; i < varmin_pairs.size(); i++)
    {
        for (pair<string, uint64_t> varmin_pair : varmin_pairs)
        {
            if (selectable_vars[varmin_pair.first] && !selected_vars[varmin_pair.first])
            {
                gao.push_back(varmin_pair.first);
                selected_vars[varmin_pair.first] = true;
                for (set<string>::iterator it = related_vars[varmin_pair.first].begin(); it != related_vars[varmin_pair.first].end(); ++it)
                {
                    selectable_vars[*it] = true;
                }
                break;
            }
        }
    }

    for (pair<string, uint64_t> varmin_pair : varmin_pairs)
    {
        if (!selected_vars[varmin_pair.first])
        {
            gao.push_back(varmin_pair.first);
        }
    }

    for (string s : single_vars)
    {
        gao.push_back(s);
    }

    return gao;
}

#endif
