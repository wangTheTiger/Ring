#ifndef LEAPFROG_OP_H
#define LEAPFROG_OP_H

#include <iostream>
#include <map>
#include <chrono>
#include "Triple.h"
#include "ring_spo.hpp"
#include "Iterator.hpp"
#include <unordered_map>
#include <algorithm>
#include "crc_arrays.hpp"

class LeapfrogOP {
public:
    vector<Triple*>* query;
    ring_spo* graph;
    map<string, vector<Iterator*>> query_iterators;
    vector<Iterator*> all_iterators;
    map<string, vector<Triple *>> triples_var;
    vector<string> single_vars;
    map<string, set<string>> related_vars;
    crc_arrays* crc;
    set<string> processed_vars;
    unordered_map<string,uint64_t> processed_vars_map;
    bool is_empty;
    int number_of_vars;
    map<string, uint64_t>* bindings;
    int* number_of_results;
    //Stores the attribute elimination order during execution.
    std::vector<std::string> gao;
    LeapfrogOP(ring_spo* graph, vector<Triple*>* query, crc_arrays* crc, map<string, uint64_t>* bindings, int * number_of_results) {
        this->graph = graph;
        this->query = query;
        this->is_empty = false;
        this->crc = crc;
        this->number_of_vars = 0;
        this->bindings = bindings;
        this->number_of_results = number_of_results;
        //Calculate the initial attribute (var) to eliminate by LTJ.
        gao.push_back(get_next_eliminated_variable());
        /*
        for (Triple* triple_pattern : *query) {
            Iterator* triple_iterator = new Iterator(triple_pattern, graph);
            if (triple_iterator->is_empty) {
                this->is_empty = true;
            }
            for (string varname : (*this->gao)) {
                if (triple_pattern->contains_variable(varname)) {
                    this->query_iterators[varname].push_back(triple_iterator);
                }
            }
            all_iterators.push_back(triple_iterator);
        }*/
    }

    ~LeapfrogOP() {
        for (Iterator* triple_iterator : all_iterators) {
            delete triple_iterator;
        }
    }

    void print_query() {
        cout << "Query: " << endl;
        for (Triple* triple : *this->query) {
            triple->serialize_as_triple_pattern();
        }
    }

    void print_gao() {
        cout << "GAO: " << endl;
        for (string varname : gao) {
            cout << varname << " ";
        }
        cout << endl;
    }

    void serialize() {
        cout << "QUERY ITERATORS: " << endl;
        for (string varname : gao) {
            cout << "VAR:" << varname << endl;
            for (Iterator* triple_iterator : this->query_iterators[varname]) {
                cout << "Triple pattern: ";
                triple_iterator->triple->serialize_as_triple_pattern();
                cout << "Index name: " << triple_iterator->index_name << endl;
            }
        }
    }
    //! TODO: Used to calculate the starting variable of the adaptive gao approach.
    /*!
    * \returns TODO:
    */
    std::unordered_map<string, uint64_t> get_num_diff_values(Triple *triple_pattern)
    {
        std::unordered_map<string, uint64_t> hash_map;
        const uint64_t nTriples = graph->get_n_triples();
        //First case: ?S ?P ?O
        if (triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable)
        {
            bwt_interval open_interval = graph->open_SPO();
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
            //we need to do graph->get_number_distinct_values_BWT_S() and for ?O we need the inverse ring, that is graph_sop.get_number_distinct_values_BWT_O().
            bwt_interval open_interval = graph->open_PSO();
            //First we get the minimum value of the range
            uint64_t cur_p = graph->min_P(open_interval);
            //Then the next value >= triple_pattern->p->constant.
            cur_p = graph->next_P(open_interval, triple_pattern->p->constant);
            if (cur_p == 0 || cur_p != triple_pattern->p->constant)
            {
                hash_map.insert({triple_pattern->s->varname, 0});
                hash_map.insert({triple_pattern->o->varname, 0});
                return hash_map;//return 0
            }
            else
            {
                bwt_interval i_p = graph->down_P(cur_p); // Range in C array pointing to the S value.
                // Ring => Going from P to S: values must be shifted back to the left, by subtracting nTriples. Remember P is between nTriples + 1 and 2*nTriples.
                uint64_t num_distinct_values_s = crc->get_number_distinct_values_spo_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples);
                // Reverse Ring => Going from P to O: values must be shifted back to the left, by subtracting nTriples. Remember P is between 2 * nTriples + 1 and 3*nTriples.
                // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
                uint64_t num_distinct_values_o = crc->get_number_distinct_values_sop_BWT_O(i_p.left() - nTriples, i_p.right() - nTriples);
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
            bwt_interval open_interval = graph->open_OSP();
            uint64_t cur_o = graph->min_O(open_interval);
            cur_o = graph->next_O(open_interval, triple_pattern->o->constant);
            if (cur_o == 0 || cur_o != triple_pattern->o->constant)
            {
                hash_map.insert({triple_pattern->p->varname, 0});
                hash_map.insert({triple_pattern->s->varname, 0});
                return hash_map;//return 0;
            }
            else
            {
                bwt_interval i_o = graph->down_O(cur_o); // Range in C array pointing to the O value.
                // Ring => Going from O to P: values must be shifted back to the left, by subtracting nTriples. Remember O is between 2* nTriples + 1 and 3*nTriples.
                uint64_t num_distinct_values_p = crc->get_number_distinct_values_spo_BWT_P(i_o.left() - nTriples, i_o.right() - nTriples);
                // Reverse Ring => Going from O to S: values must be shifted back to the left, by subtracting nTriples. Remember O is between nTriples + 1 and 2*nTriples.
                // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
                uint64_t num_distinct_values_s = crc->get_number_distinct_values_sop_BWT_S(i_o.left() - nTriples, i_o.right() - nTriples);
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
            bwt_interval open_interval = graph->open_SPO();
            uint64_t cur_s = graph->min_S(open_interval);
            cur_s = graph->next_S(open_interval, triple_pattern->s->constant);
            if (cur_s == 0 || cur_s != triple_pattern->s->constant)
            {
                hash_map.insert({triple_pattern->o->varname, 0});
                hash_map.insert({triple_pattern->p->varname, 0});
                return hash_map;//return 0;
            }
            else
            {
                bwt_interval i_s = graph->down_S(cur_s);// Range in C array pointing to the S value.
                // Ring => Going from S to O: values must be shifted to the right, by adding 2 * nTriples. Remember S is between  1 and nTriples.
                uint64_t num_distinct_values_o = crc->get_number_distinct_values_spo_BWT_O(i_s.left() + 2 * nTriples, i_s.right() + 2 * nTriples);
                // Reverse Ring => Going from S to P:  values must be shifted to the right, by adding 2 * nTriples. Remember S is between  1 and nTriples.
                // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
                uint64_t num_distinct_values_p = crc->get_number_distinct_values_sop_BWT_P(i_s.left() + 2 * nTriples, i_s.right() + 2 * nTriples);
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
            bwt_interval open_interval = graph->open_SPO();
            uint64_t cur_s = graph->min_S(open_interval);
            cur_s = graph->next_S(open_interval, triple_pattern->s->constant);
            if (cur_s == 0 || cur_s != triple_pattern->s->constant)
            {
                hash_map.insert({triple_pattern->o->varname, 0});
                return hash_map;//return 0;
            }
            else
            {
                bwt_interval i_s = graph->down_S(cur_s);
                uint64_t cur_p = graph->min_P_in_S(i_s, cur_s);
                cur_p = graph->next_P_in_S(i_s, cur_s, triple_pattern->p->constant);
                if (cur_p == 0 || cur_p != triple_pattern->p->constant)
                {
                    hash_map.insert({triple_pattern->o->varname, 0});
                    return hash_map;//return 0;
                }
                bwt_interval i_p = graph->down_S_P(i_s, cur_s, cur_p);
                // Ring => Going from P to O: values must be shifted to the right, by adding nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
                uint64_t num_distinct_values_o = crc->get_number_distinct_values_spo_BWT_O(i_p.left() + nTriples, i_p.right() + nTriples);
                //std::cout << "num_distinct_values O = " << num_distinct_values_o << " vs. interval size = " << i_p.size() << std::endl;
                hash_map.insert({triple_pattern->o->varname, num_distinct_values_o});
                return hash_map;
            }
        }
        //Sixth case S ?P O
        else if (!triple_pattern->s->isVariable && triple_pattern->p->isVariable && !triple_pattern->o->isVariable)
        {
            bwt_interval open_interval = graph->open_SOP();
            uint64_t cur_s = graph->min_S(open_interval);
            cur_s = graph->next_S(open_interval, triple_pattern->s->constant);
            if (cur_s == 0 || cur_s != triple_pattern->s->constant)
            {
                hash_map.insert({triple_pattern->p->varname, 0});
                return hash_map;//return 0;
            }
            else
            {
                bwt_interval i_s = graph->down_S(cur_s);
                uint64_t cur_o = graph->min_O_in_S(i_s);
                cur_o = graph->next_O_in_S(i_s, triple_pattern->o->constant);
                if (cur_o == 0 || cur_o != triple_pattern->o->constant)
                {
                    hash_map.insert({triple_pattern->p->varname, 0});
                    return hash_map;//return 0;
                }
                bwt_interval i_o = graph->down_S_O(i_s, cur_o);
                // Ring => Going from O to P: values must be shifted to the left, by subtracting nTriples. Remember O is between 2 * nTriples + 1 and 3 * nTriples.
                uint64_t num_distinct_values_p = crc->get_number_distinct_values_spo_BWT_P(i_o.left() - nTriples, i_o.right() - nTriples);
                //std::cout << "num_distinct_values P = " << num_distinct_values_p << " vs. interval size = " << i_o.size() << std::endl;
                hash_map.insert({triple_pattern->p->varname, num_distinct_values_p});
                return hash_map;
            }
        }
        //Seventh case ?S P O
        else if (triple_pattern->s->isVariable && !triple_pattern->p->isVariable && !triple_pattern->o->isVariable)
        {
            bwt_interval open_interval = graph->open_POS();
            uint64_t cur_p = graph->min_P(open_interval);
            cur_p = graph->next_P(open_interval, triple_pattern->p->constant);
            if (cur_p == 0 || cur_p != triple_pattern->p->constant)
            {
                hash_map.insert({triple_pattern->s->varname, 0});
                return hash_map;//return 0;
            }
            else
            {
                bwt_interval i_p = graph->down_P(cur_p);
                uint64_t cur_o = graph->min_O_in_P(i_p, cur_p);
                cur_o = graph->next_O_in_P(i_p, cur_p, triple_pattern->o->constant);
                if (cur_o == 0 || cur_o != triple_pattern->o->constant)
                {
                    hash_map.insert({triple_pattern->s->varname, 0});
                    return hash_map;//return 0;
                }
                i_p = graph->down_P_O(i_p, cur_p, cur_o);
                // Ring => Going from P to S: values must be shifted to the left, by subtracting nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
                uint64_t num_distinct_values_s = crc->get_number_distinct_values_spo_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples);
                //std::cout << "num_distinct_values S = " << num_distinct_values_s << " vs. interval size = " << i_p.size() << std::endl;
                hash_map.insert({triple_pattern->s->varname, num_distinct_values_s});
                return hash_map;
            }
        }
        return hash_map;//return 0;
    }
    //! TODO:
    /*!
    * \returns TODO:
    */
    uint64_t get_num_diff_values(string cur_var, string last_processed_var, uint64_t last_processed_val, Triple *triple_pattern){
        
        const uint64_t nTriples = graph->get_n_triples();
        bool is_second_case = false, is_third_case = false, is_fourth_case = false;
        //First case: ?S ?P ?O
        if (triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable)
        {
            //Then checking what variable is currently being processed.
            if(triple_pattern->s->varname == cur_var){
                //It is really the fourth case : S ?P ?O
                is_fourth_case = true;
            } else if(triple_pattern->p->varname == cur_var){
                //It is really the second case : ?S P ?O
                is_second_case = true;
            } else if(triple_pattern->o->varname == cur_var){
                //It is really the third case : ?S ?P O
                is_third_case = true;
            }
        }
        //Second case : ?S P ?O.
        //If processing_var's type is "?S" and has a 'last_processed_var', then is eq to Fifth case:   S P ?O.
        //Or processing_var's type is "?O" and has a 'last_processed_var', then is eq to Seventh case: ?S P O.
        //otherwise both ?S and ?O are unbounded.
        else if (is_second_case || (triple_pattern->s->isVariable && !triple_pattern->p->isVariable && triple_pattern->o->isVariable))
        {
            //TODO: falta pensar cuando venimos desde el first case. quizas no fue una buena idea usar las booleanas is_second_case, is_... por la propagación del cachito.
            if(triple_pattern->s->varname == last_processed_var){
                //Eq to fifth case: S P ?O. whereas P is constant and S is variable with an existing binding.
                //std::cout << "DEBUG -- second case with bounded S var" << std::endl;
                bwt_interval open_interval = graph->open_SPO();
                uint64_t cur_s = graph->min_S(open_interval);
                cur_s = graph->next_S(open_interval, last_processed_val);
                if (cur_s == 0 || cur_s != last_processed_val)
                {
                    return 0;
                }
                else
                {
                    bwt_interval i_s = graph->down_S(cur_s);
                    uint64_t cur_p = graph->min_P_in_S(i_s, cur_s);
                    cur_p = graph->next_P_in_S(i_s, cur_s, triple_pattern->p->constant);
                    if (cur_p == 0 || cur_p != triple_pattern->p->constant)
                    {
                        return 0;
                    }
                    bwt_interval i_p = graph->down_S_P(i_s, cur_s, cur_p);//TODO: opt if i_p.size() <= SOMETHING then just return i_p.size() :-)
                    if(triple_pattern->o->varname == cur_var){
                        // Ring => Going from P to O: values must be shifted to the right, by adding nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
                        return crc->get_number_distinct_values_spo_BWT_O(i_p.left() + nTriples, i_p.right() + nTriples);
                    }
                    return 0;
                }
            } else if(triple_pattern->o->varname == last_processed_var){
                //Eq to seventh case: ?S P O. whereas P is constant and O is variable with an existing binding.
                //std::cout << "DEBUG -- second case with bounded O var" << std::endl;
                bwt_interval open_interval = graph->open_POS();
                uint64_t cur_p = graph->min_P(open_interval);
                cur_p = graph->next_P(open_interval, triple_pattern->p->constant);
                if (cur_p == 0 || cur_p != triple_pattern->p->constant)
                {
                    return 0;
                }
                else
                {
                    bwt_interval i_p = graph->down_P(cur_p);
                    uint64_t cur_o = graph->min_O_in_P(i_p, cur_p);
                    cur_o = graph->next_O_in_P(i_p, cur_p, last_processed_val);
                    if (cur_o == 0 || cur_o != last_processed_val)
                    {
                        return 0;
                    }
                    i_p = graph->down_P_O(i_p, cur_p, cur_o);//TODO: opt if i_p.size() <= SOMETHING then just return i_p.size() :-)
                    if(triple_pattern->s->varname == cur_var){
                        // Ring => Going from P to S: values must be shifted to the left, by subtracting nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
                    return crc->get_number_distinct_values_spo_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples);
                    }
                    return 0;
                }
            }else{
                //When only P is a constant, and both ?S and ?O are unbounded variables.
                bwt_interval open_interval = graph->open_PSO();
                //First we get the minimum value of the range (for fixed P)
                uint64_t cur_p = graph->min_P(open_interval);
                //Then the next value >= triple_pattern->p->constant.
                cur_p = graph->next_P(open_interval, triple_pattern->p->constant);
                if (cur_p == 0 || cur_p != triple_pattern->p->constant)
                {
                    return 0;
                }
                else
                {
                    bwt_interval i_p = graph->down_P(cur_p); // Range in C array pointing to the S value.//TODO: opt if i_p.size() <= SOMETHING then just return i_p.size() :-)
                    if(triple_pattern->s->varname == cur_var){
                        // Ring => Going from P to S: values must be shifted back to the left, by subtracting nTriples. Remember P is between nTriples + 1 and 2*nTriples.
                        return crc->get_number_distinct_values_spo_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples);
                    } else if(triple_pattern->o->varname == cur_var){
                        // Reverse Ring => Going from P to O: values must be shifted back to the left, by subtracting nTriples. Remember P is between 2 * nTriples + 1 and 3*nTriples.
                        // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
                        return crc->get_number_distinct_values_sop_BWT_O(i_p.left() - nTriples, i_p.right() - nTriples);
                    }
                    return 0;
                }
            }
        }
        //Third case ?S ?P O
        else if (is_third_case || (triple_pattern->s->isVariable && triple_pattern->p->isVariable && !triple_pattern->o->isVariable))
        {
            bwt_interval open_interval = graph->open_OSP();
            uint64_t cur_o = graph->min_O(open_interval);
            cur_o = graph->next_O(open_interval, triple_pattern->o->constant);
            if (cur_o == 0 || cur_o != triple_pattern->o->constant)
            {
                return 0;
            }
            else
            {
                bwt_interval i_o = graph->down_O(cur_o); // Range in C array pointing to the O value.
                if(triple_pattern->p->varname == cur_var){
                    // Ring => Going from O to P: values must be shifted back to the left, by subtracting nTriples. Remember O is between 2* nTriples + 1 and 3*nTriples.
                    return crc->get_number_distinct_values_spo_BWT_P(i_o.left() - nTriples, i_o.right() - nTriples);
                } else if(triple_pattern->s->varname == cur_var){
                    // Reverse Ring => Going from O to S: values must be shifted back to the left, by subtracting nTriples. Remember O is between nTriples + 1 and 2*nTriples.
                    // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
                    return crc->get_number_distinct_values_sop_BWT_S(i_o.left() - nTriples, i_o.right() - nTriples);
                }
                return 0;
            }
        }
        //Fourth case S ?P ?O
        else if (is_fourth_case || (!triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable))
        {
            bwt_interval open_interval = graph->open_SPO();
            uint64_t cur_s = graph->min_S(open_interval);
            cur_s = graph->next_S(open_interval, triple_pattern->s->constant);
            if (cur_s == 0 || cur_s != triple_pattern->s->constant)
            {
                return 0;
            }
            else
            {
                bwt_interval i_s = graph->down_S(cur_s);// Range in C array pointing to the S value.
                if(triple_pattern->o->varname == cur_var){
                    // Ring => Going from S to O: values must be shifted to the right, by adding 2 * nTriples. Remember S is between  1 and nTriples.
                    return crc->get_number_distinct_values_spo_BWT_O(i_s.left() + 2 * nTriples, i_s.right() + 2 * nTriples);
                } else if(triple_pattern->p->varname == cur_var){
                    // Reverse Ring => Going from S to P:  values must be shifted to the right, by adding 2 * nTriples. Remember S is between  1 and nTriples.
                    // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
                    return crc->get_number_distinct_values_sop_BWT_P(i_s.left() + 2 * nTriples, i_s.right() + 2 * nTriples);
                }
                return 0;
            }
        }
        //Fifth case S P ?O
        else if (!triple_pattern->s->isVariable && !triple_pattern->p->isVariable && triple_pattern->o->isVariable)
        {
            bwt_interval open_interval = graph->open_SPO();
            uint64_t cur_s = graph->min_S(open_interval);
            cur_s = graph->next_S(open_interval, triple_pattern->s->constant);
            if (cur_s == 0 || cur_s != triple_pattern->s->constant)
            {
                return 0;
            }
            else
            {
                bwt_interval i_s = graph->down_S(cur_s);
                uint64_t cur_p = graph->min_P_in_S(i_s, cur_s);
                cur_p = graph->next_P_in_S(i_s, cur_s, triple_pattern->p->constant);
                if (cur_p == 0 || cur_p != triple_pattern->p->constant)
                {
                    return 0;
                }
                bwt_interval i_p = graph->down_S_P(i_s, cur_s, cur_p);
                if(triple_pattern->o->varname == cur_var){
                    // Ring => Going from P to O: values must be shifted to the right, by adding nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
                    return crc->get_number_distinct_values_spo_BWT_O(i_p.left() + nTriples, i_p.right() + nTriples);
                }
                return 0;
            }
        }
        //Sixth case S ?P O
        else if (!triple_pattern->s->isVariable && triple_pattern->p->isVariable && !triple_pattern->o->isVariable)
        {
            bwt_interval open_interval = graph->open_SOP();
            uint64_t cur_s = graph->min_S(open_interval);
            cur_s = graph->next_S(open_interval, triple_pattern->s->constant);
            if (cur_s == 0 || cur_s != triple_pattern->s->constant)
            {
                return 0;
            }
            else
            {
                bwt_interval i_s = graph->down_S(cur_s);
                uint64_t cur_o = graph->min_O_in_S(i_s);
                cur_o = graph->next_O_in_S(i_s, triple_pattern->o->constant);
                if (cur_o == 0 || cur_o != triple_pattern->o->constant)
                {
                    return 0;
                }
                bwt_interval i_o = graph->down_S_O(i_s, cur_o);
                if(triple_pattern->p->varname == cur_var){
                    // Ring => Going from O to P: values must be shifted to the left, by subtracting nTriples. Remember O is between 2 * nTriples + 1 and 3 * nTriples.
                    return crc->get_number_distinct_values_spo_BWT_P(i_o.left() - nTriples, i_o.right() - nTriples);
                }
                return 0;
            }
        }
        //Seventh case ?S P O
        else if (triple_pattern->s->isVariable && !triple_pattern->p->isVariable && !triple_pattern->o->isVariable)
        {
            bwt_interval open_interval = graph->open_POS();
            uint64_t cur_p = graph->min_P(open_interval);
            cur_p = graph->next_P(open_interval, triple_pattern->p->constant);
            if (cur_p == 0 || cur_p != triple_pattern->p->constant)
            {
                return 0;
            }
            else
            {
                bwt_interval i_p = graph->down_P(cur_p);
                uint64_t cur_o = graph->min_O_in_P(i_p, cur_p);
                cur_o = graph->next_O_in_P(i_p, cur_p, triple_pattern->o->constant);
                if (cur_o == 0 || cur_o != triple_pattern->o->constant)
                {
                    return 0;
                }
                i_p = graph->down_P_O(i_p, cur_p, cur_o);
                if(triple_pattern->s->varname == cur_var){
                    // Ring => Going from P to S: values must be shifted to the left, by subtracting nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
                   return crc->get_number_distinct_values_spo_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples);
                }
                return 0;
            }
        }
        return 0;
    }
    void evaluate(int level, std::chrono::steady_clock::time_point begin, std::string last_processed_var = "") {
        assert(number_of_vars > 0);
        if (this->is_empty) {
            return;
        }
        auto varname = gao[0];
        if(level > 0){
            varname = get_next_eliminated_variable();
        }
        std::cout << "DEBUG -- next variable : " << varname << std::endl;
        vector<Iterator*>* var_iterators = &this->query_iterators[varname];

        if ((*var_iterators).size() == 1 && ((*var_iterators)[0]->current_level == 1)) {
            Iterator* triple_iterator = (*var_iterators)[0];
            vector<pair<uint64_t, uint64_t>> iterator_last = triple_iterator->down_last();
            // cout << "SIZE: " << varname << " " << iterator_last.size() << endl;
            for (pair<uint64_t, uint64_t> binding_last : iterator_last) {
                // Limit
                if (*number_of_results >= 1000) {
                    break;
                }
                // TO
                std::chrono::steady_clock::time_point to_time = std::chrono::steady_clock::now();
                if (std::chrono::duration_cast<std::chrono::seconds> (to_time - begin).count() >= 600) {
                    break;
                }
                if (level >= (number_of_vars - 1)) {
                    // Print Answers
                    (*bindings)[varname] = binding_last.second;

                    for(auto it = (*bindings).cbegin(); it != (*bindings).cend(); ++it) {
                        cout << "(" << it->first << ": " << (*bindings)[it->first] << ") ";
                    }
                    cout << endl;

                    (*number_of_results) = (*number_of_results) + 1;

                } else {
                    (*bindings)[varname] = binding_last.second;
                    processed_vars.insert(varname);
                    processed_vars_map[varname] = (*bindings)[varname];
                    int new_level = level + 1;
                    evaluate(new_level, begin, varname);
                }
            }
            for (Iterator* triple_iterator : *var_iterators) {
                triple_iterator->up();
            }
        }
        else {
            for (Iterator* triple_iterator : *var_iterators) {
                triple_iterator->down();
            }
            bool search = true;
            while (search) {
                // Limit
                if (*number_of_results >= 1000) {
                    for (Iterator* triple_iterator : *var_iterators) {
                        triple_iterator->up();
                    }
                    break;
                }
                // TO
                std::chrono::steady_clock::time_point to_time = std::chrono::steady_clock::now();
                if (std::chrono::duration_cast<std::chrono::seconds> (to_time - begin).count() >= 600) {
                    for (Iterator* triple_iterator : *var_iterators) {
                        triple_iterator->up();
                    }
                    break;
                }
                sort((*var_iterators).rbegin(), (*var_iterators).rend(), compare_by_current_value);
                if ((*var_iterators)[0]->current_value() == (*var_iterators)[var_iterators->size() - 1]->current_value()) {
                    if (level >= (number_of_vars - 1)) {
                        // Print Answers
                        (*bindings)[varname] = (*var_iterators)[0]->current_value();

                        for(auto it = (*bindings).cbegin(); it != (*bindings).cend(); ++it) {
                            cout << "(" << it->first << ": " << (*bindings)[it->first] << ") ";
                        }
                        cout << endl;

                        (*number_of_results) = (*number_of_results) + 1;
                        int next_value = (*var_iterators)[0]->current_value() + 1;
                        (*var_iterators)[0]->seek(next_value);
                    } else {
                        (*bindings)[varname] = (*var_iterators)[0]->current_value();
                        processed_vars.insert(varname);
                        processed_vars_map[varname] = (*bindings)[varname];
                        int new_level = level + 1;
                        evaluate(new_level, begin, varname);
                        int next_value = (*var_iterators)[0]->current_value() + 1;
                        (*var_iterators)[0]->seek(next_value);
                    }
                }

                int seek_value = (*var_iterators)[0]->current_value();
                if (seek_value != 0) {
                    for (int i = 1; i < var_iterators->size(); i++) {
                        (*var_iterators)[i]->seek(seek_value);
                    }
                }

                for (Iterator* triple_iterator : (*var_iterators)) {
                    if (triple_iterator->current_value() == 0) {
                        for (Iterator* triple_iterator : *var_iterators) {
                            triple_iterator->up();
                        }
                        search = false;
                    }
                }

            }
        }

    }

    static bool compare_by_current_value(Iterator* a, Iterator* b) {
        return a->current_value() < b->current_value();
    }


    //! TODO: Adaptative calculation of Global Attribute Elimination Order (gao). Similar to 'get_gao' (utils.hpp) function.
    /*!
    * \returns std::string representing the next attribute to eliminate by LTJ algorithm.
    */
    std::string get_next_eliminated_variable()
    {
        //First variable.
        if(gao.size() == 0){
            map<string, vector<uint64_t>> triple_values;
            for (Triple *triple_pattern : *this->query)
            {
                auto hash_map = get_num_diff_values(triple_pattern);
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
            number_of_vars = triples_var.size();

            vector<pair<string, uint64_t>> varmin_pairs;

            //Related vars
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
            //Single vars
            for (auto it = (triple_values).cbegin(); it != (triple_values).cend(); ++it)
            {
                if (triple_values[it->first].size() == 1)
                {
                    single_vars.push_back(it->first);
                }
                else
                {
                    varmin_pairs.push_back(pair<string, uint64_t>(it->first, *min_element(triple_values[it->first].begin(), triple_values[it->first].end())));
                }
            }

            sort((varmin_pairs).begin(), (varmin_pairs).end(), [&](pair<string, int> a, pair<string, int> b) -> bool { return a.second < b.second; });
            //return the first variable that is not marked as a single variable.
            for(auto& candidate_var : varmin_pairs)
            {
                if(std::find((single_vars).begin(), (single_vars).end(), candidate_var.first) != (single_vars).end())
                {
                    return candidate_var.first;
                }
            }

            //Otherwise return the first sigle variable.
            for(auto& candidate_var : varmin_pairs)
            {
                if(std::find((single_vars).begin(), (single_vars).end(), candidate_var.first) == (single_vars).end())
                {
                    return candidate_var.first;
                }
            }

        }
        else
        {
            //Any other variable but the first. Here we have more information about the BGP.
            //Vector of pairs: variable and its minimum value.
            std::vector<std::pair<std::string, uint64_t>> varmin_pairs;
            //last processed variable and its value.
            auto last_processed_var = gao[gao.size() - 1];
            assert ( last_processed_var != "");
            auto& last_processed_value = processed_vars_map[last_processed_var];
            std::cout << "DEBUG -- evaluate level : " << gao.size() << " last bound var : " << last_processed_var << " last bound value : " << last_processed_value << std::endl;
            //1. Related vars
            auto& rel_vars = (this->related_vars)[last_processed_var];
            for(auto candidate_var : rel_vars){
                if(std::find((this->single_vars).begin(), (this->single_vars).end(), candidate_var) == (this->single_vars).end() && //(1)
                    std::find(processed_vars.begin(), processed_vars.end(), candidate_var) == processed_vars.end() ) // (2)
                {
                    //(1) variable 'candidate_var' is not marked as a single variable.
                    //(2) variable 'candidate_var' has not been processed.
                    //Until here we just have a set of variables that can be the next one to be processed.
                    //We need the triples related to this variable.
                    auto& triples = (this->triples_var)[candidate_var]; //vector<Triple *>
                    //Then I can pass those triples to: get_num_diff_values
                    //This function has to consider that our variable is only once for each triple. therefore only 1 CRC calculation has to be performed (per triple).
                    uint64_t min_num_diff_vals = std::numeric_limits<uint64_t>::max();
                    for (Triple *triple_pattern : triples)
                    {
                        uint64_t aux = get_num_diff_values(candidate_var, last_processed_var, last_processed_value, triple_pattern);
                        min_num_diff_vals = aux < min_num_diff_vals ? aux : min_num_diff_vals;
                        //std::cout << "DEBUG -- candidate_var: " << candidate_var << " last bound var : " << last_processed_var << " last bond value: "<< last_processed_value << " num of distinct values : " << aux << " for triple: ";
                        //triple_pattern->serialize_as_triple_pattern();
                    }
                    //std::cout << "DEBUG -- minimum value for candidate_var " << candidate_var << " is : " << min_num_diff_vals << std::endl;
                    if((*bindings)["?x2"] == 10216663){
                        std::cout << "llegue" <<std::endl;
                    }
                    if((*bindings)["?x2"] == 10216663 && (*bindings)["?x1"] == 16797400){//  && (*bindings)["?x4"] == 14445472
                        std::cout << "DEBUG -- candidate_var: " << candidate_var << " last bound var : " << last_processed_var << " last bound value: "<< last_processed_value << " minimum num of distinct values : " << min_num_diff_vals << std::endl;
                    }
                    //assert(min_num_diff_vals > 0);
                    varmin_pairs.push_back(std::pair<std::string, uint64_t>(candidate_var,min_num_diff_vals));//TODO: quedarme con el m'inimo y tambi'en recordar que debo setear el score.
                }
            }
            //Getting the related variable with minimum number of different values.
            if(varmin_pairs.size() > 0)
            {
                sort((varmin_pairs).begin(), (varmin_pairs).end(), [&](pair<string, int> a, pair<string, int> b) -> bool { return a.second < b.second; });
                return varmin_pairs[0].first;
            }else
            {
                //2. Single vars
                for(auto &candidate_var : this->single_vars)
                {
                    if(std::find(processed_vars.begin(), processed_vars.end(), candidate_var) == processed_vars.end() ) // (2)
                    {
                        return candidate_var;
                    }
                }
            }
        }
    }
};


#endif //LEAPFROG_ITERATOR_H
