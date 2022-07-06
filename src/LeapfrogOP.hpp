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
    vector<string>* gao;
    vector<Triple*>* query;
    ring_spo* graph;
    map<string, vector<Iterator*>> query_iterators;
    vector<Iterator*> all_iterators;
    map<string, vector<Triple *>>* triples_var;
    vector<string>* single_vars;
    map<string, set<string>>* related_vars;
    crc_arrays* crc;
    vector<string> processed_vars;
    unordered_map<string,uint64_t> processed_vars_map;
    bool is_empty;
    
    LeapfrogOP(vector<string>* gao, ring_spo* graph, vector<Triple*>* query,map<string, vector<Triple *>>* triples_var, map<string, set<string>>* related_vars, vector<string>* single_vars, crc_arrays* crc) {
        this->gao = gao;
        this->graph = graph;
        this->query = query;
        this->is_empty = false;
        this->triples_var = triples_var;
        this->related_vars = related_vars;
        this->single_vars = single_vars;
        this->crc = crc;

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
        }
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
        for (string varname : (*this->gao)) {
            cout << varname << " ";
        }
        cout << endl;
    }

    void serialize() {
        cout << "QUERY ITERATORS: " << endl;
        for (string varname : (*this->gao)) {
            cout << "VAR:" << varname << endl;
            for (Iterator* triple_iterator : this->query_iterators[varname]) {
                cout << "Triple pattern: ";
                triple_iterator->triple->serialize_as_triple_pattern();
                cout << "Index name: " << triple_iterator->index_name << endl;
            }
        }
    }

    //! TODO:
    /*!
    * \returns TODO:
    */
    uint64_t get_num_diff_values(string cur_var, string last_processed_var, uint64_t last_processed_val, Triple *triple_pattern, ring_spo &graph_spo, crc_arrays& crc_arrays){
        
        const uint64_t nTriples = graph_spo.get_n_triples();
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
        else if (is_second_case || (triple_pattern->s->isVariable && !triple_pattern->p->isVariable && triple_pattern->o->isVariable))
        {
            //TODO: falta pensar cuando venimos desde el first case. quizas no fue una buena idea usar las booleanas is_second_case, is_... por la propagación del cachito.
            /*
            OLD CODE:
            bwt_interval open_interval = graph_spo.open_PSO();
            //First we get the minimum value of the range (for fixed P)
            uint64_t cur_p = graph_spo.min_P(open_interval);
            //Then the next value >= triple_pattern->p->constant.
            cur_p = graph_spo.next_P(open_interval, triple_pattern->p->constant);
            if (cur_p == 0 || cur_p != triple_pattern->p->constant)
            {
                return 0;
            }
            else
            {
                bwt_interval i_p = graph_spo.down_P(cur_p); // Range in C array pointing to the S value.
                if(triple_pattern->s->varname == processing_var){
                    // Ring => Going from P to S: values must be shifted back to the left, by subtracting nTriples. Remember P is between nTriples + 1 and 2*nTriples.
                    return crc_arrays.get_number_distinct_values_spo_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples);
                } else if(triple_pattern->o->varname == processing_var){
                    // Reverse Ring => Going from P to O: values must be shifted back to the left, by subtracting nTriples. Remember P is between 2 * nTriples + 1 and 3*nTriples.
                    // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
                    return crc_arrays.get_number_distinct_values_sop_BWT_O(i_p.left() - nTriples, i_p.right() - nTriples);
                }
                return 0;
            }*/
            if(triple_pattern->s->varname == last_processed_var){
                //Eq to fifth case: S P ?O. whereas P is constant and S is variable with an existing binding.
                std::cout << " second case with binding S" << std::endl;
                bwt_interval open_interval = graph_spo.open_SPO();
                uint64_t cur_s = graph_spo.min_S(open_interval);
                cur_s = graph_spo.next_S(open_interval, last_processed_val);
                if (cur_s == 0 || cur_s != last_processed_val)
                {
                    return 0;
                }
                else
                {
                    bwt_interval i_s = graph_spo.down_S(cur_s);
                    uint64_t cur_p = graph_spo.min_P_in_S(i_s, cur_s);
                    cur_p = graph_spo.next_P_in_S(i_s, cur_s, triple_pattern->p->constant);
                    if (cur_p == 0 || cur_p != triple_pattern->p->constant)
                    {
                        return 0;
                    }
                    bwt_interval i_p = graph_spo.down_S_P(i_s, cur_s, cur_p);
                    if(triple_pattern->o->varname == cur_var){
                        // Ring => Going from P to O: values must be shifted to the right, by adding nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
                        return crc_arrays.get_number_distinct_values_spo_BWT_O(i_p.left() + nTriples, i_p.right() + nTriples);
                    }
                    return 0;
                }
            } else if(triple_pattern->o->varname == last_processed_var){
                //Eq to seventh case: ?S P O. whereas P is constant and O is variable with an existing binding.
                std::cout << " second case with binded O" << std::endl;
                bwt_interval open_interval = graph_spo.open_POS();
                uint64_t cur_p = graph_spo.min_P(open_interval);
                cur_p = graph_spo.next_P(open_interval, triple_pattern->p->constant);
                if (cur_p == 0 || cur_p != triple_pattern->p->constant)
                {
                    return 0;
                }
                else
                {
                    bwt_interval i_p = graph_spo.down_P(cur_p);
                    uint64_t cur_o = graph_spo.min_O_in_P(i_p, cur_p);
                    cur_o = graph_spo.next_O_in_P(i_p, cur_p, last_processed_val);
                    if (cur_o == 0 || cur_o != last_processed_val)
                    {
                        return 0;
                    }
                    i_p = graph_spo.down_P_O(i_p, cur_p, cur_o);
                    if(triple_pattern->s->varname == cur_var){
                        // Ring => Going from P to S: values must be shifted to the left, by subtracting nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
                    return crc_arrays.get_number_distinct_values_spo_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples);
                    }
                    return 0;
                }
            }
            return 0;
        }
        //Third case ?S ?P O
        else if (is_third_case || (triple_pattern->s->isVariable && triple_pattern->p->isVariable && !triple_pattern->o->isVariable))
        {
            bwt_interval open_interval = graph_spo.open_OSP();
            uint64_t cur_o = graph_spo.min_O(open_interval);
            cur_o = graph_spo.next_O(open_interval, triple_pattern->o->constant);
            if (cur_o == 0 || cur_o != triple_pattern->o->constant)
            {
                return 0;
            }
            else
            {
                bwt_interval i_o = graph_spo.down_O(cur_o); // Range in C array pointing to the O value.
                if(triple_pattern->p->varname == cur_var){
                    // Ring => Going from O to P: values must be shifted back to the left, by subtracting nTriples. Remember O is between 2* nTriples + 1 and 3*nTriples.
                    return crc_arrays.get_number_distinct_values_spo_BWT_P(i_o.left() - nTriples, i_o.right() - nTriples);
                } else if(triple_pattern->s->varname == cur_var){
                    // Reverse Ring => Going from O to S: values must be shifted back to the left, by subtracting nTriples. Remember O is between nTriples + 1 and 2*nTriples.
                    // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
                    return crc_arrays.get_number_distinct_values_sop_BWT_S(i_o.left() - nTriples, i_o.right() - nTriples);
                }
                return 0;
            }
        }
        //Fourth case S ?P ?O
        else if (is_fourth_case || (!triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable))
        {
            bwt_interval open_interval = graph_spo.open_SPO();
            uint64_t cur_s = graph_spo.min_S(open_interval);
            cur_s = graph_spo.next_S(open_interval, triple_pattern->s->constant);
            if (cur_s == 0 || cur_s != triple_pattern->s->constant)
            {
                return 0;
            }
            else
            {
                bwt_interval i_s = graph_spo.down_S(cur_s);// Range in C array pointing to the S value.
                if(triple_pattern->o->varname == cur_var){
                    // Ring => Going from S to O: values must be shifted to the right, by adding 2 * nTriples. Remember S is between  1 and nTriples.
                    return crc_arrays.get_number_distinct_values_spo_BWT_O(i_s.left() + 2 * nTriples, i_s.right() + 2 * nTriples);
                } else if(triple_pattern->p->varname == cur_var){
                    // Reverse Ring => Going from S to P:  values must be shifted to the right, by adding 2 * nTriples. Remember S is between  1 and nTriples.
                    // Important: both ring's C_s and reverse ring's C_o contains range of P's ordered lexicographically, therefore they are equivalents.
                    return crc_arrays.get_number_distinct_values_sop_BWT_P(i_s.left() + 2 * nTriples, i_s.right() + 2 * nTriples);
                }
                return 0;
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
                return 0;
            }
            else
            {
                bwt_interval i_s = graph_spo.down_S(cur_s);
                uint64_t cur_p = graph_spo.min_P_in_S(i_s, cur_s);
                cur_p = graph_spo.next_P_in_S(i_s, cur_s, triple_pattern->p->constant);
                if (cur_p == 0 || cur_p != triple_pattern->p->constant)
                {
                    return 0;
                }
                bwt_interval i_p = graph_spo.down_S_P(i_s, cur_s, cur_p);
                if(triple_pattern->o->varname == cur_var){
                    // Ring => Going from P to O: values must be shifted to the right, by adding nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
                    return crc_arrays.get_number_distinct_values_spo_BWT_O(i_p.left() + nTriples, i_p.right() + nTriples);
                }
                return 0;
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
                return 0;
            }
            else
            {
                bwt_interval i_s = graph_spo.down_S(cur_s);
                uint64_t cur_o = graph_spo.min_O_in_S(i_s);
                cur_o = graph_spo.next_O_in_S(i_s, triple_pattern->o->constant);
                if (cur_o == 0 || cur_o != triple_pattern->o->constant)
                {
                    return 0;
                }
                bwt_interval i_o = graph_spo.down_S_O(i_s, cur_o);
                if(triple_pattern->p->varname == cur_var){
                    // Ring => Going from O to P: values must be shifted to the left, by subtracting nTriples. Remember O is between 2 * nTriples + 1 and 3 * nTriples.
                    return crc_arrays.get_number_distinct_values_spo_BWT_P(i_o.left() - nTriples, i_o.right() - nTriples);
                }
                return 0;
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
                return 0;
            }
            else
            {
                bwt_interval i_p = graph_spo.down_P(cur_p);
                uint64_t cur_o = graph_spo.min_O_in_P(i_p, cur_p);
                cur_o = graph_spo.next_O_in_P(i_p, cur_p, triple_pattern->o->constant);
                if (cur_o == 0 || cur_o != triple_pattern->o->constant)
                {
                    return 0;
                }
                i_p = graph_spo.down_P_O(i_p, cur_p, cur_o);
                if(triple_pattern->s->varname == cur_var){
                    // Ring => Going from P to S: values must be shifted to the left, by subtracting nTriples. Remember P is between nTriples + 1 and 2 * nTriples.
                   return crc_arrays.get_number_distinct_values_spo_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples);
                }
                return 0;
            }
        }
        return 0;
    }
    void evaluate(int level, map<string, uint64_t>* bindings, int * number_of_results, std::chrono::steady_clock::time_point begin) {
        if (this->is_empty) {
            return;
        }
        //Using the pre-calculated variable from gao when the recursion starts. (level == 1)
        string varname = (*this->gao)[level];
        if(level > 0){
            //last processed variable and its value.
            for ( auto var : processed_vars){
                std::cout << "processed_vars : " << var << " ";
            }
            std::cout << "" << std::endl;
            std::string last_processed_var = *(processed_vars.rbegin());//TODO: processed_vars contiene varias veces las mismas variables. Analizar...
            auto& last_processed_value = processed_vars_map[last_processed_var];
            //Related vars TODO: WHAT ABOUT SINGLE_VARS?
            auto& rel_vars = (*this->related_vars)[last_processed_var];
            for(auto candidate_var : rel_vars){
                if(std::find(processed_vars.begin(), processed_vars.end(), candidate_var) == processed_vars.end() )
                {
                    //variable 'var' has not been processed.
                    //Until here we just have a set of variables that can be our next one to be processed. How can we calculate them the muthu?
                    //Possible solution: We need the triples related to this variable.
                    auto& triples = (*this->triples_var)[candidate_var]; //vector<Triple *>
                    //Then I can pass those triples to: get_num_diff_values -> requires three parameters, can I get them?
                    //auto hash_map = get_num_diff_values(triple_pattern, graph_spo, crc_arrays); I really need a variant of this function that returns only a single value for OUR variable, not for every variable in the triple pattern.
                    //This new function has to consider that our variable is only once for each triple. therefore 1 CRC calc per triple has to be performed.
                    for (Triple *triple_pattern : triples)
                    {
                        uint64_t aux = get_num_diff_values(candidate_var, last_processed_var, last_processed_value, triple_pattern, *this->graph, *this->crc);
                        std::cout << aux << std::endl;//TODO: quedarme con el m'inimo y tambi'en recordar que debo setear el score.
                    }
                }
            }
        }
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
                if (level >= (this->gao->size() - 1)) {
                    // Print Answers
                    (*bindings)[varname] = binding_last.second;

                    /*
                    for(auto it = (*bindings).cbegin(); it != (*bindings).cend(); ++it) {
                        cout << "(" << it->first << ": " << (*bindings)[it->first] << ") ";
                    }
                    cout << endl;
                    */
                    (*number_of_results) = (*number_of_results) + 1;

                } else {
                    (*bindings)[varname] = binding_last.second;
                    processed_vars.push_back(varname);
                    processed_vars_map[varname] = (*bindings)[varname];
                    int new_level = level + 1;
                    evaluate(new_level, bindings, number_of_results, begin);
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
                    if (level >= (this->gao->size() - 1)) {
                        // Print Answers
                        (*bindings)[varname] = (*var_iterators)[0]->current_value();

                        /*
                        for(auto it = (*bindings).cbegin(); it != (*bindings).cend(); ++it) {
                            cout << "(" << it->first << ": " << (*bindings)[it->first] << ") ";
                        }
                        cout << endl;
                        */

                        (*number_of_results) = (*number_of_results) + 1;
                        int next_value = (*var_iterators)[0]->current_value() + 1;
                        (*var_iterators)[0]->seek(next_value);
                    } else {
                        (*bindings)[varname] = (*var_iterators)[0]->current_value();
                        processed_vars.push_back(varname);
                        processed_vars_map[varname] = (*bindings)[varname];
                        int new_level = level + 1;
                        evaluate(new_level, bindings, number_of_results, begin);
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

};


#endif //LEAPFROG_ITERATOR_H
