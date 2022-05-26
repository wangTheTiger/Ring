#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <set>
#include "Triple.h"
#include "triple_bwt.hpp"

uint64_t get_size_interval(Triple * triple_pattern, triple_bwt & graph) {
    const uint64_t nTriples = graph.get_n_triples();
    if (triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph.open_SPO();
        return open_interval.size();//This should be faster than calculating the CRC WM + range in the whole BWT.L object.
    } else if (triple_pattern->s->isVariable && !triple_pattern->p->isVariable && triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph.open_PSO();
        uint64_t cur_p = graph.min_P(open_interval);//TODO: Is this function useful? we do the same thing in next_P.
        cur_p = graph.next_P(open_interval, triple_pattern->p->constant);
        if (cur_p == 0 || cur_p != triple_pattern->p->constant) {
            return 0;
        } else{
            bwt_interval i_p = graph.down_P(cur_p);
            uint64_t num_distinct_values = graph.calculate_gao_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples); //Going from P to S: values must be shifted back to the left, by substracting nTriples. Remember P is between nTriples + 1 and 2*nTriples.
            //TODO: Store these variables in a file. std::cout << "num_distinct_values = " << num_distinct_values << " vs. interval size = " << i_p.size() << std::endl;
            //return i_p.size();
            return num_distinct_values;
        }
    } else if (triple_pattern->s->isVariable && triple_pattern->p->isVariable && !triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph.open_OSP();
        uint64_t cur_o = graph.min_O(open_interval);
        cur_o = graph.next_O(open_interval, triple_pattern->o->constant);
        if (cur_o == 0 || cur_o != triple_pattern->o->constant) {
            return 0;
        } else{
            bwt_interval i_s = graph.down_O(cur_o);
            uint64_t num_distinct_values = graph.calculate_gao_BWT_S(i_s.left() - nTriples * 2, i_s.right() - nTriples * 2); //Going from O to S: values must be shifted back to the left, by substracting nTriples. Remember O is between 2*nTriples + 1 and 3*nTriples.
            //return i_s.size();
            return num_distinct_values;
        }
    } else if (!triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph.open_SPO();
        uint64_t cur_s = graph.min_S(open_interval);
        cur_s = graph.next_S(open_interval, triple_pattern->s->constant);
        if (cur_s == 0 || cur_s != triple_pattern->s->constant) {
            return 0;
        } else{
            bwt_interval i_s = graph.down_S(cur_s);
            uint64_t num_distinct_values = graph.calculate_gao_BWT_P(i_s.left() + nTriples, i_s.right() + nTriples); //Going from S to P: values must be shifted to the right by adding nTriples. Remember S is between 0 and nTriples and we want to go to P.
            //return i_s.size();
            return num_distinct_values;
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
            uint64_t num_distinct_values = graph.calculate_gao_BWT_O(i_p.left() + nTriples, i_p.right() + nTriples); //Going from P to O: values must be shifted to the right, by adding nTriples. Remember O is between 2* nTriples + 1 and 3*nTriples.
            //return i_p.size();
            return num_distinct_values;
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
            uint64_t num_distinct_values = graph.calculate_gao_BWT_P(i_o.left() - nTriples, i_o.right() - nTriples); //Going from O to P:values must be shifted back to the left, by substracting nTriples. Remember O is between 2*nTriples + 1 and 3*nTriples.
            //return i_o.size();
            return num_distinct_values;
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
            uint64_t num_distinct_values = graph.calculate_gao_BWT_S(i_o.left() - nTriples * 2, i_o.right() - nTriples * 2); //Going from O to S: values must be shifted back to the left, by substracting nTriples. Remember O is between 2*nTriples + 1 and 3*nTriples.
            //return i_o.size();
            return num_distinct_values;
        }
    }
    return 0;
}

//! TODO:
    /*!
        * \author Fabrizio Barisione
        * \returns TODO:
        */
uint64_t get_num_diff_values(Triple * triple_pattern, triple_bwt & graph_spo, triple_bwt & graph_sop) {
    const uint64_t nTriples = graph_spo.get_n_triples();
    if (triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph_spo.open_SPO();
        return open_interval.size();//This should be faster than calculating the CRC WM + range in the whole BWT.L object.
    } else if (triple_pattern->s->isVariable && !triple_pattern->p->isVariable && triple_pattern->o->isVariable) {
        //?S P ?O
        bwt_interval open_interval = graph_spo.open_PSO();
        uint64_t cur_p = graph_spo.min_P(open_interval);//TODO: Is this function useful? we do the same thing in next_P.
        cur_p = graph_spo.next_P(open_interval, triple_pattern->p->constant);
        if (cur_p == 0 || cur_p != triple_pattern->p->constant) {
            return 0;
        } else{
            bwt_interval i_p = graph_spo.down_P(cur_p); //Range C array pointing to S values.
            uint64_t num_distinct_values_s = graph_spo.calculate_gao_BWT_S(i_p.left() - nTriples, i_p.right() - nTriples); //Going from P to S: values must be shifted back to the left, by substracting nTriples. Remember P is between nTriples + 1 and 2*nTriples.
            //TODO: Store these variables in a file. std::cout << "num_distinct_values = " << num_distinct_values << " vs. interval size = " << i_p.size() << std::endl;
            //return i_p.size();
            //TODO: calcular lo mismo usando el otro orden.
            uint64_t num_distinct_values_o = graph_spo.calculate_gao_BWT_O(i_p.left() - nTriples, i_p.right() - nTriples);
            std::cout << "num_distinct_values S = " << num_distinct_values_s << " vs. interval size = " << i_p.size() << std::endl;
            std::cout << "num_distinct_values O = " << num_distinct_values_o << " vs. interval size = " << i_p.size() << std::endl;
            return num_distinct_values_s;
        }
    } else if (triple_pattern->s->isVariable && triple_pattern->p->isVariable && !triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph_spo.open_OSP();
        uint64_t cur_o = graph_spo.min_O(open_interval);
        cur_o = graph_spo.next_O(open_interval, triple_pattern->o->constant);
        if (cur_o == 0 || cur_o != triple_pattern->o->constant) {
            return 0;
        } else{
            bwt_interval i_s = graph_spo.down_O(cur_o);
            uint64_t num_distinct_values = graph_spo.calculate_gao_BWT_S(i_s.left() - nTriples * 2, i_s.right() - nTriples * 2); //Going from O to S: values must be shifted back to the left, by substracting nTriples. Remember O is between 2*nTriples + 1 and 3*nTriples.
            //return i_s.size();
            return num_distinct_values;
        }
    } else if (!triple_pattern->s->isVariable && triple_pattern->p->isVariable && triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph_spo.open_SPO();
        uint64_t cur_s = graph_spo.min_S(open_interval);
        cur_s = graph_spo.next_S(open_interval, triple_pattern->s->constant);
        if (cur_s == 0 || cur_s != triple_pattern->s->constant) {
            return 0;
        } else{
            bwt_interval i_s = graph_spo.down_S(cur_s);
            uint64_t num_distinct_values = graph_spo.calculate_gao_BWT_P(i_s.left() + nTriples, i_s.right() + nTriples); //Going from S to P: values must be shifted to the right by adding nTriples. Remember S is between 0 and nTriples and we want to go to P.
            //return i_s.size();
            return num_distinct_values;
        }
    } else if (!triple_pattern->s->isVariable && !triple_pattern->p->isVariable && triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph_spo.open_SPO();
        uint64_t cur_s = graph_spo.min_S(open_interval);
        cur_s = graph_spo.next_S(open_interval, triple_pattern->s->constant);
        if (cur_s == 0 || cur_s != triple_pattern->s->constant) {
            return 0;
        } else{
            bwt_interval i_s = graph_spo.down_S(cur_s);
            uint64_t cur_p = graph_spo.min_P_in_S(i_s, cur_s);
            cur_p = graph_spo.next_P_in_S(i_s, cur_s, triple_pattern->p->constant);
            if (cur_p == 0 || cur_p != triple_pattern->p->constant) {
              return 0;
            }
            bwt_interval i_p = graph_spo.down_S_P(i_s, cur_s, cur_p);
            uint64_t num_distinct_values = graph_spo.calculate_gao_BWT_O(i_p.left() + nTriples, i_p.right() + nTriples); //Going from P to O: values must be shifted to the right, by adding nTriples. Remember O is between 2* nTriples + 1 and 3*nTriples.
            //return i_p.size();
            return num_distinct_values;
        }
    } else if (!triple_pattern->s->isVariable && triple_pattern->p->isVariable && !triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph_spo.open_SOP();
        uint64_t cur_s = graph_spo.min_S(open_interval);
        cur_s = graph_spo.next_S(open_interval, triple_pattern->s->constant);
        if (cur_s == 0 || cur_s != triple_pattern->s->constant) {
            return 0;
        } else{
            bwt_interval i_s = graph_spo.down_S(cur_s);
            uint64_t cur_o = graph_spo.min_O_in_S(i_s);
            cur_o = graph_spo.next_O_in_S(i_s, triple_pattern->o->constant);
            if (cur_o == 0 || cur_o != triple_pattern->o->constant) {
              return 0;
            }
            bwt_interval i_o = graph_spo.down_S_O(i_s, cur_o);
            uint64_t num_distinct_values = graph_spo.calculate_gao_BWT_P(i_o.left() - nTriples, i_o.right() - nTriples); //Going from O to P:values must be shifted back to the left, by substracting nTriples. Remember O is between 2*nTriples + 1 and 3*nTriples.
            //return i_o.size();
            return num_distinct_values;
        }
    } else if (triple_pattern->s->isVariable && !triple_pattern->p->isVariable && !triple_pattern->o->isVariable) {
        bwt_interval open_interval = graph_spo.open_POS();
        uint64_t cur_p = graph_spo.min_P(open_interval);
        cur_p = graph_spo.next_P(open_interval, triple_pattern->p->constant);
        if (cur_p == 0 || cur_p != triple_pattern->p->constant) {
            return 0;
        } else{
            bwt_interval i_p = graph_spo.down_P(cur_p);
            uint64_t cur_o = graph_spo.min_O_in_P(i_p, cur_p);
            cur_o = graph_spo.next_O_in_P(i_p, cur_p, triple_pattern->o->constant);
            if (cur_o == 0 || cur_o != triple_pattern->o->constant) {
              return 0;
            }
            bwt_interval i_o = graph_spo.down_P_O(i_p, cur_p, cur_o);
            uint64_t num_distinct_values = graph_spo.calculate_gao_BWT_S(i_o.left() - nTriples * 2, i_o.right() - nTriples * 2); //Going from O to S: values must be shifted back to the left, by substracting nTriples. Remember O is between 2*nTriples + 1 and 3*nTriples.
            //return i_o.size();
            return num_distinct_values;
        }
    }
    return 0;
}
bool compare_by_second(pair<string, int> a, pair<string, int> b) {
    return a.second < b.second;
}

// Cambiar retorno
vector<string> get_gao_min_gen(vector<Triple*> query, triple_bwt & graph) {
    map<string, vector<uint64_t>> triple_values;
    map<string, vector<Triple*>> triples_var;
     for (Triple * triple_pattern : query) {
        uint64_t triple_size = get_size_interval(triple_pattern, graph);
        //Storing the triple_size of each variable for all triple_patterns.
        if (triple_pattern->s->isVariable) {
          triple_values[triple_pattern->s->varname].push_back(triple_size);
          triples_var[triple_pattern->s->varname].push_back(triple_pattern);
        }
        if (triple_pattern->p->isVariable) {
          triple_values[triple_pattern->p->varname].push_back(triple_size);
          triples_var[triple_pattern->p->varname].push_back(triple_pattern);
        }
        if (triple_pattern->o->isVariable) {
          triple_values[triple_pattern->o->varname].push_back(triple_size);
          triples_var[triple_pattern->o->varname].push_back(triple_pattern);
        }
    }

    vector<string> gao;
    vector<pair<string, uint64_t>> varmin_pairs;
    vector<string> single_vars;
    map<string, bool> selectable_vars;
    map<string, bool> selected_vars;
    map<string, set<string>> related_vars;

    for(auto it = (triples_var).cbegin(); it != (triples_var).cend(); ++it) {
        for (Triple * triple_pattern : triples_var[it->first]) {
            if (triple_pattern->s->isVariable && it->first.compare(triple_pattern->s->varname) != 0) {
                related_vars[it->first].insert(triple_pattern->s->varname);
            }
            if (triple_pattern->p->isVariable && it->first.compare(triple_pattern->p->varname) != 0) {
                related_vars[it->first].insert(triple_pattern->p->varname);
            }
            if (triple_pattern->o->isVariable && it->first.compare(triple_pattern->o->varname) != 0) {
                related_vars[it->first].insert(triple_pattern->o->varname);
            }
        }    
    }


    for(auto it = (triple_values).cbegin(); it != (triple_values).cend(); ++it) {
        if (triple_values[it->first].size() == 1) {
            single_vars.push_back(it->first);
        } else {
            varmin_pairs.push_back(pair<string, uint64_t>(it->first, *min_element(triple_values[it->first].begin(), triple_values[it->first].end()))); 
            selectable_vars[it->first] = false;
        }
    }
    

    sort((varmin_pairs).begin(), (varmin_pairs).end(), compare_by_second);
    if (varmin_pairs.size() > 0) {
        selectable_vars[varmin_pairs[0].first] = true;
    }
    for (int i = 0; i < varmin_pairs.size(); i++) {
        for (pair<string, uint64_t> varmin_pair : varmin_pairs) {
            if (selectable_vars[varmin_pair.first] && !selected_vars[varmin_pair.first]) {
                gao.push_back(varmin_pair.first);
                selected_vars[varmin_pair.first] = true;
                for (set<string>::iterator it=related_vars[varmin_pair.first].begin(); it!=related_vars[varmin_pair.first].end(); ++it) {
                    selectable_vars[*it] = true;
                }
                break;
            }
        }
    }

    for (pair<string, uint64_t> varmin_pair : varmin_pairs) {
        if (!selected_vars[varmin_pair.first]) {
            gao.push_back(varmin_pair.first);
        }
    }

    for (string s : single_vars) {
        gao.push_back(s);
    }

    return gao;

}

//! TODO:
    /*!
        * \author Fabrizio Barisione
        * \returns Vector of strings representing variable order.
        */
vector<string> get_gao(vector<Triple*> query, triple_bwt & graph_spo, triple_bwt & graph_sop) {
    map<string, vector<uint64_t>> triple_values;
    map<string, vector<Triple*>> triples_var;
     for (Triple * triple_pattern : query) {
        uint64_t triple_size = get_num_diff_values(triple_pattern, graph_spo, graph_sop);
        //Storing the triple_size of each variable for all triple_patterns.
        if (triple_pattern->s->isVariable) {
          triple_values[triple_pattern->s->varname].push_back(triple_size);
          triples_var[triple_pattern->s->varname].push_back(triple_pattern);
        }
        if (triple_pattern->p->isVariable) {
          triple_values[triple_pattern->p->varname].push_back(triple_size);
          triples_var[triple_pattern->p->varname].push_back(triple_pattern);
        }
        if (triple_pattern->o->isVariable) {
          triple_values[triple_pattern->o->varname].push_back(triple_size);
          triples_var[triple_pattern->o->varname].push_back(triple_pattern);
        }
    }

    vector<string> gao;
    vector<pair<string, uint64_t>> varmin_pairs;
    vector<string> single_vars;
    map<string, bool> selectable_vars;
    map<string, bool> selected_vars;
    map<string, set<string>> related_vars;

    for(auto it = (triples_var).cbegin(); it != (triples_var).cend(); ++it) {
        for (Triple * triple_pattern : triples_var[it->first]) {
            if (triple_pattern->s->isVariable && it->first.compare(triple_pattern->s->varname) != 0) {
                related_vars[it->first].insert(triple_pattern->s->varname);
            }
            if (triple_pattern->p->isVariable && it->first.compare(triple_pattern->p->varname) != 0) {
                related_vars[it->first].insert(triple_pattern->p->varname);
            }
            if (triple_pattern->o->isVariable && it->first.compare(triple_pattern->o->varname) != 0) {
                related_vars[it->first].insert(triple_pattern->o->varname);
            }
        }    
    }


    for(auto it = (triple_values).cbegin(); it != (triple_values).cend(); ++it) {
        if (triple_values[it->first].size() == 1) {
            single_vars.push_back(it->first);
        } else {
            varmin_pairs.push_back(pair<string, uint64_t>(it->first, *min_element(triple_values[it->first].begin(), triple_values[it->first].end()))); 
            selectable_vars[it->first] = false;
        }
    }
    

    sort((varmin_pairs).begin(), (varmin_pairs).end(), compare_by_second);
    if (varmin_pairs.size() > 0) {
        selectable_vars[varmin_pairs[0].first] = true;
    }
    for (int i = 0; i < varmin_pairs.size(); i++) {
        for (pair<string, uint64_t> varmin_pair : varmin_pairs) {
            if (selectable_vars[varmin_pair.first] && !selected_vars[varmin_pair.first]) {
                gao.push_back(varmin_pair.first);
                selected_vars[varmin_pair.first] = true;
                for (set<string>::iterator it=related_vars[varmin_pair.first].begin(); it!=related_vars[varmin_pair.first].end(); ++it) {
                    selectable_vars[*it] = true;
                }
                break;
            }
        }
    }

    for (pair<string, uint64_t> varmin_pair : varmin_pairs) {
        if (!selected_vars[varmin_pair.first]) {
            gao.push_back(varmin_pair.first);
        }
    }

    for (string s : single_vars) {
        gao.push_back(s);
    }

    return gao;

}

#endif 
