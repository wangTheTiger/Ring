/*
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

#include <iostream>
#include <utility>
#include <map>
#include "ring_spo.hpp"
#include "ring_sop.hpp"
#include <chrono>
#include "Term.h"
#include "Triple.h"
#include "Iterator.hpp"
#include "LeapfrogOP.hpp"
#include <boost/algorithm/string.hpp>
#include "utils.hpp"
#include<map>

using namespace boost;

//#include<chrono>
//#include<ctime>

using namespace std::chrono;

bool get_file_content(string filename, vector<string> & vector_of_strings)
{
    // Open the File
    ifstream in(filename.c_str());
    // Check if object is valid
    if(!in)
    {
        cerr << "Cannot open the File : " << filename << endl;
        return false;
    }
    string str;
    // Read the next line from File untill it reaches the end.
    while (getline(in, str))
    {
        // Line contains string of length > 0 then save it in vector
        if(str.size() > 0)
            vector_of_strings.push_back(str);
    }
    //Close The File
    in.close();
    return true;
}

bool is_number(string & s)
{
    return !s.empty() && find_if(s.begin(),
        s.end(), [](unsigned char c) { return !isdigit(c); }) == s.end();
}

Triple* get_triple(string & s, vector<Term*> & terms_created) {
    vector<string> terms_strings;
    split(terms_strings, s, is_any_of(" "), token_compress_on);

    Term* t1;
    Term* t2;
    Term* t3;

    if (is_number(terms_strings[0])) {
        uint64_t value;
        istringstream iss(terms_strings[0]);
        iss >> value;
        t1 = new Term(value);
    } else {
        t1 = new Term(terms_strings[0]);
    }

    if (is_number(terms_strings[1])) {
        uint64_t value;
        istringstream iss(terms_strings[1]);
        iss >> value;
        t2 = new Term(value);
    } else {
        t2 = new Term(terms_strings[1]);
    }

    if (is_number(terms_strings[2])) {
        uint64_t value;
        istringstream iss(terms_strings[2]);
        iss >> value;
        t3 = new Term(value);
    } else {
        t3 = new Term(terms_strings[2]);
    }
    terms_created.push_back(t1);
    terms_created.push_back(t2);
    terms_created.push_back(t3);
    return new Triple(t1, t2, t3);
}
/*TODO: remove it.
void set_scores(vector<Triple*>& query, vector<string>& gao) {
    for (Triple* triple_pattern : query) {
        triple_pattern->set_scores(gao);
    }
}
*/
int main(int argc, char* argv[])
{
    vector<string> dummy_queries;
    if(argc != 3){
        std::cout << "Missing parameters: ./query-index data_file query_file" << std::endl;
    }
    bool result = get_file_content(argv[2], dummy_queries);
    ring_spo graph_spo;
    cout << " Loading the spo index..."; fflush(stdout);
    graph_spo.load(string(argv[1]));
    cout << endl << " spo index loaded " << graph_spo.size() << " bytes" << endl;
    /*
    ring_sop graph_sop;
    cout << " Loading the sop index..."; fflush(stdout);
    graph_sop.load(string(argv[1])+"_sop");
    */
    crc_arrays crc_arrays;
    crc_arrays.load(string(argv[1]));
    cout << endl << " CRC arrays loaded " << endl;

    std::ifstream ifs;
    uint64_t nQ = 0;

    std::chrono::high_resolution_clock::time_point start, stop;
    double total_time = 0.0;
    duration<double> time_span;

    if(result)
    {
        int count = 1;
        for (string query_string : dummy_queries) {

            vector<Term*> terms_created;
            vector<Triple*> query;

            vector<string> tokens_query;
            split(tokens_query, query_string, is_any_of("."), token_compress_on);

            for (string token : tokens_query) {
                trim(token);
                Triple* triple_pattern = get_triple(token, terms_created);
                query.push_back(triple_pattern);
            }

            start = std::chrono::high_resolution_clock::now();// try with std::chrono::steady_clock
            //map<string, vector<Triple *>> triples_var;
            //vector<string> single_vars;
            //map<string, set<string>> related_vars;
            //vector<string> initial_gao = get_next_eliminated_variable(query, graph_spo, crc_arrays, triples_var, related_vars, single_vars);
            //TODO: no hay necesidad de obtener el gao completo, a menos que tome tiempo constante.
            //set_scores(query, initial_gao);

            map<string, uint64_t> bindings;
            int number_of_results = 0;
            LeapfrogOP lf(&graph_spo, &query, &crc_arrays, &bindings, &number_of_results);
            /*
            cout << "Query Details:" << endl;
            lf.print_query();
            lf.print_gao();
            lf.serialize();
            cout << "##########" << endl;
            */

            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            lf.evaluate(0, begin);
            //std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();

            stop = std::chrono::high_resolution_clock::now();
            //total_time = duration_cast<microseconds>(stop - start).count();
            time_span = duration_cast<microseconds>(stop - start);
            total_time = time_span.count();

            cout << nQ <<  ";" << number_of_results << ";" << (unsigned long long)(total_time*1000000000ULL) << std::endl;
            //cout << nQ <<  ";" << number_of_results << ";" << total_time << ";" << aux << ";" << total_time - aux << endl;
            nQ++;

            // cout << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << std::endl;

            //cout << "RESULTS QUERY " << count << ": " << number_of_results << endl;

            // Delete pointers and empty vectors
            for (Triple* triple_pattern : query) {
                delete triple_pattern;
            }

            for (Term* term : terms_created) {
                delete term;
            }
        count += 1;
        }

    }

	return 0;
}

