
#include<iostream>
#include "../ring_spo.hpp"

int main(int argc, char* argv[])
{
    ring_spo graph_spo;
    cout << " Loading the spo index..." << endl;
    graph_spo.load(string(argv[1]));

    uint64_t subject = 1;
    auto predicates = graph_spo.get_P_given_S(subject);
    for(auto &predicate: predicates){
        std::cout << "subject_id: " << subject << ", predicate_id: " << predicate << "\n";
    }

    uint64_t object = 19650523;
    auto subjects = graph_spo.get_S_given_O(object);
    for(auto &subject: subjects){
        std::cout << "object: " << object << ", subject_id: " << subject << "\n";
    }
}