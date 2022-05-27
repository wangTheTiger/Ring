
#include<iostream>
#include "../ring_sop.hpp"

int main(int argc, char* argv[])
{
    ring_sop graph_sop;
    cout << " Loading the sop index..." << endl;
    graph_sop.load(string(argv[1])+"_sop");

    uint64_t subject = 1;
    auto predicates = graph_sop.get_P_given_S(subject);
    for(auto &predicate: predicates){
        std::cout << "subject_id: " << subject << ", predicate_id: " << predicate << "\n";
    }
    /*
    uint64_t object = 19650523;
    auto subjects = graph_sop.get_S_given_O(object);
    for(auto &subject: subjects){
        std::cout << "object: " << object << ", subject_id: " << subject << "\n";
    }
    */
}