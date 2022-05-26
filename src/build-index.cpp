#include <iostream>
#include "triple_bwt.hpp"
#include <fstream>
#include <sdsl/construct.hpp>

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;


int main(int argc, char **argv)
{

    vector<spo_triple> D, E;

    std::ifstream ifs(argv[1]);
    uint64_t s, p , o;
    //Index SPO - OSP - POS (Cyclic)
    do {
        ifs >> s >> p >> o;
        D.push_back(spo_triple(s, p, o));
    } while (!ifs.eof());

    D.shrink_to_fit();
    cout << "--Indexing " << D.size() << " triples" << endl;
    memory_monitor::start();
    auto start = timer::now();

    triple_bwt ring_spo(D);
    auto stop = timer::now();
    memory_monitor::stop();
    cout << "  Index built " << ring_spo.size() << " bytes" << endl;

    ring_spo.save(string(argv[1]));
    cout << "Index saved" << endl;

    //Reverse Index SOP - PSO - OPS (Cyclic)
    for (uint64_t i = 0; i < D.size(); i++){
        auto& t = D[i];
        auto s = std::get<0>(t);
        auto p = std::get<1>(t);
        auto o = std::get<2>(t);
        //std::cout << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t) << std::endl;
        t = spo_triple(s, o, p);
        //std::cout << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t) << std::endl;
    }


    D.shrink_to_fit();
    cout << "--Indexing " << D.size() << " triples" << endl;
    memory_monitor::start();
    start = timer::now();

    triple_bwt ring_sop(D);
    D.clear(); //clear the D array since we don't need it anymore.
    stop = timer::now();
    memory_monitor::stop();
    cout << "  Index built " << ring_sop.size() << " bytes" << endl;

    ring_sop.save(string(argv[1]) + "_sop");
    cout << "Index saved" << endl;
    cout << duration_cast<seconds>(stop-start).count() << " seconds." << endl;
    cout << memory_monitor::peak() << " bytes." << endl;
    return 0;
}

