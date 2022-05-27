#include <iostream>
#include "ring_spo.hpp"
#include "ring_sop.hpp"
#include <fstream>
#include <sdsl/construct.hpp>

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;


int main(int argc, char **argv)
{
    vector<spo_triple> D;
    //1. Read the source file.
    std::ifstream ifs(argv[1]);
    uint64_t s, p , o;
    do {
        ifs >> s >> p >> o;
        D.push_back(spo_triple(s, p, o));
    } while (!ifs.eof());

    D.shrink_to_fit();

    //2. Building the Index SPO - OSP - POS (Cyclic)
    /*{
        cout << "--Indexing " << D.size() << " triples" << endl;
        memory_monitor::start();
        auto start = timer::now();

        ring_spo ring_spo(D);
        auto stop = timer::now();
        memory_monitor::stop();
        cout << "  Index built " << ring_spo.size() << " bytes" << endl;

        ring_spo.save(string(argv[1]));
        cout << "Index saved" << endl;

        cout << duration_cast<seconds>(stop-start).count() << " seconds." << endl;
        cout << memory_monitor::peak() << " bytes." << endl;
    }*/

    for (uint64_t i = 0; i < D.size(); i++){
        auto& t = D[i];
        auto s = std::get<0>(t);
        auto p = std::get<1>(t);
        auto o = std::get<2>(t);
        //Due to space efficiency spo is used to build sop ring.
        t = spo_triple(s, o, p);
    }
    //3. Building the reverse Index SOP - PSO - OPS (Cyclic)
    {
        cout << "--Indexing " << D.size() << " triples" << endl;
        memory_monitor::start();
        auto start = timer::now();

        ring_sop ring_sop(D);
        D.clear();
        auto stop = timer::now();
        memory_monitor::stop();
        cout << "  Index built " << ring_sop.size() << " bytes" << endl;

        ring_sop.save(string(argv[1]) + "_sop");
        cout << "Index saved" << endl;

        cout << duration_cast<seconds>(stop-start).count() << " seconds." << endl;
        cout << memory_monitor::peak() << " bytes." << endl;
    }

    return 0;
}

