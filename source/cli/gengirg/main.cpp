
#include <iostream>
#include <chrono>
#include <map>
#include <list>
#include <string>
#include <cstring>
#include <fstream>
#include <algorithm>

#include <omp.h>

#include <girgs/girgs-version.h>
#include <girgs/Generator.h>
#include <girgs/BitManipulation.h>


using namespace std;
using namespace chrono;

map<string, string> parseArgs(int argc, char** argv) {
    map<string, string> params;
    for (int i = 1; i < argc; i++) {
        // Get current and next argument
        if (argv[i][0] != '-')
            continue;
        std::string arg = argv[i] + 1; // +1 to skip the -
        // advance one additional position if next is used
        std::string next = (i + 1 < argc ? argv[i++ + 1] : ""); 
        params[std::move(arg)] = std::move(next);
    }
    return params;
}


template<typename T>
void logParam(T value, const string& name) {
    cout << "\t" << name << "\t=\t" << value << '\n';
}

template<typename T>
void rangeCheck(T value, T min, T max, string name, bool lex = false, bool hex = false) {
    if (value < min || value > max || (value == min && lex) || (value == max && hex)) {
        cerr << "ERROR: parameter " << name << " = " << value << " is not in range "
            << (lex ? "(" : "[") << min << "," << max << (hex ? ")" : "]") << '\n';
        exit(1);
    }
    logParam(value, name);
}

template<typename T>
void listCheck(T value, list<T> allowed_values, string name, bool lex = false, bool hex = false) {
    if (std::find(allowed_values.begin(), allowed_values.end(), value) == allowed_values.end()) {
        cerr << "ERROR: parameter " << name << " = " << value << " is not one of:\n  ";
        for (auto val : allowed_values)
            cerr << val << " ";
        cerr << '\n';
        exit(1);
    }
    logParam(value, name);
}



int main(int argc, char* argv[]) {

    // write help
    if (argc < 2 || 0 == strcmp(argv[1], "--help") || 0 == strcmp(argv[1], "-help")) {
        clog << "usage: ./gengirg\n"
            << "\t\t[-n anInt]                      // number of nodes                          default 10000\n"
            << "\t\t[-d anInt]                      // dimension of geometry    range [1,5]     default 1\n"
            << "\t\t[-ple aFloat]                   // power law exponent       range (2,3]     default 2.5\n"
            << "\t\t[-alpha aFloat]                 // model parameter          range (1,inf]   default infinity\n"
            << "\t\t[-deg aFloat]                   // average degree           range [1,n)     default 10\n"
            << "\t\t[-wseed anInt]                  // weight seed                              default 12\n"
            << "\t\t[-pseed anInt]                  // position seed                            default 130\n"
            << "\t\t[-sseed anInt]                  // sampling seed                            default 1400\n"
            << "\t\t[-threads anInt]                // number of threads to use                 default 1\n"
            << "\t\t[-file aString]                 // file name for output (w/o ext)           default \"graph\"\n"
            << "\t\t[-wMethod generate|read|fixed]  // method to obtain weights                 default generate\n"
            << "\t\t[-wFile aString]                // file name for input of weights           default \"in-weights.txt\"\n"
            << "\t\t[-wFixed aFloat]                // fixed weight value                       default 40.0\n"
            << "\t\t[-pMethod generate|read]        // method to obtain positions               default generate\n"
            << "\t\t[-pFile aString]                // file name for input of positions         default \"in-positions.txt\"\n"
            << "\t\t[-noscale 0|1]                  // disable weight scaling                   default 0\n"
            << "\t\t[-dot 0|1]                      // write result as dot (.dot)               default 0\n"
            << "\t\t[-edge 0|1]                     // write result as edgelist (.txt)          default 0\n";
        return 0;
    }

    // write version
    if(argc > 1 && 0 == strcmp(argv[1], "--version")) {
        cout << "GIRGs command line interface.\n\n"
             << GIRGS_NAME_VERSION << '\n'
             << GIRGS_PROJECT_DESCRIPTION << '\n'
             << GIRGS_AUTHOR_ORGANIZATION << '\n'
             << GIRGS_AUTHOR_DOMAIN << '\n'
             << GIRGS_AUTHOR_MAINTAINER << '\n';
        return 0;
    }

    // read params
    auto params = parseArgs(argc, argv);
    auto n      = !params["n"    ].empty()  ? stoi(params["n"    ])     : 10000;
    auto d      = !params["d"    ].empty()  ? stoi(params["d"    ])     : 1;
    auto ple    = !params["ple"  ].empty()  ? stod(params["ple"  ])     : 2.5;
    auto alpha  = !params["alpha"].empty()  ? stod(params["alpha"])     : std::numeric_limits<double>::infinity();
    auto deg    = !params["deg"  ].empty()  ? stod(params["deg"  ])     : 10.0;
    auto wseed  = !params["wseed"].empty()  ? stoi(params["wseed"])     : 12;
    auto pseed  = !params["pseed"].empty()  ? stoi(params["pseed"])     : 130;
    auto sseed  = !params["sseed"].empty()  ? stoi(params["sseed"])     : 1400;
    auto threads= !params["threads"].empty()? stoi(params["threads"])   : 1;
    auto file   = !params["file" ].empty()  ? params["file"]            : "graph";
    auto wMethod= !params["wMethod"].empty()? params["wMethod"]         : "generate";
    auto wFile  = !params["wFile"].empty()  ? params["wFile"]           :  "in-weights.txt";
    auto wFixed = !params["wFixed"].empty() ? stod(params["wFixed"])    : 40.0;
    auto pMethod= !params["pMethod"].empty()? params["pMethod"]         : "generate";
    auto pFile  = !params["pFile"].empty()  ? params["pFile"]           :  "in-positions.txt";
    auto dot    = params["dot" ] == "1";
    auto edge   = params["edge"] == "1";
    auto noscale= params["noscale"] == "1";

    // log params and range checks
    cout << "using:\n";
    logParam(n, "n");
    rangeCheck(d, 1, 5, "d");
    rangeCheck(ple, 2.0, 3.0, "ple", true, false);
    rangeCheck(alpha, 1.0, std::numeric_limits<double>::infinity(), "alpha", true);
    rangeCheck(deg, 1.0, n-1.0, "deg");
    logParam(wseed, "wseed");
    logParam(pseed, "pseed");
    logParam(sseed, "sseed");
    rangeCheck(threads, 1, omp_get_max_threads(), "threads");
    omp_set_num_threads(threads);
    logParam(file, "file");
    listCheck(wMethod, {"generate", "read", "fixed"}, "wMethod");
    if (wMethod == "read")
        logParam(wFile, "wFile");
    if (wMethod == "fixed")
        logParam(wFixed, "wFixed");
    listCheck(pMethod, {"generate", "read"}, "pMethod");
    if (pMethod == "read")
        logParam(pFile, "pFile");
    logParam(noscale, "noscale");
    logParam(dot, "dot");
    logParam(edge, "edge");
    logParam(girgs::BitManipulation<1>::name(), "morton");
    cout << "\n";

    auto t1 = high_resolution_clock::now();

    std::vector<double> weights;
    auto t2 = high_resolution_clock::now();
    if (wMethod == "generate") {
        cout << "generating weights ...\t\t" << flush;
        weights = girgs::generateWeights(n, ple, wseed, threads > 1);
    } else if (wMethod == "fixed") {
        cout << "setting fixed weights ...\t" << flush;
        weights = std::vector<double>(n, wFixed);
    } else {
        cout << "reading weights from file ...\t" << flush;
        weights = girgs::readWeights(n, wFile);
    }
    cout << "done in " << duration_cast<milliseconds>(t2 - t1).count() << "ms\tlargest = ";
    cout << *max_element(weights.begin(), weights.end()) << endl;

    std::vector<std::vector<double>> positions;
    auto t3 = high_resolution_clock::now();
    if (pMethod == "read") {
        cout << "reading positions from file ...\t" << flush;
        positions = girgs::readPositions(n, d, pFile);;
    } else {
        cout << "generating positions ...\t" << flush;
        positions = girgs::generatePositions(n, d, pseed);
    }
    cout << "done in " << duration_cast<milliseconds>(t3 - t2).count() << "ms" << endl;


    cout << "find weight scaling ...\t\t" << flush;
    auto t4 = high_resolution_clock::now();
    if (!noscale) {
        auto scaling = girgs::scaleWeights(weights, deg, d, alpha);
        cout << "done in " << duration_cast<milliseconds>(t4 - t3).count() << "ms\tscaling = " << scaling << endl;
    } else {
        cout << "skipped" << endl;
    }

    cout << "sampling edges ...\t\t" << flush;
    auto edges = girgs::generateEdges(weights, positions, alpha, sseed);
    auto t5 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t5 - t4).count() << "ms\tavg deg = " << edges.size()*2.0/n << endl;

    if (dot) {
        cout << "writing .dot file ...\t\t" << flush;
        auto t6 = high_resolution_clock::now();
        girgs::saveDot(weights, positions, edges, file+".dot");
        auto t7 = high_resolution_clock::now();
        cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    }

    if (edge) {
        cout << "writing edge list (.txt) ...\t" << flush;
        auto t6 = high_resolution_clock::now();
        ofstream f{file+".txt"};
        if(!f.is_open()) throw std::runtime_error{"Error: failed to open file \"" + file + ".txt\""};
        f << n << ' ' << edges.size() << "\n\n";
        for(auto& each : edges)
            f << each.first << ' ' << each.second << '\n';
        auto t7 = high_resolution_clock::now();
        cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    }

    return 0;
}
