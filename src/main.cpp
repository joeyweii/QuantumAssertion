#include <boost/program_options.hpp>
#include <sys/time.h> 
#include <fstream>

#include "VanQiRA.h"
#include "memMeasure.h"

extern void qasmParser(std::ifstream &inFile, std::vector<GateType> &gates, std::vector<std::vector<int> > &qubits, int &n);

int main(int argc, char **argv)
{
    // Program options
    namespace po = boost::program_options;
    po::options_description description("Options");
    description.add_options()
    ("help", "produce help message.")
    ("reorder", po::value<bool>()->default_value(true), "allow variable reordering or not.\n"
                                                             "0: disable 1: enable") 
    ("U", po::value<std::string>()->implicit_value(""), "(QASM file) circuit under assertion.")
    ("Ua", po::value<std::string>()->implicit_value(""), "(QASM file) output assertion circuit.")
	;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1)
    {
	    std::cout << description;
	    return 0;
	}


    bool isReorder = vm["reorder"].as<bool>();

    // Parse QASM files
    std::vector<GateType> gates;
    std::vector<std::vector<int>> qubits;
    int n; 

    std::ifstream inFile;

    inFile.open(vm["U"].as<std::string>()); 
    if (!inFile)
    {
        std::cerr << "Circuit file doesn't exist." << std::endl;
        return 0;
    }
    qasmParser(inFile, gates, qubits, n);
    inFile.close();

    struct timeval tStart, tFinish;
    double elapsedTime;
    double runtime;
    size_t memPeak;

    gettimeofday(&tStart, NULL);

    VanQiRA vanqira(gates, qubits, n, isReorder);

    vanqira.synthesis(vm["Ua"].as<std::string>());

    gettimeofday(&tFinish, NULL);
    elapsedTime = (tFinish.tv_sec - tStart.tv_sec) * 1000.0;
    elapsedTime += (tFinish.tv_usec - tStart.tv_usec) / 1000.0;

    runtime = elapsedTime / 1000.0;
    memPeak = getPeakRSS();
    
    vanqira.printInfo(runtime, memPeak);

    return 0;
}
