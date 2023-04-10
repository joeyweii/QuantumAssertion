#include <boost/program_options.hpp>
#include <sys/time.h> 
#include <fstream>

#include "qcSim.h"
#include "memMeasure.h"

extern Circuit* qasmParser(const std::string& filename);

int main(int argc, char **argv)
{
    // Program options
    namespace po = boost::program_options;
    po::options_description description("Options");
    description.add_options()
    ("help", "produce help message.")
    ("reorder", po::value<bool>()->default_value(true), "allow variable reordering or not.\n"
                                                             "0: disable 1: enable default: 1") 
	("bitwidth_control", po::value<int>()->default_value(0), "bitwidth control when overflowing\n"
												 "0: extendBitwidth 1: dropLSB")
	("init_bitwidth", po::value<int>()->default_value(32), "initial bitwidth r\n"
															"default: 32")
    ("circuit", po::value<std::string>()->implicit_value(""), "circuit under simulation.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1)
    {
	    std::cout << description;
	    return 0;
	}

    bool fReorder = vm["reorder"].as<bool>();
	int	 fBitWidthControl = vm["bitwidth_control"].as<int>();
	int	 fInitBitWidth = vm["init_bitwidth"].as<int>();

    Circuit *circuit = qasmParser(vm["circuit"].as<std::string>());

	int nQubits = circuit->getNumberQubits();

    struct timeval tStart, tFinish;
    double elapsedTime;
    double runtime;
    size_t memPeak;

    gettimeofday(&tStart, NULL);

    VanQiRA vanqira(nQubits, fInitBitWidth, fBitWidthControl, fReorder);

    vanqira.simulateCircuit(circuit);

    gettimeofday(&tFinish, NULL);
    elapsedTime = (tFinish.tv_sec - tStart.tv_sec) * 1000.0;
    elapsedTime += (tFinish.tv_usec - tStart.tv_usec) / 1000.0;

    runtime = elapsedTime / 1000.0;
    memPeak = getPeakRSS();
    
    std::cout << "----- Circuit Info. -----\n";
    std::cout << "#Qubits: " << nQubits << '\n';
    std::cout << "#Gates in circuit: " << circuit->getGateCount() << '\n';

    std::cout << "----- Resource Usage -----\n";
    std::cout << "Runtime: " << runtime << " seconds\n";
    std::cout << "Peak memory usage: " << memPeak << " bytes\n"; 
    std::cout << "--------------------------\n";

    return 0;
}
