#include <sys/time.h>
#include <fstream>

#include "VanQiRA.h"
#include "memMeasure.h"

int main(int argc, char **argv) {
    bool fReorder = true;
    int fInitBitWidth = 4;
    std::string UFilename;
    std::string UaFilename;
    AssertPointMode assertPointMode;
    double dp = 0.;

    for (int i = 1; i < argc;) {
        if (std::string(argv[i]) == "--U") {
            UFilename = std::string(argv[i + 1]);
            i += 2;
        } else if (std::string(argv[i]) == "--Ua") {
            UaFilename = std::string(argv[i + 1]);
            i += 2;
        } else if (std::string(argv[i]) == "--reorder") {
            if (std::string(argv[i + 1]) == "0")
                fReorder = false;
            else if (std::string(argv[i + 1]) == "1")
                fReorder = true;
            else
                assert(0 && "Parse -reorder option fail");
            i += 2;
        } else if (std::string(argv[i]) == "--init_bitwidth") {
            fInitBitWidth = std::stoi(std::string(argv[i + 1]));
            i += 2;
        } else if (std::string(argv[i]) == "--dp") {
            dp = std::stod(std::string(argv[i + 1]));
            i += 2;
        } else if (std::string(argv[i]) == "--assert_point_scenario") {
            if (std::string(argv[i + 1]) == "1/5")
                assertPointMode = AssertPointMode::TwentyPercent;
            else if (std::string(argv[i + 1]) == "2/5")
                assertPointMode = AssertPointMode::FortyPercent;
            else if (std::string(argv[i + 1]) == "3/5")
                assertPointMode = AssertPointMode::SixtyPercent;
            else if (std::string(argv[i + 1]) == "4/5")
                assertPointMode = AssertPointMode::EightyPercent;
            else if (std::string(argv[i + 1]) == "final")
                assertPointMode = AssertPointMode::Final;
            else if (std::string(argv[i + 1]) == "success_rate")
                assertPointMode = AssertPointMode::SR;
            else if (std::string(argv[i + 1]) == "expected_gate_count")
                assertPointMode = AssertPointMode::EG;

            i += 2;
        } else
            assert(0 && "Undefined options.");
    }

    Circuit *circuit = parseQASM(UFilename);

    int nQubits = circuit->getNumberQubits();

    struct timeval tStart, tFinish;
    double elapsedTime;
    double runtime;
    size_t memPeak;

    gettimeofday(&tStart, NULL);

    VanQiRA vanqira(nQubits);
    vanqira.setAutoReorder(fReorder);
    vanqira.setInitBitWidth(fInitBitWidth);

    vanqira.simUfindAssertPoint(circuit, assertPointMode, dp);
    vanqira.synUa(UaFilename);

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
