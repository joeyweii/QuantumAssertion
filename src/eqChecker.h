#ifndef _EQCHECKER_H_
#define _EQCHECKER_H_

#include "bddSystem.h"

class EquivalenceChecker : public BDDSystem
{
public:
    // Constructor and Destructor
    EquivalenceChecker
    (
        std::vector<std::vector<GateType>>& gates,
        std::vector<std::vector<std::vector<int>>>& qubits,
        int n,
        bool isReorder
    );

    ~EquivalenceChecker()  
    {
        clear();
    }

    void check();
    void printInfo(double runtime, size_t memPeak) const;

private:
    std::vector<std::vector<GateType>> _gates;              // gates in circuits. [nCircuit]*[#gate]
    std::vector<std::vector<std::vector<int>>> _qubits;     // ith qubits of gates in circuits. [nCircuit]*[#gates]*[#qubits]
    int _ratio;                                             // gate count ratio. |circuit2|/|circuit1|
    bool _isEq;                                             // if the result is equivalent or not.
    bool _isGatesSwap;                                      // if circuit1 and circuit2 are swapped.

    void invertCircuit(std::vector<GateType> &gate);
    void init();
    void applyGate(GateType type, std::vector<int> qubit, bool right);
    void calculateMiter();
    void checkFeq();
    void printResult() const;

    // Clean up EquivalenceChecker
    void clear() 
    {
        _gates[0].clear();
        _gates[1].clear();
        _gates.clear();
        _qubits[0].clear();
        _qubits[1].clear();
        _qubits.clear();
    };
};

#endif
