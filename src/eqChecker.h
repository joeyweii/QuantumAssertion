#ifndef _EQCHECKER_H_
#define _EQCHECKER_H_

#include "bddSystem.h"

class EquivalenceChecker : public BDDSystem
{
public:
    // Constructor and Destructor
    EquivalenceChecker
    (
        std::vector<GateType>& gates,
        std::vector<std::vector<int>>& qubits,
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
    std::vector<GateType> _gates;              // gates in the circuits. [#gate]
    std::vector<std::vector<int>> _qubits;     // ith qubits of gates in the circuit. [#gates]*[#qubits]
	QuantumData	*_stateVector; 

    void init();
    void initState(QuantumData *quanData);
    void applyGate(GateType type, std::vector<int> qubit, bool right);
    void calculateMiter();
    void printResult();
	double calSparsity(QuantumData *quanData);

    // Clean up EquivalenceChecker
    void clear() 
    {
        _gates.clear();
        _qubits.clear();

        for (int i = 0; i < _w; i++)
            for (int j = 0; j < _r; j++)
                Cudd_RecursiveDeref(_ddManager, _stateVector->_allBDD[i][j]);

		for (int i = 0; i < _w; i++)
			delete[] _stateVector->_allBDD[i];

        Cudd_Quit(_ddManager);
    };
};

#endif
