#ifndef _VANQIRA_H_
#define _VANQIRA_H_

#include "bddSystem.h"

class VanQiRA : public BDDSystem
{
public:
    // Constructor and Destructor
    VanQiRA
    (
        std::vector<GateType>& gates,
        std::vector<std::vector<int>>& qubits,
        int n,
        bool isReorder
    );

    ~VanQiRA()  
    {
        clear();
    }

    void synthesis(std::string pFileNameOut);
    void printInfo(double runtime, size_t memPeak) const;

private:
    std::vector<GateType> _gates;              // gates in the circuits. [#gate]
    std::vector<std::vector<int>> _qubits;     // ith qubits of gates in the circuit. [#gates]*[#qubits]
	QuantumData	*_stateVector; 
	DdNode* _S;

    void init();
    void initState();
    void applyGate(GateType type, std::vector<int> qubit, bool right);
    void simulate();
	void getVanishingEntries();
    void printResult();

    // Clean up VanQiRA
    void clear() 
    {
        _gates.clear();
        _qubits.clear();

		deleteQuantumData(_stateVector);

        Cudd_Quit(_ddManager);
    };
};

#endif
