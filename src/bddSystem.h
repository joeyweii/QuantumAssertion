#ifndef _BDDSYSTEN_H_
#define _BDDSYSTEM_H_

#include <iostream>
#include <cstdlib> 
#include <string> 
#include <vector>

#include "../cudd/cudd/cudd.h"
#include "../cudd/cudd/cuddInt.h"
#include "../cudd/util/util.h"
#include "gateType.h"

#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899

class BDDSystem
{

friend class EquivalenceChecker;

struct QuantumData
{
	DdNode ***_allBDD;
	int		_k;
};

public:
    BDDSystem(bool isReorder)
    :   _ddManager(nullptr), _zeroNode(nullptr), _identityNode(nullptr),
        _n(0), _r(32), _w(4), _inc(3), _isReorder(isReorder), _nodeCount(0)
    {}

    ~BDDSystem()  
    {
        clear();
    }

    /* gateOpe.cpp */
    void Toffoli(QuantumData *quanData, int targ, std::vector<int> cont, std::vector<int> ncont);
    void Fredkin(QuantumData *quanData, int swapA , int swapB, std::vector<int> cont);
    void Hadamard(QuantumData *quanData, int iqubit);
    void rx_pi_2(QuantumData *quanData, int iqubit, bool dagger);
    void ry_pi_2(QuantumData *quanData, int iqubit, bool tanspose);
    void Phase_shift(QuantumData *quanData, int phase, int iqubit); // phase can only be 2 to the power of an integer
    void Phase_shift_dagger(QuantumData *quanData, int phase, int iqubit);
    void PauliX(QuantumData *quanData, int iqubit);
    void PauliY(QuantumData *quanData, int iqubit, bool transpose);
    void PauliZ(QuantumData *quanData, std::vector<int> iqubit); // Z or CZ

private:
    DdManager *_ddManager;      // BDD manager.
    DdNode *_zeroNode;          // pointer to the zero node in BDD.
    DdNode *_identityNode;      // pointer to the root node of the identity BDD.
    int _n;                     // # of qubits.
    int _r;                     // resolution of integers.
    int _w;                     // # of integers = 4.
    int _inc;                   // add inc BDDs when overflow occurs, used in allocBDD.
    bool _isReorder;            // using reorder or not for BDD.
    unsigned long _nodeCount;   // node count.

    /* misc.cpp */
    void ddInitialize();
    void allocBDD(DdNode ***allBDD, bool extend);
    int overflow3(DdNode *g, DdNode *h, DdNode *crin) const;
    int overflow2(DdNode *g, DdNode *crin) const;
    void updateNodeCount();

    // Clean up BDD system
    void clear() {};
};

#endif
