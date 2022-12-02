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

public:
    BDDSystem(bool isReorder)
    :   _ddManager(nullptr), _allBDD(nullptr), _zeroNode(nullptr), _identityNode(nullptr),
        _k(0), _n(0), _r(32), _w(4), _inc(3), _isReorder(isReorder), _nodeCount(0)
    {}

    ~BDDSystem()  
    {
        clear();
    }

    /* gateOpe.cpp */
    void Toffoli(int targ, std::vector<int> cont, std::vector<int> ncont);
    void Fredkin(int swapA , int swapB, std::vector<int> cont);
    void Hadamard(int iqubit);
    void rx_pi_2(int iqubit, bool dagger);
    void ry_pi_2(int iqubit, bool tanspose);
    void Phase_shift(int phase, int iqubit); // phase can only be 2 to the power of an integer
    void Phase_shift_dagger(int phase, int iqubit);
    void PauliX(int iqubit);
    void PauliY(int iqubit, bool transpose);
    void PauliZ(std::vector<int> iqubit); // Z or CZ

private:
    DdManager *_ddManager;      // BDD manager.
    DdNode ***_allBDD;         // BDDs. [w=4][r]
    DdNode *_zeroNode;          // pointer to the zero node in BDD.
    DdNode *_identityNode;      // pointer to the root node of the identity BDD.
    int* _k;                    // k in algebraic representation.
    int _n;                     // # of qubits.
    int _r;                     // resolution of integers.
    int _w;                     // # of integers = 4.
    int _inc;                   // add inc BDDs when overflow occurs, used in allocBDD.
    bool _isReorder;            // using reorder or not for BDD.
    unsigned long _nodeCount;   // node count.

    /* misc.cpp */
    void ddInitialize();
    void initIdentity();
    void allocBDD(DdNode ***Bdd, bool extend);
    int overflow3(DdNode *g, DdNode *h, DdNode *crin) const;
    int overflow2(DdNode *g, DdNode *crin) const;
    void updateNodeCount();

    // Clean up BDD system
    void clear() 
    {
        for (int i = 0; i < _w; i++)
            for (int j = 0; j < _r; j++)
                Cudd_RecursiveDeref(_ddManager, _allBDD[i][j]);

		for (int i = 0; i < _w; i++)
			delete[] _allBDD[i];

        Cudd_Quit(_ddManager);
    };
};

#endif
