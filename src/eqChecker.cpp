#include "eqChecker.h"

// Constructor
EquivalenceChecker::EquivalenceChecker
(
    std::vector<GateType>& gates,
    std::vector<std::vector<int>>& qubits,
    int n,
    bool isReorder
)
:   
	BDDSystem
    ( 
        isReorder
    )
{
    _gates = gates;
    _qubits = qubits;
    _n = n;
	_stateVector = newQuantumData();
}

/**Function*************************************************************

  Synopsis    [Run the checking procedure.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::check()
{
    init();
    calculateMiter();
    printResult();
}

/**Function*************************************************************

  Synopsis    [Initialize checker.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::init()
{
    ddInitialize();
    initState(_stateVector);
}

/**Function*************************************************************

  Synopsis    [Initialize a basic state vector.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void EquivalenceChecker::initState(QuantumData* quanData)
{
	auto &allBDD = quanData->_allBDD;
	auto &r = quanData->_r;

	int *basicState = new int[_n];
    for (int i = 0; i < _n; i++)
        basicState[i] = 0;

	DdNode *var, *tmp;
    allBDD = new DdNode **[_w];
    for (int i = 0; i < _w; i++)
        allBDD[i] = new DdNode *[r];

    for (int i = 0; i < r; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < _w - 1; j++)
            {
                allBDD[j][i] = Cudd_Not(Cudd_ReadOne(_ddManager));
                Cudd_Ref(allBDD[j][i]);
            }
            allBDD[_w - 1][i] = Cudd_ReadOne(_ddManager);
            Cudd_Ref(allBDD[_w - 1][i]);
            for (int j = _n - 1; j >= 0; j--)
            {
                var = Cudd_bddIthVar(_ddManager, j);
                if (basicState[j] == 0)
                    tmp = Cudd_bddAnd(_ddManager, Cudd_Not(var), allBDD[_w - 1][i]);
                else
                    tmp = Cudd_bddAnd(_ddManager, var, allBDD[_w - 1][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, allBDD[_w - 1][i]);
                allBDD[_w - 1][i] = tmp;
            }
        }
        else
        {
            for (int j = 0; j < _w; j++)
            {
                allBDD[j][i] = Cudd_Not(Cudd_ReadOne(_ddManager));
                Cudd_Ref(allBDD[j][i]);
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [Apply a gate.]

  Description [Apply a gate to the left side of the matrix if right = 0, the right side otherwise.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::applyGate(GateType type, std::vector<int> qubit, bool right)
{
    if (right) for (int i = 0; i < qubit.size(); i++) qubit[i] += _n;

    if (type == GateType::X) PauliX(_stateVector, qubit[0]);
    else if (type == GateType::Y) PauliY(_stateVector, qubit[0], right);
    else if (type == GateType::Z) PauliZ(_stateVector, qubit);
    else if (type == GateType::H) Hadamard(_stateVector, qubit[0]);
    else if (type == GateType::S) Phase_shift(_stateVector, 2, qubit[0]);
    else if (type == GateType::SDG) Phase_shift_dagger(_stateVector, -2, qubit[0]);
    else if (type == GateType::T) Phase_shift(_stateVector, 4, qubit[0]);
    else if (type == GateType::TDG) Phase_shift_dagger(_stateVector, -4, qubit[0]);
    else if (type == GateType::RX_PI_2) rx_pi_2(_stateVector, qubit[0], false);
    else if (type == GateType::RX_PI_2_DG) rx_pi_2(_stateVector, qubit[0], true);
    else if (type == GateType::RY_PI_2) ry_pi_2(_stateVector, qubit[0], right^false);
    else if (type == GateType::RY_PI_2_DG) ry_pi_2(_stateVector, qubit[0], right^true);
    else if (type == GateType::CX)
    {
        std::vector<int> ncont(0);
        int targ = qubit[1];
        qubit.pop_back();
        Toffoli(_stateVector, targ, qubit, ncont);
        ncont.clear();
    }
    else if (type == GateType::CZ) PauliZ(_stateVector, qubit);
    else if (type == GateType::SWAP)
    {
        std::vector<int> cont(0);
        Fredkin(_stateVector, qubit[0], qubit[1], cont);
        cont.clear();
    }
    else if (type == GateType::CSWAP)
    {
        int swapA = qubit[1], swapB = qubit[2];
        qubit.pop_back();
        qubit.pop_back();
        Fredkin(_stateVector, swapA, swapB, qubit);
    }
    else if (type == GateType::CCX)
    {
        std::vector<int> ncont(0);
        int targ = qubit.back();
        qubit.pop_back();
        Toffoli(_stateVector, targ, qubit, ncont);
        ncont.clear();
    }

    if (_ddManager != NULL)
        updateNodeCount();
}

/**Function*************************************************************

  Synopsis    [Calculate the miter]

  Description [
               Apply gates in circuit1 and circuit2 to evolve matrx from identity interleavingly.
               The longer gateuit (gates[1]) is seen as circuit1 (applied to the left side of the matrix)
               to save computation overhead
              ]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::calculateMiter()
{
    int cntCir = 0;

    if (_isReorder) Cudd_AutodynEnable(_ddManager, CUDD_REORDER_SYMM_SIFT);

    while (cntCir < _gates.size())
    {
        applyGate(_gates[cntCir], _qubits[cntCir], false);
        cntCir++;
    }

    if (_isReorder) Cudd_AutodynDisable(_ddManager);
}

/**Function*************************************************************

  Synopsis    [Print the equivalence checking result.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::printResult()
{
    std::cout << "{\n";
    std::cout << "\t#Qubits (n): " << _n << '\n';
    std::cout << "\tGatecount of circuit: " << _gates.size() << '\n';
	std::cout << "\tr: " << _stateVector->_r << '\n';		
	std::cout << "\tSparsity: " << calSparsity(_stateVector) << std::endl;
    std::cout << "}\n";
}

/**Function*************************************************************

  Synopsis    [Print statistics.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void EquivalenceChecker::printInfo(double runtime, size_t memPeak) const
{
    std::cout << '\n';
    std::cout << "Runtime: " << runtime << " seconds\n";
    std::cout << "Peak memory usage: " << memPeak << " bytes\n"; 
}


/**Function*************************************************************

  Synopsis    [Calculate the sparsity of a state vector.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

double EquivalenceChecker::calSparsity(QuantumData *quanData)
{
    DdNode* ddS;
    ddS = Cudd_ReadLogicZero(_ddManager);
    Cudd_Ref(ddS);

    for (int i = 0; i < _w; i++)
    {
        for (int j = 0, end_j = quanData->_r; j < end_j; j++)
        {
            DdNode* tem = ddS;
            ddS = Cudd_bddOr(_ddManager, ddS, quanData->_allBDD[i][j]);
            Cudd_Ref(ddS);
            Cudd_RecursiveDeref(_ddManager, tem);
        }
    }

    double sparsity = 1 - Cudd_CountMinterm(_ddManager, ddS, _n)/pow(2, _n); 
    return sparsity;
}
