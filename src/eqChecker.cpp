#include "eqChecker.h"

// Constructor
EquivalenceChecker::EquivalenceChecker
(
    std::vector<GateType>& gates,
    std::vector<std::vector<int>>& qubits,
    int n,
    bool isReorder
)
:   BDDSystem
    ( 
        isReorder
    )
{
    _gates = gates;
    _qubits = qubits;
    _n = n;
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
    initIdentity();
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

    if (type == GateType::X) PauliX(qubit[0]);
    else if (type == GateType::Y) PauliY(qubit[0], right);
    else if (type == GateType::Z) PauliZ(qubit);
    else if (type == GateType::H) Hadamard(qubit[0]);
    else if (type == GateType::S) Phase_shift(2, qubit[0]);
    else if (type == GateType::SDG) Phase_shift_dagger(-2, qubit[0]);
    else if (type == GateType::T) Phase_shift(4, qubit[0]);
    else if (type == GateType::TDG) Phase_shift_dagger(-4, qubit[0]);
    else if (type == GateType::RX_PI_2) rx_pi_2(qubit[0], false);
    else if (type == GateType::RX_PI_2_DG) rx_pi_2(qubit[0], true);
    else if (type == GateType::RY_PI_2) ry_pi_2(qubit[0], right^false);
    else if (type == GateType::RY_PI_2_DG) ry_pi_2(qubit[0], right^true);
    else if (type == GateType::CX)
    {
        std::vector<int> ncont(0);
        int targ = qubit[1];
        qubit.pop_back();
        Toffoli(targ, qubit, ncont);
        ncont.clear();
    }
    else if (type == GateType::CZ) PauliZ(qubit);
    else if (type == GateType::SWAP)
    {
        std::vector<int> cont(0);
        Fredkin(qubit[0], qubit[1], cont);
        cont.clear();
    }
    else if (type == GateType::CSWAP)
    {
        int swapA = qubit[1], swapB = qubit[2];
        qubit.pop_back();
        qubit.pop_back();
        Fredkin(swapA, swapB, qubit);
    }
    else if (type == GateType::CCX)
    {
        std::vector<int> ncont(0);
        int targ = qubit.back();
        qubit.pop_back();
        Toffoli(targ, qubit, ncont);
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
void EquivalenceChecker::printResult() const
{
    std::cout << "{\n";
    std::cout << "\t#Qubits (n): " << _n << '\n';
    std::cout << "\tGatecount of circuit: " << _gates.size() << '\n';
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
