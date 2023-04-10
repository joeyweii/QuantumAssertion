#include "qcSim.h"


extern void synESOP(DdManager *ddManager, DdNode* ddNode, int nVars, std::vector<std::string> &ret);

// Constructor
VanQiRA::VanQiRA
(
    int nQubits,
	int fInitBitWidth,
	int fBitWidthControl,
    bool fReorder
)
:   
	BDDSystem
    ( 
		nQubits,
		fBitWidthControl,
        fReorder
    )
{
	_state = newTensor(fInitBitWidth, nQubits);
}

// Destructor
VanQiRA::~VanQiRA()
{
	deleteTensor(_state);
}

/**Function*************************************************************

  Synopsis    [Run the simulation.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void VanQiRA::simulateCircuit(const Circuit *circuit)
{
    initState(_state);

    bool fTranspose = false;
	for(int i = 0; i < circuit->getGateCount(); ++i)
		applyGate(circuit->getGate(i), _state, fTranspose);

    DdNode* S = sparsityDD(_state);
    synESOP(_ddManager, S, circuit->getNumberQubits(), _esop); 
    Cudd_RecursiveDeref(_ddManager, S);

    std::cout << "----- ESOP -----\n";
    for(const std::string &cube: _esop)
        std::cout << cube << '\n';
}

/**Function*************************************************************

  Synopsis    [Initialize a basic state vector.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void VanQiRA::initState(Tensor *tensor)
{
	int *basicState = new int[tensor->_rank];
    for (int i = 0; i < tensor->_rank; i++)
		basicState[i] = 0;

	DdNode *var, *tmp;

    for (int i = 0; i < tensor->_r; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < _w - 1; j++)
            {
                tensor->_allBDD[j][i] = Cudd_Not(Cudd_ReadOne(_ddManager));
                Cudd_Ref(tensor->_allBDD[j][i]);
            }
            tensor->_allBDD[_w - 1][i] = Cudd_ReadOne(_ddManager);
            Cudd_Ref(tensor->_allBDD[_w - 1][i]);
            for (int j = tensor->_rank - 1; j >= 0; j--)

            {
                var = Cudd_bddIthVar(_ddManager, j);
                if (basicState[j] == 0)
                    tmp = Cudd_bddAnd(_ddManager, Cudd_Not(var), tensor->_allBDD[_w - 1][i]);
                else
                    tmp = Cudd_bddAnd(_ddManager, var, tensor->_allBDD[_w - 1][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[_w - 1][i]);
                tensor->_allBDD[_w - 1][i] = tmp;
            }
        }
        else
        {
            for (int j = 0; j < _w; j++)
            {
                tensor->_allBDD[j][i] = Cudd_Not(Cudd_ReadOne(_ddManager));
                Cudd_Ref(tensor->_allBDD[j][i]);
            }
        }
    }
	delete basicState;
}
