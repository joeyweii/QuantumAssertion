#include "VanQiRA.h"

#include <limits>
#include <fstream>
#include <cmath>

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
    ),
    _nQubits(nQubits), _state(nullptr), _assertPoint(-1), _S(nullptr)
{
	_state = newTensor(fInitBitWidth, nQubits);
}

// Destructor
VanQiRA::~VanQiRA()
{
    if(_S) Cudd_RecursiveDeref(_ddManager, _S);
	deleteTensor(_state);
}

/**Function*************************************************************

  Synopsis    [Run the simulation.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void VanQiRA::simUfindAssertPoint(const Circuit *circuit, const AssertPointMode assertPointMode, const double dp)
{
    initState(_state);

    bool fTranspose = false;
    double minCost = std::numeric_limits<double>::max();
	for(int position = 0; position < circuit->getGateCount(); ++position)
    {
		applyGate(circuit->getGate(position), _state, fTranspose);

        switch(assertPointMode)
        {
            case AssertPointMode::Point1:
                if(position >= circuit->getGateCount()/5 && position < circuit->getGateCount()*2/5)
                {
                    double cost = 1 - sparsity(_state);
                    if(cost < minCost)
                    {
                        minCost = cost;
                        recordAssertPoint(position);
                    }
                }
                break;
            case AssertPointMode::Point2:
                if(position >= circuit->getGateCount()*2/5 && position < circuit->getGateCount()*3/5)
                {
                    double cost = 1 - sparsity(_state);
                    if(cost < minCost)
                    {
                        minCost = cost;
                        recordAssertPoint(position);
                    }
                }
                break;
            case AssertPointMode::Point3:
                if(position >= circuit->getGateCount()*3/5 && position < circuit->getGateCount()*4/5)
                {
                    double cost = 1 - sparsity(_state);
                    if(cost < minCost)
                    {
                        minCost = cost;
                        recordAssertPoint(position);
                    }
                }
                break;
            case AssertPointMode::Point4:
                if(position >= circuit->getGateCount()*4/5 && position < circuit->getGateCount())
                {
                    double cost = 1 - sparsity(_state);
                    if(cost < minCost)
                    {
                        minCost = cost;
                        recordAssertPoint(position);
                    }
                }
                break;
            case AssertPointMode::Final:
                if(position == circuit->getGateCount()-1)
                    recordAssertPoint(position);
                break;
            case AssertPointMode::SR:
                if(position >= circuit->getGateCount()/5)
                {
                    double d1 = position+1;
                    double d2 = circuit->getGateCount()-1-position;
                    double p = 1. - 1./(d1+d2);
                    double DR = sparsity(_state);
                    double cost = 1 - (pow(p, d1+dp+d2) / (pow(p, d1+dp) + (1-pow(p, d1+dp))*(1-DR)));

                    if(cost < minCost)
                    {
                        minCost = cost;
                        recordAssertPoint(position);
                    }
                }
                break;
            case AssertPointMode::EG:
                if(position >= circuit->getGateCount()/5)
                {
                    double d1 = position+1;
                    double d2 = circuit->getGateCount()-1-position;
                    double p = 1. - 1./(d1+d2);
                    double DR = sparsity(_state);
                    double cost = ((d1+dp)*(1-pow(p, d1+dp)) + (d1+dp+d2)*(1-DR*(1-pow(p, d1+dp))))/(pow(p, d1+dp+d2));

                    if(cost < minCost)
                    {
                        minCost = cost;
                        recordAssertPoint(position);
                    }
                }
                break;
            default:
                assert(0);
        }
    }
}

/**Function*************************************************************

  Synopsis    [Record the current state as the assertion point.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void VanQiRA::recordAssertPoint(int assertPoint)
{
    if(_S) Cudd_RecursiveDeref(_ddManager, _S);
    _S = sparsityDD(_state);
    _assertPoint = assertPoint;
}

/**Function*************************************************************

  Synopsis    [Synthesize the assertion circuit.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void VanQiRA::synUa(const std::string filename)
{
    std::vector<std::string> esop;

    assert(_S); assert(_assertPoint != -1);
	double sp = Cudd_CountMinterm(_ddManager, _S, _state->_rank) / pow(2, _state->_rank);
    std::cout << "Sparsity: " << sp << '\n';
    synESOP(_ddManager, _S, _nQubits, esop); 

    writeQASM(filename, esop);
}

/**Function*************************************************************

  Synopsis    [Convert a esop to a qasm file.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void VanQiRA::writeQASM(const std::string filename, const std::vector<std::string> &esop)
{
    std::ofstream outFile;

    outFile.open(filename, std::ios::out);
    if(!outFile) 
    {
        std::cerr << "Output file \"" << filename << "\" cannot be opened." << std::endl;
        return;
    }

    outFile << "#AssertPoint " << _assertPoint << '\n';
    outFile << "OPENQASM 2.0;\n";
    outFile << "include \"qelib1.inc\";\n";
    outFile << "qreg q[" << std::to_string(_nQubits+1) << "];\n";  

    std::vector<bool> isFlipped(_nQubits, false);
    for(const std::string& cube: esop)
    {
        assert(cube.size() == _nQubits);
        int nControls = cubeCountNumControl(cube);
        std::string tofGateStr;
        switch(nControls)
        {
            case 0: tofGateStr += "x "; break;
            case 1: tofGateStr += "cx "; break;
            case 2: tofGateStr += "ccx "; break;
            default: tofGateStr += "mcx "; break;

        }

        bool fPrintComma = false;
        for(int i = 0, end_i = _nQubits; i < end_i; ++i)
        {
            if(cube[i] == '1' || cube[i] == '0')
            {
                if
                (
                    (cube[i] == '1' && isFlipped[i]) ||
                    (cube[i] == '0' && !isFlipped[i])
                )
                {
                    outFile << "x q[" << std::to_string(i) << "];\n"; 
                    isFlipped[i] = !isFlipped[i];
                }

                if(!fPrintComma) fPrintComma = true;
                else tofGateStr += ", ";

                tofGateStr += "q[";
                tofGateStr += std::to_string(i);
                tofGateStr += "]"; 
            }
        }

        if(fPrintComma) tofGateStr += ", ";
        tofGateStr += "q[";
        tofGateStr += std::to_string(_nQubits);
        tofGateStr += "];\n";
        outFile << tofGateStr;
    }

    for(int i = 0, end_i = _nQubits; i < end_i; ++i)
    {
        if(isFlipped[i])
            outFile << "x q[" << std::to_string(i) << "]\n"; 
    }

    outFile.close();
}

int VanQiRA::cubeCountNumControl(const std::string &cube)
{
    int nControl = 0;
    for(const char literal: cube)
    {
        if(literal == '0' || literal == '1') 
            ++nControl;
    }
    return nControl;
}

/**Function*************************************************************

  Synopsis    [Initialize a basic state vector.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void VanQiRA::initState(Tensor *tensor)
{
    std::vector<int> basicState(tensor->_rank, 0);

	DdNode *var, *tmp;
    for(int i = 0; i < _w; ++i)
    {
        for(int j = 0; j < tensor->_r; ++j)
        {
            if(i == _w - 1 && j == 0)
            {
                tensor->_allBDD[_w - 1][j] = Cudd_ReadOne(_ddManager);
                Cudd_Ref(tensor->_allBDD[_w - 1][j]);
                for(int ivar = 0; ivar < tensor->_rank; ++ivar)
                {
                    var = Cudd_bddIthVar(_ddManager, ivar);
                    if (basicState[ivar] == 0)
                        tmp = Cudd_bddAnd(_ddManager, Cudd_Not(var), tensor->_allBDD[_w - 1][j]);
                    else
                        tmp = Cudd_bddAnd(_ddManager, var, tensor->_allBDD[_w - 1][j]);
                    Cudd_Ref(tmp);
                    Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[_w - 1][j]);
                    tensor->_allBDD[_w - 1][j] = tmp;
                }
            
            }
            else
            {
                tensor->_allBDD[i][j] = Cudd_Not(Cudd_ReadOne(_ddManager));
                Cudd_Ref(tensor->_allBDD[i][j]);
            }
        }
    }
}
