#include "bddSystem.h"

/**Function*************************************************************

  Synopsis    [Initialize BDD manager/_zeroNode/_identityNode.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::ddInitialize()
{
    _ddManager = Cudd_Init(_n, _n, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0); // 0~(n-1): 0-variables, n~(2n-1): 1-variables
}

/**Function*************************************************************

  Synopsis    [Alloc a new quantum data object and return.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

QuantumData* BDDSystem::newQuantumData()
{
	QuantumData* quanData = new QuantumData();
	quanData->_allBDD = nullptr;
	quanData->_r = 32;
	quanData->_k = 0;	
	return quanData;
}

/**Function*************************************************************

  Synopsis    [Delete a quantum data object.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void BDDSystem::deleteQuantumData(QuantumData* quanData)
{
	for (int i = 0; i < _w; i++)
		for (int j = 0, end_j = quanData->_r; j < end_j; j++)
			Cudd_RecursiveDeref(_ddManager, quanData->_allBDD[i][j]);

	for (int i = 0; i < _w; i++)
		delete[] quanData->_allBDD[i];

	delete quanData;
}

/**Function*************************************************************

  Synopsis    [Allocate new BDDs for each integer vector.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::allocBDD(DdNode*** allBDD, int r, bool extend)
{
    DdNode *tmp;

    DdNode ***W = new DdNode **[_w];
    for (int i = 0; i < _w; i++)
        W[i] = new DdNode *[r];

    for (int i = 0; i < r - _inc; i++)
        for (int j = 0; j < _w; j++)
            W[j][i] = allBDD[j][i];

    for (int i = 0; i < _w; i++)
        delete[] allBDD[i];

    for (int i = 0; i < _w; i++)
        allBDD[i] = W[i];

    if (extend)
    {
        for (int i = r - _inc; i < r; i++)
        {
            for (int j = 0; j < _w; j++)
            {
                allBDD[j][i] = Cudd_ReadOne(_ddManager);
                Cudd_Ref(allBDD[j][i]);
                tmp = Cudd_bddAnd(_ddManager, allBDD[j][r - _inc - 1], allBDD[j][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, allBDD[j][i]);
                allBDD[j][i] = tmp;
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [Detect overflow in integer vectors.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int BDDSystem::overflow3(DdNode *g, DdNode *h, DdNode *crin) const
{
    DdNode *tmp, *dd1, *dd2;
    int overflow;

    dd1 = Cudd_bddXor(_ddManager, g, crin);
    Cudd_Ref(dd1);

    dd2 = Cudd_bddXnor(_ddManager, g, h);
    Cudd_Ref(dd2);

    tmp = Cudd_bddAnd(_ddManager, dd1, dd2);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(_ddManager, dd1);
    Cudd_RecursiveDeref(_ddManager, dd2);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(_ddManager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [Detect overflow in integer vectors -- for the case that h is 0.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int BDDSystem::overflow2(DdNode *g, DdNode *crin) const
{
    DdNode *tmp;
    int overflow;

    tmp = Cudd_bddAnd(_ddManager, Cudd_Not(g), crin);
    Cudd_Ref(tmp);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(_ddManager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [Update max #nodes.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::updateNodeCount()
{
    _nodeCount = std::max(_nodeCount, static_cast<unsigned long>(Cudd_ReadNodeCount(_ddManager)));
}
