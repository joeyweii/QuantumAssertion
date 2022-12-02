#include "bddSystem.h"

/**Function*************************************************************

  Synopsis    [Initialize BDD manager/_zeroNode/_identityNode.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::ddInitialize()
{
    _ddManager = Cudd_Init(2*_n, 2*_n, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0); // 0~(n-1): 0-variables, n~(2n-1): 1-variables

    _zeroNode = Cudd_Not(Cudd_ReadOne(_ddManager));
    Cudd_Ref(_zeroNode);
}


/**Function*************************************************************

  Synopsis    [Initialize a basic state vector.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void BDDSystem::initState()
{
	int *basicState = new int[_n];
    for (int i = 0; i < _n; i++)
        basicState[i] = 0;

	DdNode *var, *tmp;
    _allBDD = new DdNode **[_w];
    for (int i = 0; i < _w; i++)
        _allBDD[i] = new DdNode *[_r];

    for (int i = 0; i < _r; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < _w - 1; j++)
            {
                _allBDD[j][i] = Cudd_Not(Cudd_ReadOne(_ddManager));
                Cudd_Ref(_allBDD[j][i]);
            }
            _allBDD[_w - 1][i] = Cudd_ReadOne(_ddManager);
            Cudd_Ref(_allBDD[_w - 1][i]);
            for (int j = _n - 1; j >= 0; j--)
            {
                var = Cudd_bddIthVar(_ddManager, j);
                if (basicState[j] == 0)
                    tmp = Cudd_bddAnd(_ddManager, Cudd_Not(var), _allBDD[_w - 1][i]);
                else
                    tmp = Cudd_bddAnd(_ddManager, var, _allBDD[_w - 1][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, _allBDD[_w - 1][i]);
                _allBDD[_w - 1][i] = tmp;
            }
        }
        else
        {
            for (int j = 0; j < _w; j++)
            {
                _allBDD[j][i] = Cudd_Not(Cudd_ReadOne(_ddManager));
                Cudd_Ref(_allBDD[j][i]);
            }
        }
    }

}

/**Function*************************************************************

  Synopsis    [Allocate new BDDs for each integer vector.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::allocBDD(DdNode ***Bdd, bool extend)
{
    DdNode *tmp;

    DdNode ***W = new DdNode **[_w];
    for (int i = 0; i < _w; i++)
        W[i] = new DdNode *[_r];

    for (int i = 0; i < _r - _inc; i++)
        for (int j = 0; j < _w; j++)
            W[j][i] = Bdd[j][i];

    for (int i = 0; i < _w; i++)
        delete[] Bdd[i];

    for (int i = 0; i < _w; i++)
        Bdd[i] = W[i];

    if (extend)
    {
        for (int i = _r - _inc; i < _r; i++)
        {
            for (int j = 0; j < _w; j++)
            {
                Bdd[j][i] = Cudd_ReadOne(_ddManager);
                Cudd_Ref(Bdd[j][i]);
                tmp = Cudd_bddAnd(_ddManager, Bdd[j][_r - _inc - 1], Bdd[j][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, Bdd[j][i]);
                Bdd[j][i] = tmp;
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
