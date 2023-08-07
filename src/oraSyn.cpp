#include "oraSyn.h"
#include "cudd.h"

extern void esopMinimize(std::vector<std::string> &esop, std::vector<std::string> &esop_min, int nVars);

EsopExtractionManager::EsopExtractionManager(DdManager *ddManager,
                                             DdNode *FRoot,
                                             int nVars)
    : _ddManager(ddManager)
    , _FRoot(FRoot)
    , _nVars(nVars)
    , _values(nVars, VarValue::DONTCARE) {}

EsopExtractionManager::~EsopExtractionManager() {}

void EsopExtractionManager::extract() {
    if (_FRoot == NULL) return;

    _hash.clear();
    _esop.clear();

    fullExpand(_FRoot);
    genPSDKRO(_FRoot);
}

void EsopExtractionManager::genPSDKRO(DdNode *F) {
    if (F == Cudd_ReadLogicZero(_ddManager)) return;
    if (F == Cudd_ReadOne(_ddManager)) {
        cube c;
        for (auto var : _vars) {
            if (_values[var] != VarValue::DONTCARE) c._iscare.set(var);
            if (_values[var] == VarValue::POSITIVE) c._polarity.set(var);
        }

        _esop.push_back(c);
        return;
    }

    // Find the best expansion by a cache lookup
    auto it = _hash.find(F);
    assert(it != _hash.end());
    ExpType expansion = it->second.first;

    // Determine the top-most variable
    auto varIdx = Cudd_NodeReadIndex(F);
    _vars.push_back(varIdx);

    // Calculate f0, f1, f2
    DdNode *F0, *F1, *F2;
    F0 = Cudd_NotCond(Cudd_E(F), Cudd_IsComplement(F));
    F1 = Cudd_NotCond(Cudd_T(F), Cudd_IsComplement(F));
    F2 = Cudd_bddXor(_ddManager, F0, F1);

    // Generate psdkro of the branches
    if (expansion == ExpType::pD) {
        _values[varIdx] = VarValue::DONTCARE;
        genPSDKRO(F0);
        _values[varIdx] = VarValue::POSITIVE;
        genPSDKRO(F2);
    } else if (expansion == ExpType::nD) {
        _values[varIdx] = VarValue::DONTCARE;
        genPSDKRO(F1);
        _values[varIdx] = VarValue::NEGATIVE;
        genPSDKRO(F2);
    } else {
        _values[varIdx] = VarValue::NEGATIVE;
        genPSDKRO(F0);
        _values[varIdx] = VarValue::POSITIVE;
        genPSDKRO(F1);
    }

    Cudd_RecursiveDeref(_ddManager, F2);
    _vars.pop_back();
    _values[varIdx] = VarValue::DONTCARE;
}

int EsopExtractionManager::fullExpand(DdNode *F) {
    if (F == Cudd_ReadLogicZero(_ddManager)) return 0u;
    if (F == Cudd_ReadOne(_ddManager)) return 1u;

    auto it = _hash.find(F);
    if (it != _hash.end()) return it->second.second;

    DdNode *F0, *F1, *F2;
    F0 = Cudd_NotCond(Cudd_E(F), Cudd_IsComplement(F));
    F1 = Cudd_NotCond(Cudd_T(F), Cudd_IsComplement(F));
    F2 = Cudd_bddXor(_ddManager, F0, F1);
    Cudd_Ref(F2);

    int cost0, cost1, cost2;
    cost0 = fullExpand(F0);
    cost1 = fullExpand(F1);
    cost2 = fullExpand(F2);

    int costmax = std::max(std::max(cost0, cost1), cost2);

    std::pair<ExpType, int> ret;
    if (costmax == cost0)
        ret = std::make_pair(ExpType::nD, cost1 + cost2);
    else if (costmax == cost1)
        ret = std::make_pair(ExpType::pD, cost0 + cost2);
    else
        ret = std::make_pair(ExpType::Sh, cost0 + cost1);

    _hash[F] = ret;
    return ret.second;
}

void EsopExtractionManager::getESOP(std::vector<std::string> &ret) {
    for (auto &cube : _esop)
        ret.push_back(cube.str(_nVars));
}

void synESOP(DdManager *ddManager,
             DdNode *ddNode,
             int nVars,
             std::vector<std::string> &ret) {
    EsopExtractionManager m(ddManager, ddNode, nVars);
    m.extract();
    std::vector<std::string> esop;
    m.getESOP(esop);
    esopMinimize(esop, ret, nVars);
}
