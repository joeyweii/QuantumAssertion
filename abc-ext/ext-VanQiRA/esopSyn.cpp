#include "esopSyn.h"

extern "C" Vec_Wec_t* MyExorcism( Vec_Wec_t * vEsop, int nIns);

EsopExtractionManager::EsopExtractionManager(DdManager* ddManager, DdNode* FRoot, int nVars)
: _ddManager(ddManager), _FRoot(FRoot), _nVars(nVars), _values(nVars, VarValue::DONTCARE)
{}

EsopExtractionManager::~EsopExtractionManager()
{}

void EsopExtractionManager::extract()
{
	if (_FRoot == NULL) return;

	_hash.clear();
	_esop.clear();

	fullExpand(_FRoot);
    genPSDKRO(_FRoot);
}


void EsopExtractionManager::genPSDKRO(DdNode *F)
{
	if (F == Cudd_ReadLogicZero(_ddManager))
		return;
	if (F == Cudd_ReadOne(_ddManager)) 
    {
		cube c; 
		for (auto var : _vars) {
			if(_values[var] != VarValue::DONTCARE) c._iscare.set(var);
			if(_values[var] == VarValue::POSITIVE) c._polarity.set(var);
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
	if (expansion == ExpType::pD)
    {
		_values[varIdx] = VarValue::DONTCARE;
		genPSDKRO(F0);
		_values[varIdx] = VarValue::POSITIVE;
		genPSDKRO(F2);
	} 
    else if (expansion == ExpType::nD)
    {
		_values[varIdx] = VarValue::DONTCARE;
		genPSDKRO(F1);
		_values[varIdx] = VarValue::NEGATIVE;
		genPSDKRO(F2);
	} 
    else 
    { 
		_values[varIdx] = VarValue::NEGATIVE;
		genPSDKRO(F0);
		_values[varIdx] = VarValue::POSITIVE;
		genPSDKRO(F1);
	}

	Cudd_RecursiveDeref(_ddManager, F2);
	_vars.pop_back();
	_values[varIdx] = VarValue::DONTCARE;
}

int EsopExtractionManager::fullExpand(DdNode *F)
{
	if (F == Cudd_ReadLogicZero(_ddManager))
		return 0u;
	if (F == Cudd_ReadOne(_ddManager))
		return 1u;
		
	auto it = _hash.find(F);
	if (it != _hash.end())
		return it->second.second;

    DdNode *F0, *F1, *F2;
	F0 = Cudd_NotCond(Cudd_E(F), Cudd_IsComplement(F));
	F1 = Cudd_NotCond(Cudd_T(F), Cudd_IsComplement(F));
	F2 = Cudd_bddXor(_ddManager, F0, F1); Cudd_Ref(F2);

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

Vec_Wec_t* EsopExtractionManager::getESOPWec()
{
	Vec_Wec_t *vEsop = Vec_WecAlloc(0);
	for(auto &cube: _esop)
	{
		Vec_Int_t *vCube = Vec_WecPushLevel(vEsop);
        Vec_IntGrow(vCube, _nVars + 2);
		for(int i = 0; i < _nVars; ++i)
		{
			if(cube._iscare.test(i))
			{
				if(cube._polarity.test(i))
					Vec_IntPush(vCube, 2*i);
				else
					Vec_IntPush(vCube, 2*i+1);
			}
		}

		Vec_IntPush(vCube, -1);
	}

	return vEsop;
}

int EsopExtractionManager::getNumCubes() const
{
    return _esop.size();
}


static void writeESOPWecIntoQasm(Vec_Wec_t *vEsop, int nVars, std::string &pFileName)
{
	std::fstream outFile;
	outFile.open(pFileName, std::ios::out);

	if(!outFile.is_open())
	{
		std::cerr << "[ERROR] Output QASM file cannot be opened." << std::endl;
		return;
	}

	outFile << "OPENQASM 2.0;\n";
    outFile << "include \"qelib1.inc\";\n";
	outFile << "qreg q[" << nVars+1 << "];\n";

	Vec_Int_t * vCube;	
	int c, k, Lit;

	std::vector<bool> phaseList(nVars, true);
	std::vector<int> controlList;

	Vec_WecForEachLevel( vEsop, vCube, c )
	{
		controlList.clear();

		Vec_IntForEachEntry( vCube, Lit, k )
		{
			if(Lit < 0) continue;
			int Var = Lit/2;
			if(Lit%2 == 0 && phaseList[Var] == false)
			{
				outFile << "x q[" << Var << "];\n";
				phaseList[Var] = true;
			}
			else if(Lit%2 == 1 && phaseList[Var] == true)
			{
				outFile << "x q[" << Var << "];\n";
				phaseList[Var] = false;
			}
			controlList.push_back(Var);
		}

		if(controlList.size() == 0)
			outFile << "x";
		else if(controlList.size() == 1)
			outFile << "cx";
		else if(controlList.size() == 2)
			outFile << "ccx";
		else
			outFile << "mcx";

		for(int i = 0, end_i = controlList.size(); i < end_i; ++i)
		{
			if(i != 0) outFile << ',';
			outFile << " q[" <<  controlList[i] << "]";
		}

		if(controlList.size() != 0) outFile << ',';
		outFile << " q[" << nVars << "];\n";
	}

	for(int i = 0; i < nVars; ++i)
	{
		if(phaseList[i] == false)
			outFile << "x q[" << i << "];\n";
	}

	outFile.close();
}

void synESOP(DdManager *ddManager, DdNode* ddNode, int nVars, std::string pFileNameOut)
{	
	// ESOP extraction
	EsopExtractionManager m(ddManager, ddNode, nVars);	
	m.extract();
	Vec_Wec_t *vEsop = m.getESOPWec();

	// ESOP minimization
	Vec_Wec_t *vEsop_min = MyExorcism(vEsop, nVars);

	// write ESOP into QASM
	writeESOPWecIntoQasm(vEsop_min, nVars, pFileNameOut);
}
