#include <string>
#include <vector>
#include <iostream>

#include "misc/vec/vec.h"

extern "C" Vec_Wec_t* MyExorcism( Vec_Wec_t * vEsop, int nIns);

static Vec_Wec_t* vector2wec(std::vector<std::string> &esop, int nVars)
{
	Vec_Wec_t *vEsop = Vec_WecAlloc(0);
	for(auto &cube: esop)
	{
		Vec_Int_t *vCube = Vec_WecPushLevel(vEsop);
        Vec_IntGrow(vCube, nVars + 2);
		for(int i = 0; i < nVars; ++i)
		{
			if(cube[i] == '0')
				Vec_IntPush(vCube, 2*i+1);
			else if(cube[i] == '1')
				Vec_IntPush(vCube, 2*i);
		}
		Vec_IntPush(vCube, -1);
	}
	return vEsop;
}

static void wec2vector(Vec_Wec_t *vEsop, std::vector<std::string> &esop, int nVars)
{
	Vec_Int_t * vCube;	
	int c, k, Lit;

	Vec_WecForEachLevel( vEsop, vCube, c )
	{
		std::string s(nVars, '-');
		Vec_IntForEachEntry( vCube, Lit, k )
		{
			if(Lit < 0) continue;
			int Var = Lit/2;
			if(Lit%2 == 0)
				s[Var] = '1';
			else if(Lit%2 == 1)
				s[Var] = '0';
		}
		esop.push_back(s);
	}
}

/*
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
*/

void esopMinimize(std::vector<std::string> &esop, std::vector<std::string> &esop_min, int nVars)
{
	Vec_Wec_t* vEsop = vector2wec(esop, nVars);
	Vec_Wec_t* vEsop_min = MyExorcism(vEsop, nVars);
	wec2vector(vEsop_min, esop_min, nVars);
	Vec_WecFree(vEsop); 
	Vec_WecFree(vEsop_min); 
}
