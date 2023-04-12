#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include "bddSystem.h"

enum class AssertPointMode
{
    Final,
    Middle,
    Sparsity
};

class VanQiRA : public BDDSystem
{
public:

    // Constructor and Destructor
    explicit VanQiRA
    (
        int nQubits,
		int fInitBitWidth,
		int fBitWidthMode,
        bool fReorder
    );

    ~VanQiRA();

    void simUfindAssertPoint(const Circuit *circuit, const AssertPointMode assertPointMode);
    void synUa(const std::string filename);
private:
    int _nQubits;
	Tensor *_state;
    int _assertPoint;
    DdNode *_S;

    void initState(Tensor *tensor);
    void esop2qasm(const std::vector<std::string> &esop, const std::string filename);
    int cubeCountNumControl(const std::string &cube);
    void recordAssertPoint(int assertPoint);
};

#endif
