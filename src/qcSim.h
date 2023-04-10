#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include "bddSystem.h"

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

    void simulateCircuit(const Circuit *circuit);
private:
	Tensor *_state;
    std::vector<std::string> _esop;
    void initState(Tensor *tensor);
};

#endif
