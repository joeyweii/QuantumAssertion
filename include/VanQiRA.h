#pragma once

#include "bddSystem.h"

enum class AssertPointMode {
    TwentyPercent,
    FortyPercent,
    SixtyPercent,
    EightyPercent,
    Final,
    SR,
    EG
};

class VanQiRA : public BDDSystem {
public:
    // Constructor and Destructor
    explicit VanQiRA(int nQubits);

    ~VanQiRA();

    void simUfindAssertPoint(const Circuit *circuit,
                             const AssertPointMode assertPointMode,
                             const double dp);
    void synUa(const std::string filename);

private:
    int _nQubits;
    Tensor *_state;
    int _assertPoint;
    DdNode *_S;

    void initState(Tensor *tensor);
    void writeQASM(const std::string filename,
                   const std::vector<std::string> &esop);
    int cubeCountNumControl(const std::string &cube);
    void recordAssertPoint(int assertPoint);
};
