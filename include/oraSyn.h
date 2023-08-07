#pragma once

#include "bddSystem.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <bitset>
#include <cassert>
#include <string>

namespace psdkro {

    constexpr int bitwidth = 128;

    enum class VarValue : std::uint8_t {
        POSITIVE, // var = 1
        NEGATIVE, // var = 0
        DONTCARE  // var don't care
    };

    enum class ExpType : std::uint8_t {
        pD, // Positive Davio
        nD, // Negative Davio
        Sh, // Shannon

        // for incompletely specified function
        C0, // case C0 = 0
        C1, // case C1 = 0
        F0  // case F' = 0
    };

    struct cube {
        std::bitset<bitwidth> _polarity;
        std::bitset<bitwidth> _iscare;

        cube() {
            _polarity.reset();
            _iscare.reset();
        }

        // Return this cube in string type. '0'/'1'/'-'
        std::string str(const int nVar) const {
            std::string s;
            for (auto i = 0; i < nVar; ++i) {
                if (!_iscare.test(i))
                    s.push_back('-');
                else if (_polarity.test(i))
                    s.push_back('1');
                else
                    s.push_back('0');
            }
            return s;
        }
    };
} // namespace psdkro

using namespace psdkro;

class EsopExtractionManager {

public:
    // Constructor and Destructor
    EsopExtractionManager(DdManager *ddManager, DdNode *FRoot, int nVars);
    ~EsopExtractionManager();

    // extract algorithm
    void extract();
    void getESOP(std::vector<std::string> &ret);

private:
    // First pass: dicide the best expansion and calculate the cost
    int fullExpand(DdNode *F);

    // Second pass: generate PSDKRO
    void genPSDKRO(DdNode *F);

private:
    DdManager *_ddManager;         // cudd manager
    DdNode *_FRoot;                // root node of function to be extracted
    int _nVars;                    // the number of variables
    std::vector<int> _vars;        // for generating psdkro
    std::vector<VarValue> _values; // for generating psdkro
    std::unordered_map<DdNode *, std::pair<ExpType, int>>
        _hash; // the mapping between 1) BDD node and 2) expansion type & cost
    std::vector<cube> _esop; // storing the resulting esop
};

void synESOP(DdManager *ddManager,
             DdNode *ddNode,
             int nVars,
             std::vector<std::string> &ret);
