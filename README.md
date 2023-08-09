## Introduction
This repository helps you reproduce the experiments of [1].

## Build Project
Clone the project
```
$ git clone --recurse-submodules https://github.com/joeyweii/QuantumAssertion.git <path>
$ cd <path>
```
Build binary with CMake
```
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build .
$ cd ..
```
A binary called ```VanQiRA``` under the directory ```build```will be generated.

## Files Description
- ```benchmarks/```: benchmarks used in the experiments of [1]
- ```lib/sliqsys```: bit-slicing BDD framework for manipulating quantum entities
- ```lib/abc-exorcism-api```: exorcism algorithm for ESOP minimization
- ```inlcude/```: header files of source code
- ```src/```: source codes
- ```scripts/```: python scripts for experiements

## Table I
Synthesize the assertion circuit $U_a$ given a quantum circuit $U$ with assertion point being the final of $U$
```
$ ./build/VanQiRA --U <U qasm> --Ua <Ua qasm> --assert_point_scenario final
```

Example:
Tp synthesize the assertion circuit for the circuit ```bv_100.qasm```
```
$ ./build/VanQiRA --U ./benchmarks/bv_100.qasm --Ua Ua.qasm --assert_point_scenirio final
```

## Table II
Synthesize the assertion circuit $U_a$ given a quantum circuit $U$ with assertion point determined by sparsity
```
$ ./build/VanQiRA --U <U qasm> --Ua <Ua qasm> --assert_point_scenario 2/5
```
Run the script of obtaining tp/tn/fp/fn/detection rate
```
$ python3 scripts/detect.py <U qasm> <Ua qasm> <erro probability>
```

Example:
Synthesize the assertion circuit of ```bv_7.qasm``` and obtain its tp/tn/fp/tn/detection rate with error probability 1.3%
```
$ ./build/VanQiRA --U benchmarks/bv_7.qasm --Ua Ua.qasm --assert_point_scenario 2/5
$ python3 scripts/detect.py benchmarks/bv_7.qasm Ua.qasm 0.013
```

## Figure 7
Synthesize the assertion circuit $U_a$ given a quantum circuit $U$ with assertion point scenarios success rate/ gate count.
```
$ ./build/VanQiRA --U <U qasm> --Ua <Ua qasm> --assert_point_scenario success_rate --dp <Ua size>
$ ./build/VanQiRA --U <U qasm> --Ua <Ua qasm> --assert_point_scenario expected_gate_count --dp <Ua size>
```
Run the script of obtaining success rate/expected gate count
```
$ python3 scripts/executeErr.py <U.qasm> <Ua.qasm>
```

Example:
Synthesize the assertion circuit of ```hwb7.qasm``` and obtain its success rate/expected gate count with different assertion point scenarios 
```
$ ./build/VanQiRA --U benchmarks/hwb7.qasm --Ua Ua.qasm --assert_point_scenario success_rate --dp 8
$ ./build/VanQiRA --U benchmarks/hwb7.qasm --Ua Ua.qasm --assert_point_scenario expected_gate_count --dp 8
$ python3 scripts/executeErr.py benchmarks/hwb7.qasm Ua.qasm
```
