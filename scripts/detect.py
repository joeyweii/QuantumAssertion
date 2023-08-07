# Usage: python3 detect.py <U qasm> <Ua qasm> <error probability>

import re
import sys
import numpy as np
import random as rd

from utils import *
from qiskit import *
from qiskit.quantum_info import Statevector

if __name__ == '__main__':
    trials = 10

    U = readQasm(sys.argv[1])
    Ua = readQasm(sys.argv[2])
    assertPoint = readAssertPoint(sys.argv[2])

    nQubits = U.num_qubits

    U1, U2 = splitCirc(U, assertPoint)

    err_prob = float(sys.argv[3])

    gdnState = Statevector.from_instruction(U1)

    fp = 0
    tp = 0
    fn = 0
    tn = 0
    for i in range(trials):
        print(i)
        U1_err = circInsertError(U1, err_prob)
        Ua_err = circInsertError(Ua, err_prob)

        U1Ua_err = QuantumCircuit(nQubits+1)
        addCirc(U1Ua_err, U1_err)
        addCirc(U1Ua_err, Ua_err)

        U1Ua_err_state = Statevector.from_instruction(U1Ua_err)

        U1Ua_err_state0 = []
        for i in range(pow(2, nQubits)):
            U1Ua_err_state0.append(U1Ua_err_state[i])

        U1Ua_err_state1 = []
        prob = 0

        for i in range(pow(2, nQubits), pow(2, nQubits+1)):
            U1Ua_err_state1.append(U1Ua_err_state[i])
            prob += abs(U1Ua_err_state[i])**2

        if(1-prob != 0):
            U1Ua_err_state0 = [i/(1-prob) for i in U1Ua_err_state0]
        if(prob != 0):
            U1Ua_err_state1 = [i/prob for i in U1Ua_err_state1]

        U1Ua_err_fid0 = fidelity(gdnState, U1Ua_err_state0)
        U1Ua_err_fid1 = fidelity(gdnState, U1Ua_err_state1)

        if(fidNotOne(U1Ua_err_fid1)):
            tp += prob
        else:
            fp += prob

        if(fidNotOne(U1Ua_err_fid0)):
            fn += 1-prob
        else:
            tn += 1-prob

        delete_last_line()
    
    fp = (fp/trials)*100
    tp = (tp/trials)*100
    fn = (fn/trials)*100
    tn = (tn/trials)*100

    detectionRate = tp / (tp + fn) * 100

    print("#Qubits: ", nQubits)
    print("|U|: ", len(U.data))
    print("|Ua|: ", len(Ua.data))
    print("AssertPoint: ", assertPoint+1)
    print("ErrProp: ", err_prob)
    print("TP: ", tp, " %")
    print("TN: ", tn, " %")
    print("FP: ", fp, " %")
    print("FN: ", fn, " %")
    print("Detection Rate: ", detectionRate, " %")
