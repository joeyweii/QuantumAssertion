import re
import sys
import numpy as np
import random as rd

from utils import *
from qiskit import *
from qiskit.quantum_info import Statevector


if __name__ == '__main__':
    successCount = 100

    U = readQasm(sys.argv[1])
    Ua = readQasm(sys.argv[2])
    assertPoint = readAssertPoint(sys.argv[2])

    nQubits = U.num_qubits
    err_prob = 1. / len(U.data)
    gdnState = Statevector.from_instruction(U)
    U1, U2 = splitCirc(U, assertPoint)

    attempts = 0
    accNumGates = 0
    for i in range(successCount):
        print(i)

        while(True):
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

            if(rd.random() <= prob):
                accNumGates += (len(U1.data) + len(Ua.data))
                continue
            
            U2_err = circInsertError(U2, err_prob)
            U2_err.initialize(U1Ua_err_state0)
            U1UaU2_err_state = Statevector.from_instruction(U2_err)
            U1UaU2_err_fid = fidelity(gdnState, U1UaU2_err_state)

            attempts += 1
            accNumGates += (len(U1.data) + len(Ua.data) + len(U2.data))

            if(fidNotOne(U1UaU2_err_fid) == False):
                break

        delete_last_line()

    print("|Ua|: ", len(Ua.data))
    print("AssertPoint: ", assertPoint+1)
    print("Success Rate: %.2f"%(successCount/attempts))
    print("Average #Gate: %.2f"%(accNumGates/successCount))
    
