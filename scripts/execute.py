#Usage: python3 execute.py <U qasm> <Ua qasm>

import re
import sys
import numpy as np
import random as rd

from utils import *
from qiskit import *
from qiskit.quantum_info import Statevector


if __name__ == '__main__':
    successCount = 5000

    U = readQasm(sys.argv[1])
    nQubits = U.num_qubits

    err_prob = 1. / len(U.data)

    gdnState = Statevector.from_instruction(U)

    attempts = 0
    accNumGates = 0
    for i in range(successCount):
        while(True):
            U_err = circInsertError(U, err_prob)

            U_err_state = Statevector.from_instruction(U_err)
            U_err_fid = fidelity(gdnState, U_err_state)

            attempts += 1
            accNumGates += (len(U.data))

            if(fidNotOne(U_err_fid) == False):
                break

    print("#Qubits: ", nQubits)
    print("ErrProp: ", err_prob)
    print("|U|: ", len(U.data))
    print("Success Rate: %.2f"%(successCount/attempts))
    print("Average #Gate: %.2f"%(accNumGates/successCount))
