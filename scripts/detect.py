import re
import sys
import numpy as np
import random as rd

from qiskit import *
from qiskit.quantum_info import Statevector

def readQasm(filename, fPlusOneQubit = False):
    with open(filename) as file:
        lines = file.readlines()
    circ = None
    for line in lines:
        if(line[0] == '#'):
            continue

        words = re.split(r',| ', line)
        while('' in words):
            words.remove('')

        if words[0] == 'OPENQASM':
            continue
        elif words[0] == 'include':
            continue
        elif words[0] == 'qreg':
            nQubits = int(words[1][ words[1].index('q[')+2 : words[1].index(']')])
            if(fPlusOneQubit):
                nQubits += 1
            circ = QuantumCircuit(QuantumRegister(nQubits, 'q'))
        elif words[0] == 'creg':
            continue
        elif words[0] == 'measure':
            continue
        elif words[0] in ['h', 's', 't', 'sdg', 'tdg', 'x', 'y', 'z', 'ccx', 'cx', 'mcx']:
            qubits = []
            for item in words[1:]:
                qubits.append(int(item[ item.index('q[')+2 : item.index(']') ]))
            if words[0]=='h':
                assert(len(qubits) == 1)
                assert(qubits[0] < nQubits)
                circ.h(qubits[0])
            elif words[0]=='s':
                assert(len(qubits) == 1)
                assert(qubits[0] < nQubits)
                circ.s(qubits[0])
            elif words[0]=='t':
                assert(len(qubits) == 1)
                assert(qubits[0] < nQubits)
                circ.t(qubits[0])
            elif words[0]=='sdg':
                assert(len(qubits) == 1)
                assert(qubits[0] < nQubits)
                circ.sdg(qubits[0])
            elif words[0]=='tdg':
                assert(len(qubits) == 1)
                assert(qubits[0] < nQubits)
                circ.tdg(qubits[0])
            elif words[0]=='x':
                assert(len(qubits) == 1)
                assert(qubits[0] < nQubits)
                circ.x(qubits[0])
            elif words[0]=='y':
                assert(len(qubits) == 1)
                assert(qubits[0] < nQubits)
                circ.y(qubits[0])
            elif words[0]=='z':
                assert(len(qubits) == 1)
                assert(qubits[0] < nQubits)
                circ.z(qubits[0])
            elif words[0]=='ccx':
                assert(len(qubits) == 3)
                for i in range(3):
                    assert(qubits[i] < nQubits)
                circ.ccx(qubits[0], qubits[1], qubits[2])
            elif words[0]=='cx':
                assert(len(qubits) == 2)
                for i in range(2):
                    assert(qubits[i] < nQubits)
                circ.cx(qubits[0], qubits[1])
            elif words[0]=='mcx':
                assert(len(qubits) > 3)
                for i in range(len(qubits)):
                    assert(qubits[i] < nQubits)
                circ.mcx(qubits[:-1], qubits[-1])
        else:
            print("readQasm: Line \""+line[:-1]+"\" not supported")
    return circ

def readAssertPoint(filename):
    with open(filename) as file:
        line = file.readline()
        words = re.split(r',| ', line)
        while('' in words):
            words.remove('')

        assert(words[0] == "#AssertPoint")

        return int(words[1])


def circInsertError(circ, p):
    nGates = len(circ.data)
    nQubits = circ.num_qubits

    circ_error = QuantumCircuit(QuantumRegister(nQubits, 'q'))
    assertPointOffset = 0
    for i in range(nGates):
        circ_error.data.append(circ.data[i])
        if (rd.random() > p):
            continue

        for qubitReg in circ.data[i].qubits:
            qubit, _ = QuantumCircuit.find_bit(circ, qubitReg)

            gatetype = rd.randint(0, 2)
            if gatetype == 0: # X
                circ_error.x(qubit)
            elif gatetype == 1: # Y
                circ_error.y(qubit)
            else: # Z
                circ_error.z(qubit)

    return circ_error

def splitCirc(circ, point):
    nGates = len(circ.data)
    assert(point < nGates)
    nQubits = circ.num_qubits
    circ1 = QuantumCircuit(QuantumRegister(nQubits, 'q'))
    circ2 = QuantumCircuit(QuantumRegister(nQubits, 'q'))
    
    for i in range(nGates):
        if(i <= point):
            circ1.data.append(circ.data[i])
        else:
            circ2.data.append(circ.data[i])
    return circ1, circ2

def printCirc(circ):
    for gate in circ.data:
        print(gate.operation.name, end=' ')
        for qubitReg in gate.qubits:
            qubit, _ = QuantumCircuit.find_bit(circ, qubitReg)
            print(qubit, end=' ')
        print('')

def addCirc(circ, addCirc):
    assert(circ.num_qubits >= addCirc.num_qubits)
    for gate in addCirc.data:
        gateType = gate.operation.name
        qubits = []
        for qubitReg in gate.qubits:
            qubit, _ = QuantumCircuit.find_bit(addCirc, qubitReg)
            qubits.append(qubit)
        if(gateType == 'h'):
            circ.h(qubits)
        elif(gateType == 's'):
            circ.s(qubits)
        elif(gateType == 't'):
            circ.t(qubits)
        elif(gateType == 'sdg'):
            circ.sdg(qubits)
        elif(gateType == 'tdg'):
            circ.tdg(qubits)
        elif(gateType == 'x'):
            circ.x(qubits)
        elif(gateType == 'y'):
            circ.y(qubits)
        elif(gateType == 'z'):
            circ.z(qubits)
        elif(gateType == 'cx'):
            circ.cx(qubits[0], qubits[1])
        elif(gateType == 'ccx'):
            circ.ccx(qubits[0], qubits[1], qubits[2])
        elif(gateType == 'mcx'):
            circ.mcx(qubits[:-1], qubits[-1])
        else:
            assert(False) 

def fidNotOne(fid):
    return (abs(1 - fid) > 0.05)
    
if __name__ == '__main__':
    err_prob = 0.01
    trials = 1000

    U = readQasm(sys.argv[1])
    Ua = readQasm(sys.argv[2])
    assertPoint = readAssertPoint(sys.argv[2])

    nQubits = U.num_qubits

    U1, U2 = splitCirc(U, assertPoint)
    gdnState = Statevector.from_instruction(U1)

    fp = 0
    tp = 0
    fn = 0
    tn = 0
    for i in range(trials):
        U1_err = circInsertError(U1, err_prob)
        Ua_err = circInsertError(Ua, err_prob)


        errState = Statevector.from_instruction(U1_err)

        U_err_fid = abs(np.inner(gdnState, np.conj(errState)))**2

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

        U1Ua_err_fid0 = abs(np.inner(gdnState, np.conj(U1Ua_err_state0)))**2
        U1Ua_err_fid1 = abs(np.inner(gdnState, np.conj(U1Ua_err_state1)))**2

        if(fidNotOne(U1Ua_err_fid1)):
            tp += prob
        else:
            fp += prob

        if(fidNotOne(U1Ua_err_fid0)):
            fn += 1-prob
        else:
            tn += 1-prob

    fp = (fp/trials)*100
    tp = (tp/trials)*100
    fn = (fn/trials)*100
    tn = (tn/trials)*100

    print("TP: ", tp, " %")
    print("TN: ", tn, " %")
    print("FP: ", fp, " %")
    print("FN: ", fn, " %")
