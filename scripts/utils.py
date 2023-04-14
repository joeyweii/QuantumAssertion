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
        elif words[0] in ['h', 's', 't', 'sdg', 'tdg', 'x', 'y', 'z', 'ccx', 'cx', 'mcx', 'swap', 'cswap']:
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
                assert(len(qubits) >= 2)
                for i in range(len(qubits)):
                    assert(qubits[i] < nQubits)
                circ.mcx(qubits[:-1], qubits[-1])
            elif words[0]=='swap':
                assert(len(qubits) == 2)
                for i in range(2):
                    assert(qubits[i] < nQubits)
                circ.swap(qubits[0], qubits[1])
            elif words[0]=='cswap':
                assert(len(qubits) == 3)
                for i in range(3):
                    assert(qubits[i] < nQubits)
                circ.cswap(qubits[0], qubits[1], qubits[2])
                        

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
        elif(gateType == 'mcx' or gateType == 'mcx_gray'):
            circ.mcx(qubits[:-1], qubits[-1])
        elif(gateType == 'swap'):
            circ.swap(qubits[0], qubits[1])
        elif(gateType == 'cswap'):
            circ.cswap(qubits[0], qubits[1], qubits[2])
        else:
            assert(False) 

def fidelity(state1, state2):
    return abs(np.inner(state1, np.conj(state2)))**2

def fidNotOne(fid):
    return (abs(1 - fid) > 0.05)

def delete_last_line():
    sys.stdout.write('\x1b[1A')
    sys.stdout.write('\x1b[2K')  
