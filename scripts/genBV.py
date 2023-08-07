# Usage: python3 genBV.py <#qubits> <output file>

import sys
import random
import os

def generate_astring(nqubits, prob=1.0):
    """
        generate a random binary string as a hidden bit string
    """
    answer = []
    for i in range(nqubits):
        if random.random() <= prob:
            answer.append("1")
        else:
            answer.append("0")

    return "".join(answer)

nQubits = int(sys.argv[1])

hiddenString = generate_astring(nQubits-1)

filename = sys.argv[2]
with open(filename, "wt") as fout:
    fout.write('OPENQASM 2.0;\ninclude \"qelib1.inc\";\n')
    fout.write('qreg q[' + str(nQubits) + '];\n')
    fout.write('creg c[' + str(nQubits) + '];\n')

    # Apply Hadamard gates to the first
    # (nQubits - 1) before querying the oracle
    for i in range(nQubits - 1):
        fout.write('h q[' + str(i) +'];\n')

    # Apply 1 and Hadamard gate to the last qubit
    # for storing the oracle's answer
    fout.write('x q[' + str(nQubits-1) +'];\n')
    fout.write('h q[' + str(nQubits-1) +'];\n')

    # Apply the inner-product oracle
    hiddenString = hiddenString[::-1]
    for i in range(len(hiddenString)):
        if hiddenString[i] == "1":
            fout.write('cx q[' + str(i) +'], q[' + str(nQubits-1) +'];\n')
    hiddenString = hiddenString[::-1]

    # Apply Hadamard gates after querying the oracle
    for i in range(nQubits - 1):
        fout.write('h q[' + str(i) +'];\n')
