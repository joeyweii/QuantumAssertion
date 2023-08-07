# Usage: python3 genEntanglement.py <#qubits> <output file>

import sys

if __name__ == "__main__":
     
    nQubits = int(sys.argv[1])
    filename = sys.argv[2]

    assert(nQubits >= 1)
    with open(filename, "wt") as fout:
        fout.write('OPENQASM 2.0;\ninclude \"qelib1.inc\";\n')
        fout.write('qreg q[' + str(nQubits) + '];\n')
        fout.write('creg c[' + str(nQubits) + '];\n')
        
        fout.write('h q[' + str(0) +'];\n')

        for i in range(1, nQubits):
            fout.write('cx q[0], q[{}];\n'.format(i))
