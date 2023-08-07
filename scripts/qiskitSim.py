# Usage: python3 ./qiskitSim.py <qasm file>

import numpy as np
import random as rd
import sys
from qiskit import *
from qiskit.quantum_info import Statevector

if __name__ == '__main__':
  threshold = 1e-12

  circ = QuantumCircuit.from_qasm_file(sys.argv[1])
  nQubits = circ.num_qubits 
  nGates = len(circ.data)
  
  statevector = Statevector.from_instruction(circ)
  
  print(statevector)
