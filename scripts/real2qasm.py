# Usage: python3 real_2_qasm.py <real file> <qasm file> 

import sys
import os
import re

f = open(sys.argv[1], 'r')
lines = f.readlines()

with open(sys.argv[2], 'w') as outFile:
    outFile.write('OPENQASM 2.0;\n')
    outFile.write('include "qelib1.inc";\n')
    qubitList = {}
    for line in lines:
        words = line.split()

        if words[0] == '.numvars':
            outFile.write('qreg q[' + words[1] + '];\n')
            outFile.write('creg c[' + words[1] + '];\n')
        elif words[0] == '.variables':
            for i in range(1, len(words)):
                qubitList[words[i]] = i-1
        elif words[0] == '.version' or words[0] == '.constants' or words[0] == '.garbage' or words[0] == '.begin':
            continue
        elif words[0] == 'h1':
            assert(len(words) == 2)
            assert(qubitList.get(words[1]) != None)
            outFile.write('h q[' + str(qubitList[words[1]]) + '];\n')
        elif words[0] == 't1':
            assert(len(words) == 2)
            assert(qubitList.get(words[1]) != None)
            outFile.write('x q[' + str(qubitList[words[1]]) + '];\n')
        elif words[0] == 'y1':
            assert(len(words) == 2)
            assert(qubitList.get(words[1]) != None)
            outFile.write('y q[' + str(qubitList[words[1]]) + '];\n')
        elif words[0] == 'z1':
            assert(len(words) == 2)
            assert(qubitList.get(words[1]) != None)
            outFile.write('z q[' + str(qubitList[words[1]]) + '];\n')
        elif words[0] == 's1':
            assert(len(words) == 2)
            assert(qubitList.get(words[1]) != None)
            outFile.write('s q[' + str(qubitList[words[1]]) + '];\n')
        elif words[0] == 'q1:-2':
            assert(len(words) == 2)
            assert(qubitList.get(words[1]) != None)
            outFile.write('sdg q[' + str(qubitList[words[1]]) + '];\n')
        elif words[0] == 'q1:4':
            assert(len(words) == 2)
            assert(qubitList.get(words[1]) != None)
            outFile.write('t q[' + str(qubitList[words[1]]) + '];\n')
        elif words[0] == 'q1:-4':
            assert(len(words) == 2)
            assert(qubitList.get(words[1]) != None)
            outFile.write('tdg q[' + str(qubitList[words[1]]) + '];\n')
        elif words[0] == 'rx1:2' or words[0] == 'rx(pi/2)':
            assert(len(words) == 2)
            assert(qubitList.get(words[1]) != None)
            outFile.write('rx(pi/2) q[' + str(qubitList[words[1]]) + '];\n')
        elif words[0] == 'ry1:2' or words[0] == 'ry(pi/2)':
            assert(len(words) == 2)
            assert(qubitList.get(words[1]) != None)
            outFile.write('ry(pi/2) q[' + str(qubitList[words[1]]) + '];\n')
        elif words[0] == 't2':
            assert(len(words) == 3)
            outFile.write('cx')
            for i in range(1, len(words)):
                assert(qubitList.get(words[i]) != None)
                if(i != 1):
                    outFile.write(',')
                outFile.write(' q[' + str(qubitList[words[i]]) + ']')
            outFile.write(';\n')
        elif words[0] == 'z2':
            assert(len(words) == 3)
            outFile.write('cz')
            for i in range(1, len(words)):
                assert(qubitList.get(words[i]) != None)
                if(i != 1):
                    outFile.write(',')
                outFile.write(' q[' + str(qubitList[words[i]]) + ']')
            outFile.write(';\n')
        elif words[0][0] == 't':
            try:
                val = int(words[0][1:])
                assert(val >= 0)
            except ValueError:
                assert(0)

            assert(len(words) == val+1)
            outFile.write('mcx')
            for i in range (1, val+1):
                assert(qubitList.get(words[i]) != None)
                if(i != 1):
                    outFile.write(',')
                outFile.write(' q[' + str(qubitList[words[i]]) + ']')
            outFile.write(';\n')
        elif words[0][0] == 'f':
            try:
                val = int(words[0][1:])
                assert(val >= 0)
            except ValueError:
                assert(0)

            assert(len(words) == val+1)
            outFile.write('cswap')
            for i in range (1, val+1):
                assert(qubitList.get(words[i]) != None)
                if(i != 1):
                    outFile.write(',')
                outFile.write(' q[' + str(qubitList[words[i]]) + ']')
            outFile.write(';\n')
        elif words[0] == '.end':
            outFile.close()
            break
        elif words[0] == '.inputs':
            continue
        elif words[0] == '.outputs':
            continue
        elif words[0] == '.inputbus':
            continue
        elif words[0] == '.outputbus':
            continue  
        elif words[0][0] == '#':
            continue
        else:
            sys.exit(words[0] + ' not supported')
