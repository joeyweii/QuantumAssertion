OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
ccx q[8],q[6],q[2];
ccx q[8],q[6],q[0];
ccx q[9],q[8],q[2];
ccx q[9],q[8],q[0];
ccx q[8],q[4],q[1];
ccx q[7],q[4],q[2];
ccx q[7],q[4],q[0];
cx q[6],q[2];
cx q[6],q[1];
ccx q[5],q[4],q[2];
ccx q[9],q[5],q[2];
ccx q[9],q[5],q[1];
ccx q[9],q[7],q[2];
x q[8];
ccx q[8],q[7],q[2];
ccx q[8],q[7],q[0];
mcx q[8],q[5],q[3],q[2];
mcx q[8],q[5],q[3],q[0];
mcx q[9],q[8],q[7],q[6],q[3],q[0];
mcx q[9],q[8],q[7],q[5],q[3],q[0];
mcx q[8],q[7],q[6],q[5],q[4],q[0];
mcx q[9],q[8],q[6],q[5],q[4],q[0];
x q[8];
x q[3];
ccx q[9],q[3],q[2];
ccx q[7],q[3],q[2];
ccx q[7],q[3],q[1];
ccx q[7],q[3],q[0];
mcx q[8],q[6],q[5],q[4],q[3],q[0];
x q[4];
ccx q[4],q[3],q[2];
ccx q[4],q[3],q[0];
mcx q[9],q[7],q[6],q[5],q[4],q[3],q[0];
x q[3];
mcx q[9],q[8],q[7],q[4],q[3],q[0];
x q[8];
ccx q[8],q[4],q[2];
ccx q[8],q[4],q[1];
ccx q[8],q[4],q[0];
x q[8];
x q[6];
ccx q[6],q[3],q[2];
ccx q[6],q[3],q[0];
mcx q[7],q[6],q[5],q[2];
mcx q[9],q[6],q[4],q[2];
x q[5];
x q[3];
mcx q[8],q[5],q[3],q[2];
mcx q[8],q[5],q[3],q[0];
mcx q[9],q[7],q[6],q[5],q[4],q[3],q[0];
mcx q[9],q[7],q[6],q[5],q[3],q[0];
x q[3];
mcx q[8],q[6],q[5],q[4],q[3],q[0];
mcx q[6],q[5],q[4],q[3],q[0];
x q[4];
x q[6];
x q[9];
mcx q[9],q[6],q[4],q[2];
mcx q[9],q[6],q[4],q[0];
ccx q[9],q[5],q[1];
x q[6];
mcx q[9],q[7],q[6],q[4],q[0];
x q[4];
mcx q[9],q[8],q[6],q[5],q[4],q[0];
x q[5];
x q[7];
ccx q[7],q[3],q[1];
ccx q[7],q[3],q[0];
mcx q[9],q[7],q[5],q[4],q[3],q[0];
mcx q[9],q[7],q[6],q[5],q[3],q[0];
x q[5];
mcx q[8],q[7],q[6],q[5],q[4],q[0];
x q[3];
mcx q[9],q[8],q[7],q[5],q[3],q[0];
mcx q[9],q[8],q[7],q[6],q[4],q[3],q[0];
x q[4];
x q[6];
mcx q[7],q[6],q[5],q[2];
x q[8];
mcx q[9],q[8],q[7],q[6],q[4],q[3],q[0];
