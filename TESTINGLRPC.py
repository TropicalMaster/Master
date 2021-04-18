from LRPC import *
from finiteField import *
from sympy import *
import time 
import copy
# Finite field variables
q = 2
m = 12
# Code variables
n= 12
k = 3
d = 2
s=1
num_errors = 2

word = [[40],[90]]

print("numerrors",num_errors)
start_time = time.time()
finitefield = finiteField([1,2,1,1],[0,1,2,11],q,m)


LRPCcode = LRPC(n,k,d,num_errors,finitefield)
LRPCcodeword_initial = finitefield.Codeword(copy.deepcopy(word), LRPCcode.G)
LRPCcodeword = copy.deepcopy(LRPCcodeword_initial)


for i in range(n):
    if i < num_errors:
        elem = [i]
    else:
        elem = []
    LRPCcodeword[i] = finitefield.add(LRPCcodeword[i],elem,"-")
end_init = time.time()

LRPCcode.decodingLRPC(LRPCcodeword)
print("Actual code ",LRPCcodeword_initial)
print("Wrong code  ",LRPCcodeword)
end_all = time.time()
print("init:",(end_init-start_time),"decode:",(end_all-end_init))