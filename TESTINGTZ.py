from finiteField import *
from TZ import *
from sympy import *
import copy
import time 
# Finite field variables
q = 3
# Code variables
#n = 7
n=6
k = 3
#s = 5
s=1
# 10,11 works
word = [[40],[91],[4]]
num_errors = int((2*n-k)/2)
print("numerrors",num_errors)
start_time = time.time()

FF = finiteField([2,2,2,1,2,1],[0,1,2,3,4,12],q,2*n)
TZcode = TZ(n,k,s,num_errors,FF)
TZcodeword = TZcode.encode(word.copy())
TZcodeword_initial = copy.deepcopy(TZcodeword)

error_vector = []
for i in range((2*n)):
    
    if i < num_errors:
        elem = [i]
    else:
        elem = []
    error_vector.append(elem)
    TZcodeword[i] = TZcode.FF.add(TZcodeword[i],elem,"-")
end_init = time.time()

TZcode.TZdecoding(TZcodeword)
print("TZ Codeword       ",TZcodeword_initial)
print("Wrong word        ",TZcodeword)
end_all = time.time()
print("init:",(end_init-start_time),"decode:",(end_all-end_init))


