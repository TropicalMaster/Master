from finiteField import *
from GabidulinGao import *
from sympy import *
import copy
import time 
# Finite field variables
q = 2
m = 20
# Code variables
#n = 7
n = 20
k = 3


num_errors = int((n-k)/2)
print("NUM_ERROR",num_errors)
word = [[244],[11],[12]]
start_time = time.time()
FF = finiteField([2,2,2,1,2,1],[0,1,2,3,4,12],q,n)



gaocode = GaoLike(n,k,FF)

gaocodeword_initial = FF.Codeword(word.copy(),gaocode.G)
gaocodeword = copy.deepcopy(gaocodeword_initial)

for i in range(n):
    if i < num_errors:
        elem = [i+17]
    else:
        elem = []
    gaocodeword[i] = FF.add(gaocodeword[i],elem,"-")
end_init = time.time()

gaocode.gaoDecoding(gaocodeword)
print("Actual code gao",gaocodeword_initial)
print("Wrong code  gao",gaocodeword)
end_all = time.time()
print("init:",(end_init-start_time),"decode:",(end_all-end_init))