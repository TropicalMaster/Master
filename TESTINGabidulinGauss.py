from GabidulinGauss import *
from finiteField import *
import time 
from sympy import *
import copy
# Finite field variables
q = 2
m = 20
# Code variables
#n = 7
n=20
k = 3
s=1

num_errors = int((n-k)/2)
print("NUM_ERROR",num_errors)
word = [[10],[11],[12]]
start_time = time.time()
finitefield = finiteField([1,2,1,1],[0,1,2,11],q,m)


gabidulincode = gabidulinCode(n,k,s,finitefield)

Gabidulincodeword_initial = finitefield.Codeword(word.copy(), gabidulincode.G)
Gabidulincodeword = copy.deepcopy(Gabidulincodeword_initial)



for i in range(n):
    if i < num_errors:
        elem = [i+17]
    else:
        elem = []
    Gabidulincodeword[i] = finitefield.add(Gabidulincodeword[i],elem,"-")
end_init = time.time()

gabidulincode.decodingGabidulin(Gabidulincodeword)

print("Actual code gab",Gabidulincodeword_initial)
print("Wrong code  gab",Gabidulincodeword)
end_all = time.time()
print("init:",(end_init-start_time),"decode:",(end_all-end_init))