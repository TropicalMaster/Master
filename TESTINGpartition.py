from finiteField import *
from NonAdditivePartition import *
from sympy import *
import copy
import time 
# Finite field variables
q = 3
# Code variables
n = 12
k = 3
s=7

word = [[40],[90],[67]]
num_errors = int((n-k)/2)
print("numerrors",num_errors)
start_time = time.time()
FF = finiteField([2,2,2,1,2,1],[0,1,2,3,4,12],q,n)
partition = PartitionCodes(n,k,s,FF)
partitioncodeword = partition.PartitionEncoding(word.copy())
partitioncodeword_initial = copy.deepcopy(partitioncodeword)


for i in range(n):
    if i < num_errors:
        elem = [i+1]
    else:
        elem = []
    partitioncodeword[i] = partition.FF.add(partitioncodeword[i],elem,"-")
end_init = time.time()

partition.PartitionDecoding(partitioncodeword)
print("Partition codeword",partitioncodeword_initial)
print("Wrong word        ",partitioncodeword)
end_all = time.time()
print("init:",(end_init-start_time),"decode:",(end_all-end_init))
