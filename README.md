# Master
Master thesis code

This repository contains the implementation of MRD and LRPC codes including some of their polynomial time decoding algorithms.
Three different families are included for MRD codes: Gabidulin codes, TZ codes and non-additive partition codes. Gabidulin codes have two different decoding approaches,
leading to two separate classes.

Additionally, there are different TESTING files for each of the codes. The parameters can be changed to test each implementation.
The parameters are described below:

q: prime number of F_q
m: extension degree of F_q^m
n: length of output codewords
k: length of input message 
s: automorphism for generalized codes (GG, TZ, non-additive partition codes)
d: small weight of Parity-chekc matrix (LRPC)
num_error: number of errors imposed on the codeword


To run the different codes, create an empty "__init__.py" file in the folder containing the other .py files. GitHub does not allow uploading empty files, so this has to be done
before being able to run the codes.
