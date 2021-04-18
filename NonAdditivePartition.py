from sympy import *
from math import *
import copy
from finiteField import *
from random import randint
class PartitionCodes():

    def __init__(self, n, k, s, FF):
        
        self.n = n
        self.k = k
        self.s = s
        self.FF = FF
        self.q = FF.q
        self.p = FF.p
        # Define I to be any subspace of F_q expressed in exponent representation
        # of F_q^m
        element = [0]
        element.extend([0 for i in range(n-1)])
        element = self.FF.invext(element,True)
        self.I = [element]
        self.basis = [[i] for i in range(n)]

        # Create q-Vandermonde matrix of basis (Transposed Moore matrix)
        # Shifted because [i] = q^i in this case, but still holds
        self.M = self.qvan(self.basis,self.n)
        # Matrix inversion
        self.M_inv = self.FF.inv(copy.deepcopy(self.M))


    # Create a q-vandermonde matrix of input elements
    # a: vector containing alpha element degrees of length n
    # s: the number of rows in sxn matrix, the range of q-degrees to raise the a elements to
    def qvan(self,a,s):
        # Initializing the matrix
        matrix = []
        # Going through the range s
        for i in range(s):
            row = []
            # Going through each element in the "a" vector   
            for elem in a:
                if len(elem) != 0:
                    elem = elem[0]
                    # Raising the existing element degrees to q^i (* because "a" elements are given in degrees)
                    elem = (elem * (self.q**(self.s*i))) % self.p
                    row.append([elem])
                else:
                    row.append([])
            # Append the row to the matrix
            matrix.append(row)
        return matrix


    # Calculate the norm
    def norm(self,a):
        deg = int((self.q**(2*self.n) - 1)/(self.q-1)) % self.p
        if len(a):
            norm = [(a[0]*deg) % self.p]
            return norm
        else:
            return []

    def modifiedBM(self,g):
        r = 0
        L = 0
        # Initializing polynomials with q-degree 0
        Lambda = [[0],[0]]
        B = [[0],[0]]
        while r <= len(g)-1:
            # Find delta_r as g_r + sum
            delta_r = g[r].copy()
            # Sum of Lambda coefficients multiplied with g coefficients
            for i in range(L):
                # Check whether g[r-i] is zero
                if len(g[r-(i+1)]):
                    coeff = [(Lambda[0][i+1] + g[r-(i+1)][0]*self.q**(self.s*(i+1))) % self.p]
                    delta_r = self.FF.add(delta_r,coeff,"+")
            # Condition
            if len(delta_r) == 0:
                B = self.FF.composite([0],B[0],[self.s],B[1])
            else:
                # Copy current Lambda
                Lambda_temp = copy.deepcopy(Lambda)
                # Find lambda - delta*x^[1]*B
                composite = self.FF.composite(delta_r,B[0],[self.s],B[1])
                Lambda = self.FF.addComp(Lambda[0],composite[0],Lambda[1],composite[1],"-")
                if 2*L > r:
                    B = self.FF.composite([0],B[0],[self.s],B[1])
                else:
                    # Multiply Lambda_temp with the inverse of delta_r
                    for i in range(len(Lambda_temp[0])):
                        Lambda_temp[0][i] = (Lambda_temp[0][i] - delta_r[0]) % self.p
                    # Define new B polynomial as Lambda_temp
                    B = copy.deepcopy(Lambda_temp)
                    # Increase L
                    L = r + 1 - L
            r = r + 1
            
        # Negate coefficients because of relation
        for i in range(len(Lambda[0])):
            Lambda[0][i] = self.FF.add([],[Lambda[0][i]],"-")[0]
        return [L, Lambda, B]
    
    def partialBM(self,L,Lambda,B,r,g):
        # Find delta_r as g_r + sum
        delta_r = [g[r][0]]
        # Sum of Lambda coefficients multiplied with g coefficients
        for i in range(L):
            # Check whether g[r-i] is zero
            if len(g[r-(i+1)]):
                coeff = [(Lambda[0][i+1] + g[r-(i+1)][0]*self.q**(self.s*(i+1))) % self.p]
                delta_r = self.FF.add(delta_r,coeff,"+")
        # Condition
        if len(delta_r) == 0:
            B = self.FF.composite([0],B[0],[self.s],B[1])
        else:
            # Copy current Lambda
            Lambda_temp = copy.deepcopy(Lambda)
            # Find lambda - delta*x^[1]*B
            composite = self.FF.composite(delta_r,B[0],[self.s],B[1])
            Lambda = self.FF.addComp(Lambda[0],composite[0],Lambda[1],composite[1],"-")
        return L,Lambda,B

    def PartitionEncoding(self,r):
        # Change received word according to condition
        if self.norm(r[0]) not in self.I:
            if (self.k+1) % 2 == 0:
                r.append(r[0].copy())
            else:
                r.append(self.FF.add([],r[0],"-"))
            r[0] = []
        codeword = self.FF.Codeword(r,self.M)

        return codeword
    
    def PartitionDecoding(self,f):
        
        beta = self.FF.Codeword(f,self.M_inv)
        if beta[self.k+1:self.n] == [[] for i in range(self.k+1,self.n)]:
            print("Decoded word", beta)
            return beta
        # Berlekamp-Massey algorithm on coefficients from k+1 to 2n from beta
        t0, Lambda, B = self.modifiedBM(beta[self.k+1:self.n])
        g = copy.deepcopy(beta)
        g_sols = []
        lambda_vectors = []
        if t0 == int((self.n-self.k)/2):
            t,Lambda,B = self.modifiedBM(g[self.k:self.n])
            # Find g_0
            g_0 = []
            for i in range(1,t+1):
                coeff = [(Lambda[0][i] + g[self.n-i][0]*self.q**(self.s*i)) % self.p]
                g_0 = self.FF.add(g_0,coeff,"+")
            # Check if norm(g[0]-g_0) in I
            if self.norm(self.FF.add(g[0],g_0,"-")) in self.I:
                # Add to solution set
                lambda_vectors.append(Lambda[0].copy())
                g_sols.append(copy.deepcopy(beta))
            # Create a copy of g, and add g_0 at the end to find BM(g_k+1,g_n+1)
            g_temp = copy.deepcopy(g)
            g_temp.append(g[0].copy())
            t,Lambda,B = self.modifiedBM(g_temp[self.k+1:self.n+1])
            # Find g_k
            g_k_temp = []
            for i in range(1,t):
                coeff = [(Lambda[0][i] + g[self.k+t-i][0]*self.q**(self.s*i)) % self.p]
                g_k_temp = self.FF.add(g_k_temp,coeff,"+")
            g_kt = [(g[self.k+t][0]*self.q**(self.s*(self.n+t))) % self.p]
            g_k = self.FF.add(g_kt,g_k_temp,"-")
            g_k = [(g_k[0]-Lambda[0][t]) % self.p]
            norm = self.norm(self.FF.add(g[self.k],g_k,"-"))
            # Negate norm element if nks is odd
            if (self.n*self.k*self.s) % 2 != 0:
                norm = self.FF.add([],norm,"-")
            # Check norm(g[k]-g_k) not in I
            if norm not in self.I:
                # Negate according to encoding
                if (self.k+1) % 2 != 0:
                    g_k = self.FF.add([],g_k,"-")
                # Add to solution set
                g[self.k] = g_k
                lambda_vectors.append(Lambda[0].copy())
                g_sols.append(g)
        else:
            lambda_vectors.append(Lambda[0].copy())
            g_sols.append(g)
        for l in range(len(g_sols)):
            g = g_sols[l]
            lambda_vector = lambda_vectors[l]
            # Check periodicity
            check = 0
            checklimit = 3
            for i in range(self.n+checklimit):
                g_i = []
                for j in range(1,t0+1):
                    # Subscript of g is i-j % 2n
                    k = (i-j) % (self.n)
                    if len(g[k]):
                        coeff = [(lambda_vector[j] + g[k][0]*self.q**(j*self.s)) % self.p]
                        g_i = self.FF.add(g_i, coeff,"+")
                if i < self.n:
                    g[i] = g_i
                else:
                    if g_i == g[i % self.n]:
                        check += 1 
            if check == checklimit:
                # Recover the codeword elements
                c = []
                for i in range(self.n):
                    e_i = []
                    for j in range(self.n):
                        if len(g[j]):
                            coeff = [(g[j][0] + self.basis[i][0]*self.q**(j*self.s)) % self.p]
                            e_i = self.FF.add(e_i, coeff,"+")
                    c_i = self.FF.add(f[i],e_i,"-")
                    c.append(c_i)
                print("Decoded codeword  ",c,"\n")
                return c
            else:
                if l == 2:
                    print("Decoding Failure")
                    return "Decoding Failure"

        

                

                
