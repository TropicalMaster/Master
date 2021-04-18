from sympy import *
from math import *
import copy
from finiteField import *
from random import randint
class TZ():

    def __init__(self, n, k, s, t, FF):
        self.n = n
        self.k = k
        self.s = s
        self.FF = FF
        self.q = FF.q
        self.m = FF.m
        self.p = FF.p
        self.t = t
        # Find normal basis
        self.basis = [[i] for i in range(2*n)]

        # Find an element which has a non square norm
        square = True
        while square:
            # Find the norm of a random element
            rnd_elem = randint(0,FF.p-1)
            if rnd_elem % 2 > 0:
                square = False
        # Define gamma as an element with norm 
        self.gamma = [rnd_elem]
        # Create q-Vandermonde matrix of basis (Transposed Moore matrix)
        # Shifted because [i] = q^i in this case, but still holds
        self.M = self.qvan(self.basis,2*n)
        self.M_inv = self.FF.inv(copy.deepcopy(self.M))


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
        

    # Multiply a list of elements together in F_q^m
    def mul(self,a):
        mul_all = []
        for i in range(len(a)):
            if len(a[i]):
                if len(mul_all):
                    mul_all = [(mul_all[0] + a[i][0]) % self.p]
                else:
                    mul_all = a[i].copy()
        return mul_all
    
    # Find the sum of all elements of a in F_q^m
    def sum_all(self,a):
        tot_sum = []
        for i in range(len(a)):
            if len(a[i]):
                tot_sum = self.FF.add(tot_sum,a[i],"+")
        return tot_sum


    # Encode a word
    def encode(self,r):
        
        qn = self.q**self.n
        # Normal basis element is defined as alpha
        alpha = self.FF.elements[2*self.n]
        # alpha^(q^n)
        alpha_n = [((alpha[0])*qn) % self.p]
        # r[0]^(q^n)
        r_n = [(r[0][0]*qn) % self.p]
        # Create matrix
        A = [[alpha,alpha_n,r[0]],[alpha_n,alpha,r_n]]
        # Solve for a and b
        REF = self.FF.REF(A,True,False)
        a = REF[0][0]
        b = REF[0][1]

        '''
        #Numerator
        numer = self.FF.add([(r[0][0]*qn + alpha[0]*qn) % self.p],[(r[0][0] + alpha[0])% self.p],"-")
        denom = self.FF.add([(alpha[0]*2*qn % self.p)],[(alpha[0]*2) % self.p],"-")
        a = [(numer[0]-denom[0]) % self.p]
        numer = self.FF.add(r[0],[(a[0]+alpha[0]) % self.p],"-")
        denom = [(alpha[0]*qn) % self.p]
        b = [(numer[0] - denom[0]) % self.p]
'''
        # Redefine first element of word as "a"
        r[0] = a
        # Add additional element "b"
        b = [(b[0]+self.gamma[0]) % self.p]
        r.append(b)
        codeword = self.FF.Codeword(r,self.M)
        return codeword
    
    def TZdecoding(self,f):
        # Find beta coefficients by multiplying f with M_inv
        beta = self.FF.Codeword(f,self.M_inv)
        if beta[self.k+1:2*self.n] == [[] for i in range(self.k+1,2*self.n)]:
            print("Decoded word", beta)
            return beta
        # Berlekamp-Massey algorithm on coefficients from k+1 to 2n from beta
        L, Lambda, B = self.modifiedBM(beta[self.k+1:2*self.n])
        self.t = L
        x = "No zero"
        if L == int((2*self.n-self.k)/2):
            # Find sum of Lambda coefficients multiplied by known g coefficients
            delta = []
            for i in range(1,L+1):
                coeff = [(Lambda[0][i] + beta[2*self.n-i][0]*self.q**(self.s*i)) % self.p]
                delta = self.FF.add(delta,coeff,"+")
            # Find all c coefficients
            c2 = []
            c3 = []
            # Define lambdas list to contain lambda_i values
            lambdas = []
            # Define lambda_primes to contain lambda_i' values
            lambda_primes = []
            for i in range(1,L+1):
                # g_(k+t-i)^[i]
                g_k_qs = (beta[self.k+self.t-i][0]*self.q**(self.s*i)) % self.p
                # Calculate (Lambda_i - delta*B_i^[1])
                lambda_i = [(delta[0] + B[0][i-1]*self.q**self.s) % self.p]
                lambda_i = self.FF.add([Lambda[0][i]],lambda_i,"-")
                lambdas.append(lambda_i)
                # Calculate -B_i^[1]
                lambda_prime_i = self.FF.add([],[(B[0][i-1]*self.q**self.s) % self.p],"-")
                lambda_primes.append(lambda_prime_i)
                if i < L:
                    # Calculate c2 as sum of all lambda_i*g_k_qs
                    c2 = self.FF.add(c2,[(lambda_i[0] + g_k_qs) % self.p],"+")
                    # Calculate c2 as sum of all lambda_i_prime*g_k_qs
                    c3 = self.FF.add(c3,[(lambda_prime_i[0] + g_k_qs) % self.p],"+")
            b_3 = lambda_i
            b_4 = lambda_prime_i
            
            # Test for solution x = -lambda_t/lambda_t'
            x = self.FF.add([],[(b_3[0]-b_4[0]) % self.p],"-")
            # Calculate c2 + c3*x
            temp = self.FF.add(c2,[(c3[0] + x[0]) % self.p],"+")
            # If x is not the solution try other steps
            if beta[self.k+self.t] != temp:
                # q^(n*s)
                q_ns = self.q**(self.n*self.s)
                # q^((n+t)*s)
                q_nst = self.q**((self.n+self.t)*self.s)
                # q^(t*s)
                q_t = self.q**(self.t*self.s)
                gamma_1 = [(self.gamma[0]*q_nst - self.gamma[0]*q_t) % self.p]
                # Calculate theta_1
                theta_1 = self.FF.add([(beta[0][0]*q_ns) % self.p],beta[0],"-")
                # Calculate theta_2^[t] negated
                coeff = [(gamma_1[0] + beta[self.k][0]*q_t) % self.p]
                theta_2 = self.FF.add(coeff,[(beta[self.k][0]*q_nst) % self.p],"-")
                # Calculate a's and b's
                a_4 = [(b_4[0]*q_ns) % self.p]
                # Theta_1 could be zero
                if len(theta_1):
                    a_3 = self.FF.add([(b_3[0]*q_ns) % self.p],[(a_4[0] + theta_1[0]) % self.p],"-")
                else:
                    a_3 = [(b_3[0]*q_ns) % self.p]
                a_2 = self.FF.add([],[(c3[0]*q_ns) % self.p],"-")
                temp = self.FF.add(beta[self.k+self.t],c2,"-")
                # Theta_1 could be zero
                if len(theta_1):
                    a_1 = self.FF.add([(temp[0]*q_ns) % self.p],[(theta_1[0] + (c3[0]*q_ns)) % self.p],"+")
                else:
                    a_1 = [(temp[0]*q_ns) % self.p]
                b_1 = self.FF.add([],[(gamma_1[0]+temp[0]) % self.p],"-")
                b_2 = [(gamma_1[0]+c3[0]) % self.p]

                # zeta_1 = b_4*a_2 + a_4*b_2 + a_4*b_4*theta_2
                zeta_1 = self.sum_all([self.mul([b_4,a_2]), self.mul([a_4,b_2]), self.mul([a_4,b_4,theta_2])])

                # zeta_2 = b_4*a_1 + b_3*a_2 + a_3*b_2 + a_4*b_1 + a_3*b_4*theta_2 + a_4*b_3*theta_2
                zeta_2 = self.sum_all([
                    self.mul([b_4,a_1]), self.mul([b_3,a_2]), self.mul([a_3,b_2]), self.mul([a_4,b_1]), 
                    self.mul([a_3,b_4,theta_2]), self.mul([a_4,b_3,theta_2]), 
                ]) 
                
                # zeta_3 = b_3*a_1 + a_3*b_1 + a_3*b_3*theta_2
                zeta_3 = self.sum_all([self.mul([b_3,a_1]), self.mul([a_3,b_1]), self.mul([a_3,b_3,theta_2])])

                # The case where zeta_1 != 0
                if len(zeta_1):
                    # r = zeta_2/zeta_1
                    r = (zeta_2[0] - zeta_1[0]) % self.p
                    # s = zeta_3/zeta_1
                    s = (zeta_3[0] - zeta_1[0]) % self.p
                    # Find r^2
                    r_2 = [(r*2) % self.p]
                    # Find element 4 represented in F_q
                    elem_4 = [4 % self.q]
                    # Convert to exponent representation in F_q^2n
                    elem_4.extend([0 for i in range(2*self.n-1)])
                    elem_4 = self.FF.invext(elem_4,True)
                    # Find 4*s
                    s_2 = [(s + elem_4[0]) % self.p]
                    # Find modular inverse of 2 in F_q
                    mod_inv_2 = [mod_inverse(2,self.q)]
                    mod_inv_2.extend([0 for i in range(2*self.n-1)])
                    mod_inv_2 = self.FF.invext(mod_inv_2,True)
                    # Calculate difference = (r^2-4*s)
                    difference = self.FF.add(r_2,s_2,"-")
                    if len(difference):
                        # Find modular inverse of elem_4 in F_q
                        mod_inv = [mod_inverse((4 % self.q),self.q)]
                        # Convert to exponent representation in F_q^2n
                        mod_inv.extend([0 for i in range(2*self.n-1)])
                        mod_inv = self.FF.invext(mod_inv,True)
                        # Calculate check = difference/4
                        check = (difference[0] + mod_inv[0]) % self.p
                        # Case a)
                        # Quadratic residue
                        if (check % 2) == 0:
                            # Find the square root of difference
                            difference = [int(difference[0]/2)]
                            x = self.FF.add([],[r],"-")
                            x_1 = self.FF.add(x,difference,"+")
                            x_2 = self.FF.add(x,difference,"-")
                            x_1 = [(x_1[0] + mod_inv_2[0]) % self.p]
                            x_2 = [(x_2[0] + mod_inv_2[0]) % self.p]
                            x = [x_1,x_2]
                    # Case b)
                    elif r_2 == s_2:
                        x = self.FF.add([],[r],"-")
                        x = [[(x[0] + mod_inv_2[0]) % self.p]]
                    else:
                        print("No zero: no a) or b)")
                        return "Decoding Failure"
                else:
                    # x can be uniquely determined by x = -zeta3/zeta2
                    x = [(zeta_3[0]-zeta_2[0]) % self.p]
                    x = self.FF.add([],x,"-")
                    x = [x]
        # Solutions for g coefficients
        g_sols = []
        # Having found the zeros of P(x) calculate vector of lambdas
        if x != "No zero":
            lambda_vectors = []
            for omega in x:
                # Find the lambda vector
                lambda_vector = []
                for i in range(L):
                    # Find lambda_i + lambda_i'*omega for all i
                    lambda_i = self.FF.add(lambdas[i],[(lambda_primes[i][0] + omega[0]) % self.p],"+")
                    lambda_vector.append(lambda_i)
                # Calculate g_k by 13)
                numer = self.FF.add(temp,[(c3[0] + omega[0]) % self.p],"-")
                g_k = [((numer[0] - lambda_i[0])*self.q**((2*self.n-self.t)*self.s)) % self.p]
                g = copy.deepcopy(beta)
                # g_0 = -omega
                g[0] = self.FF.add([],omega,"-")
                g[self.k] = g_k
                # Find g_1 ... g_k-1
                for i in range(1,self.k):
                    g_i = []
                    for j in range(1,self.t+1):
                        # Subscript of g is i-j % 2n
                        k = (i-j) % (2*self.n)
                        coeff = [(lambda_vector[j-1][0] + g[k][0]*self.q**(j*self.s)) % self.p]
                        g_i = self.FF.add(g_i, coeff,"+")
                    g[i] = g_i
                # Append solutions for g coefficients
                g_sols.append(g)
                lambda_vectors.append(lambda_vector)
        else:
            g = copy.deepcopy(beta)
            for i in range(self.k+1):
                g_i = []
                for j in range(1,self.t+1):
                    # Subscript of g is i-j % 2n
                    k = (i-j) % (2*self.n)
                    coeff = [(Lambda[0][j] + g[k][0]*self.q**(j*self.s)) % self.p]
                    g_i = self.FF.add(g_i, coeff,"+")
                g[i] = g_i
            g_sols.append(g)
            lambda_vectors = [[[i] for i in Lambda[0][1:]]]

        for l in range(len(g_sols)):
            g = g_sols[l]
            lambda_vector = lambda_vectors[l]
            # Check periodicity
            check = 0
            for i in range(2*self.n):
                g_i = []
                for j in range(1,self.t+1):
                    # Subscript to g taken modulo 2n
                    k = (i-j) % (2*self.n)
                    coeff = [(lambda_vector[j-1][0] + g[k][0]*self.q**(j*self.s)) % self.p]
                    g_i = self.FF.add(g_i, coeff,"+")
                # Check periodicity
                if g_i == g[i]:
                    check += 1
                else:
                    break

            if check == 2*self.n:
                # Recover the codeword elements
                c = []
                for i in range(2*self.n):
                    elem_sum = []
                    for j in range(2*self.n):
                        coeff = [(g[j][0] + self.basis[i][0]*self.q**(j*self.s)) % self.p]
                        elem_sum = self.FF.add(elem_sum, coeff,"+")
                    c_i = self.FF.add(f[i],elem_sum,"-")
                    c.append(c_i)
                print("Decoded codeword  ",c,"\n")
                return c
            else:
                if l == 1:
                    print("Decoding failure")
                    return "Decoding failure"



