from sympy import *
import copy
from itertools import product as produ
import math
class gabidulinCode():
    # Takes length, dimension and finite field as input
    # n: length of the code
    # k: dimension
    # FF: finite field in F_q^m
    def __init__(self, n, k, s, FF):
        self.k = k
        self.n = n
        self.s = s
        self.FF = FF
        self.q = FF.q
        self.m = FF.m
        self.p = FF.p
        

        basis = [[i] for i in range(n)]

        # The generator matrix is the qvan expansion of the basis vector containing normal basis
        self.G = self.FF.qvan(basis,k,s)
        # Intializing "A" matrix
        A = []
        # Raising each row to the q-degree j
        for j in range((-self.n+self.k+1),self.k):
            row = []
            for i in range(self.n):
                elem = basis[i]
                if len(elem) != 0:
                    # Raising the element to q-degree i from j
                    # Since j can have negative values we raise q to the power (j mod m)
                    coeff = (elem[0] * self.q**((s*j)%self.m)) % self.p
                    row.append([coeff])
                else:
                    # null case
                    row.append([])
            row.append([])
            # Add each row to the matrix
            # Have to use deepcopy because pointers in python
            A.append(copy.deepcopy(row))
        # Row reduction gives us a vector of h values
        self.h = self.FF.REF(A,True,False)[0]
        # Parity check matrix is the q-vandermonde expansion of h
        self.H = self.FF.qvan(self.h,(self.n-self.k),s)
        # Transpose parity check matrix
        self.H = self.FF.transpose(self.H)

    
    
    # Standard decoding algorithm for Gabidulin codes
    def decodingGabidulin(self,r):
        # Initialize syndrome vector
        s = self.FF.Codeword(r,self.H)
        print("S",s)
        # When the syndrome vector is empty, the recieved word has no errors
        if s == [[] for i in range(len(s))]:
            print("Decoded word:", r)
            return r
        # Initialize t = lower((n-k)/2) + 1
        t = int((self.n-self.k)/2) + 1
        # Initialize S_rank as 0
        S_rank = 0
        while S_rank < t:
            t -= 1
            # Initialize S matrix
            S = []
            for i in range(t,self.n-self.k):
                row = []
                for j in range(t+1):
                    if len(s[i-j]):
                        coeff = [(s[i-j][0] * self.q**(self.s*j)) % self.p]
                    else:
                        coeff = []
                    row.append(coeff.copy())
                print(row)
                row.append([])
                S.append(row)
            # keyCoeffs are the coefficients of the key equation which are the coefficients of
            # the minimal subspace polynomial
            REF = self.FF.REF(S,True,False)
            keyCoeffs = REF[0]
            S_rank = REF[1]
        # Key equations with basis as input
        keybasis = []
        # i variable is the exponent of basis elements
        for i in range(self.m):
            # keysum is the iterative sum of key equation coefficients multiplied by 
            # the basis element to q-degree j
            keysum = []
            for j in range(len(keyCoeffs)):
                # Check for null case
                if len(keyCoeffs[j]):
                    # key equation coefficient times basis element to q-degree j
                    coeff = (keyCoeffs[j][0] + i*(self.q**(self.s*j))) % self.p
                else:
                    coeff = (i*(self.q**(self.s*j))) % self.p
                keysum = self.FF.add(keysum,[coeff],"+")
            # Appends the different sums to keybasis list
            keybasis.append(keysum)
       
        # Map to vector representation matrix
        keybasis_ext = self.FF.ext(keybasis,False)
        # Transpose the matrix to get basis elements as column vectors
        keybasis_ext = self.FF.transpose(keybasis_ext)
        # Row reduced echelon form of the given ext matrix
        matrix = self.FF.RREF(keybasis_ext,False,True)[0]
        # Finding the nullspace of the RREForm matrix (kernel)
        solution = Matrix(matrix).nullspace()
        # The kernel will contain all linearly independent vectors 
        kernel = []
        # Convert the matrix to vectors in python form
        for vector in solution:
            newvector = []
            for element in vector.tolist():
                newvector.append(element[0]%self.q)
            kernel.append(newvector)
        # The kernel is the set of all vectors in the domain mapping to the zero vector of matrix
        # The kernel is equal to ext_B(a)
        # This means we need to convert the kernel to elements in F_q^m from a matrix in F_q
        # In other words find ext^-1(A)
        a = self.FF.invext(kernel,False)
        print(a)
        if len(a) == t:
            # a matrix
            a_matrix = []
            for i in range(0,-(self.n-self.k),-1):
                a_copy = copy.deepcopy(a)
                print(a_copy)
                for j in range(len(a_copy)):
                    if len(a_copy[j]):
                        # Raising a element to q-degree i
                        a_copy[j][0] = (a_copy[j][0]*self.q**((self.s*i)%self.m)) % self.p
                # Appending the syndrome element
                syndrome = s[abs(i)]
                if len(syndrome):
                    syndrome = [(syndrome[0]*self.q**((self.s*i)%self.m)) % self.p]
                a_copy.append(syndrome)
                # Add row to a_matrix
                a_matrix.append(a_copy)
            # Solving the a_matrix
            d = self.FF.REF(a_matrix,True,False)[0]
            # Convert d to vector representation
            d = self.FF.ext(d,False)
            # Finding te B matrix
            # Convert the h elements exponent representation to vector representation
            h = copy.deepcopy(self.h)
            # Convert the h elements to ground field
            hmatrix = self.FF.ext(h,False)
            # Solving for coefficients of B using hmatrix = dl
            B = Matrix()
            # Go through each dl element in d vector
            for dl in d:
                # Creating a temporary copy of hmatrix (copy cause pointers in python)
                matrix = copy.deepcopy(hmatrix)
                # Appending the dl vector to solve for
                matrix.append(dl)
                # Transpose to get the hi elements as columns
                matrix = self.FF.transpose(matrix)
                # Solve matrix
                solution = self.FF.RREF(matrix,True,True)
                # Add Brow to B matrix
                B = B.row_insert(len(B),Matrix([solution]))
            
            # Convert a vector to sympy format
            a = Matrix(kernel)
            # Product vector (should be e=a*B)
            # Transpose B to be able to multiply with a
            e = B.T*a
            # Convert r in exponent representation to r_ext vector representation matrix
            r = self.FF.ext(r,False)
            r_ext = Matrix(r)
            print("rext",r_ext)
            # Find the difference r-e
            c = r_ext - e
            print(c)
            # Apply the lambda function to reduce mod q
            c = c.applyfunc(lambda x: x % self.q)
            print("Decoded word",self.FF.invext(c.tolist(),False))
            return self.FF.invext(c.tolist(),False)
        else: 
            print("Decoding failure")
            return "Decoding failure"



            
            
            

            


