import copy
from itertools import combinations as combos
from sympy import *
import random

class LRPC():

    def __init__(self, n, k, d, r, FF):
        self.rref = True
        self.r = r
        self.n = n
        self.k = k
        self.FF = FF
        self.q = FF.q
        self.m = FF.m
        self.p = FF.p
        
        # Initialize matrix
        matrix = Matrix()
        F_basis = []
        # Will continue until the number of independent elements are equal to d
        while len(F_basis) < d:
            # Create a copy of the matrix 
            matrix_temp = matrix.copy()
            # Create a random element in F_q^m
            rndElem = rndElem = random.randint(0,self.p-1)
            # Insert the random element in extension field format as a column
            matrix_temp = matrix_temp.col_insert(len(matrix_temp),Matrix(self.FF.ext([rndElem],True)))
            # If the rank has increased add the element to the F_basis list
            if matrix_temp.rank() > len(F_basis):
                
                F_basis.append(rndElem)
                # Keep the matrix until the next iteration
                matrix = matrix_temp.copy()
        # Define the F basis as a class variable to be used later
        self.F_basis = F_basis
        # Transpose basis matrix for future matrix multiplication
        matrix = matrix.T
        # Generate a parity check matrix with random elements from all_elements
        # This matrix is completely randomized
        self.H = []
        self.H_d = []
        self.G = []
        # p is a static integer representing the dimension of parity-check columns
        p = self.n-(self.n-self.k)
        for i in range(p):
            col_d = []
            col = []
            for j in range(self.n-self.k):
                # Generate random element
                rndCoeffs,rndElem = self.rndElem(matrix,d)
                # Add coefficients of random element to column
                col_d.append(rndCoeffs)
                # Add random element to column
                col.append(rndElem)
            self.H_d.append(col_d)
            self.H.append(col)
            self.G.append(copy.deepcopy(col))
        # Create an identity matrix using rnd_diag element (rnd_diag representing 1)
        for i in range(self.n-self.k):
            rnd_diag = []
            while rnd_diag == []:
                # Create random element to be used as diagonal in H matrix 
                rnd_diag_d, rnd_diag = self.rndElem(matrix,d)
            empty_col = [[] for j in range(self.n-self.k)]
            empty_col[i] = rnd_diag.copy()
            empty_col_d = [[0 for i in range(d)] for j in range(self.n-self.k)]
            empty_col_d[i] = rnd_diag_d.copy()
            self.H_d.append(empty_col_d)
            self.H.append(empty_col)
        # Transpose to the correct format
        self.H_d = self.FF.transpose(self.H_d)
        


        # Start making the generator matrix
        # The n-k last columns of G are the negated transposed first k columns of H
        for i in range(len(self.G)):
            for j in range(len(self.G[0])):
                # Multiply by the inverse of diagonal element in H
                if len(self.G[i][j]):
                    coeff = [(self.G[i][j][0] - self.H[p+j][j][0]) % self.p]
                else:
                    coeff = []
                self.G[i][j] = self.FF.add([],coeff,"-")
        # Fill in the first k columns with the identity matrix of size k x k
        for i in range(self.k):
            for j in range(self.k):
                if i == self.k-j-1:
                    self.G[i].insert(0,[0])
                else:
                    self.G[i].insert(0,[])

        
        # Initialize empty d*r*(n-k) x nr matrix
        self.H_Fq = [[0]*self.n*r for i in range(d*r*(self.n-self.k))]
        # Iterate through rows of H keeping row iterator variable row_i
        for row_i in range(len(self.H_d)):
            # Iterate through the row keeping column iterator variable col_i
            for col_i in range(len(self.H_d[0])):
                F_element = self.H_d[row_i][col_i].copy()
                for i in range(d):
                    # Skip d rows for each new H element, skip r rows for each coefficient
                    # of F basis representation of F_element
                    row_index = i*r + row_i*r*d 
                    for j in range(r):
                        # Skip r columns for each new element while itertively moving one
                        # step along columns
                        col_index = col_i*r + j
                        self.H_Fq[row_index+j][col_index] = F_element[i]

        # Find rd x rd invertible submatrix
        if not self.rref:
            tot_len = r*d*(self.n-self.k)
            combinations = combos([i for i in range(tot_len)],n*r)
            for combination in combinations:
                H_Fq_i = []
                for row_i in combination:
                    H_Fq_i.append(self.H_Fq[row_i])
                matrix = Matrix(H_Fq_i)
                if matrix.rank() == n*r:
                    self.H_inv = self.FF.inv_q(H_Fq_i)
                    self.H_inv = Matrix(self.H_inv)
                    self.s_is = combination
                    break
        
    # Generate random element given a basis
    # basis: basis given as a matrix in F_q
    # d: number of basis elements
    def rndElem(self,basis,d):
        # Create d dimension list of random elements from q
        rndCoeffs = [random.randint(0,self.q-1) for i in range(d)]
        # Multiply with basis matrix modulo q
        rndElem_ext = Matrix(rndCoeffs).T*basis
        rndElem_ext = [i%self.q for i in rndElem_ext]
        # Convert to exponent representation
        rndElem = self.FF.invext(rndElem_ext,True)
        return rndCoeffs,rndElem
    
    def gauss(self,A):
        # Number of rows
        r = len(A)
        # Number of columns
        c = len(A[0])
        stop = False
        for i in range(int((r+(c-1))/c)):
            for j in range(c):
                if i >= r:
                    stop = True
                    break
                for k in range(i+1,r):
                    if ((A[i][j] + A[k][j]) % self.q) != 0:
                        for l in range(c):
                            A[i][l] = (A[i][l] + A[k][l]) % self.q
                for k in range(r):
                    if k != i:
                        if A[k][j] != 0:
                            for l in range(c):
                                A[k][l] = (A[k][l] + A[i][l]) % self.q
        return A
        
    
    def zassenhaus(self,A,B):
        Al = len(A)
        Bl = len(B)
        # Append every element in A to A (AA matrix)
        for i in range(Al):
            A[i].extend(A[i])
        # Initialize a zero vector of size m
        empty = [0 for i in range(self.m)]
        # Append the empty vectors to each B row
        for i in range(Bl):
            B[i].extend(empty)
            A.append(B[i])
        # Find REF of A
        gauss,t = self.FF.RREF(A,False,False)
        # Define pivot variable
        pivot = 0
        # Initialize empty basis list
        basis = []
        # True when there are no more pivots
        stop = False
        # Find all basis elements from the matrix
        for i in range(Al+Bl):
            if i == self.m:
                stop = True
            # Find all pivots
            if not stop:
                if gauss[i][pivot] == 0:
                    for j in range(i,self.m):
                        if gauss[i][j] != 0:
                            pivot = j+1
                            break
                        # Stop searching for pivots, start basis finding
                        if j == self.m-1:
                            stop = True
                            # First basis
                            candidate = gauss[i][self.m:2*self.m]
                            if candidate != empty:
                                basis.append(candidate)
                            else:
                                return []
                else:
                    pivot += 1
            else:
                # Check if zero
                candidate = gauss[i][self.m:2*self.m]
                if candidate != empty:
                    basis.append(gauss[i][self.m:2*self.m])
                else:
                    break
        # Return basis in vector representation
        return basis


    
    # Decoding algorithm for LRPC codes
    def decodingLRPC(self, y):
        d = len(self.F_basis)
        # Initialize syndrome vector
        s = self.FF.Codeword(y,self.H)
        # When the syndrome vector is empty, the recieved word has no errors
        if s == [[] for i in range(len(s))]:
            print("Decoded LRPC:", y)
            return y
        # Convert the unique syndrome elements to a sympy matrix
        s_ext = self.FF.ext(s,False)
        '''
        s_basis,t = self.FF.RREF(copy.deepcopy(s_ext),False,True)
        
        print("SBASIS_ext",s_basis)
        s_basis = self.FF.invext(s_basis,False)
        s_basisnew = []
        for i in range(len(s_basis)):
            if len(s_basis[i]) == 0:
                break
            else:
                s_basisnew.append(s_basis[i])
        print("SBASIS",s_basisnew)
        '''
        s_space = copy.deepcopy(s)
        # Multiply the syndrome space with the inverse of f_i
        all_intersect = []
        for i in range(d):
            # Find the inverse of f_i
            inv_f_i = (-self.F_basis[i]) % self.p
            # Create an empty set S_i
            S_i = []
            # Multiply inv_f_i with every element in s and add to S_i
            for element in s_space:
                if len(element):
                    S_i.append([(element[0] + inv_f_i) % self.p])
            # Find intersection of all S_i
            if i == 0:
                all_intersect = copy.deepcopy(self.FF.ext(S_i,False))
            else:
                all_intersect = self.zassenhaus(self.FF.ext(S_i,False),all_intersect)
            if len(all_intersect) == self.r or len(all_intersect) == 0:
                break
        if len(all_intersect) != self.r:
            print("Error space intersection not sufficient",all_intersect)
            return "Decoding failure"
        e_basis_ext = all_intersect
        e_basis = self.FF.invext(e_basis_ext,False)
        # Create EF basis
        EF_basis = []
        for i in range(d):
            for j in range(self.r):
                EF_basis.append([(e_basis[j][0]+self.F_basis[i])%self.p])
        # Convert all syndrome elements to our new basis EF
        # First create a matrix from the basis
        # Convert EF basis to base field
        EF_ext = self.FF.ext(EF_basis,False)

        # Find syndrome elements in EF representation
        s_EF = []
        for s_i in s_ext:
            # Copy the EF_ext elements as a matrix
            matrix = copy.deepcopy(EF_ext)
            # Append s_i element
            matrix.append(s_i)
            # Transpose the matrix
            matrix = self.FF.transpose(matrix)
            # Solve the system of equations for to give syndrome elements represented in EF
            try:
                solution = self.FF.RREF(matrix,True,False)
            except:
                print("Syndrome element not in EF space")
                return "Decoding failure"
            # Extend the whole syndrome EF element list
            s_EF.extend(solution)

        # Use inverted matrix if rref is not true
        if not self.rref:
            new_s = []
            for i in self.s_is:
                new_s.append(s_EF[i])
            mul = self.H_inv*Matrix(new_s)
            mul = mul.applyfunc(lambda x: x%self.q)
            e_coeff = []
            for coeff in mul:
                e_coeff.append(coeff)
        else:
        # Use rref if the inverted matrix is not computed
            for i in range(len(s_EF)):
                self.H_Fq[i].append(s_EF[i])
            e_coeff = self.FF.RREF(self.H_Fq,True,True)


        if e_coeff == "No solution":
            print("Decoding failure: No solution for error coefficients")
            return "Decoding failure"

        # e list represents the final error vector
        e = []
        # Initialize a zero element in the E basis
        zero_E = [0 for i in range(self.r)]
        for i in range(0,len(e_coeff),self.r):
            # Create a list e_i_coeff representing the ith e coefficient
            e_i_coeff = e_coeff[i:i+self.r]
            # In any case that e_i is not the zero element
            if e_i_coeff != zero_E:
                # Find product of e coefficient with e basis
                e_i_coeff = Matrix(e_i_coeff).T * Matrix(e_basis_ext)
                # Convert to python format
                e_i = []
                for coeff in e_i_coeff:
                    e_i.append(coeff % self.q)
                # Convert to exponent representation
                e.append(self.FF.invext(e_i,True))
            else:
                e.append([])

        decoded = []
        # Find the difference between y and e for the final decoded word
        for i in range(len(y)):
            element = self.FF.add(y[i],e[i],"-")
            decoded.append(element)
        print("Decoded LRPC", decoded)