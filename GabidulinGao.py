from sympy import *
import copy
from itertools import product as produ

class GaoLike():
    # Takes length, dimension and finite field as input
    # n: length of the code
    # k: dimension
    # FF: finite field in F_q^m
    def __init__(self, n, k, FF):
        # Specify using i-th linearized lagrange basis polynomials to define r_hat
        self.lagrange = False
        self.rref = False
        self.k = k
        self.n = n
        self.FF = FF
        self.q = FF.q
        self.m = FF.m
        self.p = FF.p
        # Use normal basis when n = m
        self.basis = [[i] for i in range(n)]
        # Calculate Mg using L function
        if n == self.m:
            self.M_coeffs = [0,self.FF.add([],[0],"-")[0]]
            self.M_degs = [0,self.m]
        else:
            self.M_coeffs,self.M_degs = self.L(-1)

        # Use i-th linearized lagrange basis polynomials to define precomputed polynomials
        if self.lagrange == True and n != self.m:
            # Find all L_i(x)/L_i(g_i)
            self.M_inv = []
            for i in range(n):
                self.M_inv.append(self.L(self.basis[i][0]))
            # The generator matrix is the qvan expansion of the basis vector 
            self.G = self.FF.qvan(self.basis,k,1)
        else:
            # The generator matrix is the qvan expansion of the basis vector 
            self.G = self.FF.qvan(self.basis,n,1)
            # Find inverse of G matrix
            self.M_inv = self.FF.inv(copy.deepcopy(self.G))

        

        # Intializing "A" matrix
        A = []
        # Raising each row to the q-degree j
        for j in range((-self.n+self.k+1),self.k):
            row = []
            for i in range(self.n):
                elem = self.basis[i]
                if len(elem) != 0:
                    # Raising the element to q-degree i from j
                    # Since j can have negative values we raise q to the power (j mod m)
                    coeff = (elem[0] * self.q**(j%self.m)) % self.p
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
        self.H = self.FF.qvan(self.h,(self.n-self.k),1)


    
    # The trace function for an element in F_q^m
    def Tr(self,a):
        # Find all the conjugates of "a" element
        conjugates = []
        for i in range(self.m):
            elem = (a[0]*self.q**i) % self.p
            if elem not in conjugates:
                conjugates.append(elem)
        # Find the sum of all conjugate elements
        elem_sum = []
        for elem in conjugates:
            elem_sum = self.FF.add(elem_sum,[elem],"+")
        if len(elem_sum) == 0:
            return 0
        else:
            return self.FF.elements[elem_sum[0]][0]
    
    

    
    # Calculating the rightDiv of two polynomials (Algorithm 2.1)
    # a_coeff: alpha exponents of the polynomial dividend
    # b_coeff: alpha exponents of the polynomial divisor
    # a_exp: q-degrees of polynomial dividend
    # b_exp: q-degrees of polynomial divisor
    def rightDiv(self, a_coeff,b_coeff,a_exp,b_exp):
        # Initializing the q-degrees of both a and b polynomial
        degA = max(a_exp)
        degB = max(b_exp)
        # Index of the largest q-degree of B (to avoid doing it every loop)
        degBindex = b_exp.index(degB)
        # The corresponding coefficient of the largest q-degree of b
        b_coeff_static = b_coeff[degBindex]
        # Initializing rturn quotient coefficients and degrees
        qCoeff = []
        qPoly = []
        while degA >= degB:
            # Calculate the current degree and coefficient of the q-polynomial
            deg = degA-degB
            # Calculate the alpha exponent of the coefficient
            coeff = a_coeff[a_exp.index(degA)]-(b_coeff_static*(self.q**(deg)))
            # Checks for inverse coefficient, and replaces it with its inverse
            if coeff < 0:
                coeff = coeff % self.p
            qCoeff.append(coeff)
            qPoly.append(deg)
            # Find the composite of q and b (q(b(x)))
            composite = self.FF.composite([coeff],b_coeff,[deg],b_exp)
            # Find the difference between a and q(b) composite
            polySum = self.FF.addComp(a_coeff,composite[0],a_exp,composite[1],"-")
            # Redefine the coefficients and degrees of a as the difference calculated above
            a_coeff = polySum[0]
            a_exp = polySum[1]
            # If the new a value is the zero element then break the loop
            if len(a_exp) == 0:
                break
            else:
                degA = max(a_exp)
        # Returns the q-polynomial list and the last remainder to be calculated
        return [[qCoeff,qPoly],[a_coeff,a_exp]]



    # Calculating the leftDiv of two polynomials (Algorithm 2.1)
    # a_coeff: alpha exponents of the polynomial dividend
    # b_coeff: alpha exponents of the polynomial divisor
    # a_exp: q-degrees of polynomial dividend
    # b_exp: q-degrees of polynomial divisor
    def leftDiv(self, a_coeff,b_coeff,a_exp,b_exp):
        # Initializing the q-degrees of both a and b polynomial
        degA = max(a_exp)
        degB = max(b_exp)
        # Index of the largest q-degree of B (to avoid doing it every loop)
        degBindex = b_exp.index(degB)
        # The corresponding coefficient of the largest q-degree of b
        b_coeff_static = b_coeff[degBindex]
        # Initializing rturn quotient coefficients and degrees
        qCoeff = []
        qPoly = []
        while degA >= degB:
            # Calculate the current degree and coefficient of the q-polynomial
            deg = degA-degB
            # Calculate the alpha exponent of the coefficient
            coeff = (a_coeff[a_exp.index(degA)]-b_coeff_static)*(self.q**(-degB % self.m))
            # Coefficient has to be an element within the extension field
            coeff = coeff % self.p
            # Calculate the iterative sum of q
            qCoeff.append(coeff)
            qPoly.append(deg)
            # Find the composite of b and q (b(q(x)))
            composite = self.FF.composite(b_coeff,[coeff],b_exp,[deg])
            # Find the difference between a and b(q) composite
            polySum = self.FF.addComp(a_coeff,composite[0],a_exp,composite[1],"-")
            # Redefine the coefficients and degrees of a as the difference calculated above
            a_coeff = polySum[0]
            a_exp = polySum[1]
            
            # If the new a value is the zero element then break the loop
            if len(a_exp) == 0:
                break
            else:
                degA = max(a_exp)
        # Returns the q-polynomial list and the last remainder to be calculated
        return [[qCoeff,qPoly],[a_coeff,a_exp]]

    # Calculating the rightLEEA of two polynomials (Algorithm 2.3)

    def rightLEEA(self,a_coeff,b_coeff,a_exp,b_exp,d):
        # Initializing remainder list
        r = [[a_coeff,a_exp],[b_coeff,b_exp]]
        # First u value after zero is -q
        u = [[[],[]],[[0],[0]]]
        # First v value after x^[0] is also x^[0]
        v = [[[0],[0]],[[],[]]]
        # Initializing degR
        degR = max(r[-1][-1])
        while degR >= d:
            # Calculating rightDivision again
            rightDiv = self.rightDiv(r[-2][0],r[-1][0],r[-2][1],r[-1][1])
            # Redefining q to the new rightDiv value
            q = rightDiv[0]
            # Appending the new remainder value to r
            r.append(rightDiv[1])
            # Finds the composite of q(u-1(x))
            QcompU = self.FF.composite(q[0],u[-1][0],q[1],u[-1][1])
            # Calculates the newest iteration of u by finding the difference between u-2(x) and q(u-1(x))
            newU = self.FF.addComp(u[-2][0],QcompU[0],u[-2][1],QcompU[1],"-")
            u.append(newU)
            # Finds the composite of q(v-1(x))
            QcompV = self.FF.composite(q[0],v[-1][0],q[1],v[-1][1])
            # Calculates the newest iteration of v by finding the difference between v-2(x) and q(v-1(x))
            newV = self.FF.addComp(v[-2][0],QcompV[0],v[-2][1],QcompV[1],"-")
            v.append(newV)

            # Redefines degR
            if len(r[-1][1])>0:
                degR = max(r[-1][1])
            else:
                break
        return [r[-1],u[-1],v[-1]]
    

    # Gao decoding algorithm
    def gaoDecoding(self,r):
        if self.lagrange and not self.rref:
            r_hat = self.FF.Codeword(r,self.M_inv)
        else:
            # Create a matrix where every row is the ith basis element to all
            # q-degrees up to n
            r_matrix = []
            for i in range(self.n):
                row = []
                for j in range(self.n):
                    row.append([(self.basis[i][0]*self.q**j) % self.p])
                row.append(r[i])
                r_matrix.append(row)
            # Solve the r_matrix resulting in r_hat
            r_hat = self.FF.REF(r_matrix,True,False)[0]
            # Convert to simple polynomial without zero
            r_coeffs = []
            r_degs = []
            for i in range(len(r_hat)):
                if len(r_hat[i]):
                    r_coeffs.append(r_hat[i][0])
                    r_degs.append(i)
            # Use Mg as defined in the constructor
            M_coeffs,M_degs = self.M_coeffs,self.M_degs

        # Define stopping degree equal to lower(n+k/2)
        d = int((self.n+self.k)/2)
        # Using rightLEEA on M_G(x) and ^r(x) 
        LEEA = self.rightLEEA(M_coeffs,r_coeffs,M_degs,r_degs,d)

        # LeftDiv should give us the decoded word or the error span polynomial
        leftDiv = self.leftDiv(LEEA[0][0],LEEA[1][0],LEEA[0][1],LEEA[1][1])

        if len(leftDiv[1][0]) == 0:
            leftDiv[0][0].reverse()
            print("Evaluation poly:",leftDiv[0][0])
            return leftDiv[0][0]
        else:
            print("Decoding failure")
            return "Decoding failure"



    def L(self, g_i):
        if g_i > self.n-1:
            return "Out of bounds"
        # Remove basis element that should not be involved
        basis_elems = [self.basis[i] for i in range(self.n) if i != g_i]
        basis_length = len(basis_elems)
        # Convert to sympy matrix of vector representation of basis elements
        basis = Matrix(self.FF.ext(basis_elems,False))
        # Create all combinations of q coefficients to be multiplied by the basis
        combinations = list(produ([i for i in range(self.q)],repeat=len(basis_elems)))
        combinations = [list(i) for i in combinations]
        # Initialize all elements list
        a = []
        # Create all elements which are a combination of the basis
        for combination in combinations:
            # Convert to sympy vector
            combination = Matrix(combination)
            combination = combination.T
            # Multiply a combination with the basis matrix
            element = combination*basis.copy()
            # Convert the sympy vector to a python list reducing modulo q
            element_ext = []
            for coeff in element:
                element_ext.append(coeff%self.q)
            # Convert to exponent representation
            element_ext = self.FF.invext(element_ext,True)
            # Add each elements to a list 
            a.append(element_ext)

        # Initialize polynomial list
        poly = []
        # Multiply all terms (x-a) together 
        for j in range(self.q**(basis_length)):
            # Create the first iteration
            if len(poly) == 0:
                poly = [self.FF.add([],a[j],"-"),[0]]
            else:
                # Insert zero (multiply expression by x)
                poly.insert(0,[])
                # Find negated a element
                neg_a = self.FF.add([],a[j],"-")[0]
                # Add negated current coefficient to all elements
                # this is the same as multiplying every element with the current
                # negated a element
                for k in range(len(poly)-1):
                    # Multiply each coefficient from k+1 to len(poly) with neg_a
                    # Add product to original polynomial 
                    if k != 0:
                        if len(poly[k+1]):
                            product = [(poly[k+1][0] + neg_a) % self.p]
                            poly[k] = self.FF.add(poly[k],product,"+")
        print(poly)
        # Convert to simple polynomial
        coeffs = []
        degs = []
        for i in range(basis_length+1):
            q_deg = self.q**i
            if len(poly[q_deg]):
                coeffs.append(poly[q_deg][0])
                degs.append(i)
        # In the case where g_i is a defined basis element we return L_i(x)/L_i(g_i)
        if g_i >= 0:
            # Find L_i(g_i)
            L_gi = []
            for i in range(len(coeffs)):
                coeff = [(coeffs[i] + g_i*self.q**degs[i]) % self.p]
                L_gi = self.FF.add(L_gi, coeff,"+")
            # Find inverse
            L_gi = (-L_gi[0] % self.p)
            # Divide polynomial by L_i(g_i)
            for i in range(len(coeffs)):
                coeffs[i] = [(coeffs[i] + L_gi) % self.p]
            return coeffs
        else:
            return coeffs,degs



            
            
            

            
