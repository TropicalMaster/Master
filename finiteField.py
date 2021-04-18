
import copy
from itertools import product as produ
from sympy import *
import random
import sys
class finiteField():
    primpolys = [
            [1,0],[2,1,0],[3,1,0],[4,1,0],[5,2,0],[6,1,0],[7,1,0],[8,4,3,2,0],[9,4,0],[10,3,0],[11,2,0],[12,6,4,1,0],[13,4,3,1,0],
            [14,5,3,1,0],[15,1,0],[16,5,3,2,0],[17,3,0],[18,5,2,1,0],[19,5,2,1,0],[20,3,0],[21,2,0],[22,1,0],[23,5,0],[24,4,3,1,0],
            [25,3,0],[26,6,2,1,0],[27,5,2,1,0],[28,3,0],[29,2,0],[30,6,4,1,0],[31,3,0],[32,7,5,3,2,1,0],[33,6,4,1,0],[34,7,6,5,2,1,0],
            [35,2,0],[36,6,5,4,2,1,0],[37,5,4,3,2,1,0],[38,6,5,1,0],[39,4,0],[40,5,4,3,0],[41,3,0],[42,5,4,3,2,1,0],[43,6,4,3,0],
            [44,6,5,2,0],[45,4,3,1,0],[46,8,5,3,2,1,0],[47,5,0],[48,7,5,4,2,1,0],[49,6,5,4,0],[50,4,3,2,0],[51,6,3,1,0],[52,3,0],
            [53,6,2,1,0],[54,6,5,4,3,2,0],[55,6,2,1,0],[56,7,4,2,0],[57,5,3,2,0],[58,6,5,1,0],[59,6,5,4,3,1,0],[60,1,0],
            [61,5,2,1,0],[62,6,5,3,0],[63,1,0],[64,4,3,1,0],[65,4,3,1,0],[66,8,6,5,3,2,0],[67,5,2,1,0],[68,7,5,1,0],[69,6,5,2,0],
            [70,5,3,1,0],[71,5,3,1,0],[72,6,4,3,2,1,0],[73,4,3,2,0],[74,7,4,3,0],[75,6,3,1,0],[76,5,4,2,0],[77,6,5,2,0],[78,7,2,1,0],
            [79,4,3,2,0],[80,7,5,3,2,1,0],[81,4,0],[82,8,7,6,4,1,0],[83,7,4,2,0],[84,8,7,5,3,1,0],[85,8,2,1,0],[86,6,5,2,0],
            [87,7,5,1,0],[88,8,5,4,3,1,0],[89,6,5,3,0],[90,5,3,2,0],[91,7,6,5,3,2,0],[92,6,5,2,0],[93,2,0],[94,6,5,1,0],
            [95,6,5,4,2,1,0],[96,7,6,4,3,2,0],[97,6,0],[98,7,4,3,2,1,0],[99,7,5,4,0],[100,8,7,2,0]
            ]
    # Constructor creates all elements of the prime field
    # coefficients: the coefficients of the primitive polynomial in increasing order (exponents icreasing order)
    # degrees: the degrees of the primitive polynomial corresponding to the coefficients in order
    # Assumes monic primitive polynomial
    # q: the prime of the prime field
    # m: the degree of the extension field
    def __init__(self, coefficients, degrees, q, m):
        self.q = q
        self.m = m
        # p is number of elements (excluding the zero element)
        self.p = self.q**self.m - 1
        
        # For q = 2, we have a hardcoded list of primitive polynomials
        if q == 2 and m <= 100:
            # Find primitive polynomial corresponding to "m" degree
            degrees = self.primpolys[m-1]
            # Sort for right format
            degrees.sort()
            # Coefficients in mod 2 are always 1
            coefficients = [1]*len(degrees)
        self.generateElements(coefficients,degrees)
    
    def generateElements(self, coefficients, degrees):
        # Finding the highest degree of the primitive polynomial 
        highestDeg = degrees.pop()
        # Removing the corresponding coefficient of the highest degree element
        coefficients.pop()
        # Converting coefficients to the primitive root
        coefficients = [-i%self.q for i in coefficients]
        # Initializing a list containing "m" number of zero elements
        primeCoeffs = [0 for i in range(highestDeg)]
        # Fills the primeCoeffs list with alpha exponents corresponding to the primitive element
        # The coefficients are placed at the index corresponding to the q-degrees (0 otherwise)
        for i in range(len(degrees)):
            primeCoeffs[degrees[i]] = coefficients[i]
        # Create a copy of primeCoeffs
        currentCoeffs = primeCoeffs[:]
        # Initializing an elements list to contain all finite field elements
        self.elements = []
        self.elements_ext = {}
        # Creates all elements of alpha exponent up to q and adds them to the coefficient list
        for i in range(highestDeg):
            zeroList = [0 for i in range(highestDeg)]
            zeroList[i] = 1
            self.elements.append(zeroList)
            self.elements_ext[str(zeroList)] = i
        # Finally adding the primitive element to the coefficient list
        self.elements.append(primeCoeffs)
        self.elements_ext[str(primeCoeffs)] = highestDeg
        # FSR
        for j in range(highestDeg+1,self.p):
            # Shifts the FSR by removing the last element and inserting a zero element at the start of list
            # Checks if the highest polynomial degree has an element with coefficient larger than 0
            # Then it adds the primitive coefficients multiplied by the last coefficient modulo q
            highestCoeff = currentCoeffs.pop()
            currentCoeffs.insert(0,0)
            if highestCoeff > 0:
                for i in range(len(degrees)):
                    currentCoeffs[degrees[i]] = (currentCoeffs[degrees[i]] + coefficients[i]*highestCoeff) % self.q
            self.elements.append(currentCoeffs.copy())
            self.elements_ext[str(currentCoeffs.copy())] = (j)
            # Break the loop if the current coefficients are the same as the first element
            if currentCoeffs == self.elements[0]:
                del self.elements[-1]
                break

        print(self.p,"?=",len(self.elements),"?=",len(self.elements_ext))


    
    # Ext mapping to ground field
    # Converts all exponent representation list elements to vector representation
    # a: list containing exponent representation elements
    # singular: boolean variable for return of singular value or vector
    def ext(self,a,singular):
        if singular:
            if len(a):
                return self.elements[a[0]].copy()
            else:
                return [0 for i in range(self.m)]
        # A is a nested list containing exponent representation elements
        A = []
        for i in range(len(a)):
            if len(a[i]) != 0:
                # Find vector representation element corresponding to exponent
                column = self.elements[a[i][0]].copy()
            else:
                column = [0 for i in range(self.m)]
            A.append(column)
        return A
    
    # Mapping to extension field
    # Converts row elements to exponent representation
    # A: matrix containing elements in F_q
    # singular: boolean variable for return of singular value or vector
    def invext(self,A,singular):
        # Null element
        Null = [0 for i in range(self.m)]
        if singular:
            if A == Null:
                return []
            else:
                return [self.elements_ext.get(str(A))]
        a = []
        for row in A:
            # Check for zero element
            if row == Null:
                a.append([])
            else:
                # Find exponent representation corresponding to vector representation
                element = self.elements_ext.get(str(row))
                a.append([element])
        return a

    # The trace function for an element in F_q^m
    def Tr(self,a):
        # Find all the conjugates of "a" element
        conjugates = []
        for i in range(self.m):
            elem = (a[0]*self.q**i) % self.p
            if elem not in conjugates:
                conjugates.append(elem)
        # Find the sum of all conjugate elements
        elemsum = []
        for elem in conjugates:
            elemsum = self.add(elemsum,[elem],"+")
        if len(elemsum) == 0:
            return 0
        else:
            return self.elements[elemsum[0]][0]

    # Add and subtract alpha values
    # This method takes in alpha exponents (x in alpha^x)
    # The alpha exponents are then reduced to the real alpha expression 
    # The expressions are then added together modulo q for each element in each expression
    # The final alpha expression is converted into an alpha exponent representation and returned
    # NB: a and b come in the form [a],[b], where [] empty list represents null element
    # [0] is then alpha^0. This is done because the null element must be represented somehow
    # a: first alpha exponent
    # b: second alpha exponent to be added or subtracted from a
    def add(self,a,b,symbol):
        # Checks for null element (empty list)
        if len(a) == 0:
            # Creates a null element in reduced alpha expression form
            reducedA = [0 for i in range(self.m)]
        else:
            # Creates a reduced alpha expression form of the element
            reducedA = self.elements[a[0]].copy()
        # Do the exact same for b element
        if len(b) == 0:
            reducedB = [0 for i in range(self.m)]
        else:
            # Check whether to add or to subtract element b
            if symbol == "+":
                reducedB = self.elements[b[0]].copy()
            else:
                reducedB = [-1*i for i in self.elements[b[0]]]
        # Checks if a is equal to negative b, which will return zero element
        if reducedA == [-1*i for i in reducedB]:
            return []
        # Initializing resulting coefficient 
        c = []
        for i in range(len(reducedA)):
            coefficient = (reducedA[i] + reducedB[i]) % self.q
            c.append(coefficient)
        # Find the alpha degree corresponding to the reduced alpha degree expression
        if len(c) != 0 and c != [0 for i in range(self.m)]:
            return self.invext(c,True)
        else:
        # null or zero element case
            return []

    # Add and subtract two polynomials with alpha coefficients
    # Pairwise adds each element in the polynomial which have the same q-degrees
    # Uses alpha representation
    # aCoeff: alpha exponents of the coefficients for each element in polynomial
    # bCoeff: alpha exponents for the polynomial to be added or subtracted
    # aDeg: the q-degree of each element in the polynomial
    # bDeg: the q-degree of each element in the polynomial to be added or subtracted
    def addComp(self, aCoeff, bCoeff, aDeg, bDeg,symbol):
        # Initializing resulting coefficient list and q-degree list
        cCoeff = []
        cDeg = []
        # Add together all coefficients of elements that have the same q-degree
        for i in range(len(aDeg)):
            # Current degree of "a" polynomial term
            degree = aDeg[i]
            if degree in bDeg:
                # Index of corresponding element in b polynomial
                b_index = bDeg.index(degree)
                # Add two elements of same degree
                coefficient = self.add([aCoeff[i]],[bCoeff[b_index]],symbol)
                # Delete the "b" element corresponding to same degree and its coefficient
                del bDeg[b_index]
                del bCoeff[b_index]
            else:
                coefficient = [aCoeff[i]]
            # If the coefficient is not the zero element, then add it to the list
            if len(coefficient) != 0:
                cCoeff.append(coefficient[0])
                cDeg.append(degree)

        # For the case of difference, we have to correct to negated coefficients
        if symbol == "-":
            # Add the rest of the polynomial elements from b
            for i in range(len(bDeg)):
                # Add b coefficient to null element, this is done so that negated values are correct
                coefficient = self.add([],[bCoeff[i]],symbol)
                cCoeff.append(coefficient[0])
                cDeg.append(bDeg[i])
        return [cCoeff,cDeg]


    # Calculating the composite between two polynomial expressions using alpha exponents
    # aCoeff: coefficients of the leftmost element in the composition (alpha exponents)
    # bCoeff: coefficients of the rightmost element in the composition (alpha exponents)
    # aDeg: degrees of the polynomial of the leftmost element in composition
    # bDeg: degrees of the polynomial of the rightmost element in composition
    def composite(self,aCoeff,bCoeff,aDeg,bDeg):
        # Zero element case (empty list)
        if len(aCoeff) == 0 or len(bCoeff) == 0:
            return [[],[]]
        # Initializing coeff and poly lists
        coeff = []
        poly = []
        # Going through each element in bDeg for each element in aDeg
        for i in range(len(aDeg)):
            for j in range(len(bDeg)):
                # Adding together the q-degree of the element to be added to the poly list
                qDeg = aDeg[i] + bDeg[j]
                # Adding the product of the coefficients to the coeff list
                coeffProduct = (aCoeff[i]+bCoeff[j]*self.q**aDeg[i]) % self.p
                # If the q-degree already exists, add its coefficient to the current coefficient
                if qDeg in poly:
                    index = poly.index(qDeg)
                    coefficient = self.add([coeff[index]],[coeffProduct],"+")
                    # Checks for null element
                    if len(coefficient) != 0:
                        coeff[index] = coefficient[0]
                    else:
                    # Deletes any element that adds up to null element
                        del coeff[index]
                        del poly[index]
                # If the q-degree does not already exist, add it to the list of q-degrees and coefficients
                else:
                    poly.append(qDeg)
                    coeff.append(coeffProduct)
        return [coeff,poly]
    

    # Create a q-vandermonde matrix of input elements
    # a: vector containing alpha element degrees of length n
    # k: the number of rows in kxn matrix, the range of q-degrees to raise the a elements to
    def qvan(self,a,k,s):
        # Initializing the matrix
        matrix = []
        # Going through the range s
        for i in range(k):
            row = []
            # Going through each element in the "a" vector   
            for elem in a:
                if len(elem) != 0:
                    elem = elem[0]
                    # Raising the existing element degrees to q^i (* because "a" elements are given in degrees)
                    elem = (elem * (self.q**(s*i))) % self.p
                    row.append([elem])
                else:
                    row.append([])
            # Append the row to the matrix
            matrix.append(row)
        return matrix

    # Vector multiplied by matrix in extension field
    # v: vector containing elements in F_q^m
    def Codeword(self,v,A):
        b = []
        # Find the product of v*A
        for i in range(len(A[0])):
            productsum = []
            for j in range(len(v)):
                if len(v[j]) and len(A[j][i]):
                    coeff = [(v[j][0] + A[j][i][0]) % self.p]
                else:
                    coeff = []
                productsum = self.add(productsum,coeff,"+")
            b.append(copy.deepcopy(productsum))
        return b
    
    # Find the transpose of a matrix
    # A: matrix to be transposed
    def transpose(self,A):
        transposed = [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))] 
        return transposed
    
    # Row reduction in F_q to rref
    # A: matrix to be reduced
    # solve: return solution or rref and rank
    def RREF(self,A,solve,rref):
        # Leading columns will be added to this list
        pivots = []
        # Copy of a for checking solutions
        copyA = copy.deepcopy(A)
        # Number of rows
        m = len(A)
        # Number of columns
        n = len(A[0])
        # Current pivot row
        h = 0
        # Current pivot column
        k = 0
        # Rank
        rank = 0
        while h < m and k < n:
            pivot = -1
            for i in range(h,m):
                if A[i][k] != 0:
                    pivot = i
                    pivots.append(k)
                    rank += 1
                    break
            if pivot == -1:
                k += 1
            else:
                # Switch rows
                A[pivot], A[h] = A[h], A[pivot]
                leadingcoeff = A[h][k]
                if rref:
                    # Divide the pivot row by leadingcoefficient
                    l = 0
                    for j in range(k,n):
                        if A[h][j] != 0:
                            A[h][j] = (A[h][j] * mod_inverse(leadingcoeff,self.q)) % self.q
                else:
                    l = h+1
                # The leading coefficient the row number above
                leadingcoeff = A[h][k]
                for i in range(l,m):
                    if i != h:
                        # If the leading coefficient is not zero ([])
                        if A[i][k] != 0:
                            # Ratio between leading coefficients of both rows 
                            ratio = (A[i][k] * mod_inverse(leadingcoeff,self.q)) % self.q
                            for j in range(k,n):
                                if A[h][j] != 0:
                                    # Product of current row element with the ratio stated above
                                    ratioproduct = (A[h][j] * ratio) % self.q
                                    # Finding the difference
                                    A[i][j] = (A[i][j] - ratioproduct) % self.q
                h += 1
                k += 1

        if not solve:
            return [A,rank]
        if n in pivots:
            return "No solution"
        # Reverse order of leadin columns
        pivots.sort(reverse = True)
        # Initialize solution list with 1 ([0]) elements
        solution = [1 for i in range(n-1)]

        # Calculating solution based of coefficients of REF matrix
        for i in range(len(pivots)):
            # Current leading column
            leading = pivots[i]
            tot_sum = 0
            # Find the sum (product) of the coefficient and solution
            for j in range(leading+1,n-1):
                if A[leading][j] != 0 and solution[j] != 0:
                    product = (A[leading][j] * solution[j]) % self.q
                else:
                    product = 0
                # Sum all elements except leading element
                tot_sum = (tot_sum + product) % self.q
            # Find b-tot_sum
            tot_sum = (A[leading][-1] - tot_sum) % self.q
            # Find solution variable as the difference (division) between the sum
            # and the leading coefficient
            if tot_sum != 0:
                leading_sol = (tot_sum * mod_inverse(A[leading][leading],self.q)) % self.q
                solution[leading] = leading_sol
            else:
                solution[leading] = 0
        return solution

    # Row echelon form
    # Row reduction to echelon form using exponent representation elements and finite field arithmetic in F_q^m
    # A: matrix to be reduced, last column is the b vector
    # solve: whether the algorithm should return a solution or the REF matrix
    # rref: whether to get rref or ref
    def REF(self, A, solve, rref):
        # Leading columns will be added to this list
        pivots = []
        # Copy of a for checking solutions
        copyA = copy.deepcopy(A)
        # Number of rows
        m = len(A)
        # Number of columns
        n = len(A[0])
        # Current pivot row
        h = 0
        # Current pivot column
        k = 0
        # Rank counter
        rank = 0
        while h < m and k < n:
            pivot = -1
            for i in range(h,m):
                if len(A[i][k]):
                    pivot = i
                    pivots.append(k)
                    rank += 1
                    break
            if pivot == -1:
                k += 1
            else:
                # Switch rows
                A[pivot], A[h] = A[h], A[pivot]
                leadingcoeff = A[h][k][0]
                # ref or rref
                if rref:
                    # Divide the pivot row by leadingcoefficient
                    l = 0
                    for j in range(k,n):
                        if len(A[h][j]):
                            A[h][j][0] = (A[h][j][0] - leadingcoeff) % self.p
                else:
                    l = h+1
                # Redefine the leading coefficient of pivot row
                leadingcoeff = A[h][k][0]
                for i in range(l,m):
                    if i != h:
                        # If the leading coefficient is not zero ([])
                        if len(A[i][k]):
                            # Ratio between leading coefficients of both rows 
                            ratio = (A[i][k][0] - leadingcoeff) % self.p
                            for j in range(k,n):
                                if len(A[h][j]):
                                    # Product of current row element with the ratio stated above
                                    ratioproduct = [(A[h][j][0]+ratio) % self.p]
                                    # Finding the difference
                                    A[i][j] = self.add(A[i][j],ratioproduct,"-")
                        

                h += 1
                k += 1

        if not solve:
            return [A,rank]
        if n in pivots:
            return "No solution"
        # Reverse order of leadin columns
        pivots.sort(reverse = True)
        # Initialize solution list with 1 ([0]) elements
        solution = [[0] for i in range(n-1)]
        # Calculating solution based of coefficients of REF matrix
        for i in range(len(pivots)):
            # Current leading column
            leading = pivots[i]
            tot_sum = []
            # Find the sum (product) of the coefficient and solution
            for j in range(leading+1,n-1):
                if len(A[leading][j]) and len(solution[j]):
                    product = [(A[leading][j][0] + solution[j][0]) % self.p]
                else:
                    product = []
                # Sum all elements except leading element
                tot_sum = self.add(tot_sum,product,"+")
            # Find b-tot_sum
            tot_sum = self.add(A[leading][-1],tot_sum,"-")
            # Find solution variable as the difference (division) between the sum
            # and the leading coefficient
            if len(tot_sum):
                leading_sol = (tot_sum[0]-A[leading][leading][0]) % self.p
                solution[leading] = [leading_sol]
            else:
                solution[leading] = []

        '''
        ############################# CHECK SOLUTIONS ############################
        for j in range(len(copyA)):
            elem_sum = []
            for k in range(len(solution)):
                if len(solution[k]) and len(copyA[j][k]):
                    elemproduct = [(solution[k][0] + copyA[j][k][0]) % self.p]
                else:
                    elemproduct = []
                elem_sum = self.add(elem_sum,elemproduct,"+")
            print(elem_sum,"=?",copyA[j][len(copyA[j])-1])
        ###########################################################################
        '''
        return [solution,len(pivots)]
    
    def normal_basis(self,n):
    # Find normal basis
        # Construct all possible candidates for normal bases
        for i in range(1,self.p):
            # Matrix to contain elements in base field
            matrix = Matrix()
            basis = []
            for j in range(n):
                # Find the corresponding alpha exponents for all j from 0 to m-1
                degree = (i*self.q**j)%self.p
                # Find base field representation of the degree 
                elem = self.elements[degree].copy()
                # Put the base field representation as a column in the matrix
                matrix = matrix.col_insert(len(matrix),Matrix(elem))
                # Append the elements to the basis vector
                basis.append([degree])
            # Find the determinant of "matrix". When it is zero the columns are dependent
            # if not, the columns are linearly independent
            independent = matrix.rref(iszerofunc=lambda x: x % self.q==0)[1]
            if len(independent) == n:
                break
        normal = basis
        return normal



    # Matrix inversion in exponent representation
    def inv(self,A):
        n = len(A)
        # Initialize a zero vetor of length n
        empty = [[] for i in range(n)]
        for i in range(n):
            # Redefine the element at i to be [0]
            identity = copy.deepcopy(empty)
            identity[i] = [0]
            # Extend the matrix with identity matrix
            A[i].extend(identity)
        # Find the RREF form of A
        inv = self.REF(A,False,True)[0]
        # Extract last n columns as inverse matrix
        inv = [row[n:2*n] for row in inv]
        return inv

    # Matrix inversion in F_q
    def inv_q(self,A):
        n = len(A)
        # Initialize a zero vetor of length n
        empty = [0 for i in range(n)]
        for i in range(n):
            # Redefine the element at i to be [0]
            identity = empty.copy()
            identity[i] = 1
            # Extend the matrix with identity matrix
            A[i].extend(identity)
        # Find the RREF form of A
        inv = self.RREF(A,False,True)[0]
        # Extract last n columns as inverse matrix
        inv = [row[n:2*n] for row in inv]
        return inv