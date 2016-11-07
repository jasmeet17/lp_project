import sys
import numpy as np
from numpy import genfromtxt
import time
from numpy.linalg import inv

class Simplex(object):

    def __init__(self):
        if len(sys.argv)!=4:
            print "Please pass all the files."
            self.dataEntered = False
            return

        ### Simplex Dictionary, Vectors and variables
        self.dict_A = genfromtxt(sys.argv[1], delimiter=',') # Dictioanry A
        self.vec_c = genfromtxt(sys.argv[2], delimiter=',') # Vector C
        self.X_b = genfromtxt(sys.argv[3], delimiter=',') # Vector b
        self.numberofConstraints = self.dict_A.shape[0] # m
        self.numberOfVaraibles = self.vec_c.shape[0]    # n

        self.B = np.identity(self.numberofConstraints)   # Identity Matrix (For Basic Vars)
        self.N = self.dict_A                        # Non-Basic columns

        self.dict_A = np.concatenate((self.N, self.B), axis=1) # [N B]

        ### track Basic-Var Columns
        self.Basic =  range(self.N.shape[1]+1, self.N.shape[1] + self.numberofConstraints+1)
        ### track Non Basic-Var Columns
        self.Non_Basic = range(1,self.N.shape[1]+1)

        ### Associated Dual Dictionary
        self.Z_N = -1 * self.vec_c

        self.dataEntered = True

    ### performs the Simplex method
    def preformSimplex(self):
        if not self.isVectorPositive(self.X_b):
            print "X_b is <= 0 "
            print "Initial solution is not primal feasible."
        else :
            print 'X_b >= 0'
            print "Initial solution is primal feasible."

        if self.isVectorPositive(self.Z_N):
            print 'Z_n >=0'
            print 'Current solution is optimal.'
        else:
            ### until theres some negative in Z_n
            while not self.isVectorPositive(self.Z_N):
                time.sleep(1)
                print "++++++++++++++++++++++++++++"
                # STEP 2
                # Get the least negative number Index( i.e. entering Index)
                j = self.Non_Basic[self.Z_N.argmin()]

                # STEP 3: Calculte delata_X_b

                # to create a unit vector, with all element zero except 1
                # np.eye(value,size_of_vector,index_of Value)
                e_j = np.eye(1, len(self.Non_Basic) , j-1)
                e_j = np.transpose(e_j)
                delata_X_b = (self.B.dot(self.N)).dot(e_j)

                # STEP 4: Calculate Primal Step Length


    ### To check all elements in the vectors are positive
    ### checks Greater than or equal to zero and returns True If so
    def isVectorPositive(self,vec):
        zero_array = np.zeros(vec.shape[0])
        compare = vec >= zero_array
        if False in compare:
            # there's a negative number
            return False
        else:
            # All numbers are Positive
            return True


    def printAllVariables(self):
        print "======== Varibales values ========"
        print "A :"
        print self.dict_A
        print self.dict_A.shape

        print "C :",
        print self.vec_c

        print "b :",
        print self.X_b

        print "m : %s" % self.numberofConstraints
        print "n : %s" % self.numberOfVaraibles

        print "B :"
        print self.B

        print "N :"
        print self.N

        print "Basic Columns : %s" % self.Basic
        print "Non Basic Columns : %s" % self.Non_Basic

        print "Z_N :",
        print self.Z_N




### create an object of Linear_Prog class
simplex = Simplex()
simplex.printAllVariables()
print "============================"


simplex.preformSimplex()

'''
# run assignment Tests
if simplex.dataEntered:
    simplex.runAssignmentTest()
'''


