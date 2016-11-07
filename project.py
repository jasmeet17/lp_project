import sys
import numpy as np
from numpy import genfromtxt

class Simplex(object):

    def __init__(self):
        if len(sys.argv)!=4:
            print "Please pass all the files."
            self.dataEntered = False
            return

        # Simplex Dictionary, Vectors and variables
        self.dict_A = genfromtxt(sys.argv[1], delimiter=',') # Dictioanry A
        self.vec_c = genfromtxt(sys.argv[2], delimiter=',') # Vector C
        self.vec_b = genfromtxt(sys.argv[3], delimiter=',') # Vector b
        self.numberofConstraints = self.dict_A.shape[0] # m
        self.numberOfVaraibles = self.vec_c.shape[0]    # n

        self.B = np.identity(self.numberofConstraints)   # Identity Matrix (For Basic Vars)
        self.N = self.dict_A                        # Non-Basic columns

        self.dict_A = np.concatenate((self.N, self.B), axis=1) # [N B]

        # track Basic-Var Columns
        self.Basic =  range(self.N.shape[1]+1, self.N.shape[1] + self.numberofConstraints+1)
        # track Non Basic-Var Columns
        self.Non_Basic = range(1,self.N.shape[1]+1)

        # Associated Dual Dictionary
        self.Z_N = -1 * self.vec_c

        self.dataEntered = True

    def preformSimplex(self):
        self.vectorPositive(self.Z_N)
        self.vectorPositive(self.vec_b)

    # To check all elements in the vectors are positive
    def vectorPositive(self,vec):
        zero_array = np.zeros(vec.shape[0])
        compare = vec > zero_array
        print compare


    def printAllVariables(self):
        print "======== Varibales values ========"
        print "A :"
        print self.dict_A

        print "C :",
        print self.vec_c
        print self.vec_c.shape

        print "b :",
        print self.vec_b

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

        print "======== ======== ======== ========"




# create an object of Linear_Prog class
simplex = Simplex()
simplex.printAllVariables()

simplex.preformSimplex()

'''
# run assignment Tests
if simplex.dataEntered:
    simplex.runAssignmentTest()
'''


