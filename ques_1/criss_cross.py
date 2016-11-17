import sys
import numpy as np
from numpy import genfromtxt
import time
from numpy.linalg import inv

class CrissCross(object):

    def __init__(self):
        if len(sys.argv)!=4:
            print "Please pass all the files."
            self.dataEntered = False
            return

        ### Simplex Dictionary, Vectors and variables
        self.dict_A = genfromtxt(sys.argv[1], delimiter=',') # Dictioanry A
        self.C_n = genfromtxt(sys.argv[2], delimiter=',') # Vector C
        self.b = genfromtxt(sys.argv[3], delimiter=',') # Vector b
        self.numberofConstraints = self.dict_A.shape[0] # m
        self.numberOfVaraibles = self.C_n.shape[0]    # n
        self.C_starts = 0
        self.B = np.identity(self.numberofConstraints)   # Identity Matrix (For Basic Vars)
        self.N = self.dict_A                        # Non-Basic columns
        self.Old_C_n = []
        self.dict_A = np.concatenate((self.N, self.B), axis=1) # [N B]

        ### track Basic-Var Columns
        self.Basic =  range(self.N.shape[1]+1, self.N.shape[1] + self.numberofConstraints+1)
        ### track Non Basic-Var Columns
        self.Non_Basic = range(1,self.N.shape[1]+1)
        # keep copy on Non_Basic for
        self.old_Non_Basic = self.Non_Basic
        # keep copy on Basic for
        self.old_Basic = self.Basic

        ### Associated Dual Dictionary
        #self.Z_n = -1 * self.C_n

        self.X_n = np.zeros(len(self.Non_Basic))
        self.X_b =  - (inv(self.B).dot(self.N)).dot(self.X_n)

        self.C_b = np.zeros(len(self.Basic))
        self.Z_n = (np.transpose(inv(self.B).dot(self.N))).dot(self.C_b) -1 * self.C_n

        self.dataEntered = True

    def preformCrissCross(self):

        # Step 1 If Z_n>=0 and X_b >=0 Stop
        if self.isVectorPositive(self.Z_n) and self.isVectorPositive(self.X_b):
            print 'Dictionary is not primal and dual infeasible'

        while self.isVectorPositive(self.Z_n) or self.isVectorPositive(self.X_b):
            # step 2
            j = np.argmax(self.Z_n<0)
            j = self.Non_Basic[j]

            # step 3
            e_j = np.eye(1, len(self.Non_Basic) , self.Non_Basic.index(j))
            e_j = np.transpose(e_j)

            delta_X_b = ((inv(self.B)).dot(self.N)).dot(e_j)

            # step 4
            i = np.argmax(delta_X_b<0)
            j = self.Non_Basic[j]

            print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            print "j : %s " % j
            print "\ne_j : \n%s " % e_j
            print "\ndelta_X_b : \n%s " % delta_X_b
            print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"


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

        print "C_n :",
        print self.C_n

        print "C_b :",
        print self.C_b

        print "b :",
        print self.b

        print "m : %s" % self.numberofConstraints
        print "n : %s" % self.numberOfVaraibles

        print "B :"
        print self.B

        print "N :"
        print self.N

        print "Basic Columns : %s" % self.Basic
        print "Non Basic Columns : %s" % self.Non_Basic

        print "X_b :",
        print self.X_b

        print "X_n :",
        print self.X_n

        print "Z_n :",
        print self.Z_n
        print



### create an object of Linear_Prog class
crissCross = CrissCross()
crissCross.printAllVariables()
crissCross.preformCrissCross()
