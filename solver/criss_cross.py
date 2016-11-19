import sys
import numpy as np
from numpy import genfromtxt
import time
from numpy.linalg import inv
import ques_1 as linear


class CrissCross(object):

    def __init__(self,dict_A,C_n,b):
        if len(sys.argv)!=4:
            # print "Please pass all the files."
            # self.dataEntered = False
            # return
            pass

        ### Simplex Dictionary, Vectors and variables
        self.dict_A = dict_A # Dictioanry A
        self.C_n = C_n # Vector C
        self.b = b # Vector b
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
        self.X_b =  inv(self.B).dot(self.b) - (inv(self.B).dot(self.N)).dot(self.X_n)

        self.C_b = np.zeros(len(self.Basic))
        self.Z_n = (np.transpose(inv(self.B).dot(self.N))).dot(self.C_b) -1 * self.C_n

        self.dataEntered = True

    def findMaxMinSubscript(self,vec,vec_NB,max_flag):
        if max_flag:
            j_z = np.where(vec>=0)
        else:
            j_z = np.where(vec<0)

        j_z = j_z[0]
        # print j_z
        if len(j_z)==1:
            j_z = j_z[0]
            return vec_NB[j_z]
        else:
            smallest = np.amax(vec_NB)
            for i in range(len(j_z)):
                if vec_NB[j_z[i]]<smallest:
                    smallest = vec_NB[j_z[i]]
            return smallest



    def callNonBasicSteps(self,j):

        e_j = np.eye(1, len(self.Non_Basic) , self.Non_Basic.index(j))
        e_j = np.transpose(e_j)
        delta_X_b = ((inv(self.B)).dot(self.N)).dot(e_j)

        i = self.findMaxMinSubscript(delta_X_b,self.Basic,True)
        leaving_index = self.Basic.index(i)


        # STEP 4: Calculate Primal Step Length
        max_val = self.primalStepLength(delta_X_b,self.X_b,leaving_index)

        t = 0
        t = 1/max_val

        e_i = np.eye(1, len(self.Basic) , self.Basic.index(i))
        e_i = np.transpose(e_i)

        delta_Z_n = - (np.transpose((inv(self.B)).dot(self.N))).dot(e_i)

        s = self.Z_n[self.Non_Basic.index(j)]/delta_Z_n[self.Non_Basic.index(j)]
        s= float(s[0])

        self.X_b = self.X_b - t * np.transpose(delta_X_b)
        self.X_b = self.X_b[0]
        self.X_b[self.Basic.index(i)]=t

        self.Z_n = self.Z_n - s * np.transpose(delta_Z_n)
        self.Z_n = self.Z_n[0]
        self.Z_n[self.Non_Basic.index(j)]=s


        # Update Basis
        self.Non_Basic[self.Non_Basic.index(j)] = i
        self.Basic[self.Basic.index(i)] = j

        b_columns = self.Basic + np.array([-1.0]*len(self.Basic))
        n_columns = self.Non_Basic + np.array([-1.0]*len(self.Non_Basic))

        self.B = self.dict_A[: , b_columns.astype(np.int64)]
        self.N = self.dict_A[: , n_columns.astype(np.int64)]



    def callBasicSteps(self,i):

        e_i = np.eye(1, len(self.Basic) , self.Basic.index(i))
        e_i = np.transpose(e_i)
        # delta_X_b = ((inv(self.B)).dot(self.N)).dot(e_i)
        delta_Z_n = - (np.transpose((inv(self.B)).dot(self.N))).dot(e_i)

        j = self.findMaxMinSubscript(delta_Z_n,self.Non_Basic,True)
        entering_index = self.Non_Basic.index(j)
        # STEP 4: Calculate Primal Step Length
        max_val = self.primalStepLength(delta_Z_n,self.Z_n,entering_index)

        s = 0
        s = 1/max_val

        e_j = np.eye(1, len(self.Non_Basic) , self.Non_Basic.index(j))
        e_j = np.transpose(e_j)

        delta_X_b = ((inv(self.B)).dot(self.N)).dot(e_j)

        t = self.X_b[self.Basic.index(i)]/delta_X_b[self.Basic.index(i)]
        t= float(t[0])

        self.X_b = self.X_b - t * np.transpose(delta_X_b)
        self.X_b = self.X_b[0]
        self.X_b[self.Basic.index(i)]=t

        self.Z_n = self.Z_n - s * np.transpose(delta_Z_n)
        self.Z_n = self.Z_n[0]
        self.Z_n[self.Non_Basic.index(j)]=s

        self.Non_Basic[self.Non_Basic.index(j)] = i
        self.Basic[self.Basic.index(i)] = j

        # Update Basis

        b_columns = self.Basic + np.array([-1.0]*len(self.Basic))
        n_columns = self.Non_Basic + np.array([-1.0]*len(self.Non_Basic))

        self.B = self.dict_A[: , b_columns.astype(np.int64)]
        self.N = self.dict_A[: , n_columns.astype(np.int64)]



    def preformCrissCross(self):

        # Step 1 If Z_n>=0 and X_b >=0 Stop
        if self.isVectorPositive(self.Z_n) and self.isVectorPositive(self.X_b):
            simplex = linear.LinearProgram(self.dict_A,self.C_n,self.b)
            flag , value =simplex.preformPrimalSimplex()

        while not self.isVectorPositive(self.Z_n) or not self.isVectorPositive(self.X_b):
            # step 2
            j_z = self.findMaxMinSubscript(self.Z_n,self.Non_Basic,False)

            i_z = self.findMaxMinSubscript(self.X_b,self.Basic,False)

            if not self.isVectorPositive(self.Z_n):
                self.callNonBasicSteps(j_z)
            else:
                self.callBasicSteps(i_z)

        value = self.getObjectiveValue([],self.C_n, self.X_b, self.Basic, self.numberOfVaraibles)
        return value

    # vec_c -> objective function
    # vec_x -> Soultion
    # basic -> index of basic vars
    # n -> number of varaibles in objective function, initially
    # returns the objective Function value
    def getObjectiveValue(self,old_c,vec_c,vec_x,basic_vec,n):

        value = 0
        if len(old_c)>0:
            vec_c=old_c

        for i in range(len(basic_vec)):
            if basic_vec[i]<=len(vec_c):
                value += vec_c[basic_vec[i]-1]*vec_x[i]


        return value

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

    def primalStepLength(self,delta_x,delta_x_i,leaving_index):

        temp_list=[]
        infinte_index = -1

        for i in range(delta_x_i.shape[0]):
            if delta_x_i[i]==0:
                if delta_x[i]==0:
                    temp_list.append(0)
                elif delta_x[i]<0:
                    temp_list.append(0)
                elif delta_x[i]>0:
                    infinte_index = i
            else:
                temp_list.append(delta_x[i]/delta_x_i[i])

        return float(temp_list[leaving_index])

'''
dict_A = genfromtxt(sys.argv[1], delimiter=',') # Dictioanry A
C_n = genfromtxt(sys.argv[2], delimiter=',') # Vector C
b = genfromtxt(sys.argv[3], delimiter=',') # Vector b

### create an object of Linear_Prog class
crissCross = CrissCross(dict_A,C_n,b)
crissCross.preformCrissCross()
'''