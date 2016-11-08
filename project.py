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
            # STEP 1: Check for optimality
            while not self.isVectorPositive(self.Z_N):
                time.sleep(1)
                print "++++++++++++++++++++++++++++"
                # STEP 2:
                # Get the least negative number Index( i.e. entering Index)
                # from Non Basic vector
                j = self.Non_Basic[self.Z_N.argmin()]

                # STEP 3: Calculte delata_X_b

                # to create a unit vector, with all element zero except 1
                # np.eye(value,size_of_vector,index_of Value)
                e_j = np.eye(1, len(self.Non_Basic) , self.Non_Basic.index(j))
                e_j = np.transpose(e_j)
                delta_X_b = (inv(self.B).dot(self.N)).dot(e_j)
                delta_X_b = np.reshape(delta_X_b,(delta_X_b.shape[0],))

                # STEP 4: Calculate Primal Step Length
                t , t_index = self.primalStepLength(delta_X_b,self.X_b)

                # Step 5: Select Leaving Variable
                # max ratio corresponds to index from Basic (Leaving Variable)
                i = self.Basic[t_index]

                # STEP 6: Compute Dual Step Direction
                # to create a unit vector, with all element zero except 1
                # np.eye(value,size_of_vector,index_of Value)

                e_i = np.eye(1, len(self.Basic) , self.Basic.index(i))
                e_i = np.transpose(e_i)
                delta_Z_n = - (np.transpose((inv(self.B)).dot(self.N))).dot(e_i)
                delta_Z_n = np.reshape(delta_Z_n,(delta_Z_n.shape[0],))

                # STEP 7: Compute Dual Step Length
                s = self.Z_N[self.Non_Basic.index(j)] / delta_Z_n[self.Non_Basic.index(j)]

                # STEP 8: Update Current Primal and Dual Solutions
                new_x = t
                self.X_b = self.X_b - t * delta_X_b
                new_z = s
                self.Z_N = self.Z_N - s * delta_Z_n

                # Step 9: Update Basis
                self.Non_Basic[self.Non_Basic.index(j)] = i
                self.Basic[self.Basic.index(i)] = j

                b_columns = self.Basic + np.array([-1.0]*len(self.Basic))
                n_columns = self.Non_Basic + np.array([-1.0]*len(self.Non_Basic))

                self.B = self.dict_A[: , b_columns.astype(np.int64)]
                self.N = self.dict_A[: , n_columns.astype(np.int64)]

                self.X_b[self.Basic.index(j)] = new_x
                self.Z_N[self.Non_Basic.index(i)] = new_z

                self.printAllVariables()
                #self.Z_N = self.Z_N * 0


    ### Calculate Primal Step Length
    ### Divide element by element (also conider 0/0 as 0)
    ### takes the max of the resulted list and return inverse and
    ### index corresponding to max
    def primalStepLength(self,delta_x,delta_x_i):
        print 'Primal Step Length'
        temp_list=[]
        for i in range(delta_x_i.shape[0]):
            temp_list.append(delta_x[i]/delta_x_i[i])

        max_val = max(temp_list)
        t = 1/max_val
        return t, temp_list.index(max_val)


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


