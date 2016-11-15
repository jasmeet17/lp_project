import sys
import numpy as np
from numpy import genfromtxt
import time
from numpy.linalg import inv
import latex_sample

old_NBasic=[]
old_Basic=[]

class LinearProgram(object):

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
        self.X_b = inv(self.B).dot(self.b) - (inv(self.B).dot(self.N)).dot(self.X_n)

        self.C_b = np.zeros(len(self.Basic))
        self.Z_n = (np.transpose(inv(self.B).dot(self.N))).dot(self.C_b) -1 * self.C_n

        self.dataEntered = True
        self.latex_text = latex_sample.getInitialMatrices(self.dict_A, self.Basic, self.Non_Basic, self.B, self.N, self.X_b, self.Z_n)

    # Since the problem in hand is not dual/primal feasible
    # applies the phase 1
    def phaseOne(self):
        print 'Phase 1'
        # keep the copy of old Objective Fucntion
        old_C_n = self.C_n

        # make new primal objective function
        self.C_n =  np.empty(self.numberOfVaraibles)
        self.C_n.fill(-1)

        # Update Z_n according to the new C_n
        self.Z_n = (np.transpose(inv(self.B).dot(self.N))).dot(self.C_b) -1 * self.C_n

        # get the initial non basic indexes array
        for x in self.old_Non_Basic:
            old_NBasic.append(x)

        # get the initial basic indexes array
        for x in self.old_Basic:
            old_Basic.append(x)

        # Now we have new primal function and the problem,
        # becomes Dual Feasible.
        print '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        self.printAllVariables()
        print "preformDualSimplex ===================================="
        self.preformDualSimplex()
        print 'Now the problem is Primal Feasible, apply Primal Simplex'

        print self.Basic
        print old_Basic
        print "@#@###@#@#@#@##@"
        print self.Non_Basic
        print old_NBasic
        print "&*&*&*&*&*&*&*&*"
        return


        # "the new A is  values are:"
        temp_A = -(inv(self.B)).dot(self.N)

        # number of terms in the Objective funtion
        n_terms = temp_A.shape[1] + 1

        # new Objective Funciton
        sum_array = np.zeros(n_terms)

        for i in range(len(old_NBasic)):
            temp_array = []
            if old_NBasic[i] in self.Basic:
                t = self.Basic.index(old_NBasic[i])
                temp_array = old_C_n[i] * np.append([self.X_b[t]],temp_A[t])
            else:
                temp_array = np.zeros(n_terms)
                temp_array[i+1] = old_C_n[i]
            # print "$$$$$"
            # print temp_array
            sum_array += temp_array

        for i in range(len(self.C_n)):
            self.C_n[i] = sum_array[i+1]

        # print '#$#$#$#$#$#$#$#$#$'
        # print self.C_n
        # return
        self.C_starts = sum_array[0]
        self.Z_n = -1 * self.C_n
        print self.C_n
        print "++++++++++++++"
        return
        self.printAllVariables()
        # Now we got the updated Objective function after
        # applying the phase 1, now apply primal simplex
        self.preformPrimalSimplex()

    ### performs the Simplex method
    def preformPrimalSimplex(self):
        iteration = 1
        if not self.isVectorPositive(self.X_b):
            print "X_b is < 0 Initial solution is not primal feasible."
            self.latex_text += latex_sample.getPrimalInitialCondition(False)
            self.phaseOne()
            return
        else :
            print 'X_b >= 0 Initial solution is primal feasible.'
            self.latex_text += latex_sample.getPrimalInitialCondition(True)

        if self.isVectorPositive(self.Z_n):
            # 'Z_n >=0 Current solution is optimal.'
            self.latex_text += latex_sample.firstStepPrimal(iteration,False)
            self.printObjectiveFunction(self.C_n)
            print "Objective Function Value : %s" % self.getObjectiveValue(self.C_n,self.X_b, self.Basic, self.numberOfVaraibles)
            return
        else:
            ### until theres some negative in Z_n
            # STEP 1: Check for optimality
            while not self.isVectorPositive(self.Z_n):
                self.latex_text += latex_sample.firstStepPrimal(iteration,True)

                # STEP 2:
                # Get the least negative number Index( i.e. entering Index)
                # from Non Basic vector
                j = self.Non_Basic[self.Z_n.argmin()]
                self.latex_text += latex_sample.secondStepPrimal(self.Z_n.argmin(),j)

                # STEP 3: Calculte delata_X_b

                # to create a unit vector, with all element zero except 1
                # np.eye(value,size_of_vector,index_of Value)
                e_j = np.eye(1, len(self.Non_Basic) , self.Non_Basic.index(j))
                self.X_n = e_j[0]
                e_j = np.transpose(e_j)
                delta_X_b = (inv(self.B).dot(self.N)).dot(e_j)
                delta_X_b = np.reshape(delta_X_b,(delta_X_b.shape[0],))
                self.latex_text += latex_sample.thirdStepPrimal(inv(self.B).dot(self.N),e_j,delta_X_b)

                # STEP 4: Calculate Primal Step Length
                t , t_index,pivot_flag = self.primalStepLength(delta_X_b,self.X_b)

                # if t less than or equal zero problem in unbounded, So STOP
                if t > 0:
                    self.latex_text += latex_sample.fourthStepPrimal(delta_X_b,self.X_b,t)
                else:
                    print delta_X_b
                    print self.X_b
                    print 'Problem in unbounded. Stop.'
                    return

                # Step 5: Select Leaving Variable
                # max ratio corresponds to index from Basic (Leaving Variable)
                i = self.Basic[t_index]
                self.latex_text += latex_sample.fifthStepPrimal(i)

                # STEP 6: Compute Dual Step Direction
                # to create a unit vector, with all element zero except 1
                # np.eye(value,size_of_vector,index_of Value)

                e_i = np.eye(1, len(self.Basic) , self.Basic.index(i))
                e_i = np.transpose(e_i)
                delta_Z_n = - (np.transpose((inv(self.B)).dot(self.N))).dot(e_i)
                delta_Z_n = np.reshape(delta_Z_n,(delta_Z_n.shape[0],))
                self.latex_text += latex_sample.sixthStepPrimal(np.transpose((inv(self.B)).dot(self.N)),e_i,delta_Z_n)

                # STEP 7: Compute Dual Step Length
                s = self.Z_n[self.Non_Basic.index(j)] / delta_Z_n[self.Non_Basic.index(j)]
                self.latex_text += latex_sample.seventhStepPrimal(s,self.Z_n[self.Non_Basic.index(j)] , delta_Z_n[self.Non_Basic.index(j)])

                # STEP 8: Update Current Primal and Dual Solutions
                new_x = t
                old_X_b = self.X_b
                self.X_b = self.X_b - t * delta_X_b
                new_z = s
                old_Z_n = self.Z_n
                self.Z_n = self.Z_n - s * delta_Z_n
                self.latex_text += latex_sample.eightStepPrimal(j,i,t,s,old_X_b,old_Z_n,delta_X_b,delta_Z_n,self.X_b,self.Z_n)

                # Step 9: Update Basis
                self.Non_Basic[self.Non_Basic.index(j)] = i
                self.Basic[self.Basic.index(i)] = j

                b_columns = self.Basic + np.array([-1.0]*len(self.Basic))
                n_columns = self.Non_Basic + np.array([-1.0]*len(self.Non_Basic))

                self.B = self.dict_A[: , b_columns.astype(np.int64)]
                self.N = self.dict_A[: , n_columns.astype(np.int64)]

                self.X_b[self.Basic.index(j)] = new_x
                self.Z_n[self.Non_Basic.index(i)] = new_z

                self.latex_text += latex_sample.ninthStepPrimal(self.Basic,self.Non_Basic,self.B,self.N,self.X_b,self.Z_n)

                #self.X_b = (inv(self.B)).dot(self.b)
                #self.X_b = (inv(self.B)).dot(self.b) - ((inv(self.B)).dot(self.N)).dot(self.X_n)
                # self.Z_n = (np.transpose(inv(self.B).dot(self.N))).dot(self.C_b) -1 * self.C_n

                # self.printAllVariables()
                iteration+=1
                self.printAllVariables()

            #print "Iteration number : %s" % iteration
            self.latex_text += latex_sample.firstStepPrimal(iteration,False)
            self.printObjectiveFunction(self.C_n)

            print "Objective Function Value : %s" % self.getObjectiveValue(self.C_n,self.X_b, self.Basic, self.numberOfVaraibles)

    ### performs Dual Simplex method
    def preformDualSimplex(self):
        if not self.isVectorPositive(self.Z_n):
            print "Z_n is <= 0 "
            print "Initial solution is not Dual feasible."
            return
        else :
            print 'Z_n >= 0'
            print "Initial solution is Dual feasible."

        if self.isVectorPositive(self.X_b):
            print 'X_b >=0'
            print 'Current solution is optimal.'
            self.printObjectiveFunction(self.C_n)
            print "Objective Function Value : %s" % self.getObjectiveValue(self.C_n,self.X_b, self.Basic, self.numberOfVaraibles)
        else:
            iteration = 1
            ### until theres some negative in Z_n
            # STEP 1: Check for optimality
            while not self.isVectorPositive(self.X_b):
                #time.sleep(1)
                # to count the iteration number
                print "Iteration number : %s" % iteration
                # STEP 2:
                # Get the least negative number Index( i.e. entering Index)
                # from Non Basic vector
                #@@@ j = self.Non_Basic[self.Z_n.argmin()]
                i = self.Basic[self.X_b.argmin()]

                # STEP 3: Calculte delata_Z_n

                # to create a unit vector, with all element zero except 1
                # np.eye(value,size_of_vector,index_of Value)
                e_i = np.eye(1, len(self.Basic) , self.Basic.index(i))
                e_i = np.transpose(e_i)
                delta_Z_n = - (np.transpose(inv(self.B).dot(self.N))).dot(e_i)
                delta_Z_n = np.reshape(delta_Z_n,(delta_Z_n.shape[0],))

                # STEP 4: Calculate Primal Step Length
                s , s_index, pivot_flag = self.primalStepLength(delta_Z_n,self.Z_n)

                # if s less than zero problem in unbounded, So STOP.
                if s >= 0:
                    pass
                else:
                    print 'Problem in unbounded. Stop.'
                    return

                # Step 5: Select Leaving Variable
                # max ratio corresponds to index from Basic (Leaving Variable)
                j = self.Non_Basic[s_index]

                # STEP 6: Compute Dual Step Direction
                # to create a unit vector, with all element zero except 1
                # np.eye(value,size_of_vector,index_of Value)

                e_j = np.eye(1, len(self.Non_Basic) , self.Non_Basic.index(j))
                e_j = np.transpose(e_j)
                delta_X_b = ((inv(self.B)).dot(self.N)).dot(e_j)
                delta_X_b = np.reshape(delta_X_b,(delta_X_b.shape[0],))

                # STEP 7: Compute Dual Step Length
                t = self.X_b[self.Basic.index(i)] / delta_X_b[self.Basic.index(i)]

                # STEP 8: Update Current Primal and Dual Solutions
                new_x = t
                self.X_b = self.X_b - t * delta_X_b
                new_z = s
                self.Z_n = self.Z_n - s * delta_Z_n

                # Step 9: Update Basis
                self.Non_Basic[self.Non_Basic.index(j)] = i
                self.Basic[self.Basic.index(i)] = j

                b_columns = self.Basic + np.array([-1.0]*len(self.Basic))
                n_columns = self.Non_Basic + np.array([-1.0]*len(self.Non_Basic))

                self.B = self.dict_A[: , b_columns.astype(np.int64)]
                self.N = self.dict_A[: , n_columns.astype(np.int64)]

                self.X_b[self.Basic.index(j)] = new_x
                self.Z_n[self.Non_Basic.index(i)] = new_z

                self.printAllVariables()
                iteration+=1
                #self.Z_n = self.Z_n * 0

            print "======= ENDS ======="
            #print "Iteration number : %s" % iteration
            self.printObjectiveFunction(self.C_n)
            print "Objective Function Value : %s" % self.getObjectiveValue(self.C_n,self.X_b, self.Basic, self.numberOfVaraibles)

    ### Calculate Primal Step Length
    ### Divide element by element (also conider 0/0 as 0)
    ### takes the max of the resulted list and return inverse and
    ### index corresponding to max
    def primalStepLength(self,delta_x,delta_x_i):
        temp_list=[]
        for i in range(delta_x_i.shape[0]):
            if delta_x_i[i]==0:
                temp_list.append(0)
            else:
                temp_list.append(delta_x[i]/delta_x_i[i])

        max_val = max(temp_list)
        if max_val == 0:
            t = 0
            pivot_flag = True
        else:
            t = 1/max_val
            pivot_flag = False

        return float(t), temp_list.index(max_val) , pivot_flag


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

    # vec_c -> objective function
    # vec_x -> Soultion
    # basic -> index of basic vars
    # n -> number of varaibles in objective function, initially
    # returns the objective Function value
    def getObjectiveValue(self,vec_c,vec_x,basic_vec,n):
        # this count just to check which has max
        print "@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@"
        print vec_c
        print vec_x
        print basic_vec
        print n
        value = 0
        for c in range(len(vec_c)):
            if basic_vec[c]<=n:
                value += vec_c[c] * vec_x[c]

        value += self.C_starts
        return value

    # Prints the Objective function in the Standard form
    # and the Value obtained by the Solution
    def printObjectiveFunction(self, vec_c):

        max_func = ""
        for c in range(vec_c.shape[0]):

            # First Check the signs, If value==Zero ,
            # Just Continue(Beause We Dont want to
            # show varaibles with zero constraints)
            if vec_c[c]==0:
                continue
            elif vec_c[c]<0:
                max_func +="  "
            elif max_func!="":
                max_func +=" + "

            max_func += str(np.round(vec_c[c], decimals=2))
            max_func += " X"+str(c+1)

        print "MAX : %s" % max_func

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
simplex = LinearProgram()
# simplex.printAllVariables()
simplex.preformPrimalSimplex()
# simplex.preformDualSimplex()


'''
from tex import latex2pdf
f = open('simplex.tex','w')
f.write(simplex.latex_text) # python will convert \n to os.linesep
f.close()
'''
print "----"












