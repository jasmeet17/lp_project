import numpy as np
import unittest

import ques_1 as linear


class Test(unittest.TestCase):
    def setUp(self):
        self.a = np.asarray([[1 , -1], [2, -1], [0, 1]]).astype(np.float)
        self.b = np.asarray([1,3,5]).astype(np.float)
        self.c = np.asarray([4,3]).astype(np.float)

        self.a1 = np.asarray([[-2 , -1], [-2, 4], [-1, 3]]).astype(np.float)
        self.b1 = np.asarray([4,-8,-7]).astype(np.float)
        self.c1 = np.asarray([-1,-1]).astype(np.float)

        self.N1 = np.asarray([[ 0 ,-1],[ 0 , 4],[ 1 , 3]]).astype(np.float)
        self.B1 = np.asarray([[ 1 ,-2 , 0], [ 0 ,-2 , 1], [ 0 ,-1 , 0]]).astype(np.float)

        self.a2 = np.asarray([[-2 , -1], [-2, 4], [-1, 3]]).astype(np.float)
        self.b2 = np.asarray([4,-8,-7]).astype(np.float)
        self.c2 = np.asarray([-1,4]).astype(np.float)

        self.a3 = np.asarray([[1 , 4, 0], [3 , -1, 1]]).astype(np.float)
        self.b3 = np.asarray([1,3]).astype(np.float)
        self.c3 = np.asarray([4, 1, 3]).astype(np.float)

        self.a4 = np.asarray([[1 , -1], [-1, -1], [2 , 1]]).astype(np.float)
        self.b4 = np.asarray([-1,-3, 4]).astype(np.float)
        self.c4 = np.asarray([3, 1]).astype(np.float)


    '''Check the Basic and non Baic Indices return are corrent'''
    def test1(self):
        print 'TEST 1---------------------------------------------- TEST 1'
        simplex = linear.LinearProgram(self.a,self.c,self.b)
        flag , value =simplex.preformPrimalSimplex()
        self.assertEqual(simplex.Non_Basic,[5,4])
        self.assertEqual(simplex.Basic,[1, 2, 3])
        print 'Basic and Non Basic Indices return'

    def test2(self):
        print 'TEST 2---------------------------------------------- TEST 2'
        simplex = linear.LinearProgram(self.a,self.c,self.b)
        flag , value =simplex.preformPrimalSimplex()
        self.assertAlmostEqual(value,31.0)
        print 'Objective value returned by Primal simplex is Correct'

    def test3(self):
        print 'TEST 3---------------------------------------------- TEST 3'
        dual = linear.LinearProgram(self.a1,self.c1,self.b1)
        flag , value =dual.preformDualSimplex()
        self.assertAlmostEqual(value,-7)
        print 'Objective value returned by Daul simplex is Correct'

    def test4(self):
        print 'TEST 4---------------------------------------------- TEST 4'
        dual = linear.LinearProgram(self.a1,self.c1,self.b1)
        flag , value = dual.preformDualSimplex()
        self.assertEqual(dual.N.all(),self.N1.all())
        self.assertEqual(dual.B.all(),self.B1.all())
        print 'The matrices B and N returned are corrent for Dual simplex'

    def test5(self):
        print 'TEST 5---------------------------------------------- TEST 5'
        simplex = linear.LinearProgram(self.a2,self.c2,self.b2)
        flag,value =simplex.preformPrimalSimplex()
        self.assertAlmostEqual(flag,-1)
        print 'The Problem is UnBounded.'

    def test6(self):
        print 'TEST 6---------------------------------------------- TEST 6'
        simplex = linear.LinearProgram(self.a3,self.c3,self.b3)
        flag,primal_value =simplex.preformPrimalSimplex()

        # new Matrices and vectors for the Dual
        a_dual = -np.transpose(self.a3)
        c_dual = -np.transpose(self.b3)
        b_dual = -np.transpose(self.c3)

        simplex = linear.LinearProgram(a_dual,c_dual,b_dual)
        flag , dual_value = simplex.preformDualSimplex()

        self.assertAlmostEqual(primal_value,-dual_value)

        print 'Duality Thoery Checked'

    def test7(self):
        print 'TEST 7---------------------------------------------- TEST 7'
        simplex = linear.LinearProgram(self.a4,self.c4,self.b4)
        flag,value =simplex.preformPrimalSimplex()

        self.assertAlmostEqual(flag,0)

        print 'Two Phase Method Problem Passed'





if __name__ == '__main__':
    unittest.main()