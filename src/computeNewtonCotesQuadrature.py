# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:48:03 2022

@author: Mitch
"""

import unittest
import basis
import numpy
import numpy as np
import math


def getNewtonCotesPts(xmin,xmax,num_points):
    q_point = np.linspace(xmin,xmax,num_points)
    if num_points == 2:
        q_weight = [1.0,1.0]
    elif num_points == 3:
        q_weight = [1/3,4/3,1/3]
    elif num_points == 4:
        q_weight = [1/4,3/4,3/4,1/4]
    elif num_points == 5:
        q_weight = [7/45, 32/45, 4/15, 32/45, 7/45]
    elif num_points == 6:
        q_weight = [19/144, 25/48, 25/72, 25/72, 25/48, 19/144]
    return q_point, q_weight
    

def computeNewtonCotesQuadrature(fun, num_points):
    q_point,q_weight = getNewtonCotesPts(-1,1,num_points)
    A_val = 0
    for i in range(0,num_points):
        A_val += fun(q_point[i]) * q_weight[i]
    return A_val

class Test_computeNewtonCotesQuadrature( unittest.TestCase ):
    def test_integrate_constant_one( self ):
        constant_one = lambda x : 1 * x**0
        for degree in range( 1, 6 ):
            num_points = degree + 1
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = constant_one, num_points = num_points ), second = 2.0, delta = 1e-12 )

    def test_exact_poly_int( self ):
        for degree in range( 1, 6 ):
            num_points = degree + 1
            poly_fun = lambda x : ( x + 1.0 ) ** degree
            indef_int = lambda x : ( ( x + 1 ) ** ( degree + 1) ) / ( degree + 1 )
            def_int = indef_int(1.0) - indef_int(-1.0)
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = poly_fun, num_points = num_points ), second = def_int, delta = 1e-12 )

    def test_integrate_sin( self ):
        sin = lambda x : math.sin(x)
        for num_points in range( 2, 7 ):
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = sin, num_points = num_points ), second = 0.0, delta = 1e-12 )

    def test_integrate_cos( self ):
        cos = lambda x : math.cos(x)
        self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = cos, num_points = 6 ), second = 2*math.sin(1), delta = 1e-4 )
        
unittest.main()        