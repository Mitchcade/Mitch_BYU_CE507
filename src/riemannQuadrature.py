# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 12:54:32 2022

@author: Mitch
"""

import unittest
import numpy 
import math

def getRiemannQuadrature(xmin,xmax,num_points):
    w = numpy.zeros(num_points)
    P = numpy.zeros(num_points)
    for i in range(0,num_points):
        w[i] = (xmax-xmin)/num_points
        P[i] = (w[i]/2) + xmin + i * w[i]
    return P,w

def riemannQuadrature(fun, num_points):
    q_point,q_weight = getRiemannQuadrature(-1,1,num_points)
    A_val = 0
    for i in range(0,num_points):
        A_val += fun(q_point[i]) * q_weight[i]
    return A_val
        

class Test_getRiemannQuadrature( unittest.TestCase ):
    def test_biunit_4_Point( self ):
        x_gold = [-3/4,-1/4,1/4,3/4]
        w_gold = [0.5,0.5,0.5,0.5]
        x_test,w_test = getRiemannQuadrature(-1,1,4)
        self.assertTrue(numpy.allclose(x_gold,x_test))
        self.assertTrue(numpy.allclose(w_gold,w_test))


class Test_computeRiemannQuadrature( unittest.TestCase ):
    def test_integrate_constant_one( self ):
        constant_one = lambda x : 1
        for num_points in range( 1, 100 ):
            self.assertAlmostEqual( first = riemannQuadrature( fun = constant_one, num_points = num_points ), second = 2.0, delta = 1e-12 )

    def test_integrate_linear( self ):
        linear = lambda x : x
        for num_points in range( 1, 100 ):
            self.assertAlmostEqual( first = riemannQuadrature( fun = linear, num_points = num_points ), second = 0.0, delta = 1e-12 )

    def test_integrate_quadratic( self ):
        linear = lambda x : x**2
        error = []
        for num_points in range( 1, 100 ):
            error.append( abs( (2.0 / 3.0) - riemannQuadrature( fun = linear, num_points = num_points ) ) )
        self.assertTrue( numpy.all( numpy.diff( error ) <= 0.0 ) )

    def test_integrate_sin( self ):
        sin = lambda x : math.sin(x)
        error = []
        for num_points in range( 1, 100 ):
            self.assertAlmostEqual( first = riemannQuadrature( fun = sin, num_points = num_points ), second = 0.0, delta = 1e-12 )

    def test_integrate_cos( self ):
        cos = lambda x : math.cos(x)
        error = []
        for num_points in range( 1, 100 ):
            error.append( abs( (2.0 / 3.0) - riemannQuadrature( fun = cos, num_points = num_points ) ) )
        self.assertTrue( numpy.all( numpy.diff( error ) <= 0.0 ) )
unittest.main()