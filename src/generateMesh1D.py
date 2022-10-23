# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 15:35:28 2022

@author: Mitch
"""

import unittest
import sympy
import math
import numpy

def generateMesh1D(xmin, xmax, num_elems, degree):
    num_nodes = (degree+1)+((num_elems-1)*degree)
    node_coords = numpy.linspace(xmin, xmax,num_nodes)
    ien_array = numpy.zeros((num_elems,degree+1))
    for i in range(num_elems):
        if i == 0:
            start_node = 0
        else:
            start_node = ien_array[i-1,-1]
        ien_array[i,:] = numpy.arange(start_node,start_node + degree+1)
    print(ien_array)
    return node_coords, ien_array

class Test_generateMesh1D( unittest.TestCase ):
    def test_make_1_linear_elem( self ):
        gold_node_coords = numpy.array( [ 0.0, 1.0 ] )
        gold_ien_array = numpy.array( [ [ 0, 1 ] ], dtype = int )
        node_coords, ien_array = generateMesh1D( xmin = 0.0, xmax = 1.0, num_elems = 1, degree = 1 )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertTrue( numpy.array_equiv( ien_array, gold_ien_array ) )

    def test_make_4_linear_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.25, 0.5, 0.75, 1.0 ] )
        gold_ien_array = numpy.array( [ [ 0, 1 ], [ 1, 2 ], [ 2, 3 ], [ 3, 4 ] ], dtype = int )
        node_coords, ien_array = generateMesh1D( xmin = 0.0, xmax = 1.0, num_elems = 4, degree = 1 )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = numpy.ndarray )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertTrue( numpy.array_equiv( ien_array, gold_ien_array ) )

    def test_make_4_quadratic_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0 ] )
        gold_ien_array = numpy.array( [ [ 0, 1, 2 ], [ 2, 3, 4 ], [ 4, 5, 6 ], [ 6, 7, 8 ] ], dtype = int )
        node_coords, ien_array = generateMesh1D( xmin = 0.0, xmax = 1.0, num_elems = 4, degree = 2 )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = numpy.ndarray )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertTrue( numpy.array_equiv( ien_array, gold_ien_array ) )


unittest.main()

