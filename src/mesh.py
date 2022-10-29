# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 23:48:52 2022

@author: Mitch
"""

import unittest
import sympy
import math
import numpy
import numpy as np


def generateMeshNonUniformDegree(xmin, xmax, degree):
    node_coords = []
    ien_array = {}
    num_elems = len(degree)
    elem_boundaries = numpy.linspace(xmin, xmax, num_elems + 1)
    for elem_idx in range(num_elems):
        elem_xmin = elem_boundaries[elem_idx]
        elem_xmax = elem_boundaries[elem_idx + 1]
        elem_nodes = numpy.linspace(elem_xmin, elem_xmax, degree[elem_idx]+1)
        if elem_idx == 0:
            node_ids = numpy.arange(0, degree[elem_idx] + 1)
            node_coords.append(elem_nodes)
        else:
            start_node = ien_array[elem_idx - 1][-1]
            node_ids = numpy.arange(start_node, start_node + degree[elem_idx] + 1)
            node_coords.append(elem_nodes[1:])
        ien_array[elem_idx] = list(node_ids)
    node_coords = numpy.concatenate(node_coords)
    return node_coords, ien_array    


def generateMeshUniformDegree(xmin, xmax, num_elems, degree):
    ien_array = {}
    num_nodes = (degree+1)+((num_elems-1)*degree)
    node_coords = numpy.linspace(xmin, xmax,num_nodes)
    for i in range(num_elems):
        if i == 0:
            start_node = 0
        else:
            start_node = ien_array[i-1][-1]
        node_ids = numpy.arange(start_node,start_node + degree+1)
        ien_array[i] = list(node_ids) 
    return node_coords, ien_array

# class Test_generateMeshUniformDegree( unittest.TestCase ):
#     def test_make_1_linear_elem( self ):
#         gold_node_coords = numpy.array( [ 0.0, 1.0 ] )
#         gold_ien_array = numpy.array( [ [ 0, 1 ] ], dtype = int )
#         node_coords, ien_array = generateMeshUniformDegree( xmin = 0.0, xmax = 1.0, num_elems = 1, degree = 1 )
#         self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
#         self.assertTrue( numpy.array_equiv( ien_array, gold_ien_array ) )

#     def test_make_4_linear_elems( self ):
#         gold_node_coords = numpy.array( [ 0.0, 0.25, 0.5, 0.75, 1.0 ] )
#         gold_ien_array = numpy.array( [ [ 0, 1 ], [ 1, 2 ], [ 2, 3 ], [ 3, 4 ] ], dtype = int )
#         node_coords, ien_array = generateMeshUniformDegree( xmin = 0.0, xmax = 1.0, num_elems = 4, degree = 1 )
#         self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
#         self.assertIsInstance( obj = ien_array, cls = numpy.ndarray )
#         self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
#         self.assertTrue( numpy.array_equiv( ien_array, gold_ien_array ) )

#     def test_make_4_quadratic_elems( self ):
#         gold_node_coords = numpy.array( [ 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0 ] )
#         gold_ien_array = numpy.array( [ [ 0, 1, 2 ], [ 2, 3, 4 ], [ 4, 5, 6 ], [ 6, 7, 8 ] ], dtype = int )
#         node_coords, ien_array = generateMeshUniformDegree( xmin = 0.0, xmax = 1.0, num_elems = 4, degree = 2 )
#         self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
#         self.assertIsInstance( obj = ien_array, cls = numpy.ndarray )
#         self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
#         self.assertTrue( numpy.array_equiv( ien_array, gold_ien_array ) )

class Test_generateMeshNonUniformDegree( unittest.TestCase ):
    def test_make_1_linear_elem( self ):
        gold_node_coords = numpy.array( [ 0.0, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1 ] }
        node_coords, ien_array = generateMeshNonUniformDegree( xmin = 0.0, xmax = 1.0, degree = [ 1 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_1_quadratic_elem( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.5, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1, 2 ] }
        node_coords, ien_array = generateMeshNonUniformDegree( xmin = 0.0, xmax = 1.0, degree = [ 2 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_2_linear_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.5, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1 ], 1: [ 1, 2 ] }
        node_coords, ien_array = generateMeshNonUniformDegree( xmin = 0.0, xmax = 1.0, degree = [ 1, 1 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_2_quadratic_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.25, 0.5, 0.75, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1, 2 ], 1: [ 2, 3, 4 ] }
        node_coords, ien_array = generateMeshNonUniformDegree( xmin = 0.0, xmax = 1.0, degree = [ 2, 2 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_4_linear_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.25, 0.5, 0.75, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1 ], 1: [ 1, 2 ], 2: [ 2, 3 ], 3: [ 3, 4 ] }
        node_coords, ien_array = generateMeshNonUniformDegree( xmin = 0.0, xmax = 1.0, degree = [ 1, 1, 1, 1 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_4_quadratic_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1, 2 ], 1: [ 2, 3, 4 ], 2: [ 4, 5, 6 ], 3: [ 6, 7, 8 ] }
        node_coords, ien_array = generateMeshNonUniformDegree( xmin = 0.0, xmax = 1.0, degree = [ 2, 2, 2, 2 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_4_p_refine_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 1.0, 1.5, 2.0, (2.0 + 1.0/3.0), (2.0 + 2.0/3.0), 3.0, 3.25, 3.5, 3.75, 4.0 ] )
        gold_ien_array = { 0: [ 0, 1 ], 1: [ 1, 2, 3 ], 2: [ 3, 4, 5, 6 ], 3: [ 6, 7, 8, 9, 10 ] }
        node_coords, ien_array = generateMeshNonUniformDegree( xmin = 0.0, xmax = 4.0, degree = [ 1, 2, 3, 4 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertEqual( first = gold_ien_array, second = ien_array )

unittest.main()
