# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 15:34:21 2022

@author: Mitch
"""

import unittest
import sympy
import math
import numpy
import basis
import mesh

# =============================================================================
# Compute Solution
# =============================================================================
def computeSolution(target_fun, domain, degree):
    node_coords, ien_array = mesh.generateMeshNonUniformDegree(domain[0], domain[1], degree)
    test_solution = list()
    n = len(node_coords)
    for i in range(0,n):
        y = target_fun(node_coords[i])
        test_solution.append(y)
    test_solution = numpy.asarray(test_solution)
    return test_solution, node_coords, ien_array


# =============================================================================
# Evaluate Solation at Point
# =============================================================================
def evaluateSolutionAt(x, coeff, node_coords, ien_array, eval_basis):
    param_domain = [-1,1] 
    elem_idx = getElementIndexContainingPoint(x, node_coords, ien_array)
    elem_nodes = getElementNodes(ien_array, elem_idx)
    elem_domain = getElemDomain(ien_array, node_coords, elem_idx)
    param_coord = spatialToParamCoords(x, elem_domain, param_domain)
    sol_at_point = 0
    for n in range(len(elem_nodes)):
        curr_node = elem_nodes[n]
        sol_at_point += coeff[curr_node] * basis.evalLagrangeBasis1D(degree = len(elem_nodes) - 1, variate = param_coord, basis_idx = n)
    return sol_at_point
    

def getElementIndexContainingPoint( x, node_coords, ien_array):
    num_elems = len(ien_array)
    for elem_idx in range( 0, num_elems ):
        elem_boundary_node_ids = [ien_array[elem_idx][0], ien_array[elem_idx][-1]]
        elem_boundary_coords = [node_coords[elem_boundary_node_ids[0]], node_coords[elem_boundary_node_ids[1]]]
        if (x >= elem_boundary_coords[0]) and ( x <= elem_boundary_coords[1]):
            return elem_idx
        
def getElemDomain(ien_array, node_coords, elem_idx):
    num_nodes = len(ien_array[elem_idx])
    coord1 = node_coords[ien_array[elem_idx][0]]
    coord2 = node_coords[ien_array[elem_idx][num_nodes -1]]
    elem_domain = numpy.array([coord1, coord2])
    return elem_domain

def spatialToParamCoords(x, spatial_domain, param_domain):
    x -= min(spatial_domain)
    jacobian = (param_domain[1] - param_domain[0])/(spatial_domain[1] - spatial_domain[0])
    x = x * jacobian
    x += min(param_domain)
    return x
        
def getElementNodes(ien_array, elem_idx):
    return ien_array[elem_idx]        
        

class Test_SpatialToParamCoords( unittest.TestCase ):
    def test_biunit_to_unit( self ):
        self.assertAlmostEqual( first = spatialToParamCoords( x = -1, spatial_domain = [-1,1], param_domain = [0,1]), second = 0)
        self.assertAlmostEqual( first = spatialToParamCoords( x = 0, spatial_domain = [-1,1], param_domain = [0,1]), second = 0.5)
        self.assertAlmostEqual( first = spatialToParamCoords( x = 1, spatial_domain = [-1,1], param_domain = [0,1]), second = 1)
    def test_unit_to_biunit( self ):
        self.assertAlmostEqual( first = spatialToParamCoords( x = 0, spatial_domain = [0,1], param_domain = [-1,1]), second = -1)
        self.assertAlmostEqual( first = spatialToParamCoords( x = 0.5, spatial_domain = [0,1], param_domain = [-1,1]), second = 0)
        self.assertAlmostEqual( first = spatialToParamCoords( x = 1, spatial_domain = [0,1], param_domain = [-1,1]), second = 1)
        
class Test_evaluateSolutionAt( unittest.TestCase ):
    def test_single_linear_element( self ):
        node_coords, ien_array = mesh.generateMeshUniformDegree(-1, 1, 1, 1 )
        coeff = numpy.array( [-1.0, 1.0 ] )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = -1.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second = -1.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x =  0.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second =  0.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = +1.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second = +1.0 )
        
    def test_two_linear_elements( self ):
        node_coords, ien_array = mesh.generateMeshUniformDegree( -1, 1, 2, 1 )
        coeff = numpy.array( [ 1.0, 0.0, 1.0 ] )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = -1.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second = +1.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x =  0.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second =  0.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = +1.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second = +1.0 )
        
    def test_single_quadratic_element( self ):
        node_coords, ien_array = mesh.generateMeshUniformDegree( -1, 1, 1, 2 )
        coeff = numpy.array( [+1.0, 0.0, 1.0 ] )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = -1.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second = +1.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x =  0.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second =  0.0 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = +1.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second = +1.0 )
        
    def test_two_quadratic_elements( self ):
        node_coords, ien_array = mesh.generateMeshUniformDegree( -2, 2, 2, 2 )
        coeff = numpy.array( [ 1.0, 0.25, 0.5, 0.25, 1.0 ] )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = -2.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second = +1.00 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = -1.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second = +0.25 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x =  0.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second = +0.50 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = +1.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second = +0.25 )
        self.assertAlmostEqual( first = evaluateSolutionAt( x = +2.0, coeff = coeff, node_coords = node_coords, ien_array = ien_array, eval_basis = basis.evalLagrangeBasis1D ), second = +1.00 )
        

class Test_computeSolution( unittest.TestCase ):
    def test_single_linear_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x, domain = [-1.0, 1.0 ], degree = [1] )
        gold_solution = numpy.array( [ -1.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )

    def test_single_quad_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], degree = [2] )
        gold_solution = numpy.array( [ 1.0, 0.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )

    def test_two_linear_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], degree = [1,1] )
        gold_solution = numpy.array( [ 1.0, 0.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )

    def test_four_quad_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], degree = [1,1,1,1] )
        gold_solution = numpy.array( [ 1.0, 0.25, 0.0, 0.25, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )
        
unittest.main()        
