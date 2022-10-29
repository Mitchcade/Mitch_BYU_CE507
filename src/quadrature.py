# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 13:17:07 2022

@author: Mitch
"""

import unittest
import basis
import numpy
import numpy as np
import math
import scipy
from scipy import optimize

# =============================================================================
# Get Quadrature Method
# =============================================================================
def quad(fun, num_points, method, spatial_domain):
    if method == "riemann":
        q_point,q_weight = getRiemannQuadrature(-1,1,num_points)
    elif method == "Newton":
        q_point,q_weight = getNewtonCotesPts(-1,1,num_points)
    elif method == "Gauss":
        q_point,q_weight = computeGaussLegendreQuadrature( num_points )
    A_val = 0
    for i in range(0,num_points):
        A_val += fun(q_point[i]) * q_weight[i]
    jacobian = (spatial_domain[1] - spatial_domain[0])/2 #The 2 is because I hard coded the domain from -1 to 1
    A_val *= jacobian
    return A_val


# =============================================================================
# RiemannQuadrature
# =============================================================================
def getRiemannQuadrature(xmin,xmax,num_points):
    w = numpy.zeros(num_points)
    P = numpy.zeros(num_points)
    for i in range(0,num_points):
        w[i] = (xmax-xmin)/num_points
        P[i] = (w[i]/2) + xmin + i * w[i]
    return P,w


# =============================================================================
# NewtonCotes Quadrature
# =============================================================================
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
    

# =============================================================================
# Gauss Legendre Quadrature
# =============================================================================
def computeGaussLegendreQuadrature( n ):
    M = numpy.zeros( 2*n, dtype = "double" )
    M[0] = 2.0
    x0 = numpy.linspace( -1, 1, n )
    sol = scipy.optimize.least_squares( lambda x : objFun( M, x ), x0, bounds = (-1, 1), ftol = 1e-14, xtol = 1e-14, gtol = 1e-14 )
    qp = sol.x
    w = solveLinearMomentFit( M, qp )
    return qp, w

def assembleLinearMomentFitSystem( degree, pts ):
    A = numpy.zeros( shape = ( degree + 1, len( pts ) ), dtype = "double" )
    ## YOUR CODE GOES HERE
    for i in range(0 , degree + 1 ):
        for j in range(0 , len(pts)):
            A[i,j] = basis.evalLegendreBasis1D(degree = i, variate = pts[j])
    return A

def solveLinearMomentFit( M, pts ):
    degree = len( M ) - 1
    A = assembleLinearMomentFitSystem( degree, pts )
    sol = scipy.optimize.lsq_linear( A, M )
    w = sol.x
    return w

def objFun( M, pts ):
    degree = len( M ) - 1
    A = assembleLinearMomentFitSystem( degree, pts )
    w = solveLinearMomentFit( M, pts )
    ## YOUR CODE GOES HERE
    obj_val = numpy.squeeze(M - A @ w) 
    return obj_val


