# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 15:46:59 2022

@author: Mitch
"""

import sympy
import math
import numpy as np

def symLinspace( num_points ):
    dx = sympy.Rational(2) / sympy.Rational(num_points - 1)
    for i in range( 0, num_points):
        if i == 0:
            xi = [sympy.Rational(-1.0)]
        else:
            xi.append( sympy.Rational(xi[i-1] + dx))
    return xi

    
def symLagrangeBasis1D(degree,basis_idx):
    x = sympy.symbols("x")
    nodes = symLinspace(degree+1)
    val = sympy.Rational( 1.0 )
    for j in range(0,degree+1):
        if j != basis_idx:    
            val = val * ((x - nodes[j])/(nodes[basis_idx]-nodes[j]))
    val = sympy.simplify( val )
    return val


x = sympy.symbols("x")
for degree in range( 0, 6 ):
    print( f"DEGREE = {degree}" )
    for basis_idx in range( 0, degree+1):
        f = symLagrangeBasis1D( degree, basis_idx )
        fi = sympy.integrate( f, (x, -1, 1) )
        print( fi )



