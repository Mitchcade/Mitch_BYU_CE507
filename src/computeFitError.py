# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 20:34:47 2022

@author: Mitch
"""

import scipy
import numpy
import math
import computeSolution as cP
import mesh
import basis
import matplotlib.pyplot as plt


def computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis ):
    num_elems = len(ien_array)
    domain = [ min( node_coords ), max( node_coords ) ]
    abs_err_fun = lambda x : abs( target_fun( x ) - cP.evaluateSolutionAt( x, coeff, node_coords, ien_array, eval_basis ) )
    fit_error, residual = scipy.integrate.quad( abs_err_fun, domain[0], domain[1], epsrel = 1e-12, limit = num_elems * 100 )
    return fit_error, residual



# linear, x**2, [0,1]
target_fun = lambda x : x**2
degree = 1
domain = [0,1]
eval_basis = basis.evalLagrangeBasis1D
elems_array = []
error_array = []
n = numpy.round(numpy.logspace(0, 3, 8), 0)
for i in range(0,8):
    num_elems = int(n[i])
    print(num_elems)
    degree_list = [degree] * num_elems
    elems_array.append(num_elems)
    coeff, node_coords, ien_array = cP.computeSolution(target_fun, domain, degree_list)
    fit_error, residual = computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis )
    error_array.append(fit_error)
elems_array = numpy.asarray(elems_array)
error_array = numpy.asarray(error_array)
ROC = numpy.log10(abs((error_array[-2] - error_array[-1]))) / numpy.log10(abs((elems_array[-2] - elems_array[-1])))
print(ROC)
# Plot
plt.figure(0)
plt.plot(elems_array, error_array, '--or', label = 'x^2 deg = 1 [0,1]')
plt.xlim(1, 1000)
plt.ylim(1e-11, 10)
plt.legend(loc="upper right")
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Number of Elements')
plt.ylabel('|Error|')
plt.savefig('fig_0', dpi=325)

   
# quadratic, x**3, [0,1]
target_fun = lambda x : x**3
degree = 2
domain = [0,1]
eval_basis = basis.evalLagrangeBasis1D
elems_array = []
error_array = []
n = numpy.round(numpy.logspace(0, 3, 8), 0)
for i in range(0,8):
    num_elems = int(n[i])
    degree_list = [degree] * num_elems
    elems_array.append(num_elems)
    coeff, node_coords, ien_array = cP.computeSolution(target_fun, domain, degree_list)
    fit_error, residual = computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis )
    error_array.append(fit_error)
elems_array = numpy.asarray(elems_array)
error_array = numpy.asarray(error_array)
ROC = numpy.log10(abs((error_array[-2] - error_array[-1]))) / numpy.log10(abs((elems_array[-2] - elems_array[-1])))
print(ROC)
# Plot
plt.figure(1)
plt.plot(elems_array, error_array, '--or', label = 'x^3 deg = 2 [0,1]')
plt.xlim(1, 1000)
plt.ylim(1e-11, 10)
plt.legend(loc="upper right")
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Number of Elements')
plt.ylabel('|Error|')
plt.savefig('fig_1', dpi=325)


# linear, sin(pi*x), [-1,1]
target_fun = lambda x : math.sin(math.pi*x)
degree = 1
domain = [-1,1]
eval_basis = basis.evalLagrangeBasis1D
elems_array = []
error_array = []
n = numpy.round(numpy.logspace(0, 3, 8), 0)
for i in range(0,8):
    num_elems = int(n[i])
    degree_list = [degree] * num_elems
    elems_array.append(num_elems)
    coeff, node_coords, ien_array = cP.computeSolution(target_fun, domain, degree_list)
    fit_error, residual = computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis )
    error_array.append(fit_error)
elems_array = numpy.asarray(elems_array)
error_array = numpy.asarray(error_array)
ROC = numpy.log10(abs((error_array[-2] - error_array[-1]))) / numpy.log10(abs((elems_array[-2] - elems_array[-1])))
print(ROC)
# Plot
plt.figure(2)
plt.plot(elems_array, error_array, '--or', label = 'sin(pi*x) deg = 1 [-1,1]')
plt.xlim(1, 1000)
plt.ylim(1e-11, 10)
plt.legend(loc="upper right")
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Number of Elements')
plt.ylabel('|Error|')
plt.savefig('fig_2', dpi=325)

# quadratic, sin(pi*x), [-1,1]
target_fun = lambda x : math.sin(math.pi*x)
degree = 2
domain = [-1,1]
eval_basis = basis.evalLagrangeBasis1D
elems_array = []
error_array = []
n = numpy.round(numpy.logspace(0, 3, 8), 0)
for i in range(0,8):
    num_elems = int(n[i])
    degree_list = [degree] * num_elems
    elems_array.append(num_elems)
    coeff, node_coords, ien_array = cP.computeSolution(target_fun, domain, degree_list)
    fit_error, residual = computeFitError( target_fun, coeff, node_coords, ien_array, eval_basis )
    error_array.append(fit_error)
elems_array = numpy.asarray(elems_array)
error_array = numpy.asarray(error_array)
ROC = numpy.log10(abs((error_array[-2] - error_array[-1]))) / numpy.log10(abs((elems_array[-2] - elems_array[-1])))
print(ROC)
# Plot
plt.figure(3)
plt.plot(elems_array, error_array, '--or', label = 'sin(pi*x) deg = 2 [-1,1]')
plt.xlim(1, 1000)
plt.ylim(1e-11, 10)
plt.legend(loc="upper right")
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Number of Elements')
plt.ylabel('|Error|')
plt.savefig('fig_3', dpi=325)