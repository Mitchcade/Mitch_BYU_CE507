# from cmath import pi, sin
import sympy
import scipy
import numpy
import matplotlib.pyplot as plt
from scipy import integrate
import math

def taylorExpansion( fun, a, order ):
    x = list( fun.atoms( sympy.Symbol ) )[0]
    t = 0
    for i in range( 0, order + 1 ):
        df = sympy.diff( fun, x, i )
        term = ( df.subs( x, a ) / sympy.factorial( i ) ) * ( x - a )**i
        t += term
    return t


# =============================================================================
# Question 32 - sin(pi*x)
# =============================================================================

def plotTayExpansion1():
    x = sympy.symbols('x', real=True)
    fun = sympy.sin(sympy.pi * x)
    plots = sympy.plot(fun,(x,-1,1), show=False)
    N = 1000
    px = numpy.linspace( -1, 1, N )
    py = numpy.zeros( N )
    fy = numpy.zeros( N )
    fig, ax = plt.subplots()
    for degree in numpy.array([0,1,3,5,7]):
        t = taylorExpansion( fun, 0, degree )
        for i in range( 0, N ):
            fy[i] = fun.subs( x, px[i] )
            py[i] = t.subs( x, px[i] )
        ax.plot( px, py, linewidth=2.0)
    ax.plot( px, fy, linewidth=2.0)
    plt.show()
    
def plotTayExpansion2():
    x = sympy.symbols('x', real=True)
    fun = sympy.sin(sympy.exp(x))
    plots = sympy.plot(fun,(x,-1,1), show=False)
    N = 1000
    px = numpy.linspace( -1, 1, N )
    py = numpy.zeros( N )
    fy = numpy.zeros( N )
    fig, ax = plt.subplots()
    for degree in numpy.array([0,1,2,3,4]):
        t = taylorExpansion( fun, 0, degree )
        for i in range( 0, N ):
            fy[i] = fun.subs( x, px[i] )
            py[i] = t.subs( x, px[i] )
        ax.plot( px, py, linewidth=2.0)
    ax.plot( px, fy, linewidth=2.0)
    plt.show()
    
def plotTayExpansion3():
    x = sympy.symbols('x', real=True)
    fun = sympy.sin(scipy.special.erfc(x))
    plots = sympy.plot(fun,(x,-1,1), show=False)
    N = 1000
    px = numpy.linspace( -1, 1, N )
    py = numpy.zeros( N )
    fy = numpy.zeros( N )
    fig, ax = plt.subplots()
    for degree in numpy.array([0,1,3,5,7]):
        t = taylorExpansion( fun, 0, degree )
        for i in range( 0, N ):
            fy[i] = fun.subs( x, px[i] )
            py[i] = t.subs( x, px[i] )
        ax.plot( px, py, linewidth=2.0)
    ax.plot( px, fy, linewidth=2.0)
    plt.show()


    
# plotTayExpansion1()
# plotTayExpansion2()
# plotTayExpansion3()

def TayExpansionError1():
    x = sympy.symbols('x', real=True)
    fun = sympy.sin(sympy.pi * x)
    error = []
    fig, ax = plt.subplots()
    degree_list = [0, 1, 2, 3, 4, 5 ,6, 7, 8, 9, 10]
    for degree in degree_list:
        # print(degree)
        t = taylorExpansion( fun, 0, degree )
        err_fun = sympy.lambdify( x, abs( t - fun ) )
        e = scipy.integrate.quad( err_fun, -1, 1 )[0]
        error.append(e)
    ax.plot( degree_list, error, linewidth=2.0)
    plt.yscale( "log" )
    plt.show()
    
def TayExpansionError2():
    x = sympy.symbols('x', real=True)
    fun = sympy.sin(sympy.exp(x))
    error = []
    fig, ax = plt.subplots()
    degree_list = [0, 1, 2,3,4,5, 6, 7, 8, 9, 10]
    for degree in degree_list:
        # print(degree)
        t = taylorExpansion( fun, 0, degree )
        err_fun = sympy.lambdify( x, abs( t - fun ) )
        e = scipy.integrate.quad( err_fun, -1, 1 )[0]
        error.append(e)
    ax.plot( degree_list, error, linewidth=2.0)
    plt.yscale( "log" )
    plt.show()
    
def TayExpansionError3():
    x = sympy.symbols('x', real=True)
    fun = sympy.sin(sympy.functions.special.error_functions.erfc(x))
    error = []
    fig, ax = plt.subplots()
    degree_list = [0, 1, 2, 3, 4, 5 , 7, 8, 9, 10]
    for degree in degree_list:
        # print(degree)
        t = taylorExpansion( fun, 0, degree )
        err_fun = sympy.lambdify( x, abs( t - fun ) )
        e = scipy.integrate.quad( err_fun, -1, 1 )[0]
        error.append(e)
    ax.plot( degree_list, error, linewidth=2.0)
    plt.yscale( "log" )
    plt.show()
        
TayExpansionError1()
TayExpansionError2()
TayExpansionError3()