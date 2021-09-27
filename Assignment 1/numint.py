"""NUI Galway CT5132/CT5148 Programming and Tools for AI (James McDermott)

Skeleton/solution for Assignment 1: Numerical Integration

By writing my name below and submitting this file, I/we declare that
all additions to the provided skeleton file are my/our own work, and
that I/we have not seen any work on this assignment by another
student/group.

Student name(s): Philip O' Connor
Student ID(s): 21249304

"""

import numpy as np
import sympy
import itertools
import math


def numint_py(f, a, b, n):
    """Numerical integration. For a function f, calculate the definite
    integral of f from a to b by approximating with n "slices" and the
    "left" scheme. This function must use pure Python, no Numpy.

    >>> round(numint_py(math.sin, 0, 1, 100), 5)
    0.45549
    >>> round(numint_py(lambda x: x, 0, 1, 100000), 3 )
    0.5
    >>> round(numint_py(lambda x: x, 0, 1, 6), 5)
    0.41667
    >>> round(numint_py(lambda x: 1, 0, 1, 100), 5)
    1.0
    >>> round(numint_py(lambda x: -1, 0, 1, 100), 5)
    -1.0
    >>> round(numint_py(math.exp, 1, 2, 100), 5)
    4.64746

    """
    A = 0
    w = (b - a) / float(n) # width of one slice
    # STUDENTS ADD CODE FROM HERE TO END OF FUNCTION
    integral = 0                                  # initialise the intergral variable
    x=a                                           # Start at with the initial position of the integral 
    for i in range(n):                            # For each of slice of the integral except the final slice (following left rectangular scheme)
            integral= (f(x)*w)+ integral          # Calculate the area of the new slice and add it to the previous slices
            x+=w                                  # Increment x by slice size w
    return(integral)



def numint(f, a, b, n, scheme="left"):
    """Numerical integration. For a function f, calculate the definite
    integral of f from a to b by approximating with n "slices" and the
    given scheme. This function should use Numpy, and no for-loop. Eg
    np.linspace() will be useful.
    
    >>> round(numint(np.sin, 0, 1, 100, 'left'), 5)
    0.45549
    >>> round(numint(lambda x: np.ones_like(x), 0, 1, 100, 'left'), 5)
    1.0
    >>> round(numint(np.exp, 1, 2, 100, 'left'), 5)
    4.64746
    >>> round(numint(np.exp, 1, 2, 100, 'midpoint'), 5)
    4.67075
    >>> round(numint(np.sin, 0, 1, 100, 'midpoint'), 5)
    0.4597
    >>> round(numint(np.exp, 1, 2, 100, 'right'), 5)
    4.69417

    """
    # STUDENTS ADD CODE FROM HERE TO END OF FUNCTION
    step =(b - a) / float(n) 
    xrange = np.linspace(a,b-step,n) # Create the linear range for values to be used for x
    #step = abs(a-xrange[1])     # Determine the step size in the linear space

    if scheme =='left':        
        return (f(xrange)*step).sum()     # Multiply f(x) by the step size to give the area per slice and sum it together excluding the final point as this is the left rectangular scheme
    
    elif scheme == 'midpoint':
        xrange2 = np.linspace(a+step,b,n) # Create a new linear space with the step size offset
        midpoint = (xrange+xrange2)/2          # Sum the two ranges together elementwise and divide by two to get the mean value
        return (f(midpoint)*step).sum()   # Multiply f(x) by the step size to give the area per slice and sum it together
    

    elif scheme == 'right': 
        xrange_right = np.linspace(a+step,b,n)
        return (f(xrange_right)*step).sum()       #  Multiply f(x) by the step size to give the area per slice and sum it together excluding the first point as this is the right rectangular schem

def true_integral(fstr, a, b):
    """Using Sympy, calculate the definite integral of f from a to b and
    return as a float. Here fstr is an expression in x, as a str. It
    should use eg "np.sin" for the sin function.

    This function is quite tricky, so you are not expected to
    understand it or change it! However, you should understand how to
    use it. See the doctest example.

    >>> true_integral("np.sin(x)", 0, 2 * np.pi)
    0.0
    >>> true_integral("x**2", 0, 1)
    0.3333333333333333
    """
    x = sympy.symbols("x")
    # make fsym, a Sympy expression in x, now using eg "sympy.sin"
    fsym = eval(fstr.replace("np", "sympy")) 
    A = sympy.integrate(fsym, (x, a, b)) # definite integral
    A = float(A.evalf()) # convert to float
    return A

def numint_err(fstr, a, b, n, scheme):
    """For a given function fstr and bounds a, b, evaluate the error
    achieved by numerical integration on n points with the given
    scheme. Return the true value (given by true_integral),
    absolute error, and relative error, as a tuple.

    Notice that the absolute error and relative error must both be
    positive.

    Notice that the relative error will be infinity when the true
    value is zero. None of the examples in our assignment will have a
    true value of zero.

    >>> print("%.4f %.4f %.4f" % numint_err("x**2", 0, 1, 10, 'left'))
    0.3333 0.0483 0.1450
    >>> print("%.4f %.4f %.4f" % numint_err("-x**2", 0, 1, 10, 'left'))
    -0.3333 0.0483 0.1450
    >>> print("%.4f %.4f %.4f" % numint_err("x**2", 0, 1, 10, 'left'))
    0.3333 0.0483 0.1450

    """
    f = eval("lambda x: " + fstr) # f is a Python function
    A = true_integral(fstr, a, b)
    # STUDENTS ADD CODE FROM HERE TO END OF FUNCTION
    integ = numint(f, a,b, n, scheme)
    
    try:
        return(A,abs(A-integ), abs((A-integ)/A ))                   # Calculations for absolute and relative error
    except ZeroDivisionError:
        print("Curve defined has 0 area - cannot calculate relative error")

def make_table(f_ab_s, ns, schemes):
    """For each function f with associated bounds (a, b), and each value
    of n and each scheme, calculate the absolute and relative error of
    numerical integration and print out one line of a table. This
    function doesn't need to return anything, just print. Each
    function and bounds will be a tuple (f, a, b), so the argument
    f_ab_s is a list of tuples.

    Hint: use print() with the format string
    "%s,%.2f,%.2f,%d,%s,%.4g,%.4g,%.4g". Hint 2: consider itertools.

    >>> make_table([("x**2", 0, 1), ("np.sin(x)", 0, 1)], [10, 100], ['left', 'midpoint'])
    x**2,0.00,1.00,10,left,0.3333,0.04833,0.145
    x**2,0.00,1.00,10,midpoint,0.3333,0.0008333,0.0025
    x**2,0.00,1.00,100,left,0.3333,0.004983,0.01495
    x**2,0.00,1.00,100,midpoint,0.3333,8.333e-06,2.5e-05
    np.sin(x),0.00,1.00,10,left,0.4597,0.04246,0.09236
    np.sin(x),0.00,1.00,10,midpoint,0.4597,0.0001916,0.0004168
    np.sin(x),0.00,1.00,100,left,0.4597,0.004211,0.009161
    np.sin(x),0.00,1.00,100,midpoint,0.4597,1.915e-06,4.167e-06
    
    """
    
    # STUDENTS ADD CODE FROM HERE TO END OF FUNCTION
    for (fstr, a, b), n, scheme in itertools.product(f_ab_s, ns, schemes):    # Unpack the inputs given and assign to their respective variables.
        error = numint_err(fstr, a, b, n, scheme)
        print(f"{fstr:s},{a:0.2f},{b:0.2f},{n:d},{scheme:s},{error[0]:.4g},{error[1]:0.4g},{error[2]:0.4g}") # Print the output in the format requested

def main():
    """Call make_table() as specified in the pdf."""
    # STUDENTS ADD CODE FROM HERE TO END OF FUNCTION
    make_table([("np.cos(x)", 0, math.pi), ("np.sin(2*x)", 0, 1), ('np.exp(x)', 0, 1)], [10, 100, 1000], ['left', 'midpoint'])

"""

INTERPRETATION: STUDENTS ADD TEXT HERE TO INTERPRET main() results.

Results from main() show that the midpoint solutions consistently gave the closest aproximations of the true integral across all f(). 
With all functions ran, and minimum relative errors compared for each function, midpoint is at least 3 orders of magnitude lower than the left method with the largest delta being 13 orders of magniude lower.
As n increases for sin(x) and exp(x) there is a non-linear decrease in observed absolute error. Cos(x) actually has an increase in absolute error at n=100 with the midpoint method, this increase in error declines at n=1000. 


"""



    
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
    main()

