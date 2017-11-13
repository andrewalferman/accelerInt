#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 10:44 2017

@author: andrewalferman
"""

import os as os
import numpy as np
import sys
sys.path.insert(0, '/scratch/Research/GRI_Mech_3/')
import pyjacob as pyjacob
# import scipy as sci
import datetime
import time as timer
import warnings
import csv

# from scipy.integrate import odeint
from scipy.integrate import ode

global functioncalls


def firstderiv(time, state, press):
    """Force the integrator to use the right arguments."""
    # Need to make sure that N2 is at the end of the state array
    dy = np.zeros_like(state)
    pyjacob.py_dydt(time, press, state, dy)
    global functioncalls
    functioncalls += 1
    return dy


def jacobval(time, state, press):
    """Force the integrator to use the right arguments."""
    # Need to get rid of N2 because PyJac doesn't compute it.
    # new = state[:-1]
    # print('Jacobian function called.')
    a = len(state)
    jacobian = np.zeros(a**2)
    # Obtain the jacobian from pyJac
    pyjacob.py_eval_jacobian(time, press, state, jacobian)
    jacobian = np.reshape(jacobian, (a, a))
    # Re-add the zeros back in
    # jacobian = np.insert(jacobian, a, np.zeros(a), axis=1)
    # jacobian = np.vstack((jacobian, np.zeros(a+1)))
    return jacobian


def jacvdp(x, y, eta):
    """Find the local Jacobian matrix of the Van der Pol equation."""
    return np.array([[0., 1.], [-1. - 2*y[0]*y[1]*eta, eta-eta*y[0]**2]])


def dydx(x, y, eta):
    """Find the local vector of the first derivative of the Van der Pol eqn."""
    # Unpack the y vector
    y1 = y[0]
    y2 = y[1]

    # Create dydx vector (y1', y2')
    f = np.array([y2, eta*y2 - y1 - eta*y2*y1**2.])
    # print(f)
    return f


def d2ydx2(x, y, eta):
    """Find the local vector of the 2nd derivative of the Van der Pol eqn."""
    # Unpack the y vector
    y1 = y[0]
    y2 = y[1]

    # Create vector of the second derivative
    y2prime = eta*y2 - y1 - eta*y2*y1**2.
    f = np.array([y2prime, eta*y2prime - y2 - 2*eta*y1*y2 - eta*y2prime*y1**2])
    return f


def derivcd4(vals, dx):
    """Take the derivative of a series using 4th order central differencing.

    Given a list of values at equally spaced points, returns the first
    derivative using the fourth order central difference formula, or forward/
    backward differencing at the boundaries.
    """
    deriv = []
    for i in range(2):
        deriv.append((-3 * vals[i] + 4 * vals[i + 1] - vals[i + 2]) / (2 * dx))
    for i in range(2, len(vals) - 2):
        deriv.append(((-1 * vals[i + 2]) + (8 * vals[i + 1]) -
                     (8 * vals[i - 1]) + vals[i - 2]) /
                     (12 * dx)
                     )
    for i in range((len(vals) - 2), len(vals)):
        deriv.append((3 * vals[i] - 4 * vals[i - 1] + vals[i - 2]) / 2 * dx)
    return deriv


def weightednorm(matrix, weights):
    """Weighted average norm function as defined in 1985 Shampine.

    Takes a matrix and 2 weights and returns the maximum value (divided by
    wi) of the sum of each value in each row multiplied by wj.  Needs to be
    passed either a matrix of m x n dimensions where m,n > 1, or a column
    vector.
    """
    # Unpack the parameters
    wi, wj = weights

    # Initialize a list that will be called later to obtain the maximum value
    colsums = np.zeros(len(matrix))

    # Try loop used because numpy didn't seem to like 1D arrays for the
    # weighted norm
    try:
        for i in range(len(matrix)):
            colsums += wj * np.abs(matrix[i])
        return np.max(colsums) / wi
    except TypeError:
        matrixcol = wj * np.abs(matrix)
        return np.sum(matrixcol) / wi


def stiffnessindex(xlist, solution, dfun, jfun, *args, **kwargs):
    """Determine the stiffness index across a solution vector.

    Function that uses stiffness parameters, the local Jacobian matrix,
    and a vector of the local function values to determine the local stiffness
    index as defined in 1985 Shampine.
    """
    SIparams = {'method': 2,
                'gamma': 1,
                'xi': 1,
                'order': 1,
                'tolerance': 1,
                'wi': 1,
                'wj': 1
                }

    for key, value in kwargs.items():
        SIparams[key] = value

    funcparams = []
    for arg in args:
        funcparams.append(arg)

    # Method 2 uses the weighted norm of the Jacobian, Method 1 uses the
    # spectral radius of the Jacobian.
    method = SIparams['method']
    # Stiffness index parameter values
    gamma = SIparams['gamma']
    xi = SIparams['xi']
    order = SIparams['order']
    tolerance = SIparams['tolerance']
    # Weighted norm parameters
    wi = SIparams['wi']
    wj = SIparams['wj']

    normweights = wi, wj

    # Obtain the derivative values for the derivative of order p
    dx = xlist[1] - xlist[0]
    dydxlist = []
    for i in range(len(solution)):
        dydxlist.append(dfun(xlist[i], solution[i, :], funcparams[0]))
    # Raise the derivative to the order we need it
    for i in range(order):
        dydxlist = derivcd4(dydxlist, dx)
    dydxlist = np.array(dydxlist)

    # Create a list to return for all the index values in a function
    indexlist = []

    # Figure out some of the values that will be multiplied many times, so that
    # each computation only needs to happen once.
    exponent = 1. / (order + 1)
    xiterm = ((np.abs(xi)**(-1 * exponent)) / np.abs(gamma))
    toleranceterm = tolerance**exponent

    # Actual computation of the stiffness index for the method specified.
    for i in range(len(solution)):
        jacobian = jfun(xlist[i], solution[i, :], funcparams[0])
        if method == 1:
            eigenvalues = np.linalg.eigvals(jacobian)
            index = toleranceterm *\
                np.max(np.abs(eigenvalues)) *\
                weightednorm(dydxlist[i, :], normweights)**(-1 * exponent) *\
                xiterm
        else:
            index = toleranceterm *\
                weightednorm(jacobian, normweights) *\
                weightednorm(dydxlist[i, :], normweights)**(-1 * exponent) *\
                xiterm
        indexlist.append(index)
    indexlist = np.array(indexlist)
    return indexlist  # , dydxlist


def stiffnessindicator(time, solution, jfun, *args):
    """
    Find the local stiffness indicator after calculating local solution.

    Given the value of the solution, find the stiffness indicator as defined by
    Soderlind 2013.
    """
    funcparams = []
    for arg in args:
        funcparams.append(arg)

    try:
        jacobian = jfun(time, solution, funcparams[0])
    except ValueError:
        jacobian = jfun(time[0], solution[0], funcparams[0])
    Hermitian = 0.5 * (jacobian + np.transpose(jacobian))
    eigvals = np.linalg.eigvals(Hermitian)
    return 0.5 * (min(eigvals) + max(eigvals))


def reftimescale(indicatorval, Tlen):
    """
    Find the local reference timescale for the stiffness indicator.

    Given the stiffness indicator values as defined by Soderlind 2013, finds
    the reference time scale.
    """
    if indicatorval >= 0:
        timescale = Tlen
    else:
        timescale = min(Tlen, -1/indicatorval)
    return timescale


def CEMA(xlist, solution, jfun, *args):
    """
    Find values for the chemical explosive mode analysis.

    Same thing as finding the maximum eigenvalue across the solution.
    """
    funcparams = []
    for arg in args:
        funcparams.append(arg)

    values = []
    try:
        for i in range(len(solution)):
            jacobian = jfun(xlist[i], solution[i], funcparams[0])
            values.append(max(np.linalg.eigvals(jacobian)))
    except TypeError:
        jacobian = jfun(xlist, solution, funcparams[0])
        values.append(max(np.linalg.eigvals(jacobian)))
    return values


def stiffnessratio(xlist, solution, jfun, *args):
    """
    Find values of the stiffness ratio.

    Ratio of the eigenvalue with the largest absolute value over the eigenvalue
    with the smallest absolute value. Ignores eigenvalues of zero.
    """
    funcparams = []
    for arg in args:
        funcparams.append(arg)

    values = []
    try:
        for i in range(len(solution)):
            jacobian = jfun(xlist[i], solution[i], funcparams[0])
            eigvals = np.array([abs(j) for j in np.linalg.eigvals(jacobian)
                                if j != 0])
            values.append(max(eigvals)/min(eigvals))
    except TypeError:
        jacobian = jfun(xlist, solution, funcparams[0])
        eigvals = np.array([abs(j) for j in np.linalg.eigvals(jacobian)
                            if j != 0])
        values.append(max(eigvals)/min(eigvals))
    return values


def stiffmetrics(xlist, solution, jfun, *args):
    """
    Find values of the stiffness ratio, stiffness indicator, and CEMA.
    """
    funcparams = []
    for arg in args:
        funcparams.append(arg)

    ratios, indicators, CEMs = [], [], []
    try:
        for i in range(len(solution)):
            jacobian = jfun(xlist[i], solution[i], funcparams[0])
            eigvals1 = np.array([abs(j) for j in np.linalg.eigvals(jacobian)
                                if j != 0])
            Hermitian = 0.5 * (jacobian + np.transpose(jacobian))
            eigvals2 = np.linalg.eigvals(Hermitian)
            ratios.append(max(eigvals1)/min(eigvals1))
            indicators.append(0.5 * (min(eigvals2) + max(eigvals2)))
            CEMs.append(max(np.linalg.eigvals(jacobian)))
    except TypeError:
        jacobian = jfun(xlist, solution, funcparams[0])
        eigvals1 = np.array([abs(j) for j in np.linalg.eigvals(jacobian)
                            if j != 0])
        Hermitian = 0.5 * (jacobian + np.transpose(jacobian))
        eigvals2 = np.linalg.eigvals(Hermitian)
        ratios = (max(eigvals1)/min(eigvals1))
        indicators = 0.5 * (min(eigvals2) + max(eigvals2))
        CEMs = max(np.linalg.eigvals(jacobian))
    return ratios, indicators, CEMs


def loadpasrdata():
    """Load the initial conditions from the full PaSR file."""
    print('Loading data...')
    filepath = os.path.join(os.getcwd(), 'ch4_full_pasr_data.npy')
    return np.load(filepath)


def rearrangepasr(Y, N2_pos):
    """Rearrange the PaSR data so it works with pyJac."""
    press_pos = 2
    temp_pos = 1
    arraylen = len(Y)

    Y_press = Y[press_pos]
    Y_temp = Y[temp_pos]
    Y_species = Y[3:arraylen]
    Ys = np.hstack((Y_temp, Y_species))

    # Put N2 to the last value of the mass species
    newarlen = len(Ys)
    Y_N2 = Ys[N2_pos]
    # Y_x = Ys[newarlen - 1]
    for i in range(N2_pos, newarlen - 1):
        Ys[i] = Ys[i + 1]
    Ys[newarlen - 1] = Y_N2
    if useN2:
        initcond = Ys
    else:
        initcond = Ys[:-1]
    return initcond, Y_press


# Finding the current time to time how long the simulation takes
starttime = datetime.datetime.now()
print('Start time: {}'.format(starttime))

ratios, indicators, CEMAvals = [], [], []
with open('../speciesdata-cvodes.csv', 'rt') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')  # , quotechar='|')
    for row in reader:
        state = row[1:-1]
        ratio, indicator, CEM = stiffmetrics(0, state, jacobval, 101325.0)
        ratios.append(ratio)
        indicators.append(indicator)
        CEMAvals.append(CEM)

print(np.shape(ratios))
