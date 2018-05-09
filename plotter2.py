#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 2017

@author: andrewalferman
"""

import numpy as np
import csv as csv
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
#import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('Agg')
import sys as sys
from matplotlib.ticker import NullFormatter

#plt.ioff()

def readingfunc(problem, solver):
    [ts, comptimes, sol] = [[] for i in range(3)]
    filename = 'timingdata-' + solver + '-' + problem + '.csv'
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            try:
                if row[-1] == '':
                    row = row[:-1]
                ts.append(float(row[1]))
                comptimes.append(float(row[2]))
                sol.append([float(row[i]) for i in range(3,len(row))])
            except:
                pass
    ts = np.array(ts)
    comptimes = np.array(comptimes)
    sol = np.array(sol)
    return ts, comptimes, sol

problems = ['vdp', 'oregonator', 'h2', 'grimech']
solvers = ['cvodes', 'exp4', 'exprb43', 'radau2a', 'rkc']

for i in range(10):
    plt.figure(i)
    plt.clf()
plt.close('all')

figformat = 'png'

for i, problem in enumerate(problems):
    plt.figure(i)

    f, axarr = plt.subplots(5, 2, sharex='col', figsize=(9.0, 7.5))
    #plt.gca().yaxis.set_minor_formatter(NullFormatter())
    f.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off',
                    right='off')

    plt.suptitle('Timing Plots for ' + problem)
    axarr[0,1].set_title('Solution')
    if problem == 'CSPtest':
        showsols = 4
        axarr[0,0].set_xscale('log')
        axarr[0,1].set_xscale('log')
        plt.xlabel('x Value')
    elif problem == 'Oregonator':
        showsols = 3
        # axarr[0,0].set_yscale('log')
        plt.xlabel('Time')
    else:
        showsols = 1
        plt.xlabel('x Value')
    for j, solver in enumerate(solvers):
        ts, comptimes, sol = readingfunc(problem, solver)
        if problem == 'CSPtest' or problem == 'Oregonator':
            axarr[j,1].legend(loc='upper left', bbox_to_anchor=(1,1),
                              fontsize='x-small')
        axarr[j,0].set_title(solver)
        axarr[j,0].plot(ts, comptimes)
        axarr[j,0].set_yscale('log')
        for k in range(showsols):
            lab = 'Y' + str(k+1)
            axarr[j,1].plot(ts, sol[:, k], label=lab)

    # Fine-tune figure; make subplots farther from each other.
    f.subplots_adjust(hspace=0.3)
    f.subplots_adjust(wspace=0.3)

    plt.savefig(problem + '.png', dpi=600)

plt.show()
