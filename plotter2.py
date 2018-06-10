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
    filename = 'AutonomousODEs/timingdata-' + solver + '-' + problem + '-ht.csv'
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

#problems = ['vdp', 'oregonator', 'csptest', 'h2', 'grimech']
problems = ['vdp', 'oregonator', 'csptest', 'h2', 'grimech']
solvers = ['cvodes', 'exp4', 'exprb43', 'radau2a', 'rkc']
printtimes = False

for i in range(10):
    plt.figure(i)
    plt.clf()
plt.close('all')

figformat = 'png'

for i, problem in enumerate(problems):
    if printtimes:
        print(problem + ':')
    plt.figure(i)

    f, axarr = plt.subplots(5, 1, sharex='col', figsize=(9.0/2, 7.5))
    #plt.gca().yaxis.set_minor_formatter(NullFormatter())
    f.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off',
                    right='off')

    #plt.suptitle('Timing Plots for ' + problem)
    #axarr[0,1].set_title('Solution')
    if problem == 'csptest':
        #showsols = 4
        axarr[0].set_xscale('log')
        #axarr[0,1].set_xscale('log')
        plt.xlabel('x Value')
        ymin = 4.5E-7
        ymax = 7.E-3
    elif problem == 'oregonator':
        #showsols = 3
        # axarr[0,0].set_yscale('log')
        plt.xlabel('Time')
        ymin = 1.E-6
        ymax = 3E-2
    elif problem == 'vdp':
        plt.xlabel('x Value')
        ymin = 9.E-7
        ymax = 1.5E-1
    elif problem == 'h2':
        plt.xlabel('Time')
        ymin = 2.E-5
        ymax = 3.E-2
    elif problem == 'grimech':
        #showsols = 1
        plt.xlabel('Time')
        ymin = 3.5E-4
        ymax = 5.E0
    for j, solver in enumerate(solvers):
        ts, comptimes, sol = readingfunc(problem, solver)
        if printtimes:
            print('Max, min times for ' + solver + ': {:.1E}, {:.1E}'.format(max(comptimes),min(comptimes)))
        else:
            axarr[j].set_title(solver)
            axarr[j].plot(ts, comptimes)
            axarr[j].set_yscale('log')
            axarr[j].set_ylim(ymin, ymax)

    if not printtimes:
        # Fine-tune figure; make subplots farther from each other.
        f.subplots_adjust(hspace=0.3)
        f.subplots_adjust(wspace=0.3)
        plt.tight_layout()

        plt.savefig(problem + '-ht.png')

if printtimes:
    sys.exit()
else:
    plt.show()
