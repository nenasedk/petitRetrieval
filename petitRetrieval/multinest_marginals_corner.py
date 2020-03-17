#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
__doc__ = """
Script that does default visualizations (marginal plots, 1-d and 2-d).
Author: Johannes Buchner (C) 2013-2019
"""
import numpy
from numpy import exp, log
import matplotlib.pyplot as plt
import sys, os
import json
import pymultinest
import corner

if len(sys.argv) != 2:
	sys.stderr.write("""SYNOPSIS: %s <output-root> 
	output-root: 	Where the output of a MultiNest run has been written to. 
	            	Example: chains/1-
%s""" % (sys.argv[0], __doc__))
	sys.exit(1)

prefix = sys.argv[1]
print('model "%s"' % prefix)
if not os.path.exists(prefix + 'params.json'):
	sys.stderr.write("""Expected the file %sparams.json with the parameter names.
For example, for a three-dimensional problem:
["Redshift $z$", "my parameter 2", "A"]
%s""" % (sys.argv[1], __doc__))
	sys.exit(2)
parameters = json.load(open(prefix + 'params.json'))
n_params = len(parameters)

a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
s = a.get_stats()

json.dump(s, open(prefix + 'stats.json', 'w'), indent=4)
#truths = [None,0.3,250.,3.4,None,None,4.0,None,-15,-3.3,-3.5,-2.5,-14,-9.,None,None]
truths = [-0.523,900.0,1.26,0.01,4.25,-2.,-7.,-3.,-2.52,-5.0,-5.0,-numpy.inf,-8.0,-numpy.inf]
widths = []
print('  marginal likelihood:')
print('    ln Z = %.1f +- %.1f' % (s['global evidence'], s['global evidence error']))
print('  parameters:')
for p, m in zip(parameters, s['marginals']):
    lo, hi = m['1sigma']
    med = m['median']
    sigma = (hi - lo) / 2
    widths.append(sigma)
    if sigma == 0:
        i = 3
    else:
        i = max(0, int(-numpy.floor(numpy.log10(sigma))) + 1)
        fmt = '%%.%df' % i
        fmts = '\t'.join(['    %-15s' + fmt + " +- " + fmt])
        print(fmts % (p, med, sigma))

print('creating marginal plot ...')
data = a.get_data()[:,2:]
weights = a.get_data()[:,0]
arange = numpy.ones(14,dtype = numpy.float64) * 0.9999995
#mask = weights.cumsum() > 1e-5
mask = weights > 1e-4

corner.corner(data[mask,:], weights=weights[mask], bins = 30,
	      labels=parameters, show_titles=True, 
              range=arange, color = 'tab:blue',
              plot_contours = True,
              truths = truths, truth_color = 'r',levels = [0.5,0.95,0.99])
plt.savefig(prefix + 'corner.pdf')
plt.savefig(prefix + 'corner.png')
plt.close()
