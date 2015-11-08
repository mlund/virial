#!/usr/bin/env python

"""
@date   november 2015, malmo
@author m. lund
"""

import matplotlib.pyplot as plt
import numpy as np
import warnings
import os
from math import sqrt, log, pi, exp, fabs, sinh
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.constants import N_A
from sys import exit

class RadialDistributionFunction:
  def __init__(self, filename):
    self.r, self.g = np.loadtxt(filename, usecols=(0,1), unpack=True)
    self.w = -np.log(self.g)

  def slice(self, rmin, rmax):
    m = (self.r>=rmin) & (self.r<=rmax)
    return self.r[m], self.g[m]

  def normalizeVolume(self,dim):
    self.g = self.g / self.r**(dim-1)
    self.w = -np.log(self.g)

def VirialCoefficient(r, w, mw):
  b2={}
  b2['hs'] = 2*pi/3*min(r)**3          # zero -> contact assuming hard spheres
  b2['tot'] = b2['hs'] + np.trapz( -2*pi*(np.exp(-w)-1)*r**2, r)
  b2['reduced'] = b2['tot'] / b2['hs']
  b2['hsrange'] = [0,min(r)]
  b2['range'] = [0,max(r)]
  if mw[0]>0 and mw[1]>0:
    b2['mlmol/g2'] = b2['tot']*N_A*1e-24/mw[0]/mw[1]
  return b2

# If run as main program

if __name__ == "__main__":
  import argparse
  from sys import stdout
  from argparse import RawTextHelpFormatter

  ps = argparse.ArgumentParser(
      description = 'Fit tail of RDFs to model pair potentials',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter )
  ps.add_argument('-z', type=float, nargs=2, default=[0,0], metavar=('z1','z2'), help='valencies')
  ps.add_argument('-a','--radii', type=float, nargs=2, default=[0,0], metavar=('a1','a2'), help='radii [angstrom]')
  ps.add_argument('-mw', type=float, nargs=2, default=[0,0], metavar=('mw1','mw2'), help='mol. weights [g/mol]')
  ps.add_argument('-lB','--bjerrum', type=float, default=7.1, metavar=('lB'), help='Bjerrum length [angstrom]')
  ps.add_argument('-p', '--plot', action='store_true', help='plot fitted w(r) using matplotlib' )
  ps.add_argument('-so','--shiftonly', action='store_true',
      help='do not replace tail w. model potential' )
  ps.add_argument('--norm', choices=['none','2d', '3d'], default='none',
      help='normalize w. volume element')
  ps.add_argument('-r', '--range', type=float, nargs=2, default=[0,0], metavar=('min','max'),
      help='fitting range [angstrom]')
  ps.add_argument('--pot', type=str,help='pair-potential -- either from list or user-defined')
  ps.add_argument('--guess', type=float, nargs='+', help='initial guess for fitting parameters')
  ps.add_argument('--potlist', action='store_true', help='show built-in potentials and quit')
  ps.add_argument('infile', type=str, help='two column input file with radial distribution function, g(r)' )
  ps.add_argument('outfile', type=str, help='three column output with manipulated r, w(r), g(r)' )
  args = ps.parse_args()

  # more convenient variable names
  lB = args.bjerrum
  a1,a2 = args.radii
  z1,z2 = args.z
  mw1,mw2 = args.mw

  # predefined potentials [ w(r)/kT, [guess parameters] ]
  potentiallist = {
      'dh'     : [ 'lB * z1 * z2 / r * np.exp(-r/a[0]) + a[1]',  [30., 0] ],
      'dhsinh' : [ 'lB * z1 * z2 * sinh(a[1]/a[0])**2 / r * np.exp(-r/a[0]) + a[2]',  [30., 10.0, 0] ],
      'zero'   : [ 'r*0 + a[0]', [0] ]
      }

  # print predefined potentials and quit
  if args.potlist==True:
    print "pre-defined pair-potentials:"
    for key, val in potentiallist.iteritems():
      print "%10s = %s" % (key,val[0])
    exit(0)

  # does the given pair-potential exist?
  if args.pot in potentiallist.keys():
    args.guess = potentiallist[args.pot][1]
    args.pot   = potentiallist[args.pot][0]

  # this is the pair-potential function used for fitting
  def pot(r, *a):
    exec 'ret='+args.pot
    return ret

  # load g(r)
  if not os.path.isfile( args.infile ):
    sys.exit( "Error: File "+args.infile+" does not exist." )
  rdf = RadialDistributionFunction( args.infile )

  # normalize with volume element?
  if args.norm=='2d':
    rdf.normalizeVolume(2)
  if args.norm=='3d':
    rdf.normalizeVolume(3)

  # cut out range to fit
  if args.range[1]<=args.range[0]:
    args.range = min(rdf.r), max(rdf.r)
  r,g = rdf.slice( *args.range )

  # fit pair-potential to data
  a = curve_fit( pot, r, -np.log(g), args.guess )[0] # fit data
  print "model potential:"
  print "  w(r)/kT =", args.pot
  print "        a =", a
  shift = a[-1]  # last element is always the shift

  # merge fitted data and calculated tail if needed
  if args.shiftonly is True:
    r, w = rdf.r, rdf.w - shift
  else:
    m = ( rdf.r<=args.range[0] )      # points below rmin 
    rd = rdf.r[m]
    wd = rdf.w[m]-shift              # w(r) from data points

    rm = np.arange( args.range[0], 4*args.range[1], rd[1]-rd[0] )
    wm = pot(rm, *a) - shift

    r, w = np.concatenate([rd,rm]), np.concatenate([wd,wm]) # final w(r)
    
  # virial coefficient
  B2 = VirialCoefficient( r, w, args.mw )
  print '\nvirial coefficient:'
  print '  B2hs = ', B2['hs'], 'A3 (', B2['hsrange'], ')'
  print '  B2   = ', B2['tot'], 'A3 =', B2.get('mlmol/g2', 'NaN'), 'mlmol/g2', ' (', B2['range'], ')'
  print '  B2*  = ', B2['reduced']

  # plot final w(r)
  if args.plot is True:
    import matplotlib.pyplot as plt
    plt.xlabel('$r$', fontsize=24)
    plt.ylabel('$\\beta w(r) = -\\ln g(r)$', fontsize=24)
    plt.plot( rdf.r, rdf.w-shift, 'r.' )
    plt.plot( r, w, 'g-' )
    plt.show()

  # save final pmf to disk
  if args.outfile:
    np.savetxt(args.outfile, np.transpose( (r,w, np.exp(-w)) ) )
