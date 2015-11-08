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

class PotentialBase(object):
  def __init__(self):
    self.popt = 0.0
    self.pcov = 0.0
    self.shift = 0.0
    self.range = []

  def info(self):
    print self.name+" fit:"
    print "  Range          = ", self.range, "A"
    print "  Fitted shift   = ", self.shift, "kT"
    self.more()

  def fit(self, r, w):
    self.range = [np.min(r), np.max(r)]
    self.popt, self.pcov = curve_fit( self.u, r, w, self.guess ) # fit data
    return self.popt

  def eval(self, rmin, rmax=1e3, dr=0.5):
    r = np.arange(rmin, rmax, dr)
    w = self.u(r, *self.popt)
    return r,w-self.shift

class Zero( PotentialBase ):
  def __init__(self):
    self.name="Zero"
    self.guess=[0]
  def u(self, r, shift):
    self.shift=shift
    return r*0 + shift
  def more(self): pass

class DebyeHuckelLimiting( PotentialBase ):

  def __init__(self, guess, fitdebye=True, fitradii=False, lB=7.1):
    self.name = "Debye-Huckel"
    self.guess = guess
    self.fitdebye = fitdebye
    self.fitradii = fitradii
    self.lB=lB

  def u( self, r, QQ, D, a, shift ):
    D=fabs(D)
    QQ = self.guess[0]
    if self.fitradii:
      a=fabs(a)
      if (a==0): a=1e-10
      QQ = (sinh(a/D) / (a/D))**2 * QQ
    if self.fitdebye is False:
      D = self.guess[1]
    self.shift = shift
    return self.lB*QQ / r * np.exp(-r/D) + shift 

  def more( self ):
    print "  Bjerrum length = ", self.lB, "A"
    print "  Debye length   = ", fabs(self.popt[1]), "A (fitted:", self.fitdebye, ")"
    print "  Charge product = ", self.popt[0]
    print "  Radius         = ", fabs(self.popt[2]), "A (fitted:", self.fitradii, ")"

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
  ps.add_argument('-D','--debye', type=float, default=1e20, metavar=('D'), help='Debye length [angstrom]')
  ps.add_argument('-m','--model', default='dh', choices=['dh','zero'], help='Model to fit')
  ps.add_argument('-p', '--plot', action='store_true', help='plot fitted w(r) using matplotlib' )
  ps.add_argument('-so','--shiftonly', action='store_true',
      help='do not replace tail w. model potential' )
  ps.add_argument('--norm', choices=['none','cylinder', 'sphere'], default='none',
      help='normalize w. volume element')
  ps.add_argument('-r', '--range', type=float, nargs=2, default=[0,0], metavar=('min','max'),
      help='fitting range [angstrom]')
  ps.add_argument('--fitradii', dest='fitradii', action='store_true', help='fit radius via sinh(ka)/ka')
  ps.add_argument('--no-fitdebye', dest='fitdebye', action='store_false', help='fit debye length')
  ps.add_argument('infile', type=str, help='two column input file with radial distribution function, g(r)' )
  ps.add_argument('outfile', type=str, help='three column output with manipulated r, w(r), g(r)' )
  args = ps.parse_args()
  
  if not os.path.isfile( args.infile ):
    sys.exit( "Error: File "+args.infile+" does not exist." )

  rdf = RadialDistributionFunction( args.infile )

  # normalize with volume element?
  if args.norm=='cylinder':
    rdf.normalizeVolume(2)
  if args.norm=='sphere':
    rdf.normalizeVolume(3)

  # cut out range to fit
  if args.range[1]<=args.range[0]:
    args.range = min(rdf.r), max(rdf.r)
  r,g = rdf.slice( *args.range )

  # set up chosen fitting model
  if args.model=='dh':
    a = (args.radii[0]+args.radii[1]) / 2.0
    guess = [ args.z[0]*args.z[1], args.debye, a, 0.0 ]
    model = DebyeHuckelLimiting(
        fitradii=args.fitradii, fitdebye=args.fitdebye, guess=guess, lB=args.bjerrum )
  if args.model=='zero':
    model = Zero()

  # fit and show info
  model.fit( r, -np.log(g) )
  model.info()

  # merge fitted data and calculated tail if needed
  if args.shiftonly is True:
    r, w = rdf.r, rdf.w - model.shift
  else:
    m = ( rdf.r<model.range[0] )                                     # points below rmin 
    rd,wd = rdf.r[m], rdf.w[m]-model.shift                           # w(r) from data points
    rm,wm = model.eval( rmin=model.range[0], rmax=2*model.range[1] ) # w(r) from model potential
    r, w = np.concatenate([rd,rm]), np.concatenate([wd,wm])

  # virial coefficient
  B2 = VirialCoefficient( r, w, args.mw )
  print
  print '  B2hs = ', B2['hs'], 'A3 (', B2['hsrange'], ')'
  print '  B2   = ', B2['tot'], 'A3 =', B2.get('mlmol/g2', 'NaN'), 'mlmol/g2', ' (', B2['range'], ')'
  print '  B2*  = ', B2['reduced']

  # plot final w(r)
  if args.plot is True:
    import matplotlib.pyplot as plt
    plt.xlabel('$r$', fontsize=24)
    plt.ylabel('$\\beta w(r) = -\\ln g(r)$', fontsize=24)
    plt.plot( rdf.r, rdf.w-model.shift, 'r.' )
    plt.plot( r, w, 'g-' )
    plt.show()

  # save final pmf to disk
  if args.outfile:
    np.savetxt(args.outfile, np.transpose( (r,w, np.exp(-w)) ) )
