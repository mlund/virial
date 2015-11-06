#!/usr/bin/env python

"""
@date   november 2015, malmo
@author m. lund
"""

import matplotlib.pyplot as plt
import numpy as np
import warnings
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

  def __init__(self, guess, fitkappa=True, fitradii=False, lB=7.1):
    self.name = "Debye-Huckel"
    self.guess = guess
    self.fitkappa = fitkappa
    self.fitradii = fitradii
    self.lB=lB

  def u( self, r, QQ, D, a, shift ):
    D=fabs(D)
    QQ = self.guess[0]
    if self.fitradii:
      a=fabs(a)
      if (a==0): a=1e-10
      QQ = (sinh(a/D) / (a/D))**2 * QQ
    if self.fitkappa is False:
      D = self.guess[1]
    self.shift = shift
    return self.lB*QQ / r * np.exp(-r/D) + shift 

  def more( self ):
    print "  Bjerrum length = ", self.lB, "A"
    print "  Debye length   = ", fabs(self.popt[1]), "A (fitted:", self.fitkappa, ")"
    print "  Charge product = ", self.popt[0]
    print "  Radius         = ", fabs(self.popt[2]), "A (fitted:", self.fitradii, ")"

class RadialDistributionFunction:
  def __init__(self, filename):
    self.r, self.g = np.loadtxt(filename, usecols=(0,1), unpack=True)
    self.w = -np.log(self.g)

  def slice(self, rmin, rmax):
    m = (self.r>=rmin) & (self.r<=rmax)
    return self.r[m], self.g[m]

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
  ps.add_argument('-nb','--nobob', action='store_true', help='do not replace tail w. model potential' )
  ps.add_argument('-r', '--range', type=float, nargs=2, default=[0,0], metavar=('min','max'),
      help='fitting range [angstrom]')
  ps.add_argument('--fitradii', dest='fitradii', action='store_true', help='fit radius via sinh(ka)/ka')
  ps.add_argument('infile', type=str, help='input rdf' )
  ps.add_argument('outfile', type=str, help='output potential of mean force' )
  args = ps.parse_args()

  rdf = RadialDistributionFunction( args.infile )

  # cut out range to fit
  if args.range[1]<=args.range[0]:
    args.range = min(rdf.r), max(rdf.r)
  r,g = rdf.slice( *args.range )

  # set up chosen fitting model
  if args.model=='dh':
    a = (args.radii[0]+args.radii[1]) / 2.0
    guess = [ args.z[0]*args.z[1], args.debye, a, 0.0 ]
    model = DebyeHuckelLimiting(
        fitradii=args.fitradii, guess=guess, lB=args.bjerrum )
  if args.model=='zero':
    model = Zero()

  # do the fitting and show info
  model.fit( r, -np.log(g) )
  model.info()

  # merge fitted data and calculated tail if needed
  if args.nobob is True:
    r, w = rdf.r, rdf.w - model.shift
  else:
    m = ( rdf.r<model.range[0] )                                   # points below rmin 
    rd,wd = rdf.r[m], rdf.w[m]-model.shift                         # w(r) from data points
    rm,wm = model.eval( rmin=model.range[0], rmax=model.range[1] ) # w(r) from model potential
    r, w = np.concatenate([rd,rm]), np.concatenate([wd,wm])

  # integrate to get virial coefficient
  b2_hs  = 2*pi/3*min(r)**3          # zero -> contact assuming hard spheres
  b2_dat = np.trapz( -2*pi*(np.exp(-w)-1)*r**2, r)
  b2 = b2_hs+b2_dat

  print "\nVirial coefficient (cubic angstrom):"
  print "  Hard sphere  [%5d:%5d] = " % ( 0, min(r) ), b2_hs
  print "  Rest         [%5d:%5d] = " % ( min(r), max(r) ), b2_dat 
  print "  TOTAL        [%5d:%5d] = " % ( 0, max(r) ), b2
  if args.mw[0]>0 and args.mw[1]>0:
    print "                             = ", b2*N_A*1e-24/args.mw[0]/args.mw[1], "ml*mol/g^2"
  print "  Reduced, B2/B2_HS          = ", b2 / b2_hs

  # plot fitted w(r)
  if args.plot is True:
    import matplotlib.pyplot as plt
    plt.plot( rdf.r, rdf.w-model.shift )
    m = (rm<max(rdf.r))
    plt.plot( rm[m], wm[m] )
    plt.show()

  # save final pmf to disk
  if args.outfile:
    np.savetxt(args.outfile, np.transpose( (r,w, np.exp(-w)) ) )
