#
# virial.py
#
# This script loads g(r) from disk, converts it to
# free energy, w(r), then shifts it so that the tail
# best fits a salt screened Debye-Huckel potential.
# The PMF is integrated to give the 2nd virial
# coefficient, B2.
#
# \author Mikael Lund
# \date Lund, December 2010
#

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, log, pi, exp
from scipy.optimize import curve_fit
from sys import exit

class PotentialBase(object):
  def __init__(self):
    self.popt = 0.0
    self.pcov = 0.0

  def fit(self, r, w):
    self.popt, self.pcov = curve_fit( self.u, r, w, self.guess ) # fit data
    return self.popt, self.pcov

class DebyeHuckelLimiting( PotentialBase ):

  def __init__(self, guess, fitcharge=False, fitkappa=True):
    self.guess = guess
    self.fitcharge = fitcharge
    self.fitkappa = fitkappa

  def u( self, r, lB, QQ, D, shift ):
    lB = self.guess[0]
    if self.fitcharge == False:
      QQ = self.guess[1]
    if self.fitkappa == False:
      D = self.guess[2]
    return lB*QQ / r * np.exp(-r/D) + shift 

class RadialDistributionFunction:
  def __init__(self, filename):
    self.r, self.g = np.loadtxt(infile, usecols=(0,1), unpack=True)
    self.w = -np.log(self.g)

  def slice(self, rmin, rmax):
    m = (self.r>=rmin) & (self.r<=rmax)
    return self.r[m], self.g[m]

infile      = "gofr.dat"
outfile     = "wofr.dat"
rmin        = 2.0                 # fitting interval (min)
rmax        = 3.00                # - / / - (max)
lB          = 7.                  # bjerrum length
Q           = [ -6.0, +6.8 ]      # charges of the g(r) particles
Mw          = [ 14400., 14400 ]   # g/mol
fit_debye   = True                # true if Debye length should be used to fit
guess_D     = 30.                 # debye length starting point
guess_shift = -2.0                # PMF shift starting point
guess       = [ guess_D, guess_shift ]

def f(r,D,shift):
  global Q,lB
  if (fit_debye==False):
    D=guess_D
  return lB*Q[0]*Q[1] / r * np.exp(-r/D) - shift


rdf = RadialDistributionFunction(infile)

r, g = rdf.slice(rmin, rmax)

dh=DebyeHuckelLimiting( fitcharge=False, guess=[7.1, 2.0, 10, 0.0] )

dh.fit( r, -np.log(g) )

exit(0)

popt, pcov = curve_fit(f, r[m], w[m], guess) # fit data
debye = popt[0]                              # get fitted values
shift = popt[1]                              # -//-
rfit = r[m]                                  # r array in fitting interval
wfit = f(rfit, *popt) + shift                # w array -//-
w    = w + shift                             # shift PMF
print "Loaded g(r) file      = ", infile
print "Saved w(r) file       = ", outfile
print "Particle charges      = ", Q
print "Particle weights      = ", Mw, "g/mol"
print "Fit range [rmin,rmax] = ", rmin, rmax
print "Fitted Debye length   = ", debye, "A"
print "Fitted ionic strength = ", (3.04/debye)**2*1000., "mM"
print "Fitted shift          = ", shift, "kT"

#
# 2ND VIRIAL COEFFICIENT
#
m      = (r<rmin)                    # slice out all data points below rmin 
r_dat  = r[m]
w_dat  = w[m]
inte   = -2*pi*( np.exp( -w_dat ) - 1 ) * r_dat**2
b2_dat = np.trapz(inte, r_dat)       # integrate data from contact -> rmin

infty=500.                           # assume pmf is ZERO after this point
r_dh  = np.arange( rmin, infty, 0.5 )# generate r array for rmin -> "infinity"
w_dh  = f(r_dh,debye,0)              # generate w array -//-
inte  = -2*pi*(np.exp(-w_dh)- 1)*r_dh**2
b2_dh = np.trapz(inte, r_dh )        # integrate using Debye-Huckel
b2_hc = 2 * pi / 3 *r[0]**3          # zero -> contact assuming hard spheres
b2_tot = b2_hc+b2_dat+b2_dh          # total B2 (A**3)

print "Virial coefficient (cubic angstrom):"
print "  Hard sphere  [%5d:%5d] = " % ( 0, r[0] ), b2_hc
print "  Loaded data  [%5d:%5d] = " % ( r[0], rmin ), b2_dat 
print "  Debye-Huckel [%5d:%5d] = " % ( rmin, infty ), b2_dh
print "  TOTAL        [%5d:%5d] = " % ( 0, infty ), b2_tot
print "                             = ", b2_tot*6.022e23*1e-23/Mw[0]/Mw[1], "ml*mol/g^2"
print "  Reduced, B2/B2_HS          = ", b2_tot / b2_hc

#
# SAVE FINAL PMF TO DISK
#
r_final = np.concatenate( (r_dat, r_dh) )  # assemble data and DH r arrays
w_final = np.concatenate( (w_dat, w_dh) )  #                - //- w arrays
np.savetxt(outfile, np.transpose( (r_final,w_final) ) )

plt.plot( r, w, "ro" )
plt.plot( rfit, wfit, linewidth=2.0 )
plt.xlabel("$\pi$")
plt.show()

