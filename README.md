# Virial
Script to fit the tail of radial distribution functions to model potential and to calculate
the resulting virial coefficients.

# Usage
All options can be viewed from the command line by typing `virial.py -h` which will give something like this,
~~~~
usage: virial.py [-h] [-z z1 z2] [-a a1 a2] [-mw mw1 mw2] [-lB lB] [-D D]
                 [-m {dh,zero}] [-p] [-nb] [--norm3d] [-r min max]
                 [--fitradii]
                 infile outfile

Fit tail of RDFs to model pair potentials

positional arguments:
  infile                input rdf
  outfile               output potential of mean force

optional arguments:
  -h, --help            show this help message and exit
  -z z1 z2              valencies (default: [0, 0])
  -a a1 a2, --radii a1 a2
                        radii [angstrom] (default: [0, 0])
  -mw mw1 mw2           mol. weights [g/mol] (default: [0, 0])
  -lB lB, --bjerrum lB  Bjerrum length [angstrom] (default: 7.1)
  -D D, --debye D       Debye length [angstrom] (default: 1e+20)
  -m {dh,zero}, --model {dh,zero}
                        Model to fit (default: dh)
  -p, --plot            plot fitted w(r) using matplotlib (default: False)
  -nb, --nobob          do not replace tail w. model potential (default:
                        False)
  --norm3d              normalize with spherical volume element (default:
                        False)
  -r min max, --range min max
                        fitting range [angstrom] (default: [0, 0])
  --fitradii            fit radius via sinh(ka)/ka (default: False)
~~~~

# Requirements
`numpy`, `scipy`, and `matplotlib`.

# Example

Take raw histogram for the distribution of distances betwee two particles, sampled in 3d
space. Firstly the data is normalized by a spherical volume element, whereafter
the tail is replaced by a fitted Yukawa potential. The resulting potential
of mean force is saved to 'wofr.dat' and plotted using matplotlib:

~~~~
virial.py --range 100 200 --model dh -z 1.0 1.0 --debye 10 --norm3d --plot gofr.dat wofr.dat
~~~~

The resulting outout will look something like this:

~~~~
Debye-Huckel fit:
  Range          =  [100.0, 200.0] A
  Fitted shift   =  -1.33174551427 kT
  Bjerrum length =  7.1 A
  Debye length   =  96.1282683426 A (fitted: True )
  Charge product =  1.0
  Radius         =  0.0 A (fitted: False )

Virial coefficient (cubic angstrom):
  Hard sphere  [    0:    5] =  261.799387799
  Rest         [    5:  199] =  229907.054727
  TOTAL        [    0:  199] =  230168.854114
  Reduced, B2/B2_HS          =  879.180261075
~~~~
