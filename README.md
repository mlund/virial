# virial
Script to fit radial distribution functions and calculate virial coefficients

# example

Take raw histogram for the distribution of distances betwee two particles, sampled in 3d
space. Firstly the data is normalized by a spherical volume element, whereafter
the tail is replaced by a fitted Yukawa potential. The resulting potential
of mean force is saved to 'wofr.dat' and plotted using matplotlib:

   $ ./virial.py --range 100 200 --model dh -z 1.0 1.0 --debye 10 --norm3d --plot gofr.dat wofr.dat


