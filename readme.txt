This package contains the files

FortranMitoNuclear.f90
FortranTree.f90
MitoNuclear.jpynb
mycolors.py
readme.txt (this file)


In order to use the jupyter notebook you need to compile the Fortran files. To do that
open a terminal in the folder containing the files and give the following commands:

f2py -c --fcompiler=gnu95 -m FortranMitoNuclear FortranMitoNuclear.f90
and
f2py -c --fcompiler=gnu95 -m FortranTree FortranTree.f90 

These commands should generate the corresponding ".so" files that are called from python.

Make sure you have python installed, together with the packages numpy, matplotlib, jupyter and jupyterlab. 

Open the notebook with the command

jupiter-lab MitoNuclear.jpynb

and run the three modules separately. The first will run the evolution program. The
second will plot the population in space, using different colours for different species.
The third will generate the phylogenetic tree and compute the gamma, alpha and Sackin indices.
The tips of the tree will have circles coloured according to the species it represents, with
the same color code used to plot the population. Colors are listed in the file mycolors.py and can
be edited by the user.

