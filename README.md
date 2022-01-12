# mito-nuclear-speciation

The simulation can be run from the file **FortranMitoNuclear.f90**. The parameters for the simulation are set in the file **input.in**. The files **seed.in** and **forbiden-sites.in** are also necessary for the simulation. The code is prepared to be run in a loop. All output files will be stored in an individual folder for each realization. The outputs are:

   - pop-new.dat  
   _saves the simulation parameters and genomes (nuclear and mitochondrial) and positions of all individuals of the population at the final generation_
   
   - ext-sizes.dat  
   _records the time step each species disappeared, the abundance when it disappeared and a marker for the class of event_
   
   - global_abundances.dat  
   _records abundances of all species at all generations_

   - global_parents.dat     
   _records the mother species of all species at all generations_
   
   - number0.dat  
   _contains the number of species at each generation_

Other Fortran files to be run are:

   - dist-intra.f90   
   _computes the nuclear and mitochondrial genetic similarities within each species_

   - dist-inter.f90   
   _calculates nuclear and mitochondrial genetic distances between pair of individuals of different species_
   
   - count-mut.f90    
   _counts the number of mutations in the nuclear and mitochondrial DNA's_
  
The plots of the manuscript were based on an ensemble of 50 runs for each value of mito-nuclear selection strenght. The data analysis made in Jupyter notebooks (Python) are contained in the files:

   - sp_hyb_ext_rates.ipynb   
   _calculation of the rates of diversification (Fig. 2)_

   - SAD_lfspan.ipynb   
   _Species Abundances Distributions and distribution of species' lifespan and fates (Fig. 3 and Fig. 4)_
   
   - gen_diversity.ipynb  and   gen_div_mito-supp.ipynb   
   _Nuclear genetic distances within and between species (Fig. 5) and mitochondrial genetic distances (Supplementary Material)_
   
   - subs_rate.ipynb    
   _Histograms of the number of substituions in the nuclear and mitochondrial genomes (Supplementary Material)_
  
  REFERENCE: _"Mito-nuclear selection induces a trade-off between species ecological dominance and evolutionary lifespan"_ Débora Princepe, Marcus A. M. de Aguiar, Joshua B. Plotkin. arXiv:2111.04866 (2021) https://arxiv.org/abs/2111.04866
