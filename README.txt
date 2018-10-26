This repository holds software for extracting coordinates from kymographs of
unfolding DNA molecules and performing a Bayesian analysis on these
in order to infer the strength of the unfolding force and friction coefficient.

By executing the dna_unfolding_analysis command, the software will analyse the 
included kymograph "p3 and l3-ZVI Export-03_molecule_1_kymograph.tif", but the routine
may be called with parameters indicating the placement of a new kymograph along with the
pixel size and frame rate of the utilised microscope setup on the form
   dna_unfolding_analysis(kymo,pixelSize,frameRate)
which will then initiate analysis of this new kymograph.

The theoretical model used for the calculation of the likelihood function is 
described in detail in the paper found at https://arxiv.org/abs/1808.02817.

The Bayesian analysis is performed using Nested Sampling implemented as in
the https://github.com/krogjens/nestedsampling repository.

NB: The code was prepared for use with Matlab 2016a or later (version 9.0).
Compatibility is not guaranteed for older version of Matlab.
 
