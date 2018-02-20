----------------
| Introduction |
----------------

This package reproduces the experiments in the paper
------------------------------------------------------------------------------
Robert Peharz, Robert Gens, Franz Pernkopf and Pedro Domingos,
"On the Latent Variable Interpretation in Sum-Product Networks",
IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI),
vol. 39(10), pp. 2030-2044, 2017.
------------------------------------------------------------------------------

In particular, it provides an implementation of the EM algorithm for 
sum-product networks with Gaussians as input distributions, and an 
implementation for MPE inference in (augmented) SPNs over discrete variables.

PLEASE NOTE THE ACCOMPANYING LICENSE FILE (modified BSD, 3-Clause). 
IF YOU USE THIS CODE FOR RESEARCH, PLEASE CITE THE PAPER ABOVE.


-------------------
| Getting Started |
-------------------

* Copy the content of this package to some folder, say foo. Within foo, you
  should then have the following folder structure.

  /CPP/src/
  /Data/
  /Matlab/
  /Models/
  /PoonDomingos/
  /Results/EM/
  /Results/MPE/

* There are C++ source files in /CPP/src/ and Matlab script files in 
  /Matlab/ -- the other folders are empty so far. The experiments are actually
  performed in C++, and Matlab is used only as script language and for 
  plotting and displaying results. So Matlab is not really necessary and can 
  be replaced with something else, e.g. Octave. In this case, probably some 
  minor changes to the scripts need to be done. We ran our experiments using 
  Matlab2015.

* To compile the C++ binaries, go to the folder /CPP/ and type 'make'. 
  This should produce /CPP/bin/ containing the binaries 'MPEsynth' and 
  'trainEM', and /CPP/obj/ containing the object files.


------------------------------
| Running the EM experiments |
------------------------------

* Download the SPN package by Hoifung Poon and Pedro Domingos 
  http://alchemy.cs.washington.edu/spn/ and unzip it into folder
  /PoonDomingos/. You should then have the paths 
  
  /PoonDomingos/code/     (not required here, but can be left there)
  /PoonDomingos/data/
  /PoonDomingos/results/
  
  You can also unzip it to some other place; in this case adapt 
  'PoonDomingosRelPath' in /Matlab/setPaths.m accordingly.

* Run 'convertData.m'. 
  This converts the datasets by Poon & Domingos into a format used here and 
  stores them in folder /Data/.
  
* Run 'runExperimentsEM.m'.
  This will run EM for 7*4*3*103 = 8652 configurations (combinations of 
  parameters weights/means/sigmas (7), original parameters and 3 random 
  initializations (4), no missing data/33% missing data/66% missing data (3), 
  and 103 datasets; see paper), and in total 8652 * 30 = 259560 EM iterations.
  If using only a single worker, this will take several weeks. However, we use
  a simple pseudo-locking mechanism to enable several workers sharing the 
  same file system: The result for every training configuration is stored in a 
  separate file in the folder /Results/EM/. When a result file is already there, 
  a worker simply skips training for the corresponding configuration. When a 
  certain result file is not yet there, the worker immediately creates a dummy
  place holder, signaling to other workers to skip this configuration.
  If something happens and some of your workers crash, they will leave these 
  dummy file in the results folder. Call /Matlab/deleteDummyFiles.m to delete 
  these.  
  
  Alternatively, to reduce running time, you might want to reduce the number 
  of training configurations. E.g. when setting 
  
  MAKEMISSING_range = [0.0];
  UPDATE_range = [0,0,1;...
    0,1,1;...
    1,0,1;...
    1,1,1];
  
  EM will only run for complete data and parameter configurations involving 
  the sum-weights. You could also reduce the number of iterations, e.g. set 
  numIter = 10, as EM converges pretty nicely within 10 iterations.

* Run 'plotLLs.m' after EM training has been completed. This reproduces Figure
  9 in the paper (EM training curves). Play around with 'updateFlags', 
  'includeNumRandomInits', 'includeOriginalInits', 'includeMAKEMISSING' in 
  order to plot training curves for various other settings.
  
* Run 'checkMonotonicityEM.m' to check monotonicity of EM on the training set
  (up to numerical errors after EM has converged).

  
-------------------------------
| Running the MPE experiments |
-------------------------------

* Run 'runExperimentsMPE.m'. This script just calls the C++ binary 'MPEsynth'
  which runs the actual experiments and writes results to /Results/MPE/.
  After that, it loads the results and displays them in tables.
