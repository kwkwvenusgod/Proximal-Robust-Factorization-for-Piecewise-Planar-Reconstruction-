PARSuMi for Practical Matrix Completion and Corruption Recovery. Version 0.2

------------------------------
Update Nov 13, 2014:

Changed PARSuMi.m and PARSuMi1.m so that they deal with fully observed matrices separately (and much faster). For this special case, One partial SVD suffices so there is no need to do matrix completion in every iteration.

1. For a demo, try out the script ``test_nomissing.m''.


------------------------------
Update Jan 12 2014:

Changed the name into into PARSuMi where P stand for ``Proximal'' and ``Partially Majorized''.
Required for the convergence guarantee.

1. PARSuMi.m is the function to call. PARSuMi1.m is the version with small additional tikhonov regularization.
2. Try out the ``testscript.m'' for a demo using synthetic data.
3. Try out ``RunDinosaur.m'' for a demo using the Oxford Dinosaur dataset.


------------------------------
ARSuMi for Practical Matrix Completion Problem.
A matlab implementation. Version 0.1
Date: June 24, 2013
Author: Yu-Xiang Wang, Choon Meng Lee.

Dependencies: 

1. CVX (used only in the ADMM vesion of huber fitting, namely in the "huber regression heuristic").

2. GRASTA by Jun He https://sites.google.com/site/hejunzz/grasta is also included in this package to run benchmarking experiments.


Usage:

1. run 'testscript.m' to figure out. It's straight forward.



Further instructions to understand the code:

1. The main ARSUMi function is in 'convergentARSuMi.m' where each alternating steps are in 'Prox_MatrixCompletion.m' and 'Prox_ErrorCorrection.m'.

2. "Huber Regression" Heuristic can be turned off by commenting the 48th line in 'convergentARSuMi.m'. (it is not needed when convex initilization is good enough)

3. The APG nuclear norm minimization for our convex initilization is implemented in folder 'ExplicitRankOutlier'. 

4. Most of the function inputs and outputs are vectorized observations instead of the full matrix, although we use capital letter to denote the variable.

