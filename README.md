# ALM_FMM_stability

This folder contains the software used in the paper:

Stability Improvements for Fast Matrix Multiplication, C. Vermeylen and M. Van Barel, ArXiv (September 2023).

In the folder 'sols', the different polyadic decompositions (PDs) and convergence results that are mentioned in the paper are given. This includes new (more stable) PDs, as well as PDs from the literature that are used in the paper to improve the stability and to compare with.

The folder 'LM+QP' contains an implementation of the well known Levenberg-Marquardt optimization method for least-squares problems. This implementation is based on the method given in:

Methods for non-linear least squares problems, K. Madsen, H.B. Nielsen, O. Tingleff. Informatics and Mathematical Modelling, Technical University of Denmark, 2nd Edition (April 2004).

The implementation enables to add certain equality constraints discussed in the paper (by using the quadratic penalty (QP) method) by changing the fields of the struct 'opts' given to the function.

Then in the folder 'ALM', our implementation of the Augmented Lagrangian (AL) method proposed in the paper is given.

The folder 'functions' contains more general functions that are used by both the LM method and the AL method.

The folder 'scripts' contains Matlab scripts that can be used to test the different functions and algorithms that are provided.
