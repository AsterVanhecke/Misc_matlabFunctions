This code demonstrates a multivariate bandwidth calculation from
[1] M. Kristan, A. Leonardis, D. Skoèaj, "Multivariate online Kernel Density
Estimation with Gaussian Kernels", Pattern Recognition, 2011.

The code includes three demos:
1. For demonstration of Multivariate KDE run: demoBW_Estimation.m  (it also compiles your code)
2. For demonstration of 1D KDE run: demoBW_Estimation1D.m
3. For demonstration of KDE with preclustering run: demoBW_with_preclustering
4. For visual demonstration run: visualizeEstimation.m

(You might have to compile the C code manually if the automatic compilation
 fails. In that case you'll have to run compileBWcalc.m, which is located in "\C_code\")

Reasons to use the bandwidth estimator from [1]:
* Reasonably fast computation
* Handles multivariate bandwidths
* Can use weighted data
* Generally produces good estimates of the bandwidths
* Can be calculated from a Gaussian mixture model, not only directly from the samples
* Avoids numerical evaluations and iterative computation -- the bandwidth is analytically computed (even form a GMM) under some approximations. 
 
Author: Matej Kristan (matej.kristan@fri.uni-lj.si)