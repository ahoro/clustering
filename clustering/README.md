This small package came out of my reading David MacKay's Information Theory book
and doing the exercises therein and at the same time, learning R. MacKay provided 
some Matlab code that implements a few of the algorithms present here, but it is 
not vectorized and apparently written to be a demonstration rather than for 
practical use. It is not completely trivial to go from an algorithms in pseudocode, 
as appear in the book, to something that someone might want to use. I decided to 
organize this code into an R package for practice writing R packages and using S3
object orientation.

Some of the formulas that the code implements come from Bill Press's "Opinionated
Lessons in Statistics", specifically, use of log space (which doesn't appear in 
MacKay) and the formula for updating the covariance matrix.

This code runs pretty fast, as it is completely vectorized, with no explicit 
loops in the clustering code itself. It also uses cholesky decomposition instead 
of matrix inversion, which is O(n^2) instead of O(n^3) in the number of data 
dimensions. It can handle fairly large data sets and high dimensional data, 
however, just because it runs does not mean that this is a reasonable application
of the algorithms used.

ClusterDemo() shows graphically the process of fitting gaussian mixture models to
synthetic data. Each run will produce new data and a new model. It shows that these
techniques do not necessarily converge to the "true" parameters, but rather ones
that fit the data.

This project has not been thoroughly tested.