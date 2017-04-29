Copyright (C) 2017 by Thomas A. Knox; this document is separate from the software in this repository and is not covered by the GNU Affero General Public License version 3 that is used for the software in this repository.

# kernelyze-base

_Numerically-optimal rank reduction of two-variable functions_

KernelyzeBase is a library, written in Fortran and with bindings to OCaml and C, that
computes low-rank approximations of two-variable functions ("kernels").

*K(x,y)* ~ *f<sub>1</sub>(x)* *g<sub>1</sub>(y)* + *f<sub>2</sub>(x)* *g<sub>2</sub>(y)* + . . . + *f<sub>n</sub>(x)* *g<sub>n</sub>(y)*

For many practically useful two-variable functions ("sign-regular kernels"), KernelyzeBase can compute new approximations which numerically achieve the smallest possible worst-case error among all rank-*n* approximations.

The worst-case error of these new approximations is much smaller than that of existing approximations such as truncated Taylor series or truncated singular function series (though KernelyzeBase also computes these more-familiar approximations, and weighting can help singular function series perform better).

This work was originally motivated by the problem of optimally approximating the behavior of bonds as interest rates move, but it has many non-financial applications as well.

This project is licensed under the GNU Affero General Public License version 3 (see the LICENSE file at the top level of the repository).  If you are interested in this software but require a different license, please contact kernelyze@gmail.com for alternative licensing options.  I authored the code in this repository (except the portions noted as third-party), but I do not hold the copyright to the code (Kernelyze LLC holds that copyright).

# Motivating Problem

Understanding the movement of a bond's price if all interest rates go up or down by the same amount (a "parallel shift" of rates)
boils down to an analysis of the kernel *exp( - s t )*, where *s* is the amount of the shift and *t* is time.  Since this is
a sign-regular kernel, KernelyzeBase can determine the (numerically best possible) functions *f<sub>1</sub>*, *g<sub>1</sub>*, . . ., *f<sub>n</sub>*,*g<sub>n</sub>*
to use in an approximation of the form *f<sub>1</sub>(s)* *g<sub>1</sub>(t)* + *f<sub>2</sub>(s)* *g<sub>2</sub>(t)* + . . . + *f<sub>n</sub>(s)* *g<sub>n</sub>(t)*.

A bond's value in shift scenarios is an integral (more specifically, a weighted sum) of the kernel *exp( - s t )*, where the sum is over the *t* variable and the weights are present values of the payments using the base (unshifted) interest rate curve.  Taking the same weighted sum over *t* of the numerically-optimal approximation of *exp( - s t )* then provides an approximation for the corresponding bond.  The numerically-optimal approximation is the best starting point for this exercise in that it minimizes (numerically) the worst-case error over shifts *s* in the user's given range for the hardest bond to approximate whose payments occur in the range of *t* given by the user.

Importantly, the approximate integrals (approximate bonds) are all linear combinations of the same *n* functions, which leads naturally via the solution of a linear system to hedging and immunization portfolios: using a
rank-*n* approximation, any bond can be approximated by at most *n - 1* other bonds.  Because the approximation works well over a range of rates, this hedging portfolio also works well over a range of rates (rather than
deteriorating in quality as duration and convexity hedging portfolios will in larger rates moves).

For sufficiently small shifts, nothing does better than a truncated Taylor series (which leads to duration and convexity approximations and their higer-order cousins).  However, the new approximation here is much better than the traditional duration or convexity approaches for even moderate movements in rates, and it adapts easily to problems in which rates curves move in different ways (steepening or flattening, yield curve twists, etc.).

Popular option pricing formulae (such as the Black (1976) formula and the normal option pricing formula discovered much earlier by Bachelier (1900)) are also sign-regular as functions of the forward and strike, so they are also subject to the new numerically-optimal approximation provided here.

Another important example of a sign-regular function of two variables is the Gaussian (or normal) probability density function, where the two variables are the mean and the usual argument of the density.  In fact, all single-parameter exponential-family densities are sign-regular as functions of their canonical parameters and sufficient statistics.  The non-central *t*-distribution also qualifies.

# Structure

The library is separated into five distinct components:
1. The foundation of the library is the set of files in the directory KernelyzeBase.  
2. The files in KernelyzeBaseTests run a battery of tests for this fundamental functionality.
    * The test_output.txt files in the x64\Release and x64\Debug directories beneath KernelyzeBaseTests show the output of this battery of tests.
3. KernelyzeBaseBindC contains the code that binds the KernelyzeBase Fortran library to C.
4. KernelyzeBaseBindCTests contains files that test the C bindings.
    * These tests simply output to stdout.
5. kernlcaml contains C stubs and OCaml code that make KernelyzeBase accessible from OCaml.
    * There is a simple set of tests for the OCaml bindings as well.

# What Are Approximation Errors Like in the New Method?

While Taylor series and, to a lesser degree, singular function series tend to have their largest errors at
the corners of a two-variable function's rectangular domain, the new approximations provided here are focused
on worst-case error.  It is intuitive that their error surfaces take on the same value at each relative
optimum of the error; otherwise, given a "nice enough" setting, parameters could be adjusted to decrease 
the magnitude of the largest error while increasing the magnitudes of smaller errors, which would 
decrease the overall worst-case error.

# Mathematical Results

The new numerically-optimal approximations provided by KernelyzeBase are tailored to minimize worst-case error.
This implies that they also provide numerically-optimal approximations to integrals and weighted sums over
one of the two variables of a sign-regular two-variable function. More precisely, they achieve a calculated lower 
bound on the worst-case error for the hardest integral to approximate, given a bound on the total variation of
the measure that is used for the integration (which translates to a bound on the sum of the absolute values of
the weights if the integral is simply a weighted sum, so that the integration measure is a linear combination 
of delta functions).

A bound (and ideally, the best possible bound) on worst-case error is important for connecting approximations
of two-variable functions to approximations of integrals over one of the two variables that may involve 
delta functions -- that is, integrals that are actually weighted sums.  

For example, a bound on the average squared error of an approximation for a two-variable function is not sufficient to bound 
the error between a weighted sum over one of the variables of the true function and a weighted sum over
one of the variables of the approximation.  This is because an approximation may achieve very good average
squared error by having very small errors for most values of its two variables, but terrible errors in a 
small area.  A weighted sum of the original two-variable function that placed weight on one of 
the small problem areas would then be far from the corresponding weighted sum for the approximate two-variable
function.  This is the oft-cited problem of "drowning in a river with an average depth of one inch."

## *n*-Widths and Approximation of Kernels

The field of mathematics which studies optimal approximation of a given (generally infinite-dimensional) space
using a space of finite dimension *n* is the study of *n*-widths.  The connection of this field to the approximation of totally positive two-variable functions is made in the 
seminal work of Charles Micchelli and Allan Pinkus.  (Total positivity is a special case of sign regularity.) This work was originally published in a series
of papers in 1977, 1978, and 1979, which are cited and nicely summarized as part of Pinkus (1985), __*n*-Widths in Approximation Theory__ 
(Springer-Verlag, New York).  Pinkus' book was a tremendously helpful inspiration in the work presented 
here, and the interested reader is advised to read it closely -- that effort will be richly repaid.

## How the New Approximation Works

Many of the insights of the Micchelli and Pinkus work are used to develop the new approximation method provided 
here, but the specific method given here does appear to be novel.  In particular, the Borsuk lower bound 
technique introduced by Tikhomirov (1960) and so nicely explained and used by Micchelli and Pinkus takes 
a new form for a kernel viewed as a mapping from the space of measures of bounded total variation to 
the space of continuous functions equipped with the sup norm.  Happily, this new form is computationally 
tractable (unlike the Micchelli and Pinkus form) using a variant of the well-known Remez (or exchange) method.

Given an *(n + 1)* vector, this modified Remez process provides a lower bound on the achievable worst-case error in reducing the rank of a sign-regular kernel.  If the vector is over *y* values let it be (*r<sub>1</sub>*, ... , *r<sub>n+1</sub>*) and let the function of two variables be *K(x,y)*.  Then take the "kernel sections" (*K(x,r<sub>1</sub>)*, ... , *K(x,r<sub>n+1</sub>)*) and execute a standard Remez process that attempts to approximate the zero function using a linear combination of these functions, where the sum of the absolute values of the coefficients in the linear combination must be one.  If the vector is over *x* values, just interchange the roles of *x* and *r* above.

The next task is to find the largest such lower bound.

The modified Remez process is used as the step of an iteration in which the Borsuk lower bound
is computed first in the *x* direction, then in the *y* direction, then in the *x* direction again, and so
on until the iteration converges (in the sense that the points at which the Remez error function achieves
its maximal magnitude converge), for a sign-regular kernel that is a function of *x* and *y*.  In the Remez
process, the points at which the error function achieves its maximal magnitude are found.  If the first Remez process uses a vector over *y* values, this results in a list of *(n + 1)* *x* values at which the Remez error function (a function of *x* in this case) achieves its greatest magnitude.  This equioscillation property is a consequence of the sign regularity of the kernel.  Take this list of *x* values as the vector (*r<sub>1</sub>*, ... , *r<sub>n+1</sub>*) of *x* values that will define the kernel sections in the next Remez process.  This next Remez process will produce a list of *(n + 1)* *y* values at which the maximal magnitude of the Remez error is achieved; use these to form the vector of *y* values that will define the kernel sections in the subsequent Remez process.

Because the coefficients in each Remez process must have absolute values that sum to one, this iteration must increase the Borsuk lower bound unless it is at a fixed point (that is, unless the iteration reaches a point at which the vectors defining the kernel sections over *x* and over *y* do not change).  To see this, suppose that the first step took sections using *y* values (*r<sub>1</sub>*, ... , *r<sub>n+1</sub>*) and obtained the *x* values (*s<sub>1</sub>*, ... , *s<sub>n+1</sub>*) where the Remez error achieved its maximal magnitude.  Now consider the Remez process in which the kernel sections are (*K(s<sub>1</sub>, y)*, ... , *K(s<sub>n+1</sub>, y)*) and the Remez process attempts to find the best approximation to the zero function (as a function of *y*) using these sections and the sum-to-one constraint on the absolute values of the coefficients.  If the *y* values obtained as error-magnitude maximizers in this Remez process are not identical to the vector (*r<sub>1</sub>*, ... , *r<sub>n+1</sub>*), then they must deliver error magnitudes greater than the Remez errors computed in the preceding step: since the coefficients on the *K(s<sub>i</sub>, y)* can
be multiplied by -1 without altering the resulting error magnitudes, take these coefficients *a<sub>i</sub>* to be signed so that the sums over *j* of *a<sub>i</sub> K(s<sub>i</sub>, r<sub>j</sub>) c<sub>j</sub>* are positive for each *j*, where the *c<sub>j</sub>* are the coefficients from the prior Remez process (in which the kernel sections were (*K(x,r<sub>1</sub>)*, ... , *K(x,r<sub>n+1</sub>)*)).  This is possible because both the coefficients and the sums over *j* of *K(s<sub>i</sub>, r<sub>j</sub>) c<sub>j</sub>* must alternate in sign (a consequence of the sign regularity of *K* and the fact that the *s<sub>i</sub>* are the *x* values at which the sum over *j* of *K(x, r<sub>j</sub>) c<sub>j</sub>* equioscillates).  Then, because the sums over *j* of *K(s<sub>i</sub>, r<sub>j</sub>) c<sub>j</sub>* equioscillate (by definition of the *s*), the sum over *i* of the sum over *j* of *a<sub>i</sub> K(s<sub>i</sub>, r<sub>j</sub>) c<sub>j</sub>* must be just the Remez error magnitude from the prior Remez process times the sum of the absolute values of the coefficients *c<sub>j</sub>*; since that sum of absolute values is one by constraint, the result is the Remez error magnitude from the prior Remez process.  This shows that it is possible to achieve the Remez error magnitude from the prior Remez process if equioscillation occurs at the *r<sub>j</sub>*; since equioscillation will occur at points chosen to maximize the magnitude of the errors given the coefficient vector, the actual error magnitude obtained will be no smaller than the prior Remez error magnitude and will only be equal to it if equioscillation occurs exactly at the *r<sub>j</sub>* (in which case the iteration has converged).

Thus, the iteration must converge to a local maximum of the Borsuk lower bound; practical experience indicates that it converges to the global maximum of the Borsuk lower bound.

Finally, the converged parameters of the iteration are used to define a rank-*(n + 1)* approximating kernel, and the rank-*n* approximating kernel is obtained from that rank-*(n + 1)* kernel by optimal approximation using the theory developed by Micchelli and Pinkus (1979) for optimal rank-reduction of a rank-*(n + 1)* kernel to a rank-*n* kernel.  To my knowledge, KernelyzeBase is the first software to actually implement that theory.

The worst-case error of this rank-*n* kernel and the converged Borsuk lower bound from the preceding iteration are within rounding error of each other for each *n* and for each sign-regular kernel that I have examined; thus the term "numerically-optimal approximation."

## Earlier Paper and Slides

Please see the directory "Papers" at the top level of the repository for a paper which studies just the exponential product kernel and looks only at square domains, along with a corresponding set of slides.  

The software here is more general than that used in the earlier paper: it is adapted to analyze an arbitrary sign-regular kernel and it accommodates general rectangular domains.  The potential asymmetry of the kernel and of the domain required a more general Borsuk lower bound iteration and a more general resulting approximation.

# Platforms Tested

The functionality here was tested with:
1. The Intel(R) Visual Fortran Compiler 16.0.4.246 on 64-bit Windows
(Visual Studio solution files and project files are included here because Intel(R) Visual Fortran has
support for Visual Studio).
2. The NAG Fortran Builder 6.1 (Build 6116) on 64-bit Windows (the project files are included in the
subdirectories suffixed with "NAG").
3. OCaml version 4.03.0 on 64-bit Windows.
3. The NAG Fortran Builder on Mac OS X (the Fortran code itself does not require any
changes, but this does require small adjustments in the test
files of KernelyzeBaseBindCTests to remove Windows-specific C elements, and the NAG
project files must be adjusted to remove backslashes in paths and point LAPACK dependencies to
the Apple-provided system libraries).

Other platforms should not present great difficulties, though older Fortran compilers may
have trouble with the modern Fortran features used by KernelyzeBase.

# Contact

Please contact kernelyze@gmail.com with any questions or comments on the software in this repository.  Kernelyze LLC holds the copyright to the code in this repository (except as noted for third-party software), which I authored.