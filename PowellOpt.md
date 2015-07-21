project: PowellOpt
project_dir: ./src
output_dir: ./doc
project_github: https://github.com/jacobwilliams/PowellOpt
summary: PowellOpt -- Optimization algorithms by M.J.D. Powell
author: Jacob Williams
github: https://github.com/jacobwilliams
website: http://degenerateconic.com
twitter: https://twitter.com/degenerateconic
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
source: true

# Brief description

This is a collection of derivative-free optimization algorithms by M.J.D. Powell.
The package contains:

* LINCOA (LINearly Constrained Optimization Algorithm)
* BOBYQA (Bound Optimization BY Quadratic Approximation)
* NEWUOA (NEW Unconstrained Optimization Algorithm)
* UOBYQA (Unconstrained Optimization BY Quadratic Approximation)
* COBYLA (Constrained Optimization BY Linear Approximations)

The original routines were written in FORTRAN 77. They have been refactored into
modern Fortran for this package. The original sourcecode was written by Powell and
released without charges or restrictions (see below). The modifications are released 
under a [BSD-style license](https://github.com/jacobwilliams/PowellOpt/blob/master/LICENSE).
