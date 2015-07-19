# PowellOpt
Optimization algorithms by M.J.D. Powell

# About

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
under a BSD-style license.

# Original Descriptions

## LINCOA

The Fortran version of LINCOA is attached. Its purpose is to seek the
least value of a function F of several variables subject to general linear
inequality constraints on the variables, when derivatives of F are not
available. The name LINCOA denotes LINearly Constrained Optimization
Algorithm. F is specified by the user through a subroutine called CALFUN.
The algorithm is intended to change the variables to values that are close
to a local constrained minimum of F. The user, however, should assume
responsibility for finding out if the calculations are adequate. It may be
helpful to employ several starting points in the space of the variables and
to try different values of the parameters ```NPT``` and ```RHOEND```. I 
intend to write
a paper that explains briefly the main features of the software.

LINCOA is not suitable for very large numbers of variables because no
attention is given to any sparsity. A few calculations with 1000 variables,
however, have been run successfully overnight, and the performance of LINCOA
is satisfactory usually for small numbers of variables. Several calculations
of the objective function may be required at points that do not satisfy the
linear constraints, especially if an equality constraint is expressed as
two inequalities.

The attachments in sequence are a suitable Makefile, followed by a main
program and a CALFUN routine for the "PtsinTet" problem, in order to provide
an example for testing. Then LINCOA and its six auxiliary routines, namely
LINCOB, GETACT, PRELIM, QMSTEP, TRSTEP and UPDATE, are given. Finally, the
output from the author's computer for the PtsinTet problem is listed.

In addition to providing CALFUN, the linear constraints, and an initial
vector of variables, the user has to set the values of RHOBEG, RHOEND and
NPT. After scaling the individual variables, so that the magnitudes of their
expected changes are similar, ```RHOBEG``` is the initial steplength for changes
to the variables, a reasonable choice being the mesh size of a coarse grid
search. Further, ```RHOEND``` should be suitable for a search on a very fine grid.
Typically, the final vector of variables is within distance ```10*RHOEND``` of
a local minimum. The parameter ```NPT``` specifies the number of interpolation
conditions on each quadratic model, the value ```NPT=2*N+1``` being recommended
for a start, where ```N``` is the number of variables.

The way of calling LINCOA should be clear from the PtsinTet example
and from the comments near the beginning of SUBROUTINE LINCOA. There are no
restrictions on or charges for the use of the software. I hope that the time
and effort I have spent on developing the package will be helpful to much
research and to many applications.

December 6th, 2013                    M.J.D. Powell (mjdp@cam.ac.uk)

## BOBYQA

The Fortran version of BOBYQA is attached. Its purpose is to seek
the least value of a function F of several variables, when derivatives
are not available, where F is specified by the user through a subroutine
called CALFUN. The name BOBYQA denotes Bound Approximation BY Quadratic
Approximation, the constraints being lower and upper bounds on every
variable, which can be set to huge values for unconstrained variables.
The algorithm is intended to change the variables to values that are close
to a local minimum of F. The user, however, should assume responsibility for
finding out if the calculations are satisfactory, by considering carefully
the values of F that occur. Details of the method of BOBYQA are given in
the report "[The BOBYQA algorithm for bound constrained optimization without
derivatives](http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf)", which 
can be reached from the "damtp.cam.ac.uk" home-page on
the web, by clicking on "Research at DAMTP", then on "Numerical Analysis"
and then on "Reports", the number of the report being 2009/NA06.

The attachments in sequence are a suitable Makefile, followed by a main
program and a CALFUN routine for the "Invdist2" problem, in order to provide
an example for testing. Then BOBYQA and its six auxiliary routines, namely
BOBYQB, ALTMOV, PRELIM, RESCUE, TRSBOX and UPDATE, are given. Finally, the
computed output that the author obtained for the Invdist2 problems is listed.

In addition to providing CALFUN, an initial vector of variables and
the lower and upper bounds, the user has to set the values of the parameters
```RHOBEG```, ```RHOEND``` and ```NPT```. After scaling the individual variables 
if necessary,
so that the magnitudes of their expected changes are similar, RHOBEG is the
initial steplength for changes to the variables, a reasonable choice being
the mesh size of a coarse grid search. Further, RHOEND should be suitable for
a search on a very fine grid. Typically, the software calculates a vector
of variables that is within distance ```10*RHOEND``` of a local minimum. Another
consideration is that every trial vector of variables is forced to satisfy
the lower and upper bounds, but there has to be room to make a search in all
directions. Therefore an error return occurs if the difference between the
bounds on any variable is less than ```2*RHOBEG```. The parameter NPT specifies
the number of interpolation conditions on each quadratic model, the value
```NPT=2*N+1``` being recommended for a start, where ```N``` is the number of variables.
It is often worthwhile to try other choices too, but much larger values tend
to be inefficient, because the amount of routine work of each iteration is
of magnitude ```NPT**2```, and because the achievement of adequate accuracy in some
matrix calculations becomes more difficult. Some excellent numerical results
have been found in the case ```NPT=N+6``` even with more than 100 variables.

The way of calling BOBYQA should be clear from the Invdist2 examples
and from the comments near the beginning of SUBROUTINE BOBYQA. There are no
restrictions on or charges for the use of the software. I hope that the time
and effort I have spent on developing the package will be helpful to much
research and to many applications.

January 5th, 2009                    M.J.D. Powell (mjdp@cam.ac.uk)

## NEWUOA

The Fortran version of NEWUOA is attached. Its purpose is to seek
the least value of a function F of several variables, when derivatives
are not available, where F is specified by the user through a subroutine
called CALFUN. The algorithm is intended to change the variables to values
that are close to a local minimum of F. The user, however, should assume
responsibility for finding out if the calculations are satisfactory, by
considering carefully the values of F that occur. The method is described
in the report "[The NEWUOA software for unconstrained optimization without
derivatives](http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2004_08.pdf)", 
which is available on the web at www.damtp.cam.ac.uk, where
you have to click on Research in DAMTP, then on Numerical Analysis and
then on Reports, the number of the report being 2004/NA08. Let ```N``` be the
number of variables. The main new feature of the method is that quadratic
models are updated using only about ```NPT=2N+1``` interpolation conditions,
the remaining freedom being taken up by minimizing the Frobenius norm of
the change to the second derivative matrix of the model.

The new software was developed from UOBYQA, which also forms quadratic
models from interpolation conditions. That method requires ```NPT=(N+1)(N+2)/2```
conditions, however, because they have to define all the parameters of the
model. The least Frobenius norm updating procedure with ```NPT=2N+1``` is usually
much more efficient when ```N``` is large, because the work of each iteration is
much less than before, and in some experiments the number of calculations
of the objective function seems to be only of magnitude ```N```.

The attachments in sequence are a suitable Makefile, followed by a main
program and a CALFUN routine for the Chebyquad problems, in order to provide
an example for testing. Then NEWUOA and its five auxiliary routines, namely
NEWUOB, BIGDEN, BIGLAG, TRSAPP and UPDATE, are given. Finally, the computed
output that the author obtained for the Chebyquad problems is listed.

The way of calling NEWUOA should be clear from the Chebyquad example
and from the comments of that subroutine. It is hoped that the software will
be helpful to much future research and to many applications. There are no
restrictions on or charges for its use. If you wish to refer to it, please
cite the published form of the DAMTP report that is mentioned above, the
full reference being "[The NEWUOA software for unconstrained minimization
without derivatives](http://link.springer.com/chapter/10.1007%2F0-387-30065-1_16)", 
in Large-Scale Nonlinear Optimization, editors G. Di
Pillo and M. Roma, Springer (2006), pages 255-297.

December 16th, 2004                    M.J.D. Powell (mjdp@cam.ac.uk)

## UOBYQA

The Fortran version of UOBYQA, written by M.J.D. Powell, is attached.
Its purpose is to seek the least value of a function F of several variables,
when derivatives are not available, where F is specified by the user through
a subroutine called CALFUN. The algorithm is intended to change the variables
to values that are close to a local minimum of F. The user, however, should
assume responsibility for finding out if the calculations are satisfactory,
by giving careful attention to the values of F that occur. The details of
the method are described in "[UOBYQA: unconstrained optimization by quadratic
approximation](http://link.springer.com/article/10.1007%2Fs101070100290)" by 
M.J.D. Powell, Mathematical Programming Series B, Volume
92, pages 555-582 (2002).

The attachments in sequence are a suitable Makefile, followed by a main
program and a CALFUN routine for the Chebyquad problems, in order to provide
an example for testing. Then UOBYQA and its three auxiliary routines, namely
UOBYQB, TRSTEP and LAGMAX, are given. Finally, the computed output that the
author obtained for the Chebyquad problems is listed.

The way of calling UOBYQA should be clear from the given example and
from the comments of that subroutine. It is hoped that the software will
be helpful to much future research and to many applications. There are no
restrictions on or charges for its use. If you wish to refer to it, please
cite the paper that is mentioned above.

## COBYLA  

Here is a single-precision Fortran implementation of the algorithm for
constrained optimization that is the subject of the report I have written on
"[A direct search optimization method that models the objective and constraint
functions by linear interpolation](http://link.springer.com/chapter/10.1007/978-94-015-8330-5_4)". 
This report has the number DAMTP 1992/NA5,
University of Cambridge, and it has been published in the proceedings of the
conference on Numerical Analysis and Optimization that was held in Oaxaca,
Mexico in January, 1992, which is the book "Advances in Optimization and
Numerical Analysis" (eds. Susana Gomez and Jean-Pierre Hennart), Kluwer
Academic Publishers (1994).

The instructions for using the Fortran code are given in the comments of
SUBROUTINE COBYLA, which is the interface between the user and the main
calculation that is done by SUBROUTINE COBYLB. There is a need for a linear
programming problem to be solved subject to a Euclidean norm trust region
constraint. Therefore SUBROUTINE TRSTLP is provided too, but you may have some
software that you prefer to use instead. These 3 subroutines are separated by
lines of hyphens below. Further, there follows the main program, the CALCFC
subroutine and the output that are appropriate to the numerical examples that
are discussed in the last section of DAMTP 1992/NA5. Please note, however,
that some cosmetic restructuring of the software has caused the given output
to differ slightly from Table 1 of the report.

There are no restrictions on the use of the software, nor do I offer any
guarantees of success. Indeed, at the time of writing this note I had applied
it only to test problems that have up to 10 variables.

Mike Powell (May 7th, 1992).

# See also
* [Original sourcecode](http://mat.uc.pt/~zhang/software.html)
