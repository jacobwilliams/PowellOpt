The code was sent by Professor Powell to Zaikun Zhang on April 6th, 2015.  
The file "email.txt" is the original email. For more information on
TOLMIN, you might contact Professor Powell (mjdp@cam.ac.uk).

April 6th, 2015                   Zaikun Zhang (www.zhangzk.net) 


Below are the remarks from Professor Powell.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

     This very long file contains the package of Fortran subroutines that is
the double length version of TOLMIN, which minimizes a general differentiable
function of several variables subject to linear constraints.  Instructions on
the use of the software are given in the comment cards of SUBROUTINE GETMIN,
which is the subroutine that has to be called by the user.  Because all the
routines are in a single file, they have to be separated before compiling and
running the Fortran.  First Makefile is given below to simplify your task if
you are using the UNIX operating system.  It is followed by the subroutines,
separated by the characters `+++'.  You should use the file names that are
adjacent to these characters for compatibility with `Makefile'.  This package
should run after you have formed the separate subroutines, except that some
compilers will give warning messages about jumping into an IF - END IF block,
and some are even so fussy that they tell you that 1.0 is not a double
precision number.  In fact the software has solved a wide range of problems
successfully.  The given calculation is the `pentagon problem', which is
studied in the report `TOLMIN: A Fortran Package for Linearly Constrained
Optimization Calculations' by M.J.D. Powell, Report Number DAMTP 1989/NA2,
University of Cambridge, the output of the calculation being given at the
end of the file.  The method of TOLMIN is described in the paper `A tolerant
algorithm for linearly constrained optimization calculations', Mathematical
Programming B, Vol. 45, pp. 547-566 (1989). Please contact the author at
mjdp@cam.ac.uk if you require further information.

April, 1990.
