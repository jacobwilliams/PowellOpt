#!/bin/bash

#
# Simple build script for PowellOpt.
#

DOCDIR='./doc/'         # build directory for documentation
SRCDIR='./src/'         # library source directory
TESTDIR='./src/tests/'  # tests source directory
BINDIR='./bin/'         # build directory for example
LIBDIR='./lib/'         # build directory for library
MODCODE='powellopt.f90' # library module file name
LIBOUT='powellopt.a'    # name of library
FORDMD='PowellOpt.md'   # FORD MD config file
FCOMPILER='gnu'         # Fortran compiler flag for FoBiS
FCOMPILERFLAGS='-c -O2' # Fortran compiler settings

# build the PowellOpt test program and documentation:

if hash FoBiS.py 2>/dev/null; then
    
    echo "build library..."
    FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${LIBDIR} -s ${SRCDIR} -dmod ./ -dobj ./obj -t ${MODCODE} -o ${LIBOUT} -mklib static -colors 
    
    echo "build test programs"
    FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${BINDIR} -s ${TESTDIR} -i ${LIBDIR} -dmod ./ -dobj ./obj -libs ${LIBDIR}${LIBOUT} -colors
    
    if hash ford 2>/dev/null; then
        ford "${FORDMD}"   
    else
	    echo "FORD not found! Install using: sudo pip install ford"
    fi 
       
else
	echo "FoBiS.py not found! Install using: sudo pip install FoBiS.py"
fi