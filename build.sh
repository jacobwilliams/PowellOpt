#!/bin/bash

# build the PowellOpt test program and documentation:

if hash FoBiS.py 2>/dev/null; then
    
    FoBiS.py build -compiler gnu -cflags '-c -O2' -s src -dbld bin
    
    if hash ford 2>/dev/null; then
        ford PowellOpt.md    
    else
	    echo "FORD not found! Install using: sudo pip install ford"
    fi 
       
else
	echo "FoBiS.py not found! Install using: sudo pip install FoBiS.py"
fi