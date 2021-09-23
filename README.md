# PowellOpt
Optimization algorithms by M.J.D. Powell

### About

This is a collection of derivative-free optimization algorithms by M.J.D. Powell.
The package contains:

* LINCOA (LINearly Constrained Optimization Algorithm)
* BOBYQA (Bound Optimization BY Quadratic Approximation)
* NEWUOA (NEW Unconstrained Optimization Algorithm)
* UOBYQA (Unconstrained Optimization BY Quadratic Approximation)
* COBYLA (Constrained Optimization BY Linear Approximations)

### Building

The [Fortran Package Manager](https://github.com/fortran-lang/fpm) (fpm) is a great package manager and build system for Fortran.
You can build using provided `fpm.toml`:
```bash
fpm build
```
To use `PowellOpt` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
PowellOpt = { git="https://github.com/jacobwilliams/PowellOpt.git" }
```

### License

The original routines were written in FORTRAN 77. They have been refactored into
modern Fortran for this package. The original sourcecode was written by Powell and
released without charges or restrictions (see below). The modifications are released
under a [BSD-style license](https://github.com/jacobwilliams/PowellOpt/blob/master/LICENSE).

### See also
* [Original sourcecode](http://mat.uc.pt/~zhang/software.html)