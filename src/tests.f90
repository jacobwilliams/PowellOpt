    program powellopt_tests
    !!
    !!  Units tests for all the optimization methods
    !!

    use bobyqa_module
    use cobyla_module
    use lincoa_module
    use newuoa_module
    use uobyqa_module

    implicit none

    call bobyqa_test()
    call cobyla_test()
    call lincoa_test()
    call newuoa_test()
    call uobyqa_test()

    end program powellopt_tests