    program powellopt_tests
    !!
    !!  Units tests for all the optimization methods
    !!

    use powellopt

    implicit none

    write(*,*) ''
    write(*,*) '===================='
    write(*,*) ' bobyqa_test'
    write(*,*) '===================='
    write(*,*) ''
    call bobyqa_test()
    
    write(*,*) ''
    write(*,*) '===================='
    write(*,*) ' cobyla_test'
    write(*,*) '===================='
    write(*,*) ''
    call cobyla_test()
    
    write(*,*) ''
    write(*,*) '===================='
    write(*,*) ' lincoa_test'
    write(*,*) '===================='
    write(*,*) ''
    call lincoa_test()
    
    write(*,*) ''
    write(*,*) '===================='
    write(*,*) ' newuoa_test'
    write(*,*) '===================='
    write(*,*) ''
    call newuoa_test()
    
    write(*,*) ''
    write(*,*) '===================='
    write(*,*) ' uobyqa_test'
    write(*,*) '===================='
    write(*,*) ''
    call uobyqa_test()

    end program powellopt_tests