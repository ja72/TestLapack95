!  TestLapack95.f90 
!****************************************************************************
!
!  PROGRAM: TestLapack95
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program TestLapack95
    use BLAS95
    use LAPACK95
    use f95_precision
    use iso_fortran_env
    implicit none
    
    integer, parameter :: float = real64
    real(real64), parameter :: prec = 2**(-53d0)
    character(len=*), parameter :: fmt0 = '(*(g0,1x))'
    character(len=*), parameter :: fmtw = '(*(g12.4,1x))'    
    ! Body of TestLapack95
    print *, '*** TESTS LAPACK95 INTERFACE FOR INTEL FORTRAN 19.1 ***'
    
    call random_seed()
    
    call test_lapack95(n=1300)
    call test_error(n_start=4, iter=20, lambda=1.67e0)
    call test_eigv(n=7)
    
    
    !call bench_lapack95(n_start=48, iter=25, lambda=0.75e0)
    ! measured speeds on large arrays > 1000
    ! MKL_SEQUENTIAL = 35,000 MFlops,   50,000 MFlops
    ! MKL_PARALLEL = 200,000 MFlops, 100,000 MFlops
    
    print *,''
    stop
    contains
    
    
        
    subroutine prepare_values(n ,A, b)
    integer, intent(in) :: n
    real(float),intent(inout) :: A(n,n), b(n)
    integer :: i
    call random_number(a)
    call random_number(b)    
    do i=1,n
        a(i,n-i+1) = a(i,i)/2.0_float + 1.0_float
        b(i) = 2*b(i)-1
    end do
    end subroutine
    
    subroutine test_eigv(n)    
    integer, intent(in) :: n
    integer :: i
    real(float), allocatable :: A(:,:), w(:), E(:,:), R(:,:), D(:,:)
    print '(*(g0,1x))', '** TEST EIGENVALUES, SIZE=', n, '**'
    allocate(A(n,n))
    allocate(w(n))
    allocate(D(n,n))
    call RANDOM_NUMBER(A)
    do i=1,n
        A(i,i) = n*abs(A(i,i))
    end do
    ! make symmetric
    A = 0.5d0*(transpose(A) + A)
    
    call report_matrix('MATRIX, A=', a)
    
    ! call syevd(a, w [,jobz] [,uplo] [,info])
    E = A ! make a copy
    call syevd(e, w)
    
    call report_vector('EIGV, w=', w)
    call report_matrix('EVEC, E=', e)
    
    D = 0d0
    do i=1,n
        D(i,i) = w(i)
    end do
        
    R = MATMUL(A-D, E)
    call report_matrix('RESIDUAL, (A-D)*E=', R)
    
    end subroutine
    
    subroutine test_lapack95(n)
    integer, intent(in) :: n
    real(float), allocatable :: A(:,:), LU(:,:)
    real(float), allocatable :: b(:), x(:), xe(:)
    integer, allocatable :: ipiv(:)
    real(float) :: ee
    
    print '(*(g0,1x))', 'NO OF EQUATIONS=', n
    allocate(A(n,n))
    allocate(b(n))
    allocate(ipiv(n))
    
    call prepare_values(n, A, b)
        
    call report_matrix('COEFFICIENT MATRIX, A=', a)
    !call report_vector('CONSTANT MATRIX, b=', b)
    
    LU = A
    x = b
    call gesv(LU,x,ipiv)
    
    !call report_vector('SOLUTION, x=', x)
    call report_vectors('CONSTANT MATRIX, b=', b, 'SOLUTION, x=', x)
    
    ! ee = maxval( ABS( b-matmul(A,x)))
    
    xe = b
    call gemv(A,x, xe, 1d0, -1d0)
    ee = maxval(abs(xe))
    
    print '(*(1x,a,g12.4))','MAX RESIDUAL, r=', ee, ', r/(n*precision)=', (ee/prec/n)
    print *,''    
    call report_matrix('LU=', LU)
    if( n<=16 )then
        print fmt0,'pivot=', ipiv
    else
        print fmt0,'pivot=', ipiv(1:8), '...', ipiv(n-7:n)
    end if
    
    end subroutine
    
    subroutine test_error(n_start, iter, lambda)
    real(8), parameter :: prec = 2**(-53d0)
    integer, intent(in) :: n_start, iter
    real(4), intent(in) :: lambda
    real(float), allocatable :: A(:,:), LU(:,:)
    real(float), allocatable :: b(:), x(:), xe(:)
    integer(kind=8) :: tic, toc, rate
    integer, allocatable :: ipiv(:)
    integer :: i, n, repeat
    real(float) :: ee, time
    if(n_start==0) then 
        n = 1000
    else
        n = n_start
    end if
    print '(1x,a8,1x,a,1x,a)', 'n', 'residual', 'r/(n*prec)'
    do i=1,iter
        allocate(A(n,n))
        allocate(b(n))
        allocate(ipiv(n))
        ee  =0
        repeat = 0
        call SYSTEM_CLOCK(tic, rate)
        do
            call prepare_values(n,A,b)
            LU = A
            x = b
            call gesv(LU,x,ipiv)
            xe = b
            call gemv(A,x, xe, 1d0, -1d0)
            ee = ee + maxval(abs(xe))
            call SYSTEM_CLOCK(toc, rate)
            time = (toc-tic)/dble(rate)
            repeat = repeat + 1
            if( time>0.25 ) then
                exit
            end if
        end do
        ee = ee/repeat
        print '(1x,i8,1x,g,1x,g0.4)', n, ee, (ee/prec/n)
        
        deallocate(A,b,ipiv)
        
        n = nint(n* (1 + lambda)**(1/3e0) ) !! next step should take +20% time
    end do    
    
    end subroutine
    
    subroutine bench_lapack95(n_start, iter, lambda)
    ! bench_lapack95(n_start, iter, lamba) - runs repeated benchmarks to LAPACK95 function GESV()
    !   n_start : the starting system size. example n_start = 1000
    !   iter    : how many size iterations to run. example iter = 8
    !   lambda  : each iteration increases in size such that the time to solve
    !             is lanbda% more. For example lambda = 0.5e0 means the size
    !             increses geometrically with a ratio r=(1+lambda)**(1/3)
    !   note:   the minimum time to solve is 1 second, so the bench repeats until
    !           it takes more than one second (see repeats variable).
    integer, intent(in) :: n_start, iter
    real(4), intent(in) :: lambda
    integer :: i, n, repeat
    integer(kind=8) :: tic, toc, rate
    real(float), allocatable :: A(:,:), LU(:,:)
    real(float), allocatable :: b(:), x(:)
    real(float) :: e
    integer, allocatable :: ipiv(:), sz(:)
    real(8), allocatable :: time(:)
    
    if(n_start==0) then 
        n = 1000
    else
        n = n_start
    end if
    allocate(time(iter))
    allocate(sz(iter))
    print '(1x,a,t15,a,t30,a,t45,a)', 'size', 'time', 'MFlop', 'Error'
    do i=1, iter   
        repeat = 0
        sz(i) = n
        ! Setup matrix
        allocate(A(n,n))
        allocate(b(n))
        allocate(ipiv(n))
        
        call prepare_values(n, A, b)
        
        ! Solve matrix
        call SYSTEM_CLOCK(tic, rate)
        do 
            LU = A
            x = b
            call gesv(LU,x,ipiv)
            call SYSTEM_CLOCK(toc, rate)
            time(i) = (toc-tic)/dble(rate)
            repeat = repeat + 1
            if( time(i)>1 ) then
                exit
            end if
        end do
        time(i) = time(i)/repeat
        e = maxval( ABS( b-matmul(A,x)))
        
        deallocate(A)
        deallocate(b)
        deallocate(ipiv)
        deallocate(LU)
        deallocate(x)
        
        !print *, 'size=',sz(i),', time=', time(i), ', size^3/time=', real(sz(i))**3/time(i)
        print '(1x,i8,t15,f11.6,t30,g0.6,t45,g0.6)', sz(i), time(i), real(sz(i))**3/(1e6*time(i)), e
        
        n = nint(n* (1 + lambda)**(1/3e0) ) !! next step should take +20% time
    end do
    
    end subroutine
    
    subroutine report_matrix(text, a)
    character(len=*), parameter :: fmtw2 = '(3(g12.4,1x),a,3(g12.4,1x))'
    character(len=*), intent(in) :: text
    real(float), intent(in) :: a(:,:)
    integer :: i, n,m    
    n = size(a,1)
    m = size(a,2)
    print *,text
    if( n<=6 ) then
        do i=1,n
            print fmtw, a(i,:)
        end do    
    else
        do i=1,3
            print fmtw2, a(i,1:3), ':', a(i, n-2:n)
        end do
        print *,' ', repeat('-', 72)
        do i=n-2,n
            print fmtw2, a(i,1:3), ':', a(i, n-2:n)
        end do
    end if
    print *,''
    end subroutine
    
    subroutine report_vector(text, b)
    character(len=*), intent(in) :: text
    real(float), intent(in) :: b(:)
    integer :: i, n
    n = size(b)
    print *,text
    if( n<=6 ) then
        do i=1,n
            print fmtw, b(i)
        end do    
    else
        do i=1,3
            print fmtw, b(i)
        end do
        print *,' ', repeat('-', 11)
        do i=n-2,n
            print fmtw, b(i)
        end do
    end if
    print *,''
    end subroutine

    subroutine report_vectors(text1, x1, text2, x2)
    character(len=*), intent(in) :: text1, text2
    real(float), intent(in) :: x1(:), x2(:)
    integer :: i, n, m, k, l
    character(len=*), parameter :: &
        fmt12 = '(1x,g22.4,1x,g22.4)', &
        fmtt = '(1x,a22,1x,a22)', &
        fmt1 = '(1x,g22.4)', &
        fmt2 = '(1x,a22,1x,g22.4)'
    n = size(x1)
    m = size(x2)
    k = min(n,m)
    l = max(n,m)
    print fmtt,text1,text2
    if( n<=6 .and. m<=6 ) then
        do i=1,k
            print fmt12, x1(i), x2(i)
        end do  
        do i=k+1,n
            print fmt1, x1(i)
        end do
        do i=k+1,m
            print fmt2, '', x2(i)
        end do  
    elseif(n==m) then
        do i=1,3
            print fmt12, x1(i), x2(i)
        end do
        print fmtt, repeat('-', 11),repeat('-', 11)
        do i=n-2,n
            print fmt12, x1(i), x2(i)
        end do
    else
        do i=1,3
            print fmt12, x1(i), x2(i)
        end do
        print fmtt, repeat('-', 11),repeat('-', 11)
        do i=1,3
            print fmt12, x1(n+i-3), x2(m+i-3)
        end do        
    end if
    print *,''
    end subroutine
    
    end program TestLapack95

