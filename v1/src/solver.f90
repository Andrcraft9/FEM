

module solver
    implicit none

    ! Maximum size of matrix and the maximum number of non-zero entries
        integer, parameter :: maxn = 100000 
    integer, parameter :: maxnz = 1000000

    ! Work arrays to keep ILU factors and 8 vectors 
    integer, parameter :: MaxWi = maxnz + 2*maxn + 1
    integer, parameter :: MaxWr = maxnz + 8*maxn
    real*8 :: rW(MaxWr)
    integer :: iW(MaxWi)

    ! BiCGStab data
    external matvec, prevec0
    integer :: ITER, INFO, NUNIT
    real*8 :: RESID

        ! ILU0 data
        integer :: ierr, ipaLU, ipjLU, ipjU, ipiw

        ! External routines from BLAS library
        real*8 :: ddot
        external ddot, dcopy    

        ! Local variables
        integer :: imatvec(1), iprevec(1), ipBCG
        real*8  :: resinit

    contains

    subroutine set_precond(N, IA, JA, A)
        implicit none

        integer*4 :: N
        integer*4, intent(in) :: IA(:), JA(:)
        real*8, intent(in) :: A(:)

        integer*4 :: nz

        print *, size(IA), size(JA), size(A)
        nz = IA(N+1)
        ipaLU = 1
        ipBCG = ipaLU + nz
        ipjU  = 1
        ipjLU = ipjU + N + 1
        ipiw  = ipjLU + nz ! work array of length n

        call ilu0(N, A, JA, IA, rW(ipaLU),iW(ipjLU),iW(ipjU),iW(ipiw),ierr)

        if (ierr .ne. 0) then
            write(*, '(A,I7)') 'initialization of ilu0 failed, ierr =', ierr
        endif

        return
    end subroutine

    subroutine solve_gmres(N, IA, JA, A, RHS, u)
        implicit none
        
        integer*4 :: N
        integer*4, intent(in) :: IA(:), JA(:)
        real*8, intent(in) :: A(:)
        real*8, intent(in) :: RHS(:)
        real*8, intent(in out) :: u(:)

!  set initial guess to 0
            call dcopy(N, 0d0, 0, u, 1)

!		print *, size(RHS)
!		print *, RHS

!  compute initial residual norm
           resinit = ddot(N, RHS, 1, RHS, 1)
          resinit = dsqrt(resinit)
             if (resinit .eq. 0d0) then
              write(*, '(A)') 'rhs=0, nothing to solve!'
            endif

! ======================================================================
!  STAGE 4: iterative solution
! ======================================================================
         ITER = 1000             ! max number of iterations	
             RESID = 1d-8 * resinit  ! threshold for \|RESID\|
            INFO  = 0               ! no troubles on input
         NUNIT = 6               ! output channel
           iprevec(1) = N          ! single entry required: system size 
            imatvec(1) = N          ! single entry required: system size 

! BiConjugate Gradient Stabilized method    
        call slpbcgs(                  &
             prevec0, iprevec, iW,rW,   &
             matvec,  imatvec, IA,JA,A, &
             rW(ipBCG), N, 8,           &
             N, RHS, u,                   &
             ITER, RESID,               &
             INFO, NUNIT)

        if (INFO .ne. 0) then
                write(*, '(A)') 'BiCGStab failed'        
           endif


        return
    end subroutine




end module
