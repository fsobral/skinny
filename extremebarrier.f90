!------------------------------------------------------------!
! SUBROUTINE EXBARR                                          !
!                                                            !
! The extreme barrier function.                              !
!                                                            !
!------------------------------------------------------------!

subroutine exbarr(x,f)

  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  real(8) :: f

  ! ARRAY ARGUMENTS
  real(8) :: x(:)

  intent(in ) :: x
  intent(out) :: f

  ! INTERFACES
  interface
     subroutine evalobjfun(n,x,f,flag)
       integer :: flag,n
       real(8) :: f
       real(8) :: x(n)

       intent(in ) :: n,x
       intent(out) :: f,flag
     end subroutine evalobjfun
  end interface

  interface
     subroutine evalconstr(n,x,ind,c,flag)
       integer :: flag,ind,n
       real(8) :: c
       real(8) :: x(n)

       intent(in ) :: ind,n,x
       intent(out) :: c,flag
     end subroutine evalconstr
  end interface

  ! LOCAL SCALARS
  integer :: flag,i,m,n
  real(8) :: c,infty,pp

  ! LOCAL ARRAYS
  real(8),pointer :: gamma(:)

  infty = engGetInfinity()
  m     =        engGetM()
  n     =        engGetN()
  gamma =>   engGetGamma()


  ! Box constraints
  ! Preliminar tests show that relaxing the box constraints
  ! increases the number of problems solved.
  c = engVerifyBounds(n,x)
  if ( c .gt. 0.0D0 ) then
     f = infty
     return
  end if

  pp = c ** 2.0D0

  do i = 1,m

     numceval = numceval + 1
     call evalconstr(n,x,i,c,flag)

     if ( flag .ne. 0 ) then
        call engSetFlag(flag)
        return
     end if

     if ( c - gamma(i) .gt. 0.0D0 ) then
        f = infty
        return
     end if

     pp = pp + max(0.0D0,c) ** 2

  end do

  numfeval = numfeval + 1
  call evalobjfun(n,x,f,flag)

!!$  OPEN(75,FILE='audetdennis.out',POSITION='APPEND')
!!$  write(75,*) numfeval,f
!!$  CLOSE(75)

  f = f + EXBARRHO * pp

  call engSetFlag(flag)

end subroutine exbarr

!------------------------------------------------------------!
! SUBROUTINE CALFUNB                                         !
!                                                            !
! This subroutine will be called by BOBYQA.                  !
!                                                            !
!------------------------------------------------------------!

subroutine calfunb(n,x,f)

  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: n
  real(8) :: f

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: n,x
  intent(out) :: f

  ! INTERFACES
  interface
     subroutine evalobjfun(n,x,f,flag)
       integer :: flag,n
       real(8) :: f
       real(8) :: x(n)

       intent(in ) :: n,x
       intent(out) :: f,flag
     end subroutine evalobjfun
  end interface

  interface
     subroutine evalconstr(n,x,ind,c,flag)
       integer :: flag,ind,n
       real(8) :: c
       real(8) :: x(n)

       intent(in ) :: ind,n,x
       intent(out) :: c,flag
     end subroutine evalconstr
  end interface

  ! LOCAL SCALARS
  integer :: flag,i,m
  real(8) :: c,infty,pp

  ! LOCAL ARRAYS
  real(8),pointer :: gamma(:)

  infty = engGetInfinity()
  gamma =>   engGetGamma()
  m     =        engGetM()

  ! Besides the methods which call this subroutine treat directly
  ! the box constraints, it seems that to penalize box infeasibility
  ! with respect to the original box-constraints makes the
  ! convergence faster.
  c  = engVerifyBounds(n,x)
  pp = c ** 2

  do i = 1,m

     numceval = numceval + 1
     call evalconstr(n,x,i,c,flag)

     if ( flag .ne. 0 ) then
        call engSetFlag(flag)
        return
     end if

     if ( c - gamma(i) .gt. 0.0D0) then
        f = infty
        return
     end if

     pp = pp + max(0.0D0,c) ** 2

  end do

  numfeval = numfeval + 1
  call evalobjfun(n,x,f,flag)

  f = f + EXBARRHO * pp

!  call engUpdateFmin(f)
!  OPEN(75,FILE='audetdennis-b.out',POSITION='APPEND')
!  write(75,*) numfeval,engGetFmin()
!  CLOSE(75)

  call engSetFlag(flag)

end subroutine calfunb

!------------------------------------------------------------!
! SUBROUTINE CALFUNSDS                                       !
!                                                            !
! This subroutine is the objective function called by SDS.   !
!                                                            !
!------------------------------------------------------------!

subroutine calfunsds(n,x,f)

  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: n
  real(8) :: f

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: n,x
  intent(out) :: f

  ! INTERFACES
  interface
     subroutine evalobjfun(n,x,f,flag)
       integer :: flag,n
       real(8) :: f
       real(8) :: x(n)

       intent(in ) :: n,x
       intent(out) :: f,flag
     end subroutine evalobjfun
  end interface

  ! LOCAL SCALARS
  integer :: flag
  
  call evalobjfun(n,x,f,flag)
  numfeval = numfeval + 1

  call engSetFlag(flag)

end subroutine calfunsds

!------------------------------------------------------------!
! SUBROUTINE CALCONSDS                                       !
!                                                            !
! This subroutine is the subroutine which evaluates if a     !
! given point 'x' is feasible or not. This subroutine is     !
! called by SDS.                                             !
!                                                            !
!------------------------------------------------------------!

subroutine calconsds(n,x,feasible)

  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: n
  logical :: feasible

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: n,x
  intent(out) :: feasible

  ! INTERFACES
  interface
     subroutine evalconstr(n,x,ind,c,flag)
       integer :: flag,ind,n
       real(8) :: c
       real(8) :: x(n)

       intent(in ) :: ind,n,x
       intent(out) :: c,flag
     end subroutine evalconstr
  end interface

  ! LOCAL SCALARS
  integer :: flag,i,m
  real(8) :: c
  
  ! LOCAL ARRAYS
  real(8),pointer :: gamma(:)

  gamma => engGetGamma()
  m     =      engGetM()

  feasible = .true.
  flag     = 0

  c = engVerifyBounds(n,x)

  if ( c .gt. 0.0D0 ) then
     feasible = .false.
     return
  end if

  do i = 1,m
     call evalconstr(n,x,i,c,flag)
     numceval = numceval + 1

     if ( flag .ne. 0 ) then
        feasible = .false.
        exit
     end if

     if ( c - gamma(i) .gt. 0.0D0 ) then
        feasible = .false.
        exit
     end if
  end do

  call engSetFlag(flag)

end subroutine calconsds
