!------------------------------------------------------------!
! HS.F90 (Hock-Schittkowski)                                 !
!                                                            !
! This problem is the Hock-Schittkowski Collection of 306    !
! nonlinear programming test problems which appear in the    !
! book (out of print):                                       !
!                                                            !
! Test Examples for Nonlinear Programming Codes,             !
! Willi Hock, Klaus Schittkowski, Springer, Lecture Notes in !
! Economics and Mathematical Systems, Vol. 187, 1981         !
!                                                            !
! and are available at                                       !
!                                                            !
! http://www.uni-bayreuth.de/departments/math/~kschittkowski/!
! home.htm                                                   !
!                                                            !
!------------------------------------------------------------!

!------------------------------------------------------------!
! SUBROUTINE INITIALIZES                                     !
!                                                            !
! This subroutine initializes the SCALAR information of the  !
! problem.                                                   !
!                                                            !
!------------------------------------------------------------!

subroutine initializes(ne,me,flag)

  use hsdata

  implicit none

  ! SCALAR ARGUMENTS
  integer,intent(out) :: flag,me,ne

  ! INTERFACES
  interface
     subroutine CONV(MODE)
       integer,intent(in) :: MODE
     end subroutine CONV
  end interface       

  ! COMMON SCALARS
  integer :: N,NILI,NINL,NELI,NENL,NTP

  ! COMMON BLOCKS
  COMMON/L1/N,NILI,NINL,NELI,NENL
  COMMON/L8/NTP

!!$  NTP = 0
!!$
!!$  write(*,FMT=0010)
!!$  read (*,*) NTP
!!$
!!$  if ( NTP .lt. 1 .or. NTP .gt. 395 .or. &
!!$       (NTP .gt. 119 .and. NTP .lt. 201) ) then
!!$     NTP = 7
!!$  end if

  call CONV(1)

  ne = N
  me = NILI + NINL + 2 * NELI + 2 * NENL

  flag = 0

  open(75,FILE='runhs.out')
  write(75,FMT=0020) NTP,N,NILI + NINL,NELI + NENL,1.0D+20,&
       1.0D+20,1.0D+20,-1
  close(75)

! NON-EXECUTABLE STATEMENTS

0010 FORMAT('Please type the number of the problem to be solved ', &
            '(1 - 119 | 201 - 395):')
0020 FORMAT(I4,1X,I4,1X,I4,1X,I4,5X,E12.5,1X,E12.5,1X,E12.5,1X,I15)


end subroutine initializes


!------------------------------------------------------------!
! SUBROUTINE INITIALIZEA                                     !
!                                                            !
! This subroutine initializes the ARRAY information of the   !
! problem.                                                   !
!                                                            !
!------------------------------------------------------------!

subroutine initializea(n,xe,le,ue,m,flag)

  use hsdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,m,n

  ! ARRAY ARGUMENTS
  real(8) :: le(n),ue(n),xe(n)

  intent(in ) :: m,n
  intent(out) :: flag,le,ue,xe

  ! LOCAL SCALARS
  integer :: error,i

  ! COMMON ARRAYS
  logical :: INDEX1(HS_MMAX),INDEX2(HS_MMAX),LXL(HS_NMAX),LXU(HS_NMAX)
  real(8) :: X(HS_NMAX),XL(HS_NMAX),XU(HS_NMAX)

  ! COMMON BLOCKS
  COMMON/L2/X
  COMMON/L11/LXL
  COMMON/L12/LXU
  COMMON/L13/XL
  COMMON/L14/XU
  COMMON/L9/INDEX1
  COMMON/L10/INDEX2

  xe(1:n) = X(1:n)

  do i = 1,n
     if ( LXL(i) ) then
        le(i) = XL(i)
     else
        le(i) = - 1.0D+20
     end if
  end do

  do i = 1,n
     if ( LXU(i) ) then
        ue(i) = XU(i)
     else
        ue(i) = 1.0D+20
     end if
  end do

  INDEX1 = .false.
  INDEX2 = .false.

  flag   = 0

end subroutine initializea

!------------------------------------------------------------!
! SUBROUTINE FINALIZE                                        !
!                                                            !
! Finalizes the problem, sometimes printing additional       !
! information related to the specific problem.               !
! This subroutine uses the following modules:                !
!                                                            !
! EPARAM: EPSFEASI                                           !
! HSDATA: HS_MMAX,HS_NMAX                                    !
!                                                            !
!------------------------------------------------------------!

subroutine finalize(ne,x,m,maxinfeas)

  use hsdata

  implicit none

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

  ! SCALAR ARGUMENTS
  integer,intent(in) :: m,ne
  real(8),intent(in) :: maxinfeas

  ! ARRAY ARGUMENTS
  real(8),intent(in) :: x(ne)

  ! LOCAL SCALARS
  character :: optm
  integer :: i,flag
  real(8) :: diff,f,reldiff

  ! COMMON SCALARS
  integer :: N,NILI,NINL,NELI,NENL,NEX,NTP
  logical :: LEX
  real(8) :: FEX

  ! COMMON ARRAYS
  real(8) :: XEX(HS_MMAX * HS_NMAX)

  ! COMMON BLOCKS
  COMMON/L1/N,NILI,NINL,NELI,NENL
  COMMON/L8/NTP
  COMMON/L20/LEX,NEX,FEX,XEX

  if ( LEX ) then
     if ( NEX .ge. 1 ) then
        diff = HS_INFTY
        do i = 1,NEX
           diff = min(diff,maxval( &
                abs(x - XEX(ne * (i - 1) + 1:ne * i))))
        end do
     end if
  else
     diff = maxval(abs(x - XEX(1:ne)))
  end if

  call evalobjfun(ne,x,f,flag)

  reldiff = (f - FEX) / max(1.0D0,abs(f),abs(FEX))

  optm = ' '
  if ( reldiff .le. 1.0D-01 .and. maxinfeas .le. 1.0D-08 ) then
     optm = '*'
  end if

  if ( NEX .eq. -1 .or. .not. LEX ) then
     write(*,FMT=1000,ADVANCE='no') NTP,N,NILI + NINL,NELI + NENL,&
          reldiff,'INF',optm,f,FEX
  else
     write(*,FMT=1000,ADVANCE='no') NTP,N,NILI + NINL,NELI + NENL,&
          reldiff,'   ',optm,f,FEX
  end if


  ! Writes output information in file (previously opened)
  write(75,FMT=1010,ADVANCE='no') NTP,N,NILI + NINL,NELI + NENL,reldiff,&
       f,maxinfeas

  ! NON-EXECUTABLE STATEMENTS

1000 FORMAT(I4,1X,I4,1X,I4,1X,I4,5X,E12.5,1X,A3,1X,A1,1X,E12.5, &
            1X,E12.5)
1010 FORMAT(I4,1X,I4,1X,I4,1X,I4,5X,E12.5,1X,E12.5,1X,E12.5)

end subroutine finalize

!------------------------------------------------------------!
! SUBROUTINE EVAOBJFUN                                       !
!                                                            !
! Defines the objective function.                            !
!                                                            !
!------------------------------------------------------------!

subroutine evalobjfun(n,xe,f,flag)

  use hsdata

  implicit none

  ! INTERFACES
  interface
     subroutine CONV(MODE)
       integer,intent(in) :: MODE
     end subroutine CONV
  end interface       

  ! SCALAR ARGUMENTS
  integer :: flag,n
  real(8) :: f

  ! ARRAY ARGUMENTS
  real(8) :: xe(n)

  intent(in ) :: n,xe
  intent(out) :: f,flag

  ! COMMON SCALARS
  real(8) :: FX

  ! COMMON ARRAYS
  real(8) :: X(HS_NMAX)

  ! COMMON BLOCKS
  COMMON/L2/X
  COMMON/L6/FX

  flag = 0

  X(1:n) = xe(1:n)

  call CONV(2)

  f = FX

end subroutine evalobjfun

!------------------------------------------------------------!
! SUBROUTINE EVALCONSTR                                      !
!                                                            !
! Defines the constraints.                                   !
!                                                            !
!------------------------------------------------------------!

subroutine evalconstr(ne,xe,ind,c,flag)

  use hsdata

  implicit none

  ! INTERFACES
  interface
     subroutine CONV(MODE)
       integer,intent(in) :: MODE
     end subroutine CONV
  end interface       

  ! SCALAR ARGUMENTS
  integer :: flag,ind,ne
  real(8) :: c

  ! ARRAY ARGUMENTS
  real(8) :: xe(ne)

  intent(in ) :: ind,ne,xe
  intent(out) :: c,flag

  ! COMMON SCALARS
  integer :: N,NILI,NINL,NELI,NENL

  ! COMMON ARRAYS
  real(8) :: G(HS_MMAX),X(HS_NMAX)
  logical :: INDEX1(HS_MMAX)

  ! COMMON BLOCKS
  COMMON/L1/N,NILI,NINL,NELI,NENL
  COMMON/L2/X
  COMMON/L3/G
  COMMON/L9/INDEX1

  ! LOCAL SCALARS
  integer :: i,mor,meq

  flag = 0

  X(1:ne) = xe(1:ne)

  mor = NILI + NINL + NELI + NENL
  meq = NELI + NENL

  if ( ind .ge. 1 .and. ind .le. mor + meq ) then
     if ( ind .ge. mor + 1) then
        i = ind - meq
     else
        i = ind
     end if
  else
     flag = - 1
     return
  end if

  INDEX1(i) = .true.

  call CONV(4)

  c = - G(i)
  if ( ind .ge. mor + 1 ) then
     c = - c
  end if

  INDEX1(i) = .false.

end subroutine evalconstr

!------------------------------------------------------------!
! SUBROUTINE EVALJACOB                                       !
!                                                            !
! Defines the Jacobian of the constraints.                   !
!                                                            !
!------------------------------------------------------------!

subroutine evaljacob(ne,xe,ind,jcvar,jcval,jcnnz,flag)

  use hsdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,jcnnz,ne
  
  ! ARRAY ARGUMENTS
  integer :: jcvar(ne)
  real(8) :: jcval(ne),xe(ne)

  intent(in ) :: ind,ne,xe
  intent(out) :: flag,jcnnz,jcval,jcvar

  ! INTERFACES
  interface
     subroutine CONV(MODE)
       integer,intent(in) :: MODE
     end subroutine CONV
  end interface    

  ! COMMON SCALARS
  integer :: N,NILI,NINL,NELI,NENL

  ! COMMON ARRAYS
  real(8) :: GG(HS_MMAX * HS_NMAX),X(HS_NMAX)
  logical :: INDEX2(HS_MMAX)

  ! COMMON BLOCKS
  COMMON/L1/N,NILI,NINL,NELI,NENL
  COMMON/L2/X
  COMMON/L5/GG
  COMMON/L10/INDEX2

  ! LOCAL SCALARS
  integer :: i,j,mor,meq

  flag = 0

  X(1:ne) = xe

  mor = NILI + NINL + NELI + NENL
  meq = NELI + NENL

  if ( ind .ge. 1 .and. ind .le. mor + meq ) then
     if ( ind .ge. mor + 1) then
        i = ind - meq
     else
        i = ind
     end if
  else
     flag = - 1
     return
  end if

  INDEX2(i) = .true.

  call CONV(5)

  jcnnz =                            ne
  jcvar =           (/ (j, j = 1,ne) /)
  do j = 1,ne
     jcval(j) = - GG((j - 1) * mor + i)
  end do

  if ( ind .ge. mor + 1 ) then
     jcval(1:jcnnz) = - jcval(1:jcnnz)
  endif

  INDEX2(i) = .false.

end subroutine evaljacob

!------------------------------------------------------------!
! SUBROUTINE DRAWSBI                                         !
!                                                            !
! This subroutine draws the state of the problem (if/when    !
! possible) at the beginning of iteration 'iter'.            !
!                                                            !
!------------------------------------------------------------!

subroutine drawsBI(n,x,m,gamma,infeas,iter)

  implicit none

  ! SCALAR ARGUMENTS
  integer :: iter,m,n
  real(8) :: infeas

  ! ARRAY ARGUMENTS
  real(8) :: gamma(m),x(n)

  intent(in) :: gamma,infeas,iter,m,n,x

end subroutine drawsBI

!------------------------------------------------------------!
! SUBROUTINE DRAWSAR                                         !
!                                                            !
! This subroutine draws the state of the problem (if/when    !
! possible) immediately after the restoration phase at       !
! iteration 'iter'.                                          !
!                                                            !
!------------------------------------------------------------!

subroutine drawsAR(n,x,m,gamma,infeas,iter)

  implicit none

  ! SCALAR ARGUMENTS
  integer :: iter,m,n
  real(8) :: infeas

  ! ARRAY ARGUMENTS
  real(8) :: gamma(m),x(n)

  intent(in) :: gamma,infeas,iter,m,n,x

end subroutine drawsAR


!------------------------------------------------------------!
! SUBROUTINE DRAWSAO                                         !
!                                                            !
! This subroutine draws the state of the problem (if/when    !
! possible) immediately after the optimization phase at      !
! iteration 'iter'.                                          !
!                                                            !
!------------------------------------------------------------!

subroutine drawsAO(ne,x,m,gamma,infeas,iter)

  use hsdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: iter,m,ne
  real(8) :: infeas

  ! ARRAY ARGUMENTS
  real(8) :: gamma(m),x(ne)

  intent(in) :: gamma,infeas,iter,m,ne,x

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
  real(8) :: f,reldiff
  integer :: flag

  ! COMMON SCALARS
  integer :: N,NEX,NILI,NINL,NELI,NENL,NTP
  logical :: LEX
  real(8) :: FEX

  ! COMMON ARRAYS
  real(8) :: XEX(HS_MMAX * HS_NMAX)

  ! COMMON BLOCKS
  COMMON/L1/N,NILI,NINL,NELI,NENL
  COMMON/L8/NTP
  COMMON/L20/LEX,NEX,FEX,XEX

  call evalobjfun(ne,x,f,flag)

  reldiff = (f - FEX) / max(1.0D0,abs(f),abs(FEX))

  open(75,FILE='runhs.out')
  write(75,FMT=0020) NTP,N,NILI + NINL,NELI + NENL,reldiff,&
       f,infeas,-1
  close(75)

  ! NON-EXECUTABLE STATEMENTS

0020 FORMAT(I4,1X,I4,1X,I4,1X,I4,5X,E12.5,1X,E12.5,1X,E12.5,1X,I15)

end subroutine drawsAO

