!------------------------------------------------------------!
! GREEDY.F90                                                 !
!                                                            !
! This file is an interface between the source of greedy     !
! problems defined in greediness.f90 and Skinny.             !
!                                                            !
!------------------------------------------------------------!

!------------------------------------------------------------!
! SUBROUTINE INITIALIZES                                     !
!                                                            !
! This subroutine initializes the SCALAR information of the  !
! problem.                                                   !
!                                                            !
!------------------------------------------------------------!

subroutine initializes(n,m,flag)

  use dftests

  implicit none

  ! SCALAR ARGUMENTS
  integer,intent(out) :: flag,m,n


!  write(*,FMT=0010)
  read (*,*) df_NP

!  if ( df_NP .lt. 1 .or. df_NP .gt. 6 ) then
!     df_NP = 1
!  end if

  flag = callp(1)

  n = df_n
  m = 2 * (df_ml + df_mn) + df_pl + df_pn

  open(75,FILE='runeng.out')
  write(75,FMT=0020) df_NP,df_n,df_pl + df_pn,df_ml + df_mn,1.0D+20,&
       1.0D+20,-1
  close(75)

! NON-EXECUTABLE STATEMENTS

0010 FORMAT('Please type the number of the problem to be solved ', &
            '(1 - 6):')
0020 FORMAT(I4,1X,I4,1X,I4,1X,I4,5X,E12.5,1X,E12.5,1X,I15,1X)


end subroutine initializes


!------------------------------------------------------------!
! SUBROUTINE INITIALIZEA                                     !
!                                                            !
! This subroutine initializes the ARRAY information of the   !
! problem.                                                   !
!                                                            !
!------------------------------------------------------------!

subroutine initializea(n,x,l,u,m,flag)

  use dftests

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,m,n

  ! ARRAY ARGUMENTS
  real(8) :: l(n),u(n),x(n)

  intent(in ) :: m,n
  intent(out) :: flag,l,u,x

  ! LOCAL SCALARS
  integer :: error,i

  x(1:n) = df_x(1:n)

  l(1:n) = df_l(1:n)
  u(1:n) = df_u(1:n)

  flag   = 0

end subroutine initializea

!------------------------------------------------------------!
! SUBROUTINE FINALIZE                                        !
!                                                            !
! Finalizes the problem, sometimes printing additional       !
! information related to the specific problem.               !
! This subroutine uses the following modules:                !
!                                                            !
!------------------------------------------------------------!

subroutine finalize(n,x,m,maxinfeas)

  use dftests

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
  integer,intent(in) :: m,n
  real(8),intent(in) :: maxinfeas

  ! ARRAY ARGUMENTS
  real(8),intent(in) :: x(n)

  ! LOCAL SCALARS
  character :: optm
  integer :: i,flag
  real(8) :: diff,f,reldiff

  call evalobjfun(n,x,f,flag)

  write(*,FMT=1000,ADVANCE='no') df_NP,df_n,df_pl + df_pn,df_ml + df_mn,&
       f,maxinfeas

  OPEN(75,FILE='runeng.out')
  write(75,FMT=1000,ADVANCE='no') df_NP,df_n,df_pl + df_pn,df_ml + df_mn,&
       f,maxinfeas
  CLOSE(75)

  flag = callp(5)

  ! NON-EXECUTABLE STATEMENTS

1000 FORMAT(I4,1X,I4,1X,I4,1X,I4,5X,E12.5,1X,E12.5,1X)

end subroutine finalize

!------------------------------------------------------------!
! SUBROUTINE EVALOBJFUN                                      !
!                                                            !
! Defines the objective function.                            !
!                                                            !
!------------------------------------------------------------!

subroutine evalobjfun(n,x,f,flag)

  use dftests

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,n
  real(8) :: f

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: n,x
  intent(out) :: f,flag

  df_x(1:n) = x(1:n)

  flag = callp(2)

  f = df_f

end subroutine evalobjfun

!------------------------------------------------------------!
! SUBROUTINE EVALCONSTR                                      !
!                                                            !
! Defines the constraints.                                   !
!                                                            !
!------------------------------------------------------------!

subroutine evalconstr(n,x,ind,c,flag)

  use dftests

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,n
  real(8) :: c

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: ind,n,x
  intent(out) :: c,flag


  ! LOCAL SCALARS
  integer :: i,mor,meq

  flag = 0

  df_x(1:n) = x(1:n)

  mor = df_ml + df_mn + df_pl + df_pn
  meq = df_ml + df_mn

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

  df_ind = i

  flag = callp(3)

  c = df_c

  if ( ind .ge. mor + 1 ) then
     c = - c
  end if

end subroutine evalconstr

!------------------------------------------------------------!
! SUBROUTINE EVALJACOB                                       !
!                                                            !
! Defines the Jacobian of the constraints.                   !
!                                                            !
!------------------------------------------------------------!

subroutine evaljacob(n,x,ind,jcvar,jcval,jcnnz,flag)

  use dftests

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,jcnnz,n
  
  ! ARRAY ARGUMENTS
  integer :: jcvar(n)
  real(8) :: jcval(n),x(n)

  intent(in ) :: ind,n,x
  intent(out) :: flag,jcnnz,jcval,jcvar


  ! LOCAL SCALARS
  integer :: i,j,mor,meq

  df_x(1:n) = x(1:n)

  mor = df_ml + df_mn + df_pl + df_pn
  meq = df_ml + df_mn

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

  df_ind = i

  flag = callp(4)

  jcnnz = df_jcnnz
  jcvar(1:jcnnz) = df_jcvar(1:jcnnz)
  jcval(1:jcnnz) = df_jcval(1:jcnnz)

  if ( ind .ge. mor + 1 ) then
     jcval(1:jcnnz) = - jcval(1:jcnnz)
  endif

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
  real(8) :: gamma,infeas

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

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
  real(8) :: gamma,infeas

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

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

subroutine drawsAO(n,x,m,gamma,infeas,iter)

  use dftests

  implicit none

  ! SCALAR ARGUMENTS
  integer :: iter,m,n
  real(8) :: gamma,infeas

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in) :: gamma,infeas,iter,m,n,x

  ! LOCAL SCALARS
  integer :: i,flag
  real(8) :: diff,f,reldiff

  call evalobjfun(n,x,f,flag)

  open(75,FILE='runeng.out')
  write(75,FMT=0020) df_NP,df_n,df_pl + df_pn,df_ml + df_mn,&
       f,infeas,-1
  close(75)

  ! NON-EXECUTABLE STATEMENTS

0020 FORMAT(I4,1X,I4,1X,I4,1X,I4,5X,E12.5,1X,E12.5,&
          1X,I15)

end subroutine drawsAO

