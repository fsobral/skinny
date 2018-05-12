!------------------------------------------------------------!
! TOYPROB.F90                                                !
!                                                            !
! This file is a simple problem for testing, taken from the  !
! Hock and Schittkowski set:                                 !
!                                                            !
! (HS6) min (1 - x_1) ** 2 s.t. 10 * (x_2 - x_1 ** 2) = 0    !
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

  implicit none

  ! SCALAR ARGUMENTS
  integer,intent(out) :: flag,m,n

  flag = 0

  n = 2
  m = 1

end subroutine initializes

!------------------------------------------------------------!
! SUBROUTINE INITIALIZEA                                     !
!                                                            !
! This subroutine initializes the ARRAY information of the   !
! problem.                                                   !
!                                                            !
!------------------------------------------------------------!

subroutine initializea(n,x,l,u,m,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,m,n

  ! ARRAY ARGUMENTS
  real(8) :: l(n),u(n),x(n)

  intent(in ) :: m,n
  intent(out) :: flag,l,u,x

  ! LOCAL SCALARS
  integer :: i

  flag = 0

  x(1) = - 1.2D0
  x(2) =   1.0D0

  do i = 1,n
     l(i) = - 1.0D+20
     u(i) =   1.0D+20
  end do

end subroutine initializea

!------------------------------------------------------------!
! SUBROUTINE FINALIZE                                        !
!                                                            !
! Finalizes the problem, sometimes printing additional       !
! information related to the specific problem.               !
!                                                            !
!------------------------------------------------------------!

subroutine finalize(ne,x,m,maxinfeas)

  implicit none

  ! SCALAR ARGUMENTS
  integer,intent(in) :: m,ne
  real(8),intent(in) :: maxinfeas

  ! ARRAY ARGUMENTS
  real(8),intent(in) :: x(ne)

end subroutine finalize

!------------------------------------------------------------!
! SUBROUTINE EVAOBJFUN                                       !
!                                                            !
! Defines the objective function.                            !
!                                                            !
!------------------------------------------------------------!

subroutine evalobjfun(n,x,f,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,n
  real(8) :: f

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: n,x
  intent(out) :: f,flag

  flag = 0

  f = (1.0D0 - x(1)) ** 2

end subroutine evalobjfun

!------------------------------------------------------------!
! SUBROUTINE EVALCONSTR                                      !
!                                                            !
! Defines the constraints.                                   !
!                                                            !
!------------------------------------------------------------!

subroutine evalconstr(n,x,ind,c,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,n
  real(8) :: c

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: ind,n,x
  intent(out) :: c,flag

  flag = 0

  if ( ind .eq. 1 ) then
     c = 10.0D0 * (x(2) - x(1) ** 2.0D0)
  else
     flag = - 1
  end if

end subroutine evalconstr

!------------------------------------------------------------!
! SUBROUTINE EVALJACOB                                       !
!                                                            !
! Defines the Jacobian of the constraints.                   !
!                                                            !
!------------------------------------------------------------!

subroutine evaljacob(n,x,ind,jcvar,jcval,jcnnz,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,jcnnz,n
  
  ! ARRAY ARGUMENTS
  integer :: jcvar(n)
  real(8) :: jcval(n),x(n)

  intent(in ) :: ind,n,x
  intent(out) :: flag,jcnnz,jcval,jcvar

  flag = 0

  if ( ind .eq. 1 ) then
     jcnnz = 2
     
     jcvar(1) = 1
     jcval(1) = - 20.0D0 * x(1)

     jcvar(2) = 2
     jcval(2) = 10.0D0
  else
     flag = - 1
  end if

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

  implicit none

  ! SCALAR ARGUMENTS
  integer :: iter,m,n
  real(8) :: gamma,infeas

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in) :: gamma,infeas,iter,m,n,x

end subroutine drawsAO
