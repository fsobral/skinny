!------------------------------------------------------------!
! DISTLINE.F90                                               !
!                                                            !
! Given a point in R^n, this problem tries to find the point !
! in a given line which has the minimum norm-2 distance to   !
! that point.                                                !
!                                                            !
!------------------------------------------------------------!


!------------------------------------------------------------!
! SUBROUTINE INIT                                            !
!                                                            !
! Initializes the problem.                                   !
!                                                            !
!------------------------------------------------------------!

subroutine init(n,x,l,u,m,flag)

  use distlinedata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,m,n

  ! ARRAY ARGUMENTS
  real(8),pointer :: l(:),u(:),x(:)

  intent(out) :: flag,m,n

  ! LOCAL SCALARS
  integer error,i

  write(*,*) 'Please enter the dimension of the problem'
  read (*,*) n

  allocate(x(n),l(n),u(n),upoint(n), STAT=error)
  if ( error .ne. 0 ) then
     flag = - 1
  end if

  write(*,*) 'Enter the coordinates of the point to be projected.'
  read(*,*) (upoint(i), i = 1,n)

  x =    upoint
  l = - 1.0D+20
  u =   1.0D+20

  m = 2

  OPEN(40,FILE='running.mp')
  WRITE(40,FMT=6000)

  ! NON-EXECUTABLE STATEMENTS

6000 FORMAT('%% This file was generated automatically by '         ,&
            'Engordamento',/,'%% ',53('-'),/,/,'verbatimtex',/     ,&
            '%&latex',/,'\documentclass{article}',/,'\usepackage[l',&
            'atin1]{inputenc}',/'\\begin{document}',/,'etex',/,/   ,&
            'vardef constr(expr x,gamma) = 2.0 - x + gamma enddef;',&
            /,/,'numeric xmin,xmax;',/,'xmin = -10;',/             ,&
            'xmax =  10;',/)

end subroutine init

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

  WRITE(40,FMT=6010)iter,gamma,gamma,-gamma,-gamma,x(1),x(2)

  ! NON-EXECUTABLE STATEMENTS

6010 FORMAT('beginfig(',I4,');',/,/,'u := 2cm;',/,'drawarrow '     ,&
            '(-2,0) * u -- (2.5,0) * u;',/,'drawarrow (0,-1) * u ' ,&
            '-- (0,2.5) * u;',/,/,'draw (xmin,constr(xmin,0)) * u ',&
            '-- (xmax,constr(xmax,0)) * u;',/,'draw (xmin,constr(x',&
            'min,',F20.15,')) * u -- (xmax,constr(xmax,',F20.15    ,&
            ')) * u dashed evenly;',/,'draw (xmin,constr(xmin,'    ,&
            F20.15,')) * u -- (xmax,constr(xmax,',F20.15,')) * u ' ,&
            'dashed evenly;',/,/,'% Last iteration',/,'draw ('     ,&
            F20.15,',',F20.15,') * u withpen pencircle scaled 4pt ',&
            'withcolor 0.8 red;',/,/,'% Current iteration')

end subroutine drawsBI

!------------------------------------------------------------!
! SUBROUTINE DRAWSAR                                         !
!                                                            !
! This subroutine draws the state of the problem (if/when    !
! possible) imediately after the restoration phase at        !
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

  WRITE(40,FMT=6020)x(1),x(2),'pensquare'

  ! NON-EXECUTABLE STATEMENTS

6020 FORMAT('draw (',F20.15,',',F20.15,') * u withpen ',A9,1X      ,&
            'scaled 4pt;')

end subroutine drawsAR


!------------------------------------------------------------!
! SUBROUTINE DRAWSAO                                         !
!                                                            !
! This subroutine draws the state of the problem (if/when    !
! possible) imediately after the optimization phase at       !
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

  WRITE(40,FMT=6020)x(1),x(2),'pensquare'
  WRITE(40,FMT=6030)

  ! NON-EXECUTABLE STATEMENTS

6020 FORMAT('draw (',F20.15,',',F20.15,') * u withpen ',A9,1X      ,&
            'scaled 4pt;')
6030 FORMAT(/,'picture all;',/,'all := currentpicture;',/          ,&
            'clearit;',/,'path box;',/,'box := (-2,-1) * u -- (-2,',&
            '2.5) * u -- (2.5,2.5) * u -- (2.5,-1) * u -- cycle;',/,&
            /,'clip all to box;',/,'draw all;',/,'endfig;',/,/)

end subroutine drawsAO

!------------------------------------------------------------!
! SUBROUTINE FINA                                            !
!                                                            !
! Finalizes the problem, sometimes printing aditional        !
! information related to the specific problem.               !
!                                                            !
!------------------------------------------------------------!

subroutine fina(n,x,m)

  implicit none

  ! SCALAR ARGUMENTS
  integer,intent(in) :: m,n

  ! ARRAY ARGUMENTS
  real(8),intent(in) :: x(n)

  WRITE(40,FMT=6040)
  CLOSE(40)

  ! NON-EXECUTABLE STATEMENTS

6040 FORMAT(/,/,'end;')

end subroutine fina

!------------------------------------------------------------!
! SUBROUTINE EVALF                                           !
!                                                            !
! Defines the objetive function.                             !
!                                                            !
!------------------------------------------------------------!

subroutine evalobjf(n,x,f,flag)

  use distlinedata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,n
  real(8) :: f

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: n,x
  intent(out) :: f,flag

  flag = 0

  f = DOT_PRODUCT(x - upoint,x - upoint)

end subroutine evalobjf

!------------------------------------------------------------!
! SUBROUTINE EVALC                                           !
!                                                            !
! Defines the contraints.                                    !
!                                                            !
!------------------------------------------------------------!

subroutine evalc(n,x,ind,c,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,n
  real(8) :: c

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: ind,n,x
  intent(out) :: c,flag

!  write(*,*) 'Entrou evalc'

  flag = 0

  if ( ind .eq. 1 ) then
     c =   SUM(x) - 2.0D0
  else if ( ind .eq. 2 ) then
     c = - SUM(x) + 2.0D0
  else
     flag = - 1
  end if

!  write(*,*) 'Saiu evalc'

end subroutine evalc

!------------------------------------------------------------!
! SUBROUTINE EVALJAC                                         !
!                                                            !
! Defines the Jacobian of the constraints.                   !
!                                                            !
!------------------------------------------------------------!

subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,jcnnz,n
  
  ! ARRAY ARGUMENTS
  integer :: jcvar(n)
  real(8) :: jcval(n),x(n)

  intent(in ) :: ind,n,x
  intent(out) :: flag,jcnnz,jcval,jcvar

  ! LOCAL SCALARS
  integer :: i

  flag = 0

!  write(*,*) 'Entrou evaljac'

  if ( ind .eq. 1 ) then
     jcnnz = n

     do i = 1,n
        jcvar(i) = i
     end do
     jcval = 1.0D0
  else if ( ind .eq. 2 ) then
     jcnnz = n

     do i = 1,n
        jcvar(i) = i
     end do
     jcval = - 1.0D0
  else
     flag = - 1
  end if

!  write(*,*) 'Saiu evaljac'

end subroutine evaljac
