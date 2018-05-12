!------------------------------------------------------------!
! SUBROUTINE INIT                                            !
!                                                            !
! Initializes the problem.                                   !
!                                                            !
!------------------------------------------------------------!

subroutine init(n,x,m,flag)

  use ellipsedata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,m,n

  ! ARRAY ARGUMENTS
  real(8),pointer :: x(:)

  intent(out) :: flag,m,n

  ! LOCAL SCALARS
  integer error

  n = 5

!!$  n = 6
!!$
  allocate(x(n), STAT=error)
  if ( error .ne. 0 ) then
     flag = - 1
  end if

  x(1:n-1) =  0.5D0
  x(n)     =  10.5D0
!!$
!!$  m = 13

  x(1) = 0.0D0
  x(2) = 0.0D0
  x(3) = 4.0D0
  x(4) = 3.0D0
  x(5) = 3.5D0

  m = 10

  OPEN(40,FILE='running.mp')
  WRITE(40,FMT=6000)

! NON-EXECUTABLE STATEMENTS

6000 FORMAT('%% This file was generated automatically by '         ,&
            'Engordamento',/,'%% ',53('-'),/,&
            /,'verbatimtex',/,'%&latex',/,'\documentclass{article}',&
            /,'\usepackage[latin1]{inputenc}',/'\\begin{document}' ,&
            /,'etex',/,/,'beginfig(1);',/,/,'u := 3cm;',/,/        ,&
            'pair udots[];')


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

  use ellipsedata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: iter,m,n
  real(8) :: gamma,infeas

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in) :: gamma,infeas,iter,m,n,x

  ! LOCAL SCALARS
  real(8) :: a,angle,b,ox,oy

  if ( MOD(iter,5) .eq. 0 ) then

     call cellipse(x,ox,oy,a,b,angle)

     WRITE(40,FMT=6010)iter,ox,oy
     WRITE(40,FMT=6020)2.0D0 * a,2.0D0 * b,angle,ox,oy,INT(iter)/100.0D0

  end if

  ! NON-EXECUTABLE STATEMENTS

6010 FORMAT('label(btex $',I4,'$ etex,(',F20.15,',',F20.15,') * u);')

6020 FORMAT('draw fullcircle xscaled ',F20.15,'u yscaled ',F20.15  ,&
            'u ','rotated ',F20.15,/,' shifted ((',F20.15,','      ,&
            F20.15,') * u) withpen pencircle scaled 0.5pt '        ,&
            'withcolor ',F20.15,'blue;',/)

end subroutine drawsAO

!------------------------------------------------------------!
! SUBROUTINE FINA                                            !
!                                                            !
! Finalizes the problem, sometimes printing aditional        !
! information related to the specific problem.               !
!                                                            !
!------------------------------------------------------------!

subroutine fina(n,x,m)

  use ellipsedata

  implicit none

  ! SCALAR ARGUMENTS
  integer,intent(in) :: m,n

  ! ARRAY ARGUMENTS
  real(8),intent(in) :: x(n)

  ! LOCAL SCALARS
  integer :: i
  real(8) :: ox,oy,a,b,angle

  call cellipse(x,ox,oy,a,b,angle)

  write(*,FMT=1010) ox,oy,a,b,angle

  WRITE(40,FMT=6030)
  WRITE(40,FMT=6050)
  CLOSE(40)

  OPEN(10,FILE='solution.mp')
  write(10,FMT=6010)
  do i = 1,4
     write(10,FMT=6020) (i - 1),p(2 * i - 1),p(2 * i)
  end do
  write(10,FMT=6030)
  write(10,FMT=6040)2.0D0 * a,2.0D0 * b,angle,ox,oy
  write(10,FMT=6050)
  CLOSE(10)

  ! NON-EXECUTABLE STATEMENTS

1010 FORMAT(/,'Ellipse data:',/,'Center: (',F20.15,',',F20.15,')',/,&
            'Major axis: ',F20.15,/,'Minor axis: ',F20.15,/,'Angle',&
            ' counterclockwise from the x-axis to the major '      ,&
            'axis: ',F20.15,/)

6010 FORMAT('%% This file was generated automatically by '         ,&
            'Engordamento',/,'%% ',53('-'),/,&
            /,'verbatimtex',/,'%&latex',/,'\documentclass{article}',&
            /,'\usepackage[latin1]{inputenc}',/'\\begin{document}' ,&
            /,'etex',/,/,'beginfig(1);',/,/,'u := 3cm;',/,/        ,&
            'pair udots[];')
6020 FORMAT('udots[',I2,'] := (',F20.15,',',F20.15,') * u;')
6030 FORMAT(/,'for i = 0 step 1 until 3:',/,'draw udots[i] withpen',&
            ' pencircle scaled 7pt withcolor red;',/,'endfor',/,/)
6040 FORMAT('draw fullcircle xscaled ',F20.15,'u yscaled ',F20.15  ,&
            'u ','rotated ',F20.15,/,' shifted ((',F20.15,','      ,&
            F20.15,') * u) withpen pencircle scaled 2pt;',/)
6050 FORMAT('endfig;',/,'end;')

end subroutine fina

!------------------------------------------------------------!
! SUBROUTINE EVALF                                           !
!                                                            !
! Defines the objetive function.                             !
!                                                            !
!------------------------------------------------------------!

subroutine evalobjf(n,x,f,flag)

  use ellipsedata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,n
  real(8) :: f

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: n,x
  intent(out) :: f,flag

  ! LOCAL SCALARS
  integer :: i,nin,ntot
  real(8) :: area

  ! LOCAL ARRAYS
  real(8) :: c(2),point(2)

  flag = 0

  c(1) = (x(1) + x(3)) / 2.0D0
  c(2) = (x(2) + x(4)) / 2.0D0

  area = 4.0D0 * x(5) ** 2

  ntot = 10000
  nin  =     0

 ! write(*,*) '->',c,x(5),evalellipse(x,c),SQRT((x(1) - x(3)) ** 2 + (x(2) - x(4)) ** 2)

  call RANDOM_SEED
  do i = 1,ntot
     call RANDOM_NUMBER(point)
     point = (point - 0.5D0) * 2.0D0
     point = c + point * x(5)
!     write(*,*) point,evalellipse(x,point)
     if ( evalellipse(x,point) .le. 0.0D0 ) then
        nin = nin + 1
     end if
  end do

  f = area * (DBLE(nin) / DBLE(ntot))

end subroutine evalobjf

!------------------------------------------------------------!
! SUBROUTINE EVALC                                           !
!                                                            !
! Defines the contraints.                                    !
!                                                            !
!------------------------------------------------------------!

subroutine evalc(n,x,ind,c,flag)

  use ellipsedata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,n
  real(8) :: c

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: ind,n,x
  intent(out) :: c,flag

  ! LOCAL SCALARS
  integer :: i

  flag = 0

  if ( ind .ge. 1 .and. ind .le. 4 ) then
     i = 2 * ind - 1
     c = evalellipse(x,p(i:i+1))
  else if ( ind .ge. 5 .and. ind .le. 8 ) then
     i = 2 * (ind - 4) - 1
     c = - evalellipse(x,p(i:i+1))
  else if ( ind .eq. 9 ) then
     c = (2.0D0 * x(5)) ** 2 - (x(1) - x(3)) ** 2 - (x(2) - x(4)) ** 2
  else if ( ind .eq. 10 ) then
     c = - x(5)
!!$  if ( ind .ge. 1 .and. ind .le. 5 ) then
!!$     i = 2 * ind - 1
!!$     c = evalellipse(x,p(i:i + 1))
!!$  else if ( ind .ge. 6 .and. ind .le. 10 ) then
!!$     i = 2 * (ind - 5) - 1
!!$     c = - evalellipse(x,p(i:i + 1))
!!$  else if ( ind .eq. 11 ) then
!!$     c =   x(1) * x(3) * x(6) + 2.0D0 * x(2) * x(4) * x(5) - &
!!$           x(3) * x(4) ** 2 - x(2) ** 2 * x(6) - x(1) * x(5) ** 2 - ZERO
!!$  else if ( ind .eq. 12 ) then
!!$     c = - x(1) * x(3) * x(6) - 2.0D0 * x(2) * x(4) * x(5) + &
!!$           x(3) * x(4) ** 2 + x(2) ** 2 * x(6) + x(1) * x(5) ** 2 + ZERO
!!$  else if ( ind .eq. 13 ) then
!!$     c = x(2) ** 2 - x(1) * x(3) - ZERO
  else
     flag = - 1
  end if

end subroutine evalc

!------------------------------------------------------------!
! SUBROUTINE EVALJAC                                         !
!                                                            !
! Defines the Jacobian of the constraints.                   !
!                                                            !
!------------------------------------------------------------!

subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

  use ellipsedata

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

!  write(*,*) 'Chamou evalc',ind

  if ( ind .ge. 1 .and. ind .le. 4 ) then
     i = 2 * ind - 1
     call evalgellipse(x,p(i:i + 1),jcvar,jcval,jcnnz)
  else if ( ind .ge. 5 .and. ind .le. 8 ) then
     i = 2 * (ind - 4) - 1
     call evalgellipse(x,p(i:i + 1),jcvar,jcval,jcnnz)
     jcval = - jcval
  else if ( ind .eq. 9 ) then
     jcnnz = 5
     
     jcvar(1:jcnnz) = (/ 1,2,3,4,5 /)
     jcval(1) = - 2.0D0 * (x(1) - x(3))
     jcval(2) = - 2.0D0 * (x(2) - x(4))
     jcval(3) =   2.0D0 * (x(1) - x(3))
     jcval(4) =   2.0D0 * (x(2) - x(4))
     jcval(5) =   8.0D0 * x(5)
  else if ( ind .eq. 10 ) then
     jcnnz = 1
     jcvar(1) = 5
     jcval(1) = - 1.0D0

!!$  if ( ind .ge. 1 .and. ind .le. 5 ) then
!!$     i = 2 * ind - 1
!!$     call evalgellipse(x,p(i:i + 1),jcvar,jcval,jcnnz)
!!$  else if ( ind .ge. 6 .and. ind .le. 10 ) then
!!$     i = 2 * (ind - 5) - 1
!!$     call evalgellipse(x,p(i:i + 1),jcvar,jcval,jcnnz)
!!$     jcval = - jcval
!!$  else if ( ind .eq. 11 ) then
!!$     jcnnz = 6
!!$
!!$     do i = 1,jcnnz
!!$        jcvar(i) = i
!!$     end do
!!$
!!$     jcval(1) = x(3) * x(6) - x(5) ** 2
!!$     jcval(2) = 2.0D0 * (x(4) * x(5) - x(2) * x(6))
!!$     jcval(3) = x(1) * x(6) - x(4) ** 2
!!$     jcval(4) = 2.0D0 * (x(2) * x(5) - x(3) * x(4))
!!$     jcval(5) = 2.0D0 * (x(2) * x(4) - x(1) * x(5))
!!$     jcval(6) = x(1) * x(3) - x(2) ** 2
!!$  else if ( ind .eq. 12 ) then
!!$     jcnnz = 6
!!$
!!$     do i = 1,jcnnz
!!$        jcvar(i) = i
!!$     end do
!!$
!!$     jcval(1) = - x(3) * x(6) + x(5) ** 2
!!$     jcval(2) = - 2.0D0 * (x(4) * x(5) - x(2) * x(6))
!!$     jcval(3) = - x(1) * x(6) + x(4) ** 2
!!$     jcval(4) = - 2.0D0 * (x(2) * x(5) - x(3) * x(4))
!!$     jcval(5) = - 2.0D0 * (x(2) * x(4) - x(1) * x(5))
!!$     jcval(6) = - x(1) * x(3) + x(2) ** 2
!!$  else if ( ind .eq. 13 ) then
!!$     jcnnz = 3
!!$
!!$     jcvar(1) = 1
!!$     jcval(1) = - x(3)
!!$
!!$     jcvar(2) = 2
!!$     jcval(2) = 2.0D0 * x(2)
!!$
!!$     jcvar(3) = 3
!!$     jcval(3) = - x(1)
  else
     flag = - 1
  end if

!  write(*,*) 'Saiu evalc'

end subroutine evaljac
