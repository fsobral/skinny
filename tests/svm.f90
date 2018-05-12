!------------------------------------------------------------!
! SVM.F90                                                    !
!                                                            !
! Given three types of points: IN, BORDER and OTHER in a     !
! quantity of 'npin', 'npbo' and 'npot',                     !
! respectively, and 'nhyper' hyperplanes, the objective is   !
! to find the position of these hyperplanes such that        !
!                                                            !
! * The BORDER points are in the border of any hyperplane    !
!   and inside the generated polytope;                       !
! * The IN points are inside the polytope;                   !
!                                                            !
! * The number of OTHER points which are inside the polytope !
!   must be minimized.                                       !
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

  use svmdata

  implicit none

  ! SCALAR ARGUMENTS
  integer,intent(out) :: flag,m,n

  ! LOCAL SCALARS
  integer error

  call load()

  ! Number of variables

  n = nhyper * (dim + 1)

  ! Number of constraints

!  m = 2 * npbo + (npbo + npin + 2) * nhyper
  m = 2 * npbo + (npbo + 2 * npin + 1) * nhyper

  flag   = 0

end subroutine initializes

!------------------------------------------------------------!
! SUBROUTINE INITIALIZEA                                     !
!                                                            !
! This subroutine initializes the ARRAY information of the   !
! problem.                                                   !
!                                                            !
!------------------------------------------------------------!

subroutine initializea(n,x,l,u,m,flag)
  
  use svmdata

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

  ! Lower and upper bounds

  l(1:n) = - 1.0D+20
  do i = 0,nhyper - 1
     u(i * (dim + 1) + 1:i * (dim + 1) + dim) = 0.0D0!1.0D+20
     u((i + 1) * (dim + 1))                   = 1.0D+20
  end do

  ! Starting point

!!$  call RANDOM_NUMBER(x(1:n))
!!$  x(1:n) = - 10.0D0 + x(1:n) * 20.0D0
!!$  do i = 0,nhyper - 1
!!$     x((i + 1) * (dim + 1)) = 200.0D0
!!$  end do
  

  do i = 1,n
     read(*,*) x(i)
  end do

  OPEN(98,file='SVMPOINTS.IN',position='APPEND')
  do i = 1,n
     write(98,FMT='(E23.16)') x(i)
  end do
  CLOSE(98)

  OPEN(40,FILE='running.mp')
  WRITE(40,FMT=6000)xmin,xmax,ymin,ymax,scale,xmin,ymin,xmin,ymax,&
       xmax,ymax,xmax,ymin
  
  if ( dim .eq. 2 ) then
     
     WRITE(40,FMT=6010)
     do i = 1,npbo
        WRITE(40,FMT=6011) i,pntbo(1,i),pntbo(2,i)
     end do
     
     WRITE(40,FMT=6020)
     do i = 1,npin
        WRITE(40,FMT=6021) i,pntin(1,i),pntin(2,i)
     end do
     
     WRITE(40,FMT=6030)
     do i = 1,min(SVMMAXDOTS,npot)
        WRITE(40,FMT=6031) i,pntot(1,i),pntot(2,i)
     end do

     WRITE(40,*)

  end if

  CLOSE(40)

  ! NON-EXECUTABLE STATEMENTS

6000 FORMAT('%% This file was generated automatically by '         ,&
            'Engordamento',/,'%% ',53('-'),/,/,'verbatimtex',/     ,&
            '%&latex',/,'\documentclass{article}',/,'\usepackage[l',&
            'atin1]{inputenc}',/'\\begin{document}',/,'etex',/,/   ,&
            'vardef hyper(expr x,b) = b - x enddef;'               ,&
            /,/,'numeric xmin,xmax,ymax,ymin;',/,'xmin = ',F20.13  ,&
            ';',/,'xmax = ',F20.13,';',/,'ymin = ',F20.13,';',/    ,&
            'ymax = ',F20.13,';',/,'u := ',F20.13,'cm;',/,/,'path ',&
            'box;',/,'box := ',4('(',F20.13,',',F20.13,') * u -- '),&
            ' cycle;' )
6010 FORMAT(/,'pair pntbo[];')
6011 FORMAT('pntbo[',I6,'] := (',F20.15,',',F20.15,');')
6020 FORMAT(/,'pair pntin[];')
6021 FORMAT('pntin[',I6,'] := (',F20.15,',',F20.15,');')
6030 FORMAT(/,'pair pntot[];')
6031 FORMAT('pntot[',I6,'] := (',F20.15,',',F20.15,');')

end subroutine initializea

!------------------------------------------------------------!
! SUBROUTINE FINALIZE                                        !
!                                                            !
! Finalizes the problem, sometimes printing additional       !
! information related to the specific problem.               !
! This subroutine uses the following modules:                !
!                                                            !
! SVMDATA                                                    !
!                                                            !
!------------------------------------------------------------!

subroutine finalize(n,x,m,maxinfeas)

  use svmdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: m,n
  real(8) :: maxinfeas

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in) :: m,maxinfeas,n,x

  ! LOCAL SCALARS
  integer :: j
  real(8) :: a,b,c,l,u

  OPEN(40,FILE='running.mp',POSITION='APPEND')

  if ( dim .eq. 2 ) then

     WRITE(40,FMT=6040)1000
     WRITE(40,FMT=6080)min(SVMMAXDOTS,npot),npin,npbo

     do j = 0,nhyper - 1
        
        a = x(j * (dim + 1) + 1)
        b = x(j * (dim + 1) + 2)
        c = x((j + 1) * (dim + 1))
        
        if ( abs(b) .le. 1.0D-06 ) then
           if ( abs(a) .gt. 1.0D-06 ) then
              WRITE(40,FMT=6050)(c / a),(c / a),- a / abs(a),(c / a)
           endif
        else
           if ( abs(a) .le. 1.0D-06 ) then
              WRITE(40,FMT=6060) (c / b),(c / b),- b / abs(b),(c / b)
           else
              if ( (a / b) .lt. 0.0D0 ) then
                 l = max(xmin,(- b * ymin + c) / a)
                 u = min(xmax,(- b * ymax + c) / a)
              else
                 l = max(xmin,(- b * ymax + c) / a)
                 u = min(xmax,(- b * ymin + c) / a)
              end if
                
              WRITE(40,FMT=6070) l,(a / b) * l,(c / b),u,(a / b) * u,(c / b),&
                                 - (a / SQRT(a ** 2 + b ** 2)),&
                                 - (b / SQRT(a ** 2 + b ** 2)),&
                                 (l + u) / 2.0D0,(a / b) * (l + u) / 2.0D0,&
                                 (c / b)
           end if
        end if

     end do
     
     WRITE(40,FMT=6090)

  end if

  call unload()

  WRITE(40,FMT=6090)
  CLOSE(40)

  ! NON-EXECUTABLE STATEMENTS

6040 FORMAT('beginfig(',I6,');',/)
6050 FORMAT(/,'draw (',F20.13,',ymax) * u -- (',F20.13,',ymin) * ' ,&
            'u;',/,'drawarrow ((0,0) * u -- (',F20.13,',0) * u) '  ,&
            'shifted ((',F20.13,',(ymin + ymax) / 2) * u);')
6060 FORMAT(/,'draw (xmax,',F20.13,') * u -- (xmin,',F20.13,') * ' ,&
            'u;',/,'drawarrow ((0,0) * u -- (0,',F20.13,') * u) '  ,&
            'shifted (((xmin + xmax) / 2,',F20.13,') * u);')
6070 FORMAT(/,'draw (',F20.13,',hyper(',F20.13,','    ,&
            F20.13,')) * u','-- (',F20.13,',hyper(',&
            F20.13,',',F20.13,')) * u;',/,'drawarrow ((0,0) * u --',&
            '(',F20.13,',',F20.13,') * u) shifted ((',&
            F20.13,',hyper(',F20.13,&
            ',',F20.13,')) * u);')
6080 FORMAT(/,/,'for i = 1 step 1 until ',I6,':',/,'draw pntot[i] ',&
            '* u withpen pencircle scaled 2pt withcolor 0.6 white;',/,&
            'endfor',/,/,'for i = 1 step 1 until ',I6,':',/,'draw ',&
            'pntin[i] * u withpen pencircle scaled 4pt withcolor ' ,&
            '0.6 green;',/,'endfor',/,/,'for i = 1 step ',&
            '1 until ',I6,':',/,'draw pntbo[i] * u withpen '       ,&
            'pensquare scaled 4pt withcolor 0.9 blue;',/,'endfor')
6090 FORMAT(/,'picture all;',/,'all := currentpicture;',/          ,&
            'clearit;',/,'clip all to box;',/,'draw all;',/,/      ,&
            'endfig;',/,'end;')

end subroutine finalize

!------------------------------------------------------------!
! SUBROUTINE EVAOBJFUN                                       !
!                                                            !
! Defines the objective function.                            !
!                                                            !
!------------------------------------------------------------!

subroutine evalobjfun(n,x,f,flag)

  use svmdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,n
  real(8) :: f

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: n,x
  intent(out) :: f,flag

  ! LOCAL SCALARS
  logical :: isin
  integer :: i,j
  real(8) :: tmp

  integer :: nin,ntot
  real(8),allocatable :: point(:)

  flag = 0

  f = 0.0D0

!!$  if ( dim .eq. 2 ) then
!!$
!!$     ntot = 100000
!!$     nin  = 0
!!$
!!$     allocate(point(2))
!!$
!!$     call RANDOM_SEED
!!$     do i = 1,ntot
!!$        call RANDOM_NUMBER(point)
!!$        
!!$        point(1) = xmin + point(1) * (xmax - xmin)
!!$        point(2) = ymin + point(2) * (ymax - ymin)
!!$        
!!$        isin = .true.
!!$        
!!$        do j = 0,nhyper - 1
!!$           tmp = DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim), &
!!$                 point) - x((j + 1) * (dim + 1))
!!$           if ( tmp .gt. 0.0D0 ) then
!!$              isin = .false.
!!$              exit
!!$           end if
!!$        end do
!!$
!!$        if ( isin ) then
!!$           nin =  nin + 1
!!$        end if
!!$
!!$     end do
!!$
!!$     f = (xmax - xmin) * (ymax - ymin) * (DBLE(nin) / DBLE(ntot))
!!$
!!$     deallocate(point)
!!$
!!$  else
!!$     flag = - 1
!!$  end if

  do i = 1,npot
     isin = .true.

     do j = 0,nhyper - 1
        tmp = DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim), &
              pntot(:,i)) - x((j + 1) * (dim + 1))
        if ( tmp .gt. 0.0D0 ) then
           isin = .false.
           exit
        end if
     end do

     if ( isin ) then
        f = f + 1.0D0
     end if
  end do

end subroutine evalobjfun

!------------------------------------------------------------!
! SUBROUTINE EVALCONSTR                                      !
!                                                            !
! Defines the contraints.                                    !
!                                                            !
!------------------------------------------------------------!

subroutine evalconstr(n,x,ind,c,flag)

  use svmdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,n
  real(8) :: c

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: ind,n,x
  intent(out) :: c,flag

  ! LOCAL SCALARS
  integer :: i,j

  flag = 0

  if ( ind .ge. 1 .and. ind .le. npbo ) then

     i =   ind
     c = 1.0D0
     do j = 0,nhyper - 1
        c = c * ( DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim), &
                  pntbo(:,i)) - x((j + 1) * (dim + 1)) )
     end do

  else if ( ind .ge. npbo + 1 .and. ind .le. 2 * npbo ) then

     i = ind - npbo
     c =   - 1.0D0
     do j = 0,nhyper - 1
        c = c * ( DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim), &
                  pntbo(:,i)) - x((j + 1) * (dim + 1)) )
     end do

  else if ( ind .ge. 2 * npbo + 1 .and. &
            ind .le. 2 * npbo + nhyper * npin ) then

     i = ((ind - 2 * npbo - 1) / nhyper) + 1
     j =    MOD((ind - 2 * npbo - 1),nhyper)

     c = x((j + 1) * (dim + 1)) - &
         DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim),pntin(:,i))
     c = SVMEPSIN ** 2 * &
         DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim), &
         x(j * (dim + 1) + 1:j * (dim + 1) + dim)) - c ** 2

  else if ( ind .ge. 2 * npbo + nhyper * npin .and. &
            ind .le. 2 * npbo + 2 * nhyper * npin ) then

     i = ((ind - (2 * npbo + nhyper * npin) - 1) / nhyper) + 1
     j =    MOD((ind - (2 * npbo + nhyper * npin) - 1),nhyper)
     c = DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim), &
         pntin(:,i)) - x((j + 1) * (dim + 1))

  else if ( ind .ge. 2 * npbo + 2 * nhyper * npin + 1 .and. &
            ind .le. 2 * npbo + nhyper * (npbo + 2 * npin) ) then

     i = ((ind - (2 * npbo + 2 * nhyper * npin) - 1) / nhyper) + 1
     j =      MOD(ind - (2 * npbo + 2 * nhyper * npin) - 1,nhyper)
     c = DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim), &
         pntbo(:,i)) - x((j + 1) * (dim + 1))

  else if ( ind .ge. 2 * npbo + (npbo + 2 * npin) * nhyper + 1 .and. &
            ind .le. 2 * npbo + (npbo + 2 * npin + 1) * nhyper ) then

     j = ind - (2 * npbo + (npbo + 2 * npin) * nhyper) - 1

     c =  1.0D0 - DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim),&
         x(j * (dim + 1) + 1:j * (dim + 1) + dim))

!!$  else if ( ind .ge. 2 * npbo + (npbo + npin + 1) * nhyper + 1 .and. &
!!$            ind .le. 2 * npbo + (npbo + npin + 2) * nhyper ) then
!!$
!!$     j = ind - 2 * npbo - (npbo + npin + 1) * nhyper - 1
!!$
!!$     c = 1.0D0 - DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim),&
!!$         x(j * (dim + 1) + 1:j * (dim + 1) + dim))

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

  use svmdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,jcnnz,n
  
  ! ARRAY ARGUMENTS
  integer :: jcvar(n)
  real(8) :: jcval(n),x(n)

  intent(in ) :: ind,n,x
  intent(out) :: flag,jcnnz,jcval,jcvar

  ! LOCAL SCALARS
  integer :: i,j,k
  real(8) :: tmpprod
  real(8),allocatable :: tmp(:)

  flag = 0

  if ( ind .ge. 1 .and. ind .le. npbo ) then
     allocate(tmp(nhyper))

     i = ind
     do j = 0,nhyper - 1
        tmp(j + 1) = DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim), &
                     pntbo(:,i)) - x((j + 1) * (dim + 1))
     end do

     jcnnz = n

     do j = 0,nhyper - 1
        tmpprod = PRODUCT(tmp(1:j)) * PRODUCT(tmp(j + 2:nhyper))

        do k = j * (dim + 1) + 1,(j + 1) * (dim + 1)
           jcvar(k) = k
        end do

        jcval(j * (dim + 1) + 1:(j + 1) * (dim + 1)) = pntbo(:,i) * tmpprod
        jcval((j + 1) * (dim + 1)                  ) =            - tmpprod
     end do

     deallocate(tmp)

  else if ( ind .ge. npbo + 1 .and. ind .le. 2 * npbo ) then
     allocate(tmp(nhyper))

     i = ind - npbo
     do j = 0,nhyper - 1
        tmp(j + 1) = DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim), &
                     pntbo(:,i)) - x((j + 1) * (dim + 1))
     end do

     jcnnz = n

     do j = 0,nhyper - 1
        tmpprod = - PRODUCT(tmp(1:j)) * PRODUCT(tmp(j + 2:nhyper))

        do k = j * (dim + 1) + 1,(j + 1) * (dim + 1)
           jcvar(k) = k
        end do

        jcval(j * (dim + 1) + 1:(j + 1) * (dim + 1)) = pntbo(:,i) * tmpprod
        jcval((j + 1) * (dim + 1)                  ) =            - tmpprod
     end do

     deallocate(tmp)

  else if ( ind .ge. 2 * npbo + 1 .and. &
            ind .le. 2 * npbo + nhyper * npin ) then
     i = ((ind - 2 * npbo - 1) / nhyper) + 1
     j =    MOD((ind - 2 * npbo - 1),nhyper)

     jcnnz = dim + 1

     tmpprod = x((j + 1) * (dim + 1)) - &
     DOT_PRODUCT(x(j * (dim + 1) + 1:j * (dim + 1) + dim),pntin(:,i))

     jcvar(1:dim + 1) = (/ (j * (dim + 1) + k, k = 1,dim + 1) /)
     jcval(1:    dim) = 2.0D0 * SVMEPSIN ** 2 *      &
          x(j * (dim + 1) + 1:j * (dim + 1) + dim) + &
          2.0D0 * tmpprod * pntin(:,i)
     jcval(  dim + 1) = - 2.0D0 * tmpprod

  else if ( ind .ge. 2 * npbo + nhyper * npin + 1 .and. &
            ind .le. 2 * npbo + 2 * nhyper * npin ) then
     i = ((ind - (2 * npbo + nhyper * npin) - 1) / nhyper) + 1
     j =    MOD((ind - (2 * npbo + nhyper * npin) - 1),nhyper)

     jcnnz = dim + 1

     jcvar(1:dim + 1) = (/ (j * (dim + 1) + k, k = 1,dim + 1) /)
     jcval(1:    dim) =                               pntin(:,i)
     jcval(  dim + 1) =                                  - 1.0D0

  else if ( ind .ge. 2 * npbo + 2 * nhyper * npin + 1 .and. &
            ind .le. 2 * npbo + nhyper * (npbo + 2 * npin) ) then
     i = ((ind - (2 * npbo + 2 * nhyper * npin) - 1) / nhyper) + 1
     j =      MOD(ind - (2 * npbo + 2 * nhyper * npin) - 1,nhyper)

     jcnnz = dim + 1

     jcvar(1:dim + 1) = (/ (j * (dim + 1) + k, k = 1,dim + 1) /)
     jcval(1:    dim) =                               pntbo(:,i)
     jcval(  dim + 1) =                                  - 1.0D0

  else if ( ind .ge. 2 * npbo + (npbo + 2 * npin) * nhyper + 1 .and. &
            ind .le. 2 * npbo + (npbo + 2 * npin + 1) * nhyper ) then

     j = ind - (2 * npbo + (npbo + 2 * npin) * nhyper) - 1

     jcnnz = dim

     jcvar(1:jcnnz) = (/ (j * (dim + 1) + k, k = 1,dim) /)
     jcval(1:jcnnz) = - 2.0D0 * x(j * (dim + 1) + 1:j * (dim + 1) + dim)

!!$  else if ( ind .ge. 2 * npbo + (npbo + npin + 1) * nhyper + 1 .and. &
!!$            ind .le. 2 * npbo + (npbo + npin + 2) * nhyper ) then
!!$
!!$     j = ind - 2 * npbo - (npbo + npin + 1) * nhyper - 1
!!$
!!$     jcnnz = dim
!!$
!!$     jcvar(1:jcnnz) = (/ (j * (dim + 1) + k, k = 1,dim) /)
!!$     jcval(1:jcnnz) = - 2.0D0 * x(j * (dim + 1) + 1:j * (dim + 1) + dim)

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

  use svmdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: iter,m,n
  real(8) :: gamma,infeas

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in) :: gamma,infeas,iter,m,n,x


  ! LOCAL SCALARS
  integer :: j
  real(8) :: a,b,c,l,u

  OPEN(40,FILE='running.mp',POSITION='APPEND')

  if ( dim .eq. 2 ) then

     WRITE(40,FMT=6040)iter
     WRITE(40,FMT=6080)min(SVMMAXDOTS,npot),npin,npbo

     do j = 0,nhyper - 1
        
        a = x(j * (dim + 1) + 1)
        b = x(j * (dim + 1) + 2)
        c = x((j + 1) * (dim + 1))
        
        if ( abs(b) .le. 1.0D-06 ) then
           if ( abs(a) .gt. 1.0D-06 ) then
              WRITE(40,FMT=6050)(c / a),(c / a),- a / abs(a),(c / a)
           endif
        else
           if ( abs(a) .le. 1.0D-06 ) then
              WRITE(40,FMT=6060) (c / b),(c / b),- b / abs(b),(c / b)
           else
              if ( (a / b) .lt. 0.0D0 ) then
                 l = max(xmin,(- b * ymin + c) / a)
                 u = min(xmax,(- b * ymax + c) / a)
              else
                 l = max(xmin,(- b * ymax + c) / a)
                 u = min(xmax,(- b * ymin + c) / a)
              end if
                
              WRITE(40,FMT=6070) l,(a / b) * l,(c / b),u,(a / b) * u,(c / b),&
                                 - (a / SQRT(a ** 2 + b ** 2)),&
                                 - (b / SQRT(a ** 2 + b ** 2)),&
                                 (l + u) / 2.0D0,(a / b) * (l + u) / 2.0D0,&
                                 (c / b)
           end if
        end if

     end do
     
     WRITE(40,FMT=6090)
     CLOSE(40)

  end if

  ! NON-EXECUTABLE STATEMENTS

6040 FORMAT('beginfig(',I6,');',/)
6050 FORMAT(/,'draw (',F20.13,',ymax) * u -- (',F20.13,',ymin) * ' ,&
            'u;',/,'drawarrow ((0,0) * u -- (',F20.13,',0) * u) '  ,&
            'shifted ((',F20.13,',(ymin + ymax) / 2) * u);')
6060 FORMAT(/,'draw (xmax,',F20.13,') * u -- (xmin,',F20.13,') * ' ,&
            'u;',/,'drawarrow ((0,0) * u -- (0,',F20.13,') * u) '  ,&
            'shifted (((xmin + xmax) / 2,',F20.13,') * u);')
6070 FORMAT(/,'draw (',F20.13,',hyper(',F20.13,','    ,&
            F20.13,')) * u','-- (',F20.13,',hyper(',&
            F20.13,',',F20.13,')) * u;',/,'drawarrow ((0,0) * u --',&
            '(',F20.13,',',F20.13,') * u) shifted ((',&
            F20.13,',hyper(',F20.13,&
            ',',F20.13,')) * u);')
6080 FORMAT(/,/,'for i = 1 step 1 until ',I6,':',/,'draw pntot[i] ',&
            '* u withpen pencircle scaled 2pt withcolor 0.6 white;',/,&
            'endfor',/,/,'for i = 1 step 1 until ',I6,':',/,'draw ',&
            'pntin[i] * u withpen pencircle scaled 4pt withcolor ' ,&
            '0.6 green;',/,'endfor',/,/,'for i = 1 step ',&
            '1 until ',I6,':',/,'draw pntbo[i] * u withpen '       ,&
            'pensquare scaled 4pt withcolor 0.9 blue;',/,'endfor')
6090 FORMAT(/,'picture all;',/,'all := currentpicture;',/          ,&
            'clearit;',/,'clip all to box;',/,'draw all;',/,/      ,&
            'endfig;',/)

end subroutine drawsAO
