!------------------------------------------------------------!
! SUBROUTINE RESTORATION                                     !
!                                                            !
! This subroutine returns a feasible point 'x'. Up to now it !
! uses ALGENCAN to solve the 'projection'                    !
!                                                            !
!                min  ||x - x_k||_2^2                        !
!                s.t. g(x) <= gamma_k                        !
!                     l <= x <= u                            !
!                                                            !
! ARGUMENTS                                                  !
!                                                            !
! NOR: integer,scalar,INPUT                                  !
!      Number of variables                                   !
!                                                            !
! X(NOR) : real(8),array,INPUT/OUTPUT                        !
!          On INPUT is the current point, on OUTPUT is the   !
!          (probably) feasible point                         !
!                                                            !
! L(NOR) : real(8),array,INPUT                               !
!          Lower bounds for the variables                    !
!                                                            !
! U(NOR) : real(8),array,INPUT                               !
!          Upper bounds for the variables                    !
!                                                            !
! MOR: integer,scalar,INPUT                                  !
!      Number of constraints                                 !
!                                                            !
! INFEAS: real(8),scalar,OUTPUT                              !
!         2 norm of infeasibility                            !
!                                                            !
! FLAG: integer,scalar,OUTPUT                                !
!       Status returned by the method:                       !
!       = 0 - solution found                                 !
!       > 0 - something happend                              !
!                                                            !
!                                                            !
! This subroutine uses the following modules                 !
!                                                            !
! - skinny (SK_PRINTE,SK_RESTTYPE,compass)                   !
! - engdata (engXPrev,ENG_FORRES)                            !
!                                                            !
!------------------------------------------------------------!

subroutine restoration(nor,x,l,u,mor,epsfeas,verbose,infeas,flag)

  use skinny, only: SK_PRINTE, SK_RESTTYPE, compass
  use engdata, only: engXPrev, ENG_FORRES

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,mor,nor
  real(8) :: epsfeas,infeas
  logical :: verbose

  ! ARRAY ARGUMENTS
  real(8) :: l(nor),u(nor),x(nor)

  intent(in   ) :: epsfeas,l,mor,nor,u,verbose
  intent(out  ) :: flag,infeas
  intent(inout) :: x

  ! INTERFACES
  interface
     subroutine algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda, &
          equatn,linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)

       logical :: checkder
       integer :: inform,iprint,m,n,ncomp
       real(8) :: cnorm,epsfeas,epsopt,f,nlpsupn,snorm
       logical :: coded(10),equatn(m),linear(m)
       real(8) :: l(n),lambda(m),u(n),x(n)

       intent(in   ) :: checkder,coded,epsfeas,epsopt,equatn,inform,&
            iprint,l,linear,m,n,ncomp,u
       intent(out  ) :: cnorm,f,lambda,nlpsupn,snorm
       intent(inout) :: x
     end subroutine algencan
     subroutine callfrest(n,x,f)
       integer :: n
       real(8) :: f
       real(8) :: x(n)

       intent(in ) :: n,x
       intent(out) :: f
     end subroutine callfrest
  end interface

  ! Uncomment for using Algencan

!!$  ! LOCAL SCALARS
!!$  logical :: checkder
!!$  integer :: iprint,m,n,ncomp
!!$  real(8) :: cnorm,epsopt,nlpsupn,snorm
!!$
!!$  ! LOCAL ARRAYS
!!$  logical :: coded(10),equatn(mor),linear(mor)
!!$  real(8) :: lambda(mor)

  select case(SK_RESTTYPE)
  case(1) ! ALGENCAN

     ! Remove these 3 lines and uncomment the code bellow
     write(*,*) 'Uncomment restoration.f90 for using Algencan.'
     flag = - 1
     return 

!!$     epsopt =  1.0D-08
!!$     iprint =       11
!!$     ncomp  = SK_PRINTE
!!$
!!$     checkder = .false.
!!$
!!$     n = nor
!!$     m = mor
!!$
!!$     coded(1:10) = .false.
!!$     coded(1: 2) =  .true.
!!$     coded(4: 5) =  .true.
!!$
!!$     lambda(1:m) =   0.0D0
!!$     equatn(1:m) = .false.
!!$     linear(1:m) = .false.
!!$
!!$     ! 'Global' variable for the objective function in restoration
!!$     ! phase. This variable is located in 'engdata'.
!!$     engXPrev = x
!!$
!!$     ! Tells that the functions used by ALGENCAN (evalf,evalc,...) have
!!$     ! to be changed to solve RESTORATION's problem.
!!$     ENG_FORRES = .true.
!!$
!!$     call algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda, &
!!$          equatn,linear,coded,checkder,infeas,cnorm,snorm,       &
!!$          nlpsupn,flag)

  case(2) ! COMPASS SEARCH

     call compass(nor,x,l,u,callfrest,epsfeas,1000,infeas)

     flag = 0

  case default
     flag = - 1
     return
  end select

end subroutine restoration

!------------------------------------------------------------!
! SUBROUTINE CALLFREST                                       !
!                                                            !
! This subroutine is used in the 'Restoration Phase'. The    !
! main reason of its existence is to find feasible points    !
! in a problem that does not have differentiable             !
! constraints.                                               !
!                                                            !
! This function is defined as the sum of the squares of the  !
! infeasibilities.                                           !
!                                                            !
!------------------------------------------------------------!

subroutine callfrest(n,x,f)

  use engdata

  ! SCALAR ARGUMENTS
  integer :: n
  real(8) :: f

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: n,x
  intent(out) :: f

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
  real(8),pointer :: gamma(:)

  flag  =              0
  gamma => engGetGamma()
  m     =      engGetM()

  f = 0.0D0
  do i = 1,m
     call evalconstr(n,x,i,c,flag)
     numceval = numceval + 1

     if ( flag .ne. 0 ) then
        call engSetFlag(flag)
        return
     end if

     f = f + max(0.0D0,c - gamma(i)) ** 2.0D0
  end do

end subroutine callfrest

!------------------------------------------------------------!
! SUBROUTINE ALGENCANB                                       !
!                                                            !
! This subroutine prepares ALGENCAN to solve a quadratic     !
! function subject to nonlinear general constraints defined  !
! by the user.                                               !
!                                                            !
! The quadratic function is constructed by BOBYQA (Powell,09)!
! as an approximation to the original function (also defined !
! by the user).                                              !
!                                                            !
! ALGENCANB is called by BOBYQA instead of subroutine        !
! TRSBOX, provided with the method, which minimizes the      !
! quadratic approximation subject to box constraints.        !
!                                                            !
! This subroutine uses the following modules                 !
!                                                            !
! - bobyqadata (GQ,HQ,PQ,XPTS,NPT,DELTA,FQ)                  !
! - engdata (engGetM,ENG_FORRES)                             !
!                                                            !
!------------------------------------------------------------!


! Uncomment this subroutine if you intend to use BOBYQA with Algencan.
! Ask for me the modified code in: fsobral at ime.unicamp.br.


!!$subroutine algencanb(n,npoints,points,xbase,x,l,u,h,hextra,g,fixf,&
!!$     trradius,xnew,d,gnew,dsq,crvmin)
!!$
!!$  use bobyqadata
!!$  use engdata
!!$
!!$  implicit none
!!$
!!$  ! SCALAR ARGUMENTS
!!$  integer :: n,npoints
!!$  real(8) :: crvmin,trradius,dsq,fixf
!!$
!!$  ! ARRAY ARGUMENTS
!!$  real(8) :: points(npoints,n),d(n),g(n),gnew(n),h(n * (n + 1) / 2), &
!!$       hextra(npoints),l(n),u(n),x(n),xnew(n),xbase(n)
!!$
!!$  intent(in ) :: fixf,h,hextra,l,n,npoints,points,trradius,u,x,xbase
!!$  intent(out) :: crvmin,d,dsq,xnew
!!$
!!$  ! INTERFACES
!!$  interface
!!$     subroutine algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda, &
!!$          equatn,linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)
!!$
!!$       logical :: checkder
!!$       integer :: inform,iprint,m,n,ncomp
!!$       real(8) :: cnorm,epsfeas,epsopt,f,nlpsupn,snorm
!!$       logical :: coded(10),equatn(m),linear(m)
!!$       real(8) :: l(n),lambda(m),u(n),x(n)
!!$
!!$       intent(in   ) :: checkder,coded,epsfeas,epsopt,equatn,inform,&
!!$            iprint,l,linear,m,n,ncomp,u
!!$       intent(out  ) :: cnorm,f,lambda,nlpsupn,snorm
!!$       intent(inout) :: x
!!$     end subroutine algencan
!!$  end interface
!!$
!!$  ! LOCAL SCALARS
!!$  integer :: error,i,ih,inform,iprint,j,m,ncomp
!!$  real(8) :: cnorm,epsfeas,epsopt,f,snorm,nlpsupn
!!$  logical :: checkder
!!$
!!$  ! LOCAL ARRAYS
!!$  logical coded(10)
!!$
!!$  !      write(*,*) 'Chamou Algencan'
!!$
!!$  ! Adds the trust region constraint
!!$  m = engGetM() + 1
!!$
!!$  ! -------------------------- !
!!$  ! Setting common information !
!!$  ! -------------------------- !
!!$
!!$  if ( (.not. allocated(GQ)) .or. (.not. allocated(PQ)) .or. &
!!$       (.not. allocated(XPTS)) .or. (.not. allocated(equatn)) .or. &
!!$       (.not. allocated(linear)) .or. (.not. allocated(lambda)) ) then
!!$
!!$     allocate(GQ(n),HQ(n * (n + 1) / 2),PQ(npoints),XPTS(npoints,n), &
!!$          equatn(m),linear(m),lambda(m), STAT=error)
!!$
!!$     if ( error .ne. 0 ) then
!!$        write(*,*) 'Error in memory allocation (FOR BOBYQA).'
!!$        stop
!!$     end if
!!$  end if
!!$
!!$  NPT   = npoints
!!$  FQ    = fixf
!!$  DELTA = trradius
!!$
!!$ !!$  write(*,*) 'NPT=',NPT,'FOPT=',fixf,'delta=',delta
!!$ !!$  write(*,*) (g(i),i = 1,n)
!!$ !!$  write(*,*) 'x=',(x(i) + xbase(i),i = 1,n)
!!$ !!$  write(*,*) 'l=',(l(i),i = 1,n)
!!$ !!$  write(*,*) 'u=',(u(i),i = 1,n)
!!$
!!$  do i = 1,n * (n + 1) / 2
!!$     HQ(i) = h(i)
!!$  end do
!!$  do i = 1,n
!!$     GQ(i) = g(i)
!!$  end do
!!$  do i = 1,npoints
!!$     do j = 1,n
!!$        XPTS(i,j) = points(i,j)
!!$     end do
!!$     !     write(*,*)'XPTS',i,(XPTS(i,j),j = 1,n)
!!$     PQ(i) = hextra(i)
!!$  end do
!!$
!!$ !!$  ih = 0
!!$ !!$  write(*,*) 'HQ'
!!$ !!$  do i = 1,n
!!$ !!$     write(*,*) 'IH',ih,(HQ(j), j = ih + 1,ih + i)
!!$ !!$     ih = ih + i
!!$ !!$  end do
!!$
!!$  ! ---------------------------- !
!!$  ! Setting ALGENCAN information !
!!$  ! ---------------------------- !
!!$
!!$  epsfeas  = 1.0D-08
!!$  epsopt   = 1.0D-08
!!$  iprint   =      11
!!$  ncomp    =       6
!!$  checkder =  .false.
!!$
!!$  equatn(1:m) = .false.
!!$  linear(1:m) = .false.
!!$  lambda(1:m) =   0.0D0
!!$  d(1:n)      =   0.0D0
!!$
!!$  coded(1: 5) =  .true.
!!$  coded(6:10) = .false.
!!$
!!$  ! -------------- !     
!!$  ! Calls ALGENCAN !
!!$  ! -------------- !     
!!$
!!$  ! To be used by the constraints.
!!$  engXPrev = x(1:n) + xbase(1:n)
!!$
!!$  ! Tells that the functions used by ALGENCAN (evalf,evalc,...) have
!!$  ! to be changed to solve BOBYQA's problem.
!!$  ENG_FORRES = .false.
!!$
!!$  call algencan(epsfeas,epsopt,iprint,ncomp,n,d,l,u,m,lambda,&
!!$       equatn,linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)
!!$
!!$  !      write(*,*) 'Depois de ALGENCAN'
!!$
!!$  ! ----------------- !     
!!$  ! Final information !
!!$  ! ----------------- !     
!!$
!!$  dsq = 0.0D0
!!$  do i = 1,n
!!$     xnew(i) = x(i) + d(i)
!!$     dsq     = dsq + d(i) ** 2
!!$  end do
!!$  !
!!$  if ( SQRT(dsq) .eq. DELTA ) then
!!$     crvmin =   0.0D0
!!$  else
!!$     crvmin = - 1.0D0
!!$  end if
!!$  !
!!$  call evalg(n,d,gnew,inform)
!!$
!!$ !!$  ! --------------- !
!!$ !!$  ! Memory cleaning !
!!$ !!$  ! --------------- !
!!$ !!$
!!$ !!$  deallocate(GQ,HQ,PQ,XPTS,equatn,linear,lambda)
!!$
!!$end subroutine algencanb
!!$

subroutine evalf(n,x,f,flag)

  use bobyqadata
  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,n
  real(8) :: f

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: n,x
  intent(out) :: f,flag

  ! LOCAL SCALARS
  integer :: i,ih,j,k
  real(8) :: temp


  !  write(*,*) 'Entrou evalf',x

  if ( ENG_FORRES ) then
     f = 5.0D-01 * DOT_PRODUCT(x - engXPrev,x - engXPrev)
  else
     f = 0.0D0
     ih = 0

     do j = 1,n
        do i = 1,j
           ih = ih + 1
           if ( i .lt. j ) then
              f = f + x(j) * HQ(ih) * x(i)
           end if
           f = f + x(i) * HQ(ih) * x(j)
        end do
     end do

     do k = 1,NPT
        if ( PQ(k) .ne. 0.0D0 ) then
           temp = 0.0D0
           do j = 1,n
              temp = temp + XPTS(k,j) * x(j)
           end do
           temp = temp * PQ(k)
           do i = 1,n
              f = f + temp * XPTS(k,i) * x(i)
           end do
        end if
     end do

     f = 0.5D0 * f

     do i = 1,n
        f = f + GQ(i) * x(i)
     end do

     f = f + FQ

  end if

  flag = 0

end subroutine evalf

subroutine evalg(n,x,g,flag)

  use bobyqadata
  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,n

  ! ARRAY ARGUMENTS
  real(8) :: g(n),x(n)

  intent(in ) :: n,x
  intent(out) :: flag,g

  ! LOCAL SCALARS
  integer :: i,ih,j,k,m
  real(8) :: temp

  if ( ENG_FORRES ) then
     g = (x - engXPrev)
  else
     ih = 0

     do j = 1,n
        g(j) = GQ(j)
        do i = 1,j
           ih = ih + 1
           if ( i .lt. j ) then
              g(j) = g(j) + HQ(ih) * x(i)
           end if
           g(i) = g(i) + HQ(ih) * x(j)
        end do
     end do

     do k = 1,NPT
        if ( PQ(k) .ne. 0.0D0 ) then
           temp = 0.0D0
           do j = 1,n
              temp = temp + XPTS(k,j) * x(j)
           end do
           temp = temp * PQ(k)
           do i = 1,n
              g(i) = g(i) + temp * XPTS(k,i)
           end do
        end if
     end do
  end if

  flag = 0

end subroutine evalg

subroutine evalh(n,x,hlin,hcol,hval,hnnz,flag)

  use bobyqadata
  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,n,hnnz

  ! ARRAY ARGUMENTS
  integer :: hcol(*),hlin(*)
  real(8) :: hval(*),x(n)

  intent(in ) :: n,x
  intent(out) :: flag,hcol,hlin,hnnz,hval


  ! LOCAL SCALARS
  integer :: i,j,k
  real(8) :: temp

  if ( ENG_FORRES ) then
     flag = - 1
  else

     ! Algencan uses the lower part of the Hessian, so we have to
     ! invert the indices.

     hnnz = 0

     do j = 1,n
        do i = 1,j
           hnnz = hnnz + 1

           hlin(hnnz) = j
           hcol(hnnz) = i
           hval(hnnz) = HQ(hnnz)

           temp = 0.0D0
           do k = 1,NPT
              temp = temp + XPTS(k,i) * (PQ(k) * XPTS(k,j))
           end do

           hval(hnnz) = hval(hnnz) + temp
        end do
     end do

     flag = 0
  end if

end subroutine evalh

subroutine evalc(n,x,ind,c,flag)

  use bobyqadata
  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,n
  real(8) :: c

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: ind,n,x
  intent(out) :: c,flag

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
  integer :: i,m
  real(8),pointer :: gamma(:)

  !  write(*,*) 'Entrou evalc'

  gamma => engGetGamma()
  m     =      engGetM()
  flag  =              0

  if ( ENG_FORRES ) then
     call evalconstr(n,x,ind,c,flag)
     numceval = numceval + 1

     if ( flag .ne. 0 ) then
        return
     end if

     c = c - gamma(ind)
  else 
     if ( ind .ge. 1 .and. ind .le. m ) then

        call evalconstr(n,engXPrev + x,ind,c,flag)
        numceval = numceval + 1

        if ( flag .ne. 0 ) then
           return
        end if

        c = c - gamma(ind)
     else if ( ind .eq. m + 1 ) then        
        c = - delta ** 2
        do i = 1,n
           c = c + x(i) ** 2
        end do
     else
        flag = - 1
     end if
  end if

  !  write(*,*) 'Saiu evalc'

end subroutine evalc

subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,ind,jcnnz,n

  ! ARRAY ARGUMENTS
  integer :: jcvar(n)
  real(8) :: jcval(n),x(n)

  intent(in ) :: ind,n,x
  intent(out) :: flag,jcnnz,jcval,jcvar

  ! INTERFACES
  interface
     subroutine evaljacob(n,x,ind,jcvar,jcval,jcnnz,flag)
       integer :: flag,ind,jcnnz,jcvar(n),n
       real(8) :: jcval(n),x(n)

       intent(in ) :: ind,n,x
       intent(out) :: flag,jcnnz,jcval,jcvar
     end subroutine evaljacob
  end interface

  ! LOCAL SCALARS
  integer :: i,m

  flag = 0
  m    = engGetM()

  !  write(*,*) 'Entrou evaljac'

  if ( ENG_FORRES ) then

     call evaljacob(n,x,ind,jcvar,jcval,jcnnz,flag)

  else 
     if ( ind .ge. 1 .and. ind .le. m ) then

        call evaljacob(n,engXPrev + x,ind,jcvar,jcval,jcnnz,flag)

     else if ( ind .eq. m + 1 ) then
        jcnnz = n
        do i = 1,n
           jcvar(i) = i
           jcval(i) = 2.0D0 * x(i)
        end do
     else
        flag = - 1
     end if
  end if

  !  write(*,*) 'Saiu evaljac'

end subroutine evaljac

subroutine evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag)

  implicit none

  !     SCALAR ARGUMENTS
  integer :: flag,hcnnz,ind,n

  !     ARRAY ARGUMENTS
  integer :: hccol(:),hclin(:)
  real(8) :: hcval(:),x(n)

  flag = - 1

end subroutine evalhc

!!$subroutine evalfc(n,x,f,m,c,flag)
!!$
!!$  implicit none
!!$
!!$  !     SCALAR ARGUMENTS
!!$  integer :: flag,m,n
!!$  real(8) :: f
!!$
!!$  !     ARRAY ARGUMENTS
!!$  real(8) :: c(m),x(n)
!!$
!!$  flag = - 1
!!$
!!$end subroutine evalfc
!!$
!!$subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag)
!!$
!!$  implicit none
!!$
!!$  !     SCALAR ARGUMENTS
!!$  integer :: flag,jcnnz,m,n
!!$
!!$  !     ARRAY ARGUMENTS
!!$  integer :: jcfun(:),jcvar(:)
!!$  real(8) :: g(n),jcval(:),x(n)
!!$
!!$  flag = - 1
!!$
!!$end subroutine evalgjac
!!$
!!$subroutine evalhl(n,x,m,lambda,scalef,scalec,hllin,hlcol,hlval, &
!!$     hlnnz,flag)
!!$
!!$  implicit none
!!$
!!$  !     SCALAR ARGUMENTS
!!$  integer :: flag,hlnnz,m,n
!!$  real(8) :: scalef
!!$
!!$  !     ARRAY ARGUMENTS
!!$  integer :: hlcol(:),hllin(:)
!!$  real(8) :: hlval(:),lambda(m),scalec(m),x(n)
!!$
!!$  flag = - 1
!!$
!!$end subroutine evalhl

subroutine evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

  implicit none

  !     SCALAR ARGUMENTS
  logical ::goth
  integer :: flag,m,n
  real(8) :: sf

  !     ARRAY ARGUMENTS
  real(8) :: hp(n),lambda(m),p(n),sc(m),x(n)

  flag = - 1

end subroutine evalhlp
