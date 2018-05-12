!------------------------------------------------------------!
! SUBROUTINE FUN                                             !
!                                                            !
!------------------------------------------------------------!

subroutine fun(n,m,x,val,c,iferr)
  
  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: n,m
  logical :: iferr
  real(8) :: val

  ! ARRAY ARGUMENTS
  real(8) :: c(m),x(n)
  
  intent(in ) :: n,m,x
  intent(out) :: val,c,iferr

  ! INTERFACES
  interface
     subroutine evalobjfun(n,x,f,flag)
       integer :: flag,n
       real(8) :: f
       real(8) :: x(n)

       intent(in ) :: n,x
       intent(out) :: f,flag
     end subroutine evalobjfun

     subroutine evalconstr(n,x,ind,c,flag)
       integer :: flag,ind,n
       real(8) :: c
       real(8) :: x(n)

       intent(in ) :: ind,n,x
       intent(out) :: c,flag
     end subroutine evalconstr
  end interface

  ! LOCAL SCALARS
  integer :: i,easy,flag
  real(8) :: cc,pp

  iferr = .false.
  flag = 0
  easy = engGetM()

  numfeval = numfeval + 1
  call evalobjfun(n,x,val,flag)

  if ( flag .ne. 0 ) then
     call engSetFlag(flag)
     return
  end if

  ! We do not deal with the difficult constraints.

  pp = 0.0D0

  do i = 1,easy
     numceval = numceval + 1
     call evalconstr(n,x,i,cc,flag)

     if ( flag .ne. 0 ) then
        call engSetFlag(flag)
        return
     end if
     
     pp = pp + max(0.0D0,cc) ** 2
  end do

  val = val + EXBARRHO * pp

  call engSetFlag(flag)

end subroutine fun

!------------------------------------------------------------!
! SUBROUTINE EASYCON                                         !
!                                                            !
!------------------------------------------------------------!

subroutine easycon(n,x,ncnln,c)

  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: n,ncnln

  ! ARRAY ARGUMENTS
  real(8) :: c(ncnln),x(n)

  intent(in ) :: n,ncnln,x
  intent(out) :: c

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
  integer :: i,flag,m

  m = engGetM()
  flag = 0

  do i = 1,m
     numceval = numceval + 1
     call evalconstr(n,x,i,c(i),flag)

     if ( flag .ne. 0 ) then
        call engSetFlag(flag)
        return
     end if
  end do

  call engSetFlag(flag)
  
end subroutine easycon

!------------------------------------------------------------!
! SUBROUTINE EASYJAC                                         !
!                                                            !
!------------------------------------------------------------!

subroutine easyjac(n,x,ncnln,nrowc,cjac)

  use engdata

  implicit none

  ! SCALAR ARGUMENTS
  integer :: ncnln,n,nrowc
  
  ! ARRAY ARGUMENTS
  real(8) :: x(n),cjac(nrowc,n)

  intent(in ) :: n,ncnln,nrowc,x
  intent(out) :: cjac

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
  integer i,j,m,flag,jcnnz

  ! LOCAL ARRAYS
  integer,allocatable :: jcvar(:)
  real(8),allocatable :: jcval(:)

  m = engGetM()
  flag = 0

  allocate(jcvar(n),jcval(n))

  do i = 1,m
     cjac(i,1:n) = 0.0D0
     call evaljacob(n,x,i,jcvar,jcval,jcnnz,flag)
 
     if ( flag .ne. 0 ) then
        call engSetFlag(flag)
        return
     end if
    
     cjac(i,jcvar(1:jcnnz)) = jcval(1:jcnnz)
!     do j = 1,jcnnz
!        cjac(i,jcvar(j)) = jcval(j)
!     end do
  end do

  deallocate(jcvar,jcval)

  call engSetFlag(flag)

end subroutine easyjac

!------------------------------------------------------------!
! SUBROUTINE EASYHESS                                        !
!                                                            !
!------------------------------------------------------------!

subroutine easyhess(k,n,x,ncnln,nrowc,chess,lambda)
  
  implicit none

  ! SCALAR ARGUMENTS
  integer :: ncnln,n,nrowc,k
  real(8) :: lambda

  ! ARRAY ARGUMENTS
  real(8) :: x(n),chess(nrowc,n)

  intent(in   ) :: k,n,x,ncnln,nrowc,lambda
  intent(inout) :: chess

end subroutine easyhess
