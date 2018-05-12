program runhstests

  use skinny
  use hsdata

  implicit none

  ! INTERFACES
  interface
     subroutine initializes(n,m,flag)
       integer,intent(out) :: flag,m,n
     end subroutine initializes
  end interface

  interface
     subroutine initializea(n,x,l,u,m,flag)
       integer :: flag,m,n
       real(8) :: l(n),u(n),x(n)
       
       intent(in ) :: m,n
       intent(out) :: flag,l,u,x
     end subroutine initializea
  end interface

  interface
     subroutine finalize(n,x,m,maxinfeas)
       integer :: m,n
       real(8) :: maxinfeas
       real(8) :: x(n)
       
       intent(in) :: m,maxinfeas,n,x
     end subroutine finalize
  end interface
    
  ! LOCAL SCALARS
  integer :: error,fcnt,flag,gcnt,i,m,n,it,sfactor
  real(8) :: f,gammadec,radius,reffort,maxinfeas,feas,opti

  ! LOCAL ARRAYS
  real(8),allocatable :: l(:),u(:),x(:)

  ! COMMON SCALARS
  integer :: NTP

  ! COMMON ARRAYS
  COMMON/L8/NTP

!  write(*,FMT=2000)'PROB','N','INEQ','EQ','REL. DIF. F','NoS','?','F',&
!       'FO','MAXINFEAS','FEVAL','CEVAL'

  do it = 1,1!395

     read(*,*) NTP

     if ( NTP .eq. 82 .or. NTP .eq. 94 .or. NTP .eq. 115 .or. &
!          NTP .eq. 67 .or. & ! Wrong derivatives
          (NTP .gt. 119 .and. NTP .lt. 200) ) then
        cycle
     end if


     ! ------------------------- !
     ! Initializes the workspace !
     ! ------------------------- !

     call initializes(n,m,flag)

     allocate(l(n),u(n),x(n), STAT=error)

     if ( error .ne. 0 ) then
        write(*,*) 'Memory problem.'
        stop
     end if

     call initializea(n,x,l,u,m,flag)


     ! ------------ !
     ! Calls Skinny !
     ! ------------ !

     gammadec = 1.0D-1
     reffort = 1.0D-1
     feas = 1.0D-8
     opti = 1.0D-8
     sfactor = 2
     radius = 1.0D0

     call sksolve(n,x,l,u,m,feas,opti,gammadec,reffort,radius,sfactor, &
          .false.,6,f,maxinfeas,fcnt,gcnt,flag)
 
     open(75,FILE='runhs.out')

     call finalize(n,x,m,maxinfeas)

     write( *,FMT=3000)maxinfeas,fcnt,gcnt,flag
     write(75,FMT=3010)fcnt
     close(75)

     ! -------------------- !
     ! Cleans the workspace !
     ! -------------------- !

     deallocate(l,u,x)

  end do

  ! NON-EXECUTABLE STATEMENTS

2000 FORMAT(/,A4,1X,A4,1X,A4,1X,A4,5X,A12,1X,A3,1X,A1,1X,A12,1X,A12,1X,A9,&
            1X,A15,1X,A15)
3000 FORMAT(1X,E9.2,1X,I15,1X,I15,1X,I3)
3010 FORMAT(1X,I15)

end program runhstests
