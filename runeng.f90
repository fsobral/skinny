program runeng

  use skinny

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

  interface
     subroutine evalobjfun(n,x,f,flag)
       integer :: flag,n
       real(8) :: f
       real(8) :: x(n)

       intent(in ) :: n,x
       intent(out) :: f,flag
     end subroutine evalobjfun
  end interface
    
  ! PARAMETERS
  integer,parameter :: MAXTRIALS = 1

  ! LOCAL SCALARS
  integer :: error,fcnt,flag,ftot,gcnt,gtot,i,m,n,rsize,tbest,trial
  real(8) :: f,fbest,maxinfeas

  ! LOCAL ARRAYS
  real(8),allocatable :: l(:),u(:),x(:)
  real(8),allocatable :: xbest(:)
  integer,allocatable :: seed(:)

  fbest = SK_INFINITY
  ftot  = 0
  gtot  = 0

  ! ----------------------------- !
  ! Random numbers initialization !
  ! ----------------------------- !

  call RANDOM_SEED(SIZE=rsize)
  allocate(seed(rsize))

  ! ------------------------------- !
  ! Starts the multistart technique !
  ! ------------------------------- !

  do trial = 1,MAXTRIALS

     do i = 1,rsize
        seed(i) = 12345678 + trial * 100
     end do

     call RANDOM_SEED(PUT=seed)

     ! ------------------------- !
     ! Initializes the workspace !
     ! ------------------------- !

     !write(*,*) 'theta and tau:'
     !read(*,*) SK_GAMDEC,SK_REFFORT

     call initializes(n,m,flag)

     allocate(l(n),u(n),x(n), STAT=error)
     if ( error .ne. 0 ) then
        write(*,*) 'Memory problem.'
        stop
     end if

     call initializea(n,x,l,u,m,flag)

     if ( .not. ALLOCATED(xbest) ) then
        allocate(xbest(n), STAT=error)
        if ( error .ne. 0 ) then
           write(*,*) 'There is a memory problem...'
           stop
        end if
     end if

     ! ------------ !
     ! Calls Skinny !
     ! ------------ !

     call sksolve(n,x,l,u,m,1.0D-08,1.0D-03,1.0D-01,1.0D-01,1.0D0,2, &
          .true.,-1,f,maxinfeas,fcnt,gcnt,flag)

     ftot = ftot + fcnt
     gtot = gtot + gcnt
     
      ! Best solution found
     
     if ( maxinfeas .lt. 1.0D-04 ) then
        
        call evalobjfun(n,x,f,flag)
        
        if ( f .lt. fbest ) then
           tbest = trial
           fbest =     f
           xbest = x(1:n)
        end if
     end if

     write(*,FMT=2000)trial,f,fbest

!     call finalize(n,xbest,m,maxinfeas)
     call finalize(n,x,m,maxinfeas)

     ! -------------------- !
     ! Cleans the workspace !
     ! -------------------- !

     deallocate(l,u,x)
  end do

  write(*,FMT=2010)fbest,tbest,MAXTRIALS,ftot,gtot,xbest

  OPEN(75,FILE='runeng.out',POSITION='APPEND')
  write(75,FMT=3000)ftot
  CLOSE(75)

  ! Clear memory

  deallocate(seed,xbest)


  ! NON-EXECUTABLE STATEMENTS

2000 FORMAT(/,'Trial',48X,I10,/                           ,&
            'Objective function value found=',12X,E20.10,/,&
            'Best function value=',23X,E20.10,/)

2010 FORMAT(/,/,63('*'),/,'The process has terminated!',&
            /,63('*'),/,/,&
            'Best objective function value=',13X,E20.10,/,&
            'Trials=',25X,I15,'/',I15,/                  ,&
            'Number of f evaluations= ',23X,I15,/        ,&
            'Number of c evaluations= ',23X,I15,/,/      ,&
            'Final point:',/,3(1X,F20.10))

3000 FORMAT(I15)

end program runeng

