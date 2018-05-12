module skinny

  ! ATTRIBUTES
  
  ! Feasibility tolerance
  real(8),private :: SK_EPSFEA  = 1.0D-08
  ! Controls the restoration effort. If 0 means that restoration will be
  ! made on the original set.
  real(8),private :: SK_REFFORT = 1.0D-01

  ! PARAMETERS

  real(8),parameter :: SK_INFINITY = 1.0D+99
  real(8),parameter :: SK_ZERO   = 1.0D-12

  ! Initial domain relaxation
  real(8),parameter :: SK_INIGAM = 1.0D+01

  ! Decrease factor for the sufficient decrease parameter
  real(8),parameter :: SK_ETADEC = 1.0D-1

  ! Selects the type of the restoration
  ! 1 - ALGENCAN (derivatives of constraints are available)
  ! 2 - COMPASS SEARCH (derivatives of constr. are NOT available)
  integer,parameter :: SK_RESTTYPE = 2

  ! Defines multiple or single valued gamma
  logical,parameter :: SK_MULTGAMMA = .false.

  integer,parameter :: SK_PRINTE = 6
  integer,parameter :: SK_MAXITE = 1000
  integer,parameter :: SK_MAXFIT = 1000
  integer,parameter :: SK_MAXFCNT = 1000 * SK_MAXITE

  ! INTERFACES

  interface
     SUBROUTINE DFO( N,NX,X,LDX,FX,CONX,IFINIV,M,C,NCLIN,NCNLN,LB,    &
          UB,A,LDA,XNAMES,PNAME,CNAMES,IT,NF,INFO,MAXIT,MAXNF,STPCRTR,&
          DELMIN,STPTHR,CNSTOLP,DELTA,PP,SCALE,IOUTP,IPRINTP)

       integer :: N,M,NX,LDX,NCLIN,NCNLN,LDA,IT,NF,INFO,MAXIT,    &
            MAXNF,STPCRTR,SCALE,IOUTP,IPRINTP
       real(8) :: X(LDX * NX),FX(NX),LB(*),UB(*),A(LDA * N),C(*), &
            CONX(*),DELMIN,STPTHR,CNSTOLP,DELTA,PP
       logical :: IFINIV
       character * 256 :: XNAMES(N),PNAME,CNAMES(*)

       intent(in   ) :: N,NX,LDX,CONX,IFINIV,M,NCLIN,NCNLN,LB,UB,A,&
            LDA,XNAMES,PNAME,CNAMES,MAXIT,MAXNF,STPCRTR,DELMIN,STPTHR,&
            CNSTOLP,DELTA,PP,SCALE,IOUTP,IPRINTP
       intent(out  ) :: C,IT,NF,INFO
       intent(inout) :: X,FX
     END SUBROUTINE DFO
     SUBROUTINE bobyqa(n,npt,x,l,u,rhobeg,rhoend,iprint,maxfun,w,calfunb)
       integer :: iprint,maxfun,n,npt
       real(8) :: rhobeg,rhoend
       real(8) :: x(n),l(n),u(n),&
            w(((n + 1) * (n + 2) / 2 + 5) * ((n + 1) * (n + 2) / 2 + n) &
            + 3 * n * (n + 5) / 2 + ((n + 1) * (n + 2) / 2 + 5))
       external :: calfunb
       
       intent(in   ) :: iprint,l,maxfun,n,npt,rhobeg,rhoend,u
       intent(out  ) :: w
       intent(inout) :: x
     END SUBROUTINE bobyqa
     SUBROUTINE sds(n,x,delta,epsopt,maxsteps,maxit,objfun,constr,f,flag)
       integer :: flag,maxit,maxsteps,n
       real(8) :: delta,epsopt,f
       real(8) :: x(n)
       external :: constr,objfun
       
       intent(in   ) :: epsopt,delta,maxit,maxsteps,n
       intent(out  ) :: f,flag
       intent(inout) :: x
     END SUBROUTINE sds
  end interface

  interface
     SUBROUTINE restoration(nor,x,l,u,mor,epsfeas,verbose,infeas,flag)

       integer :: flag,mor,nor
       real(8) :: epsfeas,infeas,l(nor),u(nor),x(nor)
       logical :: verbose

       intent(in   ) :: epsfeas,l,mor,nor,u,verbose
       intent(out  ) :: flag,infeas
       intent(inout) :: x
     END SUBROUTINE restoration
  end interface

contains

  !------------------------------------------------------------!
  ! SUBROUTINE SKSOLVE                                         !
  !                                                            !
  ! This subroutine solves the following optimization problem: !
  !                                                            !
  !                min  f(x)                                   !
  !                s.t. g(x) <= 0                              !
  !                                                            !
  ! supposing that the function f(x) does not have             !
  ! derivatives, the derivatives of the constraints are        !
  ! available and the feasible set is skinny                   !
  !                                                            !
  ! At each step of the method, it solves the following 'fat'  !
  ! problem:                                                   !
  !                                                            !
  !                min  f(x)                                   !
  !                s.t. g(x) <= \gamma_k,                      !
  !                                                            !
  ! where \gamma_k goes to zero as k goes to infinity.         !
  !                                                            !
  ! ARGUMENTS                                                  !
  !                                                            !
  ! N: integer,scalar,INPUT                                    !
  !    Number of variables                                     !
  !                                                            !
  ! X(N) : real(8),array,INPUT/OUTPUT                          !
  !        On INPUT is the initial point, on OUTPUT is the     !
  !        (probable) solution                                 !
  !                                                            !
  ! L(N) : real(8),array,INPUT                                 !
  !        Lower bounds for the variables                      !
  !                                                            !
  ! U(N) : real(8),array,INPUT                                 !
  !        Upper bounds for the variables                      !
  !                                                            !
  ! M : integer,scalar,INPUT                                   !
  !     Number of constraints                                  !
  !                                                            !
  ! FEAS : real(8),scalar,INPUT                                !
  !        Feasibility tolerance                               !
  !                                                            !
  ! OPTI : real(8),scalar,INPUT                                !
  !        Optimality tolerance (used in the improvement step) !
  !                                                            !
  ! GAMMADEC : real(8),scalar,INPUT                            !
  !            Decrease factor of the relaxation parameter     !
  !                                                            !
  ! REFFORT : real(8),scalar,INPUT                             !
  !           Restoration effort in the restoration step       !
  !                                                            !
  ! DELTA : real(8),scalar,INPUT                               !
  !         Size of the region used in the sampling step       !
  !                                                            !
  ! SFACTOR : integer,scalar,INPUT                             !
  !           Multiplier related to the number of sampling     !
  !           points in the sampling step. In the implementa-  !
  !           tion we use SFACTOR * N                          !
  !                                                            !
  ! VERBOSE : logical,scalar,INPUT                             !
  !           If true, then prints information on screen       !
  !                                                            !
  ! UOPTMTYPE : integer,scalar,INPUT                           !
  !             The type of algorithm to be used in the impro- !
  !             vement step. Use -1 for all options            !
  !                                                            !
  ! F : real(8),scalar,OUTPUT                                  !
  !     Objective function value at the solution               !
  !                                                            !
  ! MAXINFEAS : real(8),scalar,OUTPUT                          !
  !             Sup-norm of the infeasibility                  !
  !                                                            !
  ! FCNT : integer,scalar,OUTPUT                               !
  !        Number of the objective function evaluations        !
  !                                                            !
  ! GCNT : integer,scalar,OUTPUT                               !
  !        Number of evaluations of all the constraints        !
  !                                                            !
  ! FLAG: integer,scalar,OUTPUT                                !
  !       Status returned by the method:                       !
  !                                                            !
  !       = 2 - reached the maximum number of function evals.  !
  !       = 1 - the relaxation parameter is too small          !
  !       = 0 - solution found                                 !
  !       < 0 - something happened when calling subroutines    !
  !                                                            !
  !------------------------------------------------------------!

  subroutine sksolve(n,x,l,u,m,feas,opti,gammadec,reffort,delta, &
       sfactor,verbose,uoptmtype,f,maxinfeas,fcnt,gcnt,flag)
    
    use engdata
    !use Nelder_Mead ! Uncomment this for Nelder-Mead
    
    implicit none
    
    ! SCALAR ARGUMENTS
    integer :: fcnt,flag,gcnt,m,n,sfactor,uoptmtype
    logical :: verbose
    real(8) :: delta,f,feas,gammadec,maxinfeas,opti,reffort

    ! ARRAY ARGUMENTS
    real(8) :: l(n),u(n),x(n)

    intent(in   ) :: delta,feas,gammadec,l,m,n,opti,reffort, &
         sfactor,u,verbose
    intent(out  ) :: f,fcnt,flag,gcnt,maxinfeas
    intent(inout) :: x

    ! INTERFACES
    interface
       subroutine drawsBI(n,x,m,gamma,infeas,iter)
         integer :: iter,m,n
         real(8) :: infeas
         real(8) :: gamma(m),x(n)
         
         intent(in) :: gamma,infeas,iter,m,n,x
       end subroutine drawsBI
    end interface

    interface
       subroutine drawsAR(n,x,m,gamma,infeas,iter)
         integer :: iter,m,n
         real(8) :: infeas
         real(8) :: gamma(m),x(n)
         
         intent(in) :: gamma,infeas,iter,m,n,x
       end subroutine drawsAR
    end interface

    interface
       subroutine drawsAO(n,x,m,gamma,infeas,iter)
         integer :: iter,m,n
         real(8) :: infeas
         real(8) :: gamma(m),x(n)
         
         intent(in) :: gamma,infeas,iter,m,n,x
       end subroutine drawsAO
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

    interface
       subroutine evalconstr(n,x,ind,c,flag)
         integer :: flag,ind,n
         real(8) :: c
         real(8) :: x(n)
         
         intent(in ) :: ind,n,x
         intent(out) :: c,flag
       end subroutine evalconstr
    end interface

    interface
       subroutine exbarr(x,f)
         real(8),intent(out) :: f
         real(8),intent(in ) :: x(:)
       end subroutine exbarr
    end interface

    interface
       subroutine calfunb(n,x,f)
         integer :: n
         real(8) :: f
         real(8) :: x(n)

         intent(in ) :: n,x
         intent(out) :: f
       end subroutine calfunb
    end interface

    interface
       subroutine calfunsds(n,x,f)
         integer :: n
         real(8) :: f
         real(8) :: x(n)

         intent(in ) :: n,x
         intent(out) :: f
       end subroutine calfunsds
    end interface

    interface
       subroutine calconsds(n,x,feasible)
         integer :: n
         logical :: feasible
         real(8) :: x(n)

         intent(in ) :: n,x
         intent(out) :: feasible
       end subroutine calconsds
    end interface

    interface
    end interface
    
    ! LOCAL SCALARS
    integer :: error,i,j,k,maxnpt,npt,optmtype,outflag,fav,it,nx
    real(8) :: c,epsopt,eta,etas,gammasn,gammin,mshref,norm2infeas,&
               rhobeg,seed,tmpscalar,pp
    character * 256 :: pname

    ! LOCAL ARRAYS
    real(8),allocatable :: fx(:),gamma(:),step(:),var(:),xnew(:),w(:),&
         lb(:),ub(:),constr(:)
    character * 256,allocatable :: cname(:),vname(:)

    if ( verbose ) write(*,FMT=0010)
    
    if ( uoptmtype .lt. 0 ) then
       write(*,FMT=0020)
       read (*,       *) optmtype
    else
       optmtype = uoptmtype
    end if

    ! Checks the derivatives of the constraints

    if ( engJccheck() ) then
       
       call jctest(n,m)
       write(*,FMT=0030)
       read(*,*)
       
    end if
    
    ! Stores the original n, m  and bounds to be used in the extreme
    ! barrier function
    call engSetM(m)
    call engSetN(n)

    call engSetBounds(l,u)

    if ( verbose ) write(*,FMT=1000)n,m

    ! ---------------------------------- !
    ! Initialization of global variables !
    ! ---------------------------------- !

    numfeval = 0
    numceval = 0

    ! -------------- !
    ! Initialization !
    ! -------------- !

    SK_EPSFEA = feas
    gammin = feas
    SK_REFFORT = reffort

    maxnpt = (n + 1) * (n + 2) / 2
    rhobeg = min(2.0D0,minval(u - l) / 2.0D0) / 2.0D0

    allocate(step(n),var(n),xnew(n),engXPrev(n),gamma(m), &
    w((maxnpt + 5) * (maxnpt + n) + 3 * n * (n + 5) / 2), &
    fx(1),vname(n),cname(0),lb(n + m),ub(n + m),constr(m),STAT=error)

    if ( error .ne. 0 ) then
       flag = - 1
       write(*,FMT=0040)
       return
    end if
    
    ! Just for Nelder-Mead
    step = min(1.0D0,(u - l) / 2.0D0 )

    ! Just for DFO
    nx    = 1
    pname = 'PROB'
    vname = 'var'
    pp = 1.0D+3
    lb(1:n) = l(1:n)
    ub(1:n) = u(1:n)

    ! Adjusts fixed penalty parameter
    !
    call evalobjfun(n,x,f,flag)
    !
    tmpscalar = - SK_INFINITY
    do i = 1,m
       call evalconstr(n,x,i,c,flag)
       tmpscalar = max(tmpscalar,abs(c))

       if ( SK_MULTGAMMA ) gamma(i) = max(abs(c),SK_INIGAM)
    end do       
    !
    EXBARRHO = max(1.0D0,abs(f)) / max(1.0D0,tmpscalar)
    !
    ! End of the adjustment

    ! Adjusts initial value for sufficient decrease
    etas = 1.0D-2 * max(1.0D0,abs(f))
    eta  = 1.0D-8

    ! A prime number for pseudo-random number generation
    seed = 12345701.0D0

    if ( .not. SK_MULTGAMMA ) gamma = max(tmpscalar,SK_INIGAM)

    gammasn = maxval(abs(gamma))
    mshref  = min(gammasn,gammasn ** 2.0D0)

    maxinfeas = verifyFeas(n,x,m,flag)

    if ( flag .ne. 0 ) then
       flag = - 1
       return
    end if

    if ( verbose ) write(*,FMT=1010)maxinfeas,gammasn

    call engSetGamma(m,gamma)

    ! --------------------- !
    ! End of initialization !
    ! --------------------- !


    call exbarr(x,f)

    flag = engGetFlag()
    if ( flag .ne. 0 ) return
    

    !-----------!
    ! Main loop !
    !-----------!

    do k = 1,SK_MAXITE


       call drawsBI(n,x,m,gamma,maxinfeas,k)
       
       
       if ( verbose ) then
          write(*,FMT=1020)k,f,min(n,SK_PRINTE),x(1:min(n,SK_PRINTE))
          write(*,FMT=1030)gammasn,DBLE(numfeval)
          if ( SK_MULTGAMMA ) write(*,FMT=1031) min(m,SK_PRINTE),&
                              gamma(1:min(m,SK_PRINTE))
       end if
       
       ! ---------------- !
       ! Restoration step !
       ! ---------------- !

       select case(SK_RESTTYPE)
       case(1) ! ALGENCAN
          if ( verbose ) then
             write(*,FMT=1040) '(ALGENCAN).'
          end if
       case(2) ! COMPASS SEARCH
          if ( verbose ) then
             write(*,FMT=1040) '(COMPASS).'
          end if
       end select

       call engSetGamma(m,(gamma - SK_EPSFEA) * SK_REFFORT)
       call restoration(n,x,l,u,m,SK_EPSFEA,verbose,norm2infeas,&
            outflag)
       call engSetGamma(m,gamma)
     
       maxinfeas = verifyFeas(n,x,m,flag)

       if ( flag .ne. 0 ) then
          flag = - 1
          return
       end if

       call exbarr(x,f)
       
       if ( verbose ) then
          flag = engGetFlag()
          if ( flag .ne. 0 ) then
             return
          end if

          write(*,FMT=1050)outflag,f,norm2infeas,maxinfeas
          write(*,FMT=2000)x(1:min(n,SK_PRINTE))
          write(*,FMT=1060)
       end if


       call drawsAR(n,x,m,gamma,maxinfeas,k)


       !-------------------!
       ! Optimization step !
       !-------------------!


100    continue

       !---------------------!
       ! 1. Improvement step !
       !---------------------!


       if ( optmtype .eq. 1 ) then ! NELDER-MEAD

          ! Remove these 3 lines for Nelder-Mead
          write(*,FMT=9000) 'Nelder-Mead'
          flag = - 1
          return

          ! Uncomment for Nelder-Mead

!!$          if ( verbose ) then
!!$             write(*,FMT=1070)
!!$          end if
!!$
!!$          epsopt = 1.0D-08
!!$          step = min(step,gamma)
!!$
!!$          call minim(x,step,n,f,n * 1000,-1,epsopt,1,0,0.0D0,var,exbarr,&
!!$               outflag)
!!$          
!!$          if ( verbose ) then
!!$             write(*,FMT=1080)outflag,f
!!$          end if


       else if ( optmtype .eq. 2 ) then ! BOBYQA

          ! Remove these 3 lines for BOBYQA
          write(*,FMT=9000) 'BOBYQA'
          flag = - 1
          return

          ! Uncomment for BOBYQA

!!$          if ( verbose ) write(*,FMT=2040)
!!$          
!!$          npt = min(maxnpt,N + 2)
!!$          epsopt = opti
!!$          
!!$          call bobyqa(n,npt,x,l,u,min(rhobeg,mshref),epsopt,0,&
!!$          n * SK_MAXFIT,w,calfunb)
!!$
!!$          call exbarr(x,f)
!!$
!!$          flag = engGetFlag()
!!$          if ( flag .ne. 0 ) return
!!$    
!!$          if ( verbose ) write(*,FMT=2050) f
          
       else if ( optmtype .eq. 3 ) then ! NELDER-MEAD + BOBYQA

          ! Remove these 3 lines for use this strategy
          write(*,FMT=9000) 'Nelder-mead and BOBYQA'
          flag = - 1
          return

          ! Uncomment for use this strategy

!!$          if ( verbose ) then
!!$             write(*,FMT=2080)
!!$          end if
!!$          
!!$          epsopt = max(1.0D-12,opti / DBLE(k))
!!$          
!!$          call minim(x,step,n,f,n * SK_MAXFIT,-1,epsopt,1,0,0.0D0,var,&
!!$               exbarr,outflag)
!!$          if ( verbose ) then
!!$             write(*,FMT=2090)f
!!$          end if
!!$          
!!$          npt = min(maxnpt,n + 2)
!!$          
!!$          call bobyqa(n,npt,x,l,u,min(rhobeg,mshref),epsopt,0,&
!!$          n * SK_MAXFIT,w,calfunb)
!!$          call exbarr(x,f)
!!$
!!$          if ( verbose ) then
!!$             write(*,FMT=2091)f
!!$          end if
          
       else if ( optmtype .eq. 4 ) then ! COMPASS
          
          if ( verbose ) write(*,FMT=3020)

          epsopt = opti
          call compass(n,x,l,u,calfunb,epsopt,n * SK_MAXFIT,f)

          if ( verbose ) write(*,FMT=3030) f

       else if ( optmtype .eq. 5 ) then ! SDS
          
          if ( verbose ) write(*,FMT=3040)

          call sds(n,x,rhobeg,opti,2,n * SK_MAXFIT,calfunsds,&
          calconsds,f,flag)

          if ( verbose ) write(*,FMT=3050) flag,f

       else if ( optmtype .eq. 6 ) then ! DFO

          ! Remove these 3 lines for DFO
          write(*,FMT=9000) 'DFO'
          flag = - 1
          return

          ! Uncomment for DFO

!!$          if ( verbose ) write(*,FMT=3060)
!!$
!!$          lb(n + 1:n + m) = - SK_INFINITY
!!$          ub(n + 1:n + m) =   gamma(1:m) - SK_EPSFEA
!!$
!!$          fx(1) = f
!!$
!!$          epsopt = max(opti,min(1.0D-1,mshref))
!!$
!!$          CALL DFO(n,nx,x,n,fx,var,.true.,0,constr,0,m,lb,ub,var,1, &
!!$               vname,pname,cname,it,fav,flag,SK_MAXITE,n * SK_MAXFIT, &
!!$               1,epsopt,1.0D-3,SK_EPSFEA,1.0D0,pp,0,6,-1)
!!$
!!$          call exbarr(x,f)
!!$
!!$          if ( verbose ) write(*,FMT=3070) flag,f

       end if

       ! ---------------- !
       ! 2. Sampling step !
       ! ---------------- !

       maxinfeas = verifyFeas(n,x,m,flag)

       call sampling(n,x,l,u,m,delta,gamma,f,                       &
       max(eta,min(maxinfeas,etas,1.0D+16 * eta)),seed,sfactor * n, &
       verbose,outflag)

       if ( outflag .eq. 0 ) then
          goto 100
       else if ( outflag .ne. 1 ) then
          flag = outflag
       end if

       call drawsAO(n,x,m,gamma,maxinfeas,k)


       if ( verbose ) then
          write(*,FMT=2000)x(1:min(n,SK_PRINTE))
          write(*,FMT=1090)
       end if



       ! ----------------- !
       ! STOPPING CRITERIA !
       ! ----------------- !

       if ( flag .ne. 0 ) then
          flag = - 1
          return
       end if

       ! ------------------------------------------- !
       ! 1. Stops if the optimized point is feasible !
       ! ------------------------------------------- !

       if ( maxinfeas .le. SK_EPSFEA ) then
          flag = 0
          exit
       end if

       ! -------------------------------------------- !
       ! 2. Stops if the set if very skinny (failure) !
       ! -------------------------------------------- !

       if ( gammasn .le. gammin &
            .or. gammasn - SK_EPSFEA .lt. 0.0D0 ) then
          flag = 1
          exit
       end if

       ! ------------------------------------------------------------------ !
       ! 3. Stops if has reached the maximum number of function evaluations !
       ! ------------------------------------------------------------------ !

       if ( numfeval .gt. SK_MAXFCNT ) then
          flag = 2
          exit
       end if


       ! -------------- !
       ! NEXT ITERATION !
       ! -------------- !

       if ( SK_MULTGAMMA ) then
          do i = 1,m
             call evalconstr(n,x,i,c,flag)
             if ( flag .ne. 0 ) return
             
             gamma(i) = min(max(0.0D0,c),gamma(i)) * gammadec
             gamma(i) = max(gammin,gamma(i))
          end do
          gammasn = maxval(abs(gamma))
       else
          gammasn = max(gammin,min(maxinfeas,gammasn) * gammadec)
          gamma   = gammasn
       end if

       mshref = max(gammin,gammasn ** 2.0D0)
       etas   = min(etas * SK_ETADEC,gammasn)
       eta    = eta / DBLE(k + 1)

       call engSetGamma(m,gamma)
       
    end do ! End main loop


    ! --------------------------------------------------- !
    ! Printing additional information and cleaning memory !
    ! --------------------------------------------------- !
       
    fcnt = numfeval
    gcnt = numceval

    if ( verbose ) then
       write(*,FMT=2010)f,maxinfeas,numfeval,numceval,flag,x
    end if

    deallocate(engXPrev,gamma,step,var,xnew,w)

    return
    
! NON-EXECUTABLE STATEMENTS

0010 FORMAT(/,'Welcome to Skinny!',/,'------------------',/)
0020 FORMAT('Please choose optimization type:',/,&
            '1. Nelder-Mead 2. BOBYQA 3. Nelder-Mead + BOBYQA',/,&
            '4. COMPASS Search 5. SDS 6. DFO',/)
0030 FORMAT(/,'Please, type any key to continue...',/)
0040 FORMAT(/,'Error: memory problems.',/)

1000 FORMAT('Number of variables:  ',1X,I3,/,&
            'Number of constraints:',1X,I3,/)
1010 FORMAT('Max infeasibility =',1X,E10.3,'. Setting ||GAMMA|| to',&
            E10.3,/)
1020 FORMAT(/,63('-'),/,'Iteration',1X,':',42X,I10,/,63('-'),/,/,&
            'Objective function value=',18X,E20.10,/&
            'Current point (first',1X,I31,1X,'elements):',/,&
            3(1X,E20.10))
1030 FORMAT(/,'||GAMMA||=',1X,1PE10.3,25X,'FEVAL=',1X,1PE10.3,/)
1031 FORMAT('Current GAMMA (first',1X,I31,1X,'elements):',/,&
            3X,6(1X,1PE9.2))

1040 FORMAT(/,/,'Entering Restoration phase',1X,A20)
1050 FORMAT(/,'Flag:                    ',28X,I10,/,&
              'Objective function value=',18X,1PE20.10,/&
              '2-norm (dist or infeas)= ',28X,1PE10.3,/,&
              'Max infeasibility=       ',28X,1PE10.3)
1060 FORMAT(/,'Leaving Restoration phase.',/)
1070 FORMAT(/,'Entering Optimization phase (NELDER MEAD).')
1080 FORMAT(/,'Flag of MINIM:           ',28X,I10,/,&
              'Objective function value=',18X,E20.10)
1090 FORMAT(/,'Leaving Optimization phase.',/)
2000 FORMAT(/,'Point:',/,3(1X,E20.10))
2010 FORMAT(/,/,63('*'),/,'Skinny has terminated!',&
            /,63('*'),/,/,&
            'Objective function value=',18X,E20.10,/,&
            'Max infeasibility=',25X,E20.10,/,&
            'Number of f evaluations= ',28X,I10,/,&
            'Number of c evaluations= ',28X,I10,/,&
            'Flag=',53X,I5,/&
            /,'Final point:',/,3(1X,F20.10))

!!$2020 FORMAT(/,'Entering Optimization phase (ORTHOMADS).')
!!$2030 FORMAT(/,'Flag of ORTHOMADS:       ',28X,I10,/,&
!!$              'Objective function value=',18X,E20.10)
2040 FORMAT(/,'Entering Optimization phase (BOBYQA).')
2050 FORMAT(/,'Flag of BOBYQA:          ',28X,/,&
              'Objective function value=',18X,E20.10)

!!$2060 FORMAT(/,'Entering Optimization phase (NM + ORTHOMADS).')
!!$2070 FORMAT(/,'Objective function value (NM)=',13X,E20.10)
!!$2071 FORMAT('Objective function value (OM)=',13X,E20.10)

2080 FORMAT(/,'Entering Optimization phase (NM + BOBYQA).')
2090 FORMAT(/,'Objective function value (NM)=',13X,E20.10)
2091 FORMAT('Objective function value (BO)=',13X,E20.10)

!!$3000 FORMAT(/,'Entering Optimization phase (BOBYQA + ORTHOMADS).')
!!$3010 FORMAT(/,'Objective function value (BO)=',13X,E20.10)
!!$3011 FORMAT('Objective function value (OM)=',13X,E20.10)

3020 FORMAT(/,'Entering Optimization phase (COMPASS).')
3030 FORMAT(/,'Objective function value=',13X,E20.10)

3040 FORMAT(/,'Entering Optimization phase (SDS).')
3050 FORMAT(/,'Flag of SDS:             ',28X,I10,/,&
              'Objective function value=',18X,E20.10)

3060 FORMAT(/,'Entering Optimization phase (DFO).')
3070 FORMAT(/,'Flag of DFO:             ',28X,I10,/,&
              'Objective function value=',18X,E20.10)

9000 FORMAT(/,42('!'),/,'You need to download specific code to use',&
            /,'the algorithm: ',A26,/,42('!'),/)

  end subroutine sksolve

  !------------------------------------------------------------!
  ! SUBROUTINE COMPASS                                         !
  !                                                            !
  ! This subroutine applies the compass search.                !
  !                                                            !
  !------------------------------------------------------------!
  
  subroutine compass(n,x,l,u,objfnctn,epsopt,maxfeval,f)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: maxfeval,n
    real(8) :: epsopt,f

    ! ARRAY ARGUMENTS
    real(8) :: l(n),u(n),x(n)

    intent(in   ) :: epsopt,l,maxfeval,n,u
    intent(out  ) :: f
    intent(inout) :: x

    ! INTERFACES
    interface
       subroutine objfnctn(n,x,f)
         integer,intent(in ) :: n
         real(8),intent(out) :: f
         real(8),intent(in ) :: x(n)
       end subroutine objfnctn
    end interface

    ! LOCAL SCALARS
    integer :: i,fcnt
    real(8) :: alpha,fk,old

    alpha = 1.0D0
    fcnt  =     0

    call objfnctn(n,x,f)
    fcnt = fcnt + 1

    do while ( alpha .gt. epsopt )

       do i = 1,n

          old = x(i)

          if ( old - alpha .ge. l(i) ) then
             x(i) = old - alpha
             call objfnctn(n,x,fk)
             fcnt = fcnt + 1

             if ( fk .lt. f ) then
                exit
             end if
          end if

          if ( old + alpha .le. u(i) ) then
             x(i) = old + alpha
             call objfnctn(n,x,fk)
             fcnt = fcnt + 1

             if ( fk .lt. f ) then
                exit
             end if
          end if

          x(i) = old

       end do

       if ( fk .lt. f ) then
          f     = fk
       else
          alpha = alpha / 2.0D0
       end if

       ! Reached the maximum number of function evaluations

       if ( fcnt .gt. maxfeval ) then
          return
       end if

    end do

  end subroutine compass

  !------------------------------------------------------------!
  ! SUBROUTINE JCTEST                                          !
  !                                                            !
  ! This subroutine uses finite differences to test the        !
  ! evaluation of the Jacobian of the constraints.             !
  !                                                            !
  !------------------------------------------------------------!

  subroutine jctest(n,m)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: m,n

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

    interface
       subroutine evaljacob(n,x,ind,jcvar,jcval,jcnnz,flag)
         integer :: flag,ind,jcnnz,n
         integer :: jcvar(n)
         real(8) :: jcval(n),x(n)

         intent(in ) :: ind,n,x
         intent(out) :: flag,jcnnz,jcval,jcvar
       end subroutine evaljacob
    end interface

    ! LOCAL SCALARS
    integer :: flag,i,j,jcnnz
    real(8) :: c,cm,cp,diff,h,maxdiff,xor

    ! LOCAL ARRAYS
    integer,dimension(n) :: jcvar
    real(8),dimension(n) :: jacobian,jcval,x

    h       = 1.0D-08
    maxdiff =   0.0D0

    call RANDOM_SEED
    call RANDOM_NUMBER(x)

    write(*,FMT=9000)'CNSTR','VRBLE','JCBN','FDIF','DIFF'

    do i = 1,m

       call evaljacob(n,x,i,jcvar,jcval,jcnnz,flag)

       jacobian                 =          0.0D0
       jacobian(jcvar(1:jcnnz)) = jcval(1:jcnnz)

       do j = 1,n
          xor = x(j)

          x(j) = xor + h
          call evalconstr(n,x,i,cp,flag)
          x(j) = xor - h
          call evalconstr(n,x,i,cm,flag)
          x(j) = xor

          diff = ABS(((cp - cm) / (2.0D0 * h)) - jacobian(j))
          maxdiff = max(maxdiff,diff)

          write(*,FMT=9010)i,j,jacobian(j),(cp - cm) / (2.0D0 * h),diff
               
       end do
       write(*,*)

    end do
    
    write(*,FMT=9020)maxdiff
    
9000 FORMAT(A5,1X,A5,1X,A12,1X,A12,1X,A12)
9010 FORMAT(I5,1X,I5,1X,D12.5,1X,D12.5,1X,D12.5)
9020 FORMAT('Max. diff.=',27X,D12.5,/)

  end subroutine jctest

  !------------------------------------------------------------!
  ! FUNCTION VERIFYFEAS                                        !
  !                                                            !
  ! This function returns the maximum infeasibility associated !
  ! to the vector 'x'. In case something goes wrong it returns !
  ! a negative number.                                         !
  !                                                            !
  !------------------------------------------------------------!

  function verifyFeas(n,x,m,flag)

    use engdata
    
    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,m,n

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    intent(in ) :: m,n,x
    intent(out) :: flag

    ! RETURN VALUE
    real(8) :: verifyFeas

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
    integer :: i
    real(8) :: c
    
    verifyFeas = max(0.0D0,engVerifyBounds(n,x))

    do i = 1,m
       call evalconstr(n,x,i,c,flag)
       verifyFeas = max(verifyFeas,c)
       
       if ( flag .ne. 0 ) then
          return
       end if
    end do

    flag     = 0
    numceval = numceval + m

  end function verifyFeas

  !------------------------------------------------------------!
  ! SUBROUTINE SAMPLING                                        !
  !                                                            !
  ! The main idea of this subroutine is to sample in a set of  !
  ! feasible random directions.                                !
  !                                                            !
  ! flag = 0 means that a smaller objective function was found !
  ! flag = 1 means that nothing was found                      !
  ! flag < 0 means error.                                      !
  !                                                            !
  !------------------------------------------------------------!

  subroutine sampling(n,x,l,u,m,delta,gamma,f,eta,seed,nsamples,&
  verbose,flag)

    use engdata

    implicit none

    integer :: flag,m,n,nsamples
    logical :: verbose
    real(8) :: gamma(m),l(n),u(n),x(n)
    real(8) :: delta,eta,f,seed

    interface
       function drandsc(ix)
         real(8) :: ix
         real(8) :: drandsc
       end function drandsc
    end interface

    interface
       subroutine exbarr(x,f)
         real(8),intent(out) :: f
         real(8),intent(in ) :: x(:)
       end subroutine exbarr
    end interface

    integer :: i,ns,outflag
    real(8) :: d(n),dnorm,fnew,norm2infeas,pi,rad,u1,u2,y(n)

    ns   = 0
    flag = 1
    pi   = acos(-1.0D0)


    if ( verbose ) write(*,FMT=010) nsamples,eta
    if ( nsamples .le. 0 ) return

!!$    ! Generating uniformly in the box [-delta,delta]^n
!!$
!!$500 dnorm = 0.0D0
!!$    do i = 1,n
!!$       d(i) = - delta + drandsc(seed) * 2.0D0 * delta
!!$       dnorm = dnorm + d(i) ** 2
!!$    end do


    ! This generation follows Muller(59) 'A note on a method
    ! for generating points uniformly on N-dimensional spheres' to
    ! uniformly random generation on the unitary hypersphere. For
    ! generation of the normal(0,1) distribution it uses the
    ! Box-Muller(58) algorithm. See 'A note on the generation of
    ! random normal deviates'.

500 dnorm = 0.0D0
    do i = 1,n,2
       u1 = drandsc(seed)
       u2 = drandsc(seed)

       d(i)  = sqrt(- 2.0D0 * log(u1)) * cos(2.0D0 * pi * u2)
       dnorm = dnorm + d(i) ** 2
       if ( i + 1 .le. n ) then
          d(i + 1) = sqrt(- 2.0D0 * log(u1)) * sin(2.0D0 * pi * u2)
          dnorm    = dnorm + d(i + 1) ** 2
       end if
    end do

    dnorm = sqrt(dnorm)
    rad = drandsc(seed)

    ! For unifom distribution inside the ball.
    rad = rad ** (1.0D0 / DBLE(n))
    rad = rad * (delta / dnorm)

    dnorm = 0.0D0
    do i = 1,n
       d(i) = d(i) * rad
       dnorm = dnorm + d(i) ** 2
    end do

    y = x + d
    
    call engSetGamma(m,gamma - SK_EPSFEA)
    call restoration(n,y,l,u,m,SK_EPSFEA,verbose,norm2infeas,outflag)
    call engSetGamma(m,gamma)

    call exbarr(y,fnew)
    flag = engGetFlag()

    if ( verbose ) write(*,FMT=011)ns + 1,fnew,sqrt(dnorm), &
         d(1:min(n,SK_PRINTE))

    if ( fnew .lt. f - eta ) then
       f = fnew
       x = y
       flag = 0

       return
    end if

    ns = ns + 1
    if ( ns .ge. nsamples) then
       flag = 1
       return
    end if

    goto 500

    ! NON-EXECUTABLE STATEMENTS
010 FORMAT(/,'Sampling',/,'Number of samples:',I5,1X,'Eta=',1X,1PE10.2,/)
011 FORMAT(I5,1X,'Functional Value=',1PE10.2,2X,'||d||=',1PE10.2, &
           /,3(1X,1PE20.10))

  end subroutine sampling

end module skinny

