C     MINTR_ALGENCAN
C
C     This file contains the interface between Algencan[1] and DFO.  It
C     is expected that Algencan will solve the quadratic model subject
C     to the box, easy and the quadratic model of the hard constraints.
C
C     It is important to rename the subroutine called SCL because it
C     conflicts with a different function used in Algencan, which has
C     the same name.
C
C     Reference
C
C     [1] R. Andreani, E. G. Birgin, J. M. Marti'nez and
C     M. L. Schuverdt, "On Augmented Lagrangian methods with general
C     lower-level constraints", SIAM Journal on Optimization 18,
C     pp. 1286-1309, 2007.
C
C     The author of this interface was Francisco N. C. Sobral
C     fsobral at ime dot unicamp dot br
C     June, 13th 2011.

      
      SUBROUTINE MINTR( N   , X0   , MVAL  , DELTA, LWRBND, UPRBND,
     *                  A   , LDA  , NCLIN , NCNLN, WRK   , LWRK  ,
     *                  IWRK, LIWRK, INFORM, METHOD )
    

C
C  *******************************************************************
C  THIS SUBROUTINE FINDS THE MINIMUM OF THE QUADRATIC MODEL WITHIN THE
C  GIVEN REGION. THE REGION IS DEFINED BY THE INTERSECTION OF THE
C  TRUST REGION OF RADIUS DELTA AND THE ANALYTICALLY DEFINED FEASIBLE 
C  SET OF THE ORIGINAL PROBLEM.
C
C                         T       T
C                MIN [GMOD X+0.5*X HMOD X]
C         
C                  X0-DELTA <=  X   <= XO+DELTA
C        S.T.
C                               / X  \
C                        LB <= ( AX   ) <= UB
C                               \C(X)/
C  PARAMETERS:
C
C   N      (INPUT)  DIMENTION OF THE PROBLEM
C
C   X0     (INPUT)  ARRAY OF LENGTH N CONTAINING THE CENTER OF THE TRUST
C                   REGION
C          (OUTPUT) CONTAINS THE OPTIMAL POINT FOR THE MODEL
C
C   MVAL   (OUTPUT) THE VALUE OF THE MODEL AT THE OPTIMAL POINT
C
C   DELTA  (INPUT)  TRUST REGION RADIUS
C
C   LWRBND (INPUT)  ARRAY OF LENGHT N+NCLIN+NCNLN OF LOWER BOUNDS 
C
C   UPRBND (INPUT)     ''       ''         ''        UPPER   ''
C
C   NCLIN  (INPUT)  NUMBER OF LINEAR ANALYTIC CONSTRAINTS
C
C   A      (INPUT)  (LDA X N) MATRIX OF LINEAR ANALYTIC CONSTRAINTS
C  
C   NCNLN  (INPUT)  NUMBER OF NOLINEAR INEQUALITIES (DIFFICULT AND EASY)
C
C   WRK             REAL SPACE WORKING ARRAY
C
C   IWRK            INTEGER SPACE WORKING ARRAY
C
C   INFORM (OUTPUT) INFORMATION ON EXIT
C              0    SUCCESSFUL MINIMIZATION
C              1    THE DERIVATIVES OF THE CONSTRAINT OR SOME PARAMETER
C                   SET BY THE USER IS INCORRECT
C              2    MINIMIZATION FAILED FOR SOME REASON
C
C  
C   METHOD (INPUT)  METHOD FOR HANDLING CONSTRAINTS
C              1    MINIMIZE MODEL OF OBJECTIVE S.T. MODELS OF CONSTRAINTS
C              2    MINIMIZE MERIT FUNCTION OF THE MODELS OF CON AND OBJ
C             3,4   MINIMIZE MODEL OF A MERIT FUNCTION (EXACT OR QUAD)
C
C  **********************************************************************
C

      implicit none


      INTEGER           N ,  NCLIN , NCNLN, LIWRK, LWRK, IWRK(LIWRK),
     +                  LDA, INFORM, METHOD

 
      DOUBLE PRECISION  X0(N), MVAL, DELTA, LWRBND(N+NCLIN+NCNLN),
     *                  UPRBND(N+NCLIN+NCNLN), WRK(LWRK), A(LDA*N) 

C
C  COMMON VARIABLES
C

C
C  PRINTOUT PARAMETERS
C
      INTEGER          IOUT  , IPRINT
      DOUBLE PRECISION MCHEPS, CNSTOL

      COMMON / DFOCM /  IOUT, IPRINT, MCHEPS, CNSTOL
      SAVE / DFOCM /
C
C  EXTERNAL SUBROUTINES
C

      EXTERNAL          FUNOBJ, FUNCON

      DOUBLE PRECISION DDOT

      EXTERNAL         DDOT

      INCLUDE 'dfo_model_inc.f'

C     PARAMETERS
      double precision INF
      parameter(INF = 1.0D+20)

C     LOCAL SCALARS
      logical checkder
      integer aprint,i,j,k,m,ncomp,ilb,iub,ila
      double precision cnorm,epsfeas,epsopt,nlpsupn,snorm

C     LOCAL ARRAYS
      logical coded(10),equatn(NCONMX),linear(NCONMX)

C     COMMON ARRAYS
      integer cmap(NCONMX)
      double precision cb(2 * NCONMX),cmul(2 * NCONMX)

C     COMMON BLOCKS
      common/algencancm/ cmul,cb,cmap

      useipopt = 0

      epsfeas = 1.0D-8
      epsopt = 1.0D-8
      aprint = 10
      ncomp = 50

C     m + 1 is the position where the information about the extra
C     constraints will be. At the end of the loop below, m contains the
C     exact number of constraints that Algencan will work.
      m = nlin + nnln + ncon

C     Adapts the constraints to Algencan's format.

      do i = 1,nlin + nnln + ncon
         j = n + i

         equatn(i) = .false.
         linear(i) = .false.

         if ( i .le. nlin .or.
     +        i .gt. nlin + nnln ) linear(i) = .true.


         if ( lwrbnd(j) .gt. -INF ) then
            cmap(i) = i
            cmul(i) = - 1.0D0
            cb(i) = - lwrbnd(j)

            if ( lwrbnd(j) .eq. uprbnd(j) ) equatn(i) = .true.
         end if
         if ( uprbnd(j) .lt. INF .and. (.not. equatn(i)) ) then

            if ( lwrbnd(j) .gt. -INF ) then
C              Creates the extra constraints
               m = m + 1
               
               cmap(m) = i
               cmul(m) = 1.0D0
               cb(m) = uprbnd(j)

               equatn(m) = .false.
               if ( i .le. nlin .or.
     +              i .gt. nlin + nnln ) linear(m) = .true.

            else
               cmap(i) = i
               cmul(i) = 1.0D0
               cb(i) = uprbnd(j)
            end if
         end if
      end do

C      write(*,*) 'm=',m,'lin=',nlin,'easy=',nnln,'hard=',ncon

C     Stores the linear constraints.

      do j = 1,n
         do i = 1,nlin
            amat(i,j) = a((j - 1) * LDA + i)
         end do
      end do

      ilb = 1
      iub = ilb + n
      ila = iub + n

C      write(*,*) 'espaco necessario=',n + n + m

      if ( lwrk .lt. 2 * n + m ) then
         inform = 2
         return
      end if

C     Lower and upper bounds.

      do i = 1,n
         wrk((ilb - 1) + i) = max(lwrbnd(i),x0(i) - delta)
         wrk((iub - 1) + i) = min(uprbnd(i),x0(i) + delta)
      end do

C     Initial Lagrange multipliers.

      do i = 1,m
         wrk((ila - 1) + i) = 0.0D0
      end do

      checkder = .false.

      coded( 1) = .false. ! evalf
      coded( 2) = .false. ! evalg
      coded( 3) = .false. ! evalh
      coded( 4) = .false. ! evalc
      coded( 5) = .false. ! evaljac
      coded( 6) = .false. ! evalhc
      coded( 7) = .true.  ! evalfc
      coded( 8) = .true.  ! evalgjac
      coded( 9) = .false. ! evalhl
      coded(10) = .false. ! evalhlp

      
c$$$      write(*,*)
c$$$      write(*,*)
c$$$      do k = 1,ncon
c$$$         do i = 1,n
c$$$            write(*,*) (QCON((k - 1) * n * n + (i - 1) * n + j),j = 1,n)
c$$$         end do
c$$$         write(*,*)
c$$$         write(*,*) (LCON((k - 1) * n + j), j = 1,n)
c$$$         write(*,*)
c$$$         write(*,*) lwrbnd(n + nlin + nnln + k),CCON(k),
c$$$     +              uprbnd(n + nlin + nnln + k)
c$$$         write(*,*)
c$$$      end do
c$$$      write(*,*)
c$$$      write(*,*) 'Obj. F'
c$$$      do i = 1,n
c$$$         write(*,*) (HMOD(i,j),j = 1,n)
c$$$      end do
c$$$      write(*,*)
c$$$      write(*,*) (GMOD(i),i = 1,n)
c$$$      write(*,*)
c$$$      write(*,*) 'x=',(X0(i),i = 1,n)
c$$$      write(*,*) 'l=',(wrk(ilb - 1 + i),i = 1,n)
c$$$      write(*,*) 'u=',(wrk(iub - 1 + i),i = 1,n)
c$$$      write(*,*)
c$$$      write(*,*)

      call algencan(epsfeas,epsopt,aprint,ncomp,n,x0,wrk(ilb),wrk(iub),
     +     m,wrk(ila),equatn,linear,coded,checkder,mval,cnorm,snorm,
     +     nlpsupn,inform)

      if ( cnorm .gt. epsfeas ) inform = 2

C      write(*,*) 'Flag of Algencan:',inform
C      write(*,*) '-->',x0(1),x0(2)

      end


c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$      subroutine evalf(n,x,f,flag)
c$$$
c$$$      implicit none
c$$$
c$$$      integer n,flag
c$$$      double precision f,x(n)
c$$$
c$$$      flag = - 1
c$$$
c$$$      end 
c$$$
c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$      subroutine evalg(n,x,g,flag)
c$$$
c$$$      implicit none
c$$$
c$$$C     SCALAR ARGUMENTS
c$$$      integer flag,n
c$$$
c$$$C     ARRAY ARGUMENTS
c$$$      double precision g(n),x(n)
c$$$
c$$$      flag = - 1
c$$$
c$$$      end
c$$$
c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$      subroutine evalh(n,x,hlin,hcol,hval,hnnz,flag)
c$$$
c$$$      implicit none
c$$$
c$$$C     SCALAR ARGUMENTS
c$$$      integer flag,n,hnnz
c$$$
c$$$C     ARRAY ARGUMENTS
c$$$      integer hcol(*),hlin(*)
c$$$      double precision hval(*),x(n)
c$$$
c$$$      flag = - 1
c$$$
c$$$      end
c$$$
c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$      subroutine evalc(n,x,ind,c,flag)
c$$$
c$$$      integer n,ind,flag
c$$$      double precision c,x(n)
c$$$
c$$$      flag = - 1
c$$$
c$$$      end
c$$$
c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$      subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)
c$$$
c$$$      implicit none
c$$$
c$$$C     SCALAR ARGUMENTS
c$$$      integer flag,ind,jcnnz,n
c$$$
c$$$C     ARRAY ARGUMENTS
c$$$      integer jcvar(n)
c$$$      double precision x(n),jcval(n)
c$$$
c$$$      flag = - 1
c$$$
c$$$      end
c$$$
c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$      subroutine evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag)
c$$$
c$$$      implicit none
c$$$
c$$$C     SCALAR ARGUMENTS
c$$$      integer flag,hcnnz,ind,n
c$$$
c$$$C     ARRAY ARGUMENTS
c$$$      integer hccol(*),hclin(*)
c$$$      double precision hcval(*),x(n)
c$$$
c$$$      flag = - 1
c$$$
c$$$      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalfc(n,x,f,m,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

C     COMMON VARIABLES AND PARAMETERS
      include '/home/fsobral/programs/Dfo-2.0.0/dfo_model_inc.f'

C     COMMON ARRAYS
      integer cmap(NCONMX)
      double precision cb(2 * NCONMX),cmul(2 * NCONMX)

C     COMMON BLOCKS
      common/algencancm/ cmul,cb,cmap

C     LOCAL SCALARS
      integer i,j,k,l

      flag = 0

      f = 0.0D0
      do i = 1,n
         f = f + GMOD(i) * x(i)
         do j = 1,n
            f = f + 5.0D-1 * x(i) * x(j) * HMOD(i,j)
         end do
      end do      

C     Linear constraints

      do i = 1,nlin
         c(i) = 0.0D0
         do j = 1,n
            c(i) = c(i) + amat(cmap(i),j) * x(j)
         end do
         c(i) = cmul(i) * c(i) - cb(i)
      end do

C     Easy constraints

      if ( nnln .gt. 0 ) then
         call easycon(n,x,nnln,c(nlin + 1))
      end if
      do i = nlin + 1,nlin + nnln
         c(i) = cmul(i) * c(i) - cb(i)
      end do


C     Quadratic model of the hard constraints

      do k = 1,ncon
         l = k + nlin + nnln
         c(l) = ccon(k)
         do i = 1,n
            c(l) = c(l) + lcon((k - 1) * n + i) * x(i)
            do j = 1,n
               c(l) = c(l) + 5.0D-1 * qcon((k - 1) * n * n +
     +               (i - 1) * n + j) * x(i) * x(j)
            end do
         end do
         c(l) = cmul(l) * c(l) - cb(l)
      end do

C     Extra constraints

      do i = nlin + nnln + ncon + 1,m
         c(i) = - c(cmap(i))
      end do      

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

C     COMMON VARIABLES AND PARAMETERS
      include '/home/fsobral/programs/Dfo-2.0.0/dfo_model_inc.f'

C     COMMON ARRAYS
      integer cmap(NCONMX)
      double precision cb(2 * NCONMX),cmul(2 * NCONMX)

C     COMMON BLOCKS
      common/algencancm/ cmul,cb,cmap

C     LOCAL SCALARS
      integer i,j,k,l

      flag = 0

      do i = 1,n
         g(i) = GMOD(i)
         do j = 1,n
            g(i) = g(i) + HMOD(i,j) * x(j)
         end do
      end do      

C     Linear constraints (column oriented)

      jcnnz = 0

      do j = 1,n         
         do i = 1,nlin
            jcnnz = jcnnz + 1
         
            jcfun(jcnnz) = i
            jcvar(jcnnz) = j
            jcval(jcnnz) = cmul(i) * amat(cmap(i),j)
         end do
      end do

C     Easy constraints (column oriented)

      if ( nnln .gt. 0 ) then
         call easyjac(n,x,nnln,nnln,jcval(jcnnz + 1))

         do j = 1,n
            do i = 1,nnln
               l = nlin + i
               jcnnz = jcnnz + 1

               jcfun(jcnnz) = l
               jcvar(jcnnz) = j
               jcval(jcnnz) = cmul(l) * jcval(jcnnz)
            end do
         end do
      end if

C     Quadratic model of the hard constraints
C     (column oriented)

      do j = 1,n
         do i = 1,ncon
            l = nlin + nnln + i
            jcnnz = jcnnz + 1
            jcfun(jcnnz) = l
            jcvar(jcnnz) = j
            jcval(jcnnz) = lcon((i - 1) * n + j)

            do k = 1,n
               jcval(jcnnz) = jcval(jcnnz) + qcon((i - 1) * n * n +
     +               (j - 1) * n + k) * x(k)
            end do
            jcval(jcnnz) = cmul(l) * jcval(jcnnz)
         end do
      end do

C

      do j = 1,n
         do i = nlin + nnln + ncon + 1,m
            jcnnz =  jcnnz + 1

            if ( cmap(i) .le. nlin ) then
               l = (cmap(i) - 1) * n + j
            else if ( cmap(i) .le. nlin + nnln ) then
               l = nlin * n + (cmap(i) - 1) * n + j
            else
               l = nlin * n + nnln * n + (cmap(i) - 1) * n + j
            end if

            jcfun(jcnnz) = i
            jcvar(jcnnz) = j
            jcval(jcnnz) = - jcval(l)
         end do
      end do      

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhl(n,x,m,lambda,scalef,scalec,hllin,hlcol,hlval,
     +hlnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,hlnnz,m,n
      double precision scalef

C     ARRAY ARGUMENTS
      integer hlcol(*),hllin(*)
      double precision hlval(*),lambda(m),scalec(m),x(n)

C     COMMON VARIABLES AND PARAMETERS
      include '/home/fsobral/programs/Dfo-2.0.0/dfo_model_inc.f'

C     COMMON ARRAYS
      integer cmap(NCONMX)
      double precision cb(2 * NCONMX),cmul(2 * NCONMX)

C     COMMON BLOCKS
      common/algencancm/ cmul,cb,cmap

C     LOCAL SCALARS
      integer easy,i,j,k,l

      flag = 0

      hlnnz = 0
      easy = nlin + nnln

C     Objective function and quadratic model of the hard
C     constraints

      do j = 1,n
         do i = 1,n
            hlnnz = hlnnz + 1
            hllin(hlnnz) = i
            hlcol(hlnnz) = j
            hlval(hlnnz) = scalef * hmod(i,j)

            do k = 1,ncon
               hlval(hlnnz) = hlval(hlnnz) + cmul(k) *
     +              scalec(easy + k) * lambda(easy + k) *
     +              qcon((k - 1) * n * n + (i - 1) * n + j)
            end do

            do k = nlin + nnln + ncon + 1,m
               if ( cmap(k) .ge. nlin + nnln + 1 ) then
                  hlval(hlnnz) = hlval(hlnnz) + cmul(k) * 
     +                 scalec(easy + cmap(k)) * lambda(easy + cmap(k)) *
     +                 qcon((cmap(k) - 1) * n * n + (i - 1) * n + j)
               end if
            end do

         end do
      end do      

C      Easy constraints

      do i = nlin + 1,m
         lambda(i) = cmul(i) * lambda(i) * scalec(i)
      end do

      do i = 1,nnln
         call easyhess(i,n,x,nnln,nnln,hlval,lambda(nlin + i))
      end do

      do i = nlin + nnln + ncon + 1,m
         if ( cmap(i) .ge. nlin + 1 .and.
     +        cmap(i) .le. nlin + nnln ) then
            call easyhess(cmap(i),n,x,nnln,nnln,hlval,lambda(i))
         end if
      end do

      do i = nlin + 1,m
         lambda(i) = lambda(i) / (scalec(i) * cmul(i))
      end do

      end 

c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$      subroutine evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)
c$$$
c$$$      implicit none
c$$$
c$$$C     SCALAR ARGUMENTS
c$$$      logical goth
c$$$      integer flag,m,n
c$$$      double precision sf
c$$$
c$$$C     ARRAY ARGUMENTS
c$$$      double precision hp(n),lambda(m),p(n),sc(m),x(n)
c$$$
c$$$      flag = - 1
c$$$
c$$$      end
