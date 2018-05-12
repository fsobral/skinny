c  sds.for
c  Strong Directional Search method for minimizing a continuous function
c  f
c  subject to (preferently fat) constraints. 
c  Coded by J. M. Martinez in 25/1/2011 under the Tosco Programming Project.
c
c  n is the number of variables
c  x is the initial point as entry and best solution found as output.
c  The initial point  x  must be feasible.
c  The user must code the subroutine calfun(n, x, f)
c  to compute the objective function 
c  and the subroutine constraints(n, x, feasible).
c  "feasible" must be a logical variable that needs to be set .true.
c  if the point satisfies the constraints and .false. otherwise.
c  eps is a user-given stopping criterion (in the domain space)
c  delta is a user-given initial step (in the domain space)
c  delta should be bigger than eps.
c  The initial approximation must satisfy xl leq x leq xu.
c
c  See
c  J. M. Martinez (2011) "The Tosco Programming Project", to appear in
c                         Forgettable Research. 
c

c  Sample main code

c$$$      implicit double precision (a-h, o-z)
c$$$      dimension x(100)
c$$$      external calfun,constraints
c$$$
c$$$      n = 3 
c$$$      do i = 1, n
c$$$      x(i) = 0.d0
c$$$      end do
c$$$
c$$$
c$$$      delta = 0.8
c$$$      eps = 1.d-4
c$$$      maxit = 10000
c$$$C      write(*, *)' maxsteps = '
c$$$C      write(*, *)' between 1 and infinity'
c$$$      read(*, *) delta,maxsteps
c$$$
c$$$      call sds(n, x, delta, eps, maxsteps, maxit, calfun, constraints,
c$$$     +     ier)
c$$$      write(*, *)' ier = ', ier
c$$$      stop
c$$$      end
c$$$
c$$$      subroutine calfun(n, x, f)
c$$$      implicit double precision (a-h, o-z)
c$$$      dimension x(n)
c$$$      f = 10.*(x(2)-x(1)**2)**2 + (x(1)-1.)**2 + dsin(x(1)*x(3))
c$$$     *  - 100. * x(2)**3
c$$$
c$$$      f = 0.d0
c$$$      do i = 1, n
c$$$      f = f + x(i)
c$$$      end do
c$$$
c$$$      return
c$$$      end
c$$$
c$$$      subroutine constraints(n, x, feasible)
c$$$      implicit double precision (a-h, o-z)
c$$$      dimension x(n)
c$$$      logical feasible
c$$$      feasible = .true.
c$$$      do i = 1, n
c$$$      if(x(i).lt.-1.d0.or.x(i).gt.1.d0) then
c$$$      feasible = .false.
c$$$      return
c$$$      endif
c$$$      end do
c$$$      z = 0.d0
c$$$      do i = 1, n
c$$$      z = z + x(i)
c$$$      end do
c$$$      if(z.gt.1.d0) then
c$$$      feasible = .false.
c$$$      return
c$$$      endif
c$$$      return
c$$$      end


      subroutine sds (n, x, delta, eps, maxsteps, maxit, objfun, constr,
     +                f, ier)
      implicit double precision (a-h, o-z)
      double precision drandsc
      dimension x(n)
      dimension dir(100, 100), xn(100)

      logical change, feasible

      external objfun,constr

      call constr(n, x, feasible)

      if(.not.feasible) then
      write(*, *)' The initial point is not feasible'
      ier = 2
      return
      else
C      write(*, *)' The initial point is feasible'
      endif


c  ndir is the number of directions used in the search, in addition
c  to the coordinate directions
c     We set ndir = n by default.
      ndir = n
      seed = 2251948.d0
      do i = 1, 10
      z = drandsc(seed)
      end do

      kon = 0
      nef = 0
      step = 0.9*delta

      call objfun (n, x, f)
      nef = nef + 1
1     continue

      write(*, *)' SDS iteration ', kon
      write(*, *)' X = ', (x(i),i=1,n)
      write(*, *)' f(X) = ', f
      write(*, *)' step = ', step 
      write(*, *)' Function evaluations:', nef
      write(*, *)

c  Stopping criterion
      if(2.*step.le.eps) then
      ier = 0
C      write(*, *)' Convergence with precision ', eps
C      write(*,*) f,nef
      return
      endif
      if(kon.ge.maxit) then
      ier = 1
C      write(*, *)' Number ', maxit,' of iterations exhausted'
C      write(*,*) f,nef
      return
      endif

c  Recompute search directions every 10 iterations

      if(mod(kon, 10).eq.0.and.ndir.gt.0) then
      jdir = 1
2     do i = 1, n
      dir(i, jdir) = 2.d0*drandsc(seed)-1.d0
      end do
      z = 0.d0
      do i = 1, n
      z = z + dir(i, jdir)**2
      end do
      if(z.gt.1.d0) go to 2
      z = dsqrt(z)
      do i = 1, n
      dir(i, jdir) = dir(i, jdir)/z
      end do
      if(jdir.lt.ndir) then
      jdir = jdir + 1
      go to 2
      endif
      endif



c  Cycle with size = step
      naux = delta/step + 1.
      nsteps = min0 ( maxsteps, naux)
      change = .false.     

c   Coordinate search
 
      do i = 1, n
      pasobest = 0
      do i1 = 1, nsteps
      do isi = -1, 1, 2        
      paso = dfloat(i1)*step * dfloat(isi)

      if(dabs(paso).le.delta) then 

      save = x(i)
      x(i) = x(i) + paso       

      call constr (n, x, feasible)

      if(feasible) then
         call objfun(n, x, fn)
         nef = nef + 1
         if(fn.lt.f) then
            change = .true. 
            f = fn
            pasobest = paso
         endif
      endif

      x(i) = save 
      endif
      
      end do
      end do

      x(i) = x(i) + pasobest

      end do 


c  Random directional search

3     continue
  
      do j = 1, ndir
         pasobest = 0.d0
         do i1 = 1, nsteps
            do isi = -1, 1, 2        
               paso = dfloat(i1)*step * dfloat(isi)
               
               if(dabs(paso).le.delta) then 
                  do i = 1, n
                     xn(i) = x(i) + paso * dir(i, j)
                  end do
                  call constr (n, xn, feasible)
                  
                  
                  if(feasible) then
                     call objfun(n, xn, fn)
                     nef = nef + 1
                     if(fn.lt.f) then
                        pasobest = paso 
                        change = .true. 
                        f = fn
                     endif
                  endif
               endif
               
            end do
            
         end do
         
         do i = 1, n
            x(i) = x(i) + pasobest * dir(i, j)
         end do      
      end do 
 
       
      if(.not.change) step = step/2.d0
      kon = kon + 1
      go to 1

      end
 
c$$$      
c$$$ 
c$$$      double precision function drand(ix)
c$$$
c$$$C     This is the random number generator of Schrage:
c$$$
c$$$C
c$$$
c$$$C     L. Schrage, A more portable Fortran random number generator, ACM
c$$$
c$$$C     Transactions on Mathematical Software 5 (1979), 132-138.
c$$$
c$$$      double precision ix
c$$$      double precision a,p,b15,b16,xhi,xalo,leftlo,fhi,k
c$$$      data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/
c$$$
c$$$      xhi= ix/b16
c$$$
c$$$      xhi= xhi - dmod(xhi,1.d0)
c$$$
c$$$      xalo= (ix-xhi*b16)*a
c$$$
c$$$      leftlo= xalo/b16
c$$$
c$$$      leftlo= leftlo - dmod(leftlo,1.d0)
c$$$
c$$$      fhi= xhi*a + leftlo
c$$$
c$$$      k= fhi/b15
c$$$
c$$$      k= k - dmod(k,1.d0)
c$$$
c$$$      ix= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
c$$$
c$$$      if (ix.lt.0) ix= ix + p
c$$$
c$$$      drand= ix*4.656612875d-10
c$$$
c$$$      return
c$$$
c$$$      end                   
