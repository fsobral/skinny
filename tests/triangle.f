c$$$      program teste
c$$$
c$$$      implicit none
c$$$
c$$$      integer ITEMAX,FMAX,NMAX
c$$$      parameter (ITEMAX = 1000)
c$$$      parameter (FMAX   =  100)
c$$$      parameter (NMAX   = 3 * ITEMAX)
c$$$
c$$$      integer i,n
c$$$      double precision v(2),x(NMAX)
c$$$      double precision over
c$$$
c$$$      integer nite
c$$$      integer fite(ITEMAX)
c$$$      double precision lite(ITEMAX)
c$$$      common /mctdata/lite,fite,nite
c$$$
c$$$      nite = 3
c$$$
c$$$      n = 3 * nite
c$$$c$$$
c$$$c$$$      x(1) = 2.0D0
c$$$c$$$      x(2) = 2.0D0
c$$$c$$$      x(3) = 3.1415D0 / 4.0D0
c$$$c$$$      lite(1) = 1.0D0
c$$$c$$$      fite(1) = 5
c$$$c$$$
c$$$c$$$      x(4) = 1.0D0
c$$$c$$$      x(5) = 1.0D0
c$$$c$$$      x(6) = 3.0D0 * 3.1415D0 / 4.0D0
c$$$c$$$      lite(2) = 2.0D0
c$$$c$$$      fite(2) = 4
c$$$
c$$$      read(*,*) (x(i), i = 1,n)
c$$$      do i = 1,nite
c$$$         lite(i) = 1.0D0
c$$$         fite(i) = i + 2
c$$$      end do
c$$$
c$$$      call printp(n,x)
c$$$
c$$$C      call overlap(n,x,over)
c$$$
c$$$C      do i = 1,3
c$$$C         call getvert(x(1),x(3),lite(1),i,v)
c$$$C         call cregion(v,1,over)
c$$$C         write(*,*) over
c$$$C      end do
c$$$
c$$$
c$$$C      call mct(x(1),x(3),lite(1),fite(1),x(4),x(6),lite(2),fite(2),over)
c$$$
c$$$C      call feasib(9,x,over)
c$$$
c$$$      end


C     This file provides subroutines for evaluating general
C     polygon overlapping

      subroutine overlap(n,x,f)

      implicit none

      integer ITEMAX
      parameter (ITEMAX = 1000)

      integer n
      double precision f
      double precision x(n)

      integer nite
      integer fite(ITEMAX)
      double precision lite(ITEMAX)
      common /mctdata/lite,fite,nite

      double precision over
      integer i,it1,it2,j

      f = 0.0D0

      do i = 1,nite - 1
         do j = i + 1,nite

            it1 = 3 * (i - 1) + 1
            it2 = 3 * (j - 1) + 1

            call mct(x(it1),x(it1 + 2),lite(i),fite(i),
     +           x(it2),x(it2 + 2),lite(j),fite(i),over)

            f = f + over
         end do
      end do

      end

      subroutine feasib(n,x,c)

      implicit none

      integer n
      double precision c
      double precision x(n)

      integer ITEMAX
      parameter (ITEMAX = 1000)

      integer nite
      double precision lite(ITEMAX)
      common /mctdata/lite,nite

      double precision over
      integer i,it1

      c = 0.0D0

      do i = 1,nite
         call mctf(x(3 * (i - 1) + 1),x(3 * i),lite(i),over)
         c = c + over
      end do

      end

C
C     Computes the overlapping between two polygons
C

      subroutine mct(ct1,theta1,l1,nf1,ct2,theta2,l2,nf2,over)

      implicit none

      integer FMAX
      parameter (FMAX   =  100)

      integer nf1,nf2
      double precision ct1(2),ct2(2),l1,l2,over,theta1,theta2

      double precision seed,ss,tmp
      double precision llvs(2),v(3,2 * FMAX),pnt(2)
      integer i,j,ncnt,npts
      logical inside

      double precision drandsc

C     Constructs triangles

      call gettri(ct1,theta1,l1,nf1,v(1,1),llvs,ss)
      call gettri(ct2,theta2,l2,nf2,v(1,nf1 + 1),llvs,ss)

C
C     Begins Monte Carlo method
C

      seed = 12345701.0D0
      npts = 100000
      ncnt = 0

      do i = 1,npts

         pnt(1) = llvs(1) + ss * drandsc(seed)
         pnt(2) = llvs(2) + ss * drandsc(seed)

         inside = .true.

         do j = 1,nf1
            tmp = v(1,j) * pnt(1) + v(2,j) * pnt(2) - v(3,j)
            if ( tmp .gt. 0.0D0 ) then
               inside = .false.
               goto 0010
            end if
         end do
         do j = nf1 + 1,nf1 + nf2
            tmp = v(1,j) * pnt(1) + v(2,j) * pnt(2) - v(3,j)
            if ( tmp .gt. 0.0D0 ) then
               inside = .false.
               goto 0010
            end if
         end do

 0010    if ( inside ) then
            ncnt = ncnt + 1
C            write(*,FMT=1200) pnt(1),pnt(2)
         end if

      end do
      
C     Estimates intersection

      over = (DBLE(ncnt) / DBLE(npts)) * (ss ** 2.0D0)

C     NON-EXECUTABLE STATEMENTS

 1200 FORMAT('fill fullcircle scaled 1pt shifted ((',F20.12,',',
     +     F20.12,')*u) withcolor 0.5red;')
      end


C
C     Constructs one triangle
C

      subroutine gettri(ct,theta,l,nf,v,llvs,ss)

      implicit none

      integer nf
      double precision ct(2),llvs(2),theta,v(3,nf)
      double precision l,ss

      integer i,j
      double precision tmpa,tmpb,radt,sqrt3,pi,r,tmp
      double precision ctr(2),vert(2,3),tmpv(2)


      r     = l * sqrt(3.0D0) / 3.0D0      
      sqrt3 = sqrt(3.0D0) / 2.0D0
      tmpa  = r * sin(theta)
      tmpb  = r * cos(theta)

C     Constructs vertices

C     Triangle scheme
C     
C              VERT(:,1)
C              /      \
C             /        \
C            /          \
C     VERT(:,2) -------- VERT(:,3)


c$$$      vert(1,1) = ct(1) - tmpb / 2.0D0 - tmpa * sqrt3
c$$$      vert(2,1) = ct(2) + tmpb * sqrt3 - tmpa / 2.0D0
c$$$
c$$$      vert(1,2) = ct(1) - tmpb / 2.0D0 + tmpa * sqrt3
c$$$      vert(2,2) = ct(2) - tmpb * sqrt3 - tmpa / 2.0D0
c$$$
c$$$      vert(1,3) = ct(1) + tmpb
c$$$      vert(2,3) = ct(2) + tmpa

     
      pi = 2.0D0 * asin(1.0D0)
      do j = 1,nf
         tmpv(1) = r * cos(2.0D0 * pi * j / DBLE(nf))
         tmpv(2) = r * sin(2.0D0 * pi * j / DBLE(nf))

C         vert(1,j) = ct(1) + cos(theta) * tmpv(1) - sin(theta) * tmpv(2)
C         vert(2,j) = ct(2) + sin(theta) * tmpv(1) + cos(theta) * tmpv(2)

         v(1,j) = ct(1) + cos(theta) * tmpv(1) - sin(theta) * tmpv(2)
         v(2,j) = ct(2) + sin(theta) * tmpv(1) + cos(theta) * tmpv(2)

      end do

C     Constructs normal vectors
C     V(1,i) x + V(2,i) y = V(3,i)

c$$$      v(1,1) = vert(2,1) - vert(2,nf)
c$$$      v(2,1) = vert(1,nf) - vert(1,1)
c$$$      v(3,1) = v(1,1) * vert(1,1) + v(2,1) * vert(2,1)
c$$$      do i = 2,nf
c$$$         v(1,i) = vert(2,i) - vert(2,i - 1)
c$$$         v(2,i) = vert(1,i - 1) - vert(1,i) 
c$$$         v(3,i) = v(1,i) * vert(1,i) + v(2,i) * vert(2,i)
c$$$      end do

      tmpv(1) = v(2,1) - v(2,nf)
      tmpv(2) = v(1,nf) - v(1,1)
      v(3,nf) = tmpv(1) * v(1,1) + tmpv(2) * v(2,1)
      do i = 2,nf
         tmp = v(1,i - 1)
         v(1,i - 1) = v(2,i) - v(2,i - 1)
         v(2,i - 1) = tmp - v(1,i) 
         v(3,i - 1) = v(1,i - 1) * v(1,i) + v(2,i - 1) * v(2,i)
      end do
      v(1,nf) = tmpv(1)
      v(2,nf) = tmpv(2)
      
      
c$$$C     Constructs normal vectors
c$$$C     V(1,i) x + V(2,i) y = V(3,i)
c$$$
c$$$      do i = 1,3
c$$$         v(1,i) = vert(2,i) - vert(2,MOD(i + 1,3) + 1)
c$$$         v(2,i) = vert(1,MOD(i + 1,3) + 1) - vert(1,i) 
c$$$         v(3,i) = v(1,i) * vert(1,i) + v(2,i) * vert(2,i)
c$$$      end do

C     Side of the containing square

      ss  = 2.0D0 * r

C     Lower left vertex of the containing square

      llvs(1) = ct(1) - r
      llvs(2) = ct(2) - r


c$$$C     Printing
c$$$
c$$$      write(*,FMT=1000)vert(1,1),vert(2,1),vert(1,2),vert(2,2),
c$$$     +     vert(1,3),vert(2,3)
C      write(*,FMT=1100)llvs(1),llvs(2),llvs(1) + ss,llvs(2),
C     +     llvs(1) + ss,llvs(2) + ss,llvs(1),llvs(2) + ss


C     NON-EXECUTABLE STATEMENTS

 1000 FORMAT('draw (',F20.12,',',F20.12,')*u -- (',F20.12,',',F20.12,
     +     ')*u -- (',F20.12,',',F20.12,')*u -- cycle;',/)

 1100 FORMAT('draw (',F20.12,',',F20.12,')*u -- (',F20.12,',',F20.12,
     +     ')*u -- (',F20.12,',',F20.12,')*u -- (',F20.12,',',F20.12,
     +     ')*u -- cycle;')


      end

C
C     Prints the polygon
C

      subroutine printp(n,x)

      implicit none

      integer ITEMAX
      parameter (ITEMAX = 1000)

      integer n
      double precision x(n)

      integer nite
      integer fite(ITEMAX)
      double precision lite(ITEMAX)
      common /mctdata/lite,fite,nite

      integer i,it,j
      double precision v(2)

      do i = 1,nite
         it = 3 * (i - 1) + 1

         write(*,FMT=1010,ADVANCE='NO')i

         do j = 1,fite(i)
            call getvert(x(it),x(it + 2),lite(i),fite(i),j,v)
            write(*,FMT = 1000,ADVANCE='NO') v(1),v(2)
         end do

         write(*,FMT=1020)
      end do

      write(*,FMT=1050)nite

C     NON-EXECUTABLE STATEMENTS

 1000 FORMAT('(',F20.12,',',F20.12,')*u -- ')
 1010 FORMAT('p[',I3,'] :=',1X)
 1020 FORMAT(1X,'cycle;')
 1050 FORMAT(/,/,'for i=1 step 1 until ',I3,':'/,
     +     'fill p[i] withcolor 0.6white;',/,
     +     'draw p[i] withpen pencircle scaled 1pt;',/,
     +     'endfor')

      end

C
C     Gets the i-th vertex of the polygon
C

      subroutine getvert(ct,theta,l,nf,vi,v)

      implicit none

      double precision ct(2),v(2)
      double precision l,theta
      integer nf,vi

      double precision r,pi
      double precision tmpv(2)

C     Constructs vertices

C     Triangle scheme
C     
C              VERT(:,1)
C              /      \
C             /        \
C            /          \
C     VERT(:,2)   ...    VERT(:,nf)


      r = l * sqrt(3.0D0) / 3.0D0      

      pi = 2.0D0 * asin(1.0D0)
      tmpv(1) = r * cos(2.0D0 * pi * vi / DBLE(nf))
      tmpv(2) = r * sin(2.0D0 * pi * vi / DBLE(nf))

      v(1) = ct(1) + cos(theta) * tmpv(1) - sin(theta) * tmpv(2)
      v(2) = ct(2) + sin(theta) * tmpv(1) + cos(theta) * tmpv(2)

c$$$      sqrt3 = sqrt(3.0D0) / 2.0D0
c$$$      tmpa  = r * sin(theta)
c$$$      tmpb  = r * cos(theta)
c$$$
c$$$      if ( vi .eq. 1 ) then
c$$$         v(1) = ct(1) - tmpb / 2.0D0 - tmpa * sqrt3
c$$$         v(2) = ct(2) + tmpb * sqrt3 - tmpa / 2.0D0
c$$$      else if ( vi .eq. 2 ) then
c$$$         v(1) = ct(1) - tmpb / 2.0D0 + tmpa * sqrt3
c$$$         v(2) = ct(2) - tmpb * sqrt3 - tmpa / 2.0D0
c$$$      else if ( vi .eq. 3 ) then
c$$$         v(1) = ct(1) + tmpb
c$$$         v(2) = ct(2) + tmpa
c$$$      end if

      end


C
C     Gets the derivatives associated with the
C     i-th vertex of the triangle. We know that it has the
C     following aspect
C      __                  __
C     |    |                 |
C     | I  | \partial \theta |
C     |__  |               __|
C
C     so, we only return the last column (full)

      subroutine getjvert(ct,theta,l,nf,vi,jv)

      implicit none

      double precision ct(2),jv(2)
      double precision l,theta
      integer nf,vi

      double precision tmpa,tmpb,sqrt3,pi,r
      double precision tmpv(2)

      sqrt3 = sqrt(3.0D0) / 2.0D0
      tmpa  = sin(theta) * l
      tmpb  = cos(theta) * l



      r = l * sqrt(3.0D0) / 3.0D0      

      pi = 2.0D0 * asin(1.0D0)
      tmpv(1) = r * cos(2.0D0 * pi * vi / DBLE(nf))
      tmpv(2) = r * sin(2.0D0 * pi * vi / DBLE(nf))

      jv(1) = - sin(theta) * tmpv(1) - cos(theta) * tmpv(2)
      jv(2) =   cos(theta) * tmpv(1) - sin(theta) * tmpv(2)


c$$$C     Jeito do Hector
c$$$
c$$$      r     = l * sqrt(3.0D0) / 3.0D0      
c$$$      sqrt3 = sqrt(3.0D0) / 2.0D0
c$$$      tmpa  = r * sin(theta)
c$$$      tmpb  = r * cos(theta)
c$$$
c$$$      if ( vi .eq. 1 ) then
c$$$         jv(1) = tmpa / 2.0D0 - tmpb * sqrt3
c$$$         jv(2) = - tmpa * sqrt3 - tmpb / 2.0D0
c$$$      else if ( vi .eq. 2 ) then
c$$$         jv(1) = tmpa / 2.0D0 + tmpb * sqrt3
c$$$         jv(2) = tmpa * sqrt3 - tmpb / 2.0D0
c$$$      else if ( vi .eq. 3 ) then
c$$$         jv(1) = - tmpa
c$$$         jv(2) =   tmpb
c$$$      end if

      end


C
C     Verifies if the polygon is inside the generic container
C

      subroutine mctf(ct,theta,l,nf,over)

      implicit none

      integer FMAX
      parameter (FMAX   =  100)

      integer nf
      double precision ct(2),l,over,theta

      double precision seed,ss,tmp
      double precision llvs(2),v(3,FMAX),pnt(2)
      integer i,j,ncnt,npts
      logical inside

      double precision drandsc

C     Constructs polygon

      call gettri(ct,theta,l,nf,v(1,1),llvs,ss)

C
C     Begins Monte Carlo method
C

      seed = 12345701.0D0
      npts = 100000
      ncnt = 0

      do i = 1,npts

         pnt(1) = llvs(1) + ss * drandsc(seed)
         pnt(2) = llvs(2) + ss * drandsc(seed)

         inside = .true.

         do j = 1,nf
            tmp = v(1,j) * pnt(1) + v(2,j) * pnt(2) - v(3,j)
            if ( tmp .gt. 0.0D0 ) then
               inside = .false.
               goto 0010
            end if
         end do

         if ( pnt(1) ** 3.0D0 - pnt(2) .le. 0.0D0 .and.
     +        pnt(2) + pnt(1) ** 2.0D0 - 10.0D0 .le. 0.0D0 ) then
            inside = .false.
            goto 0010
         end if

 0010    if ( inside ) then
            ncnt = ncnt + 1
C            write(*,FMT=1200) pnt(1),pnt(2)
         end if

      end do
      
C     Estimates intersection

      over = (DBLE(ncnt) / DBLE(npts)) * (ss ** 2.0D0)

C     NON-EXECUTABLE STATEMENTS

 1200 FORMAT('fill fullcircle scaled 1pt shifted ((',F20.12,',',
     +     F20.12,')*u) withcolor 0.4green;')
      end      
