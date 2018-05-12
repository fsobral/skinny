! --------------------------------------------------------------- !
! COLLECTIONSET.F90                                               !
!                                                                 !
! This is a collection of interesting derivative-free problems.   !
! The problems have the same form defined in Hock-Schittkowski    !
! set of problems                                                 !
!                                                                 !
! min f(x)                                                        !
! s.t. a_i^T x - b_i <= 0  (i = 1,...,df_pl)                      !
!      g_i(x)        <= 0  (i = df_pl + 1,...,df_pl + df_pn)      !
!      a_i^T x - b_i  = 0 (i = df_pl + df_pn + 1,...,             !
!                              df_pl + df_pn + df_ml)             !
!      h_i(x)         = 0 (i = df_pl + df_pn + df_ml + 1,...,df_m)!
!                                                                 !
!                                                                 !
! --------------------------------------------------------------- !


module dftests

  implicit none

  ! COMMON SCALARS
  integer :: df_ind,df_flag,df_jcnnz,df_ml,df_mn,df_n,df_pl,df_pn
  real(8) :: df_c,df_f
  ! Defines the problem to be used
  integer :: df_NP

  ! COMMON ARRAYS
  integer,allocatable :: df_jcvar(:)
  real(8),allocatable :: df_l(:),df_jcval(:),df_u(:),df_x(:)

contains

  function callp(wtd)

    implicit none

    ! SCALAR ARGUMENTS
    integer,intent(in) :: wtd

    ! RETURN VALUE
    integer :: callp

    ! LOCAL SCALARS
    integer :: flag

    select case(df_NP)
    case(1)
       call P1(wtd,flag)
    case(2)
       call P2(wtd,flag)
    case(3)
       call P3(wtd,flag)
    case(4)
       call P4(wtd,flag)
    case(5)
       call P5(wtd,flag)
    case(6)
       call P6(wtd,flag)
    case(7)
       call P7(wtd,flag)
    case(8)
       call P8(wtd,flag)
    case(9)
       call P9(wtd,flag)
    case default
       flag = - 1
    end select

    callp = flag

  end function callp

  ! -------- !
  ! PROBLEMS !
  ! -------- !


  ! --------------------------------------------------------------- !
  ! Problem 1 from                                                  !
  !                                                                 !
  ! V. Gautschi, "Optimally scaled and optimally conditioned        !
  ! Vandermonde and Vandermonde-like matrices", BIT Numerical       !
  ! Mathematics, DOI: 10.1007/s10543-010-0293-1                     !
  !                                                                 !
  ! The problem consists of finding the optimally conditioned       !
  ! Vandermonde matrix with dimension n.                            !
  !                                                                 !
  ! --------------------------------------------------------------- !

  subroutine P1(wtd,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,wtd

    intent(in ) :: wtd
    intent(out) :: flag


    ! LOCAL SCALARS
    integer :: error,i,j
    real(8) :: maxp,maxs,p,s

    flag = 0

    select case(wtd)

    case(1) ! Initialization

       df_n  = 10 ! >= 2
       df_pl =  9
       df_pn =  0
       df_ml =  0
       df_mn =  1
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),     STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       do i = 1,df_n
          df_x(i) = df_n - DBLE(i) + 1.0D0
       end do
       
       df_l(1:df_n) = 0.0D+00
       df_u(1:df_n) = 1.0D+20

    case(2) ! Objective function

       maxs = - 1.0D+99
       maxp = - 1.0D+99

       do i = 1,df_n
          s = 0.0D0
          do j = 1,df_n
             s = s + df_x(j) ** (i - 1)
          end do
          maxs = max(maxs,s)
       end do

       do i = 1,df_n
          p = 1.0D0
          do j = 1,df_n
             if ( i .ne. j ) then
                p = p * (1.0D0 + df_x(j)) / abs(df_x(i) - df_x(j))
             end if
          end do
          maxp = max(maxp,p)
       end do
       write(*,*) maxs,maxp
       df_f = maxs * maxp

    case(3) ! Constraints

       if ( df_ind .ge. 1 .and. df_ind .le. df_n - 1 ) then
          df_c = - df_x(df_ind) + df_x(df_ind + 1)
       else if ( df_ind .eq. df_n ) then
          df_c = - DBLE(df_n)
          do i = 1,df_n
             df_c = df_c + df_x(i) ** (df_n - 1)
          end do
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .ge. 1 .and. df_ind .le. df_n - 1) then
          df_jcnnz = 2

          df_jcvar(1) =  df_ind
          df_jcval(1) = - 1.0D0

          df_jcvar(2) = df_ind + 1
          df_jcval(2) =      1.0D0
       else if ( df_ind .eq. df_n ) then
          df_jcnnz = df_n
          do i = 1,df_n
             df_jcvar(i) = i
             df_jcval(i) = DBLE(df_n - 1) * df_x(i) ** (df_n - 2)
          end do
       else 
          flag = - 1
       end if

    case(5) ! Finalization

       deallocate(df_l,df_u,df_x,df_jcval,df_jcvar)

    case default ! Incorrect mode
       flag = - 1
    end select

  end subroutine P1

  ! --------------------------------------------------------------- !
  ! Problem 4.1 from                                                !
  !                                                                 !
  ! Audet, C. and Dennis Jr., J. E. "A progressive barrier for      !
  ! derivative-free nonlinear programming", SIAM Journal on         !
  ! Optimization, Vol. 20, No 1, pp 445--472, 2009.                 !
  !                                                                 !
  ! "Optimization on a thin curved domain".                         !
  !                                                                 !
  ! This problem has 2 starting points: one feasible and one        !
  ! infeasible. Please uncomment the lines to use the desirable     !
  ! point. The standard is the infeasible.                          !
  !                                                                 !
  ! --------------------------------------------------------------- !

  subroutine P2(wtd,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,wtd

    intent(in ) :: wtd
    intent(out) :: flag


    ! LOCAL SCALARS
    integer :: error
    real(8) :: thinfactor

    flag = 0
    thinfactor = 1.0D-1
    
    select case(wtd)

    case(1) ! Initialization

       df_n  = 2
       df_pl = 0
       df_pn = 2
       df_ml = 0
       df_mn = 0
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       !df_x(1:df_n) = (/ 0.0D0, 0.0D0 /) ! Feasible
       df_x(1:df_n) = (/ 0.0D0, 10.0D0 /) ! Infeasible
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = sqrt((df_x(1) - 20.0D0) ** 2 + (df_x(2) - 1.0D0) ** 2)

    case(3) ! Constraints

       if ( df_ind .eq. 1 ) then
          df_c = sin(df_x(1)) - df_x(2) - thinfactor
       else if ( df_ind .eq. 2 ) then
          df_c = df_x(2) - sin(df_x(1))
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .eq. 1 ) then
          df_jcnnz = 2

          df_jcvar(1) = 1
          df_jcval(1) = cos(df_x(1))

          df_jcvar(2) = 2
          df_jcval(2) = - 1.0D0
       else if ( df_ind .eq. 2 ) then
          df_jcnnz = 2

          df_jcvar(1) = 1
          df_jcval(1) = - cos(df_x(1))

          df_jcvar(2) = 2
          df_jcval(2) = 1.0D0
       else 
          flag = - 1
       end if

    case(5) ! Finalization

       deallocate(df_l,df_u,df_x,df_jcval,df_jcvar)

    case default ! Incorrect mode
       flag = - 1
    end select

  end subroutine P2

  ! --------------------------------------------------------------- !
  ! Problem 4.2 from                                                !
  !                                                                 !
  ! Audet, C. and Dennis Jr., J. E. "A progressive barrier for      !
  ! derivative-free nonlinear programming", SIAM Journal on         !
  ! Optimization, Vol. 20, No 1, pp 445--472, 2009.                 !
  !                                                                 !
  ! "Linear optimization on an hypersphere".                        !
  !                                                                 !
  ! This problem has a variable number of variables and two         !
  ! different starting points: one feasible and one infeasible.     !
  ! Uncomment the respective lines to use each point.               !
  !                                                                 !
  ! --------------------------------------------------------------- !

  subroutine P3(wtd,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,wtd

    intent(in ) :: wtd
    intent(out) :: flag


    ! LOCAL SCALARS
    integer :: error,i

    flag = 0

    select case(wtd)

    case(1) ! Initialization

       df_n  = 50 ! This number is variable
       df_pl = 0
       df_pn = 1
       df_ml = 0
       df_mn = 0
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       !df_x(1:df_n) = 0.0D0 ! Feasible
       df_x(1:df_n) = 3.0D0 ! Infeasible
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = 0.0D0
       do i = 1,df_n
          df_f = df_f + df_x(i)
       end do

    case(3) ! Constraints

       if ( df_ind .eq. 1 ) then
          df_c = - 3.0D0 * df_n
          do i = 1,df_n
             df_c = df_c + df_x(i) ** 2
          end do
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .eq. 1 ) then
          df_jcnnz = df_n

          do i = 1,df_n
             df_jcvar(i) = i
             df_jcval(i) = 2.0D0 * df_x(i)
          end do
       else 
          flag = - 1
       end if

    case(5) ! Finalization

       deallocate(df_l,df_u,df_x,df_jcval,df_jcvar)

    case default ! Incorrect mode
       flag = - 1
    end select

  end subroutine P3

  ! --------------------------------------------------------------- !
  ! Problem 4.3 from                                                !
  !                                                                 !
  ! Audet, C. and Dennis Jr., J. E. "A progressive barrier for      !
  ! derivative-free nonlinear programming", SIAM Journal on         !
  ! Optimization, Vol. 20, No 1, pp 445--472, 2009.                 !
  !                                                                 !
  ! "Linear optimization over a nonconvex set".                     !
  !                                                                 !
  ! This problem has a customizable number of variables and two     !
  ! starting points: one feasible and one infeasible. Uncomment     !
  ! each point for the desirable test.                              !
  !                                                                 !
  ! --------------------------------------------------------------- !

  subroutine P4(wtd,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,wtd

    intent(in ) :: wtd
    intent(out) :: flag


    ! LOCAL SCALARS
    integer :: error,i

    flag = 0

    select case(wtd)

    case(1) ! Initialization

       df_n  = 100 ! This number is variable
       df_pl = 0
       df_pn = 2
       df_ml = 0
       df_mn = 0
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       ! Feasible
       !df_x(1:df_n) = (/ DBLE(n), (0.0D0, i = 2,df_n) /)

       ! Infeasible
       df_x(1:df_n) = (/ DBLE(df_n), (0.0D0, i = 2,df_n - 1), - DBLE(df_n) /)
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = df_x(df_n)

    case(3) ! Constraints

       if ( df_ind .eq. 1 ) then
          df_c = - df_n ** 2
          do i = 1,df_n
             df_c = df_c + (df_x(i) - 1.0D0) ** 2
          end do
       else if ( df_ind .eq. 2 ) then
          df_c = df_n ** 2
          do i = 1,df_n
             df_c = df_c - (df_x(i) + 1.0D0) ** 2
          end do
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .eq. 1 ) then
          df_jcnnz = df_n

          do i = 1,df_n
             df_jcvar(i) = i
             df_jcval(i) = 2.0D0 * (df_x(i) - 1.0D0)
          end do
       else if ( df_ind .eq. 2 ) then
          df_jcnnz = df_n

          do i = 1,df_n
             df_jcvar(i) = i
             df_jcval(i) = - 2.0D0 * (df_x(i) + 1.0D0)
          end do
       else 
          flag = - 1
       end if

    case(5) ! Finalization

       deallocate(df_l,df_u,df_x,df_jcval,df_jcvar)

    case default ! Incorrect mode
       flag = - 1
    end select

  end subroutine P4

  ! --------------------------------------------------------------- !
  ! Problem Disconnect Domain 1                                     !
  !                                                                 !
  ! The domain of this problem is disconnected, it can be a good    !
  ! idea to turn the thin domain into a fat domain to make the      !
  ! disconnections disappear.                                       !
  !                                                                 !
  ! --------------------------------------------------------------- !

  subroutine P5(wtd,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,wtd

    intent(in ) :: wtd
    intent(out) :: flag


    ! LOCAL SCALARS
    integer :: error,i,j,ngrade
    real(8) :: tmp

    flag = 0

    ngrade = 5

    select case(wtd)

    case(1) ! Initialization

       df_n  = 1
       df_pl = 0
       df_pn = 0
       df_ml = 0
       df_mn = 1
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       df_x(1:df_n) = (/ 10.0D0 /)
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = df_x(1) ** 2

    case(3) ! Constraints

       if ( df_ind .eq. 1 ) then
          df_c = 1.0D0
          do i = 1,ngrade
             df_c = df_c * (df_x(1) - DBLE(i))
          end do
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .eq. 1 ) then
          df_jcnnz    = 1
          df_jcvar(1) = 1
          df_jcval(1) = 0.0D0
          do i = 1,ngrade
             tmp = 1.0D0
             do j = 1,ngrade
                if ( j .ne. i ) tmp = tmp * (df_x(1) - DBLE(j))
             end do
             df_jcval(1) = df_jcval(1) + tmp
          end do
       else 
          flag = - 1
       end if

    case(5) ! Finalization

       deallocate(df_l,df_u,df_x,df_jcval,df_jcvar)

    case default ! Incorrect mode
       flag = - 1
    end select

  end subroutine P5

  ! --------------------------------------------------------------- !
  ! Problem Disconnect Domain 2 - Binary knapsack problem           !
  !                                                                 !
  ! The domain of this problem is disconnected, it can be a good    !
  ! idea to turn the thin domain into a fat domain to make the      !
  ! disconnections disappear.                                       !
  !                                                                 !
  ! --------------------------------------------------------------- !

  subroutine P6(wtd,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,wtd

    intent(in ) :: wtd
    intent(out) :: flag


    ! LOCAL SCALARS
    integer :: error,i,nit,cap
    real(8) :: tmp

    ! LOCAL ARRAYS
    integer :: c(10),p(10)

    flag = 0

    nit = 5
    cap = 13

    c(1:nit) = (/ 8, 3,  8, 5, 5 /)
    p(1:nit) = (/ 1, 5, 10, 1, 7 /)

    select case(wtd)

    case(1) ! Initialization

       df_n  = nit
       df_pl = 0
       df_pn = 1
       df_ml = 0
       df_mn = nit
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       df_x(1:df_n) = 0.0D0
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = 0.0D0
       do i = 1,nit
          df_f = df_f - df_x(i) * p(i)
       end do

    case(3) ! Constraints

       if ( df_ind .eq. 1 ) then
          df_c = - DBLE(cap)
          do i = 1,nit
             df_c = df_c + df_x(i) * c(i)
          end do
       else if ( df_ind .ge. 2 .and. df_ind .le. nit + 1 ) then
          df_c = df_x(df_ind - 1) * (df_x(df_ind - 1) - 1.0D0)
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints


       if ( df_ind .eq. 1 ) then
          df_jcnnz = nit

          do i = 1,nit
             df_jcvar(i) = i
             df_jcval(i) = c(i)
          end do
       else if ( df_ind .ge. 2 .and. df_ind .le. nit + 1 ) then
          df_jcnnz = 1
          df_jcvar(1) = df_ind - 1
          df_jcval(1) = 2.0D0 * df_x(df_ind - 1) - 1.0D0
       else 
          flag = - 1
       end if

    case(5) ! Finalization

       deallocate(df_l,df_u,df_x,df_jcval,df_jcvar)

    case default ! Incorrect mode
       flag = - 1
    end select

  end subroutine P6

  ! --------------------------------------------------------------- !
  ! Packing triangles                                               !
  !                                                                 !
  !                                                                 !
  ! --------------------------------------------------------------- !

  subroutine P7(wtd,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,wtd

    intent(in ) :: wtd
    intent(out) :: flag

    ! INTERFACES
    interface
       subroutine mct(ct1,theta1,l1,nf1,ct2,theta2,l2,nf2,over)
         integer :: nf1,nf2
         real(8) :: ct1(2),ct2(2),l1,l2,over,theta1,theta2

         intent(in ) :: ct1,ct2,l1,l2,nf1,nf2,theta1,theta2
         intent(out) :: over
       end subroutine mct

       subroutine getvert(ct,theta,l,nf,vi,v)
         real(8) :: ct(2),v(2)
         real(8) :: l,theta
         integer :: nf,vi

         intent(in ) :: ct,l,nf,vi,theta
         intent(out) :: v
       end subroutine getvert

       subroutine getjvert(ct,theta,l,nf,vi,jv)
         real(8) :: ct(2),jv(2)
         real(8) :: l,theta
         integer :: nf,vi

         intent(in ) :: ct,l,nf,vi,theta
         intent(out) :: jv
       end subroutine getjvert

       subroutine cregion(x,ri,c)
         integer :: ri
         real(8) :: c
         real(8) :: x(2)
         
         intent(in ) :: ri,x
         intent(out) :: c
       end subroutine cregion

       subroutine jregion(x,ri,jcnnz,jcvar,jcval)
         integer :: ri,jcnnz
         integer :: jcvar(:)
         real(8) :: jcval(:),x(2)

         intent(in ) :: ri,x
         intent(out) :: jcnnz,jcval,jcvar
       end subroutine jregion
    end interface


    ! LOCAL SCALARS
    integer :: error,i,j,nite,it,it1,it2
    real(8) :: tmp

    ! LOCAL ARRAYS
    integer :: fite(10)
    real(8) :: lite(10),v(2)

    flag = 0

    nite = 9
    lite(1:nite) = 1.0D0
    fite(1:nite) = 3

    select case(wtd)

    case(1) ! Initialization

       df_n  = 3 * nite
       df_pl = 0
       df_pn = SUM(fite(1:nite))
       df_ml = 0
       df_mn = 0
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       call RANDOM_NUMBER(df_x(1:df_n))

       df_x(1:df_n) = - 1.0D0 + 2.0D0 * df_x(1:df_n)

!       df_x(1:df_n) = (/ 2.0D0, 2.0D0, 3.1415D0 / 4.0D0, &
!            1.0D0, 2.0D0, 3.0D0 * 3.1415D0 / 4.0D0, &
!            2.0D0, 1.5D0, 3.1415D0, -1.0D0, -1.0D0, 3.1415D0 /)
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = 0.0D0

       do i = 1,nite - 1
          it1 = 3 * (i - 1) + 1
          do j = i + 1,nite

             it2 = 3 * (j - 1) + 1
             
             call mct(df_x(it1),df_x(it1 + 2),lite(i),fite(i), &
                  df_x(it2),df_x(it2 + 2),lite(j),fite(j),tmp)
             
             df_f = df_f + tmp
          end do
       end do

    case(3) ! Constraints

       if ( df_ind .ge. 1 .and. df_ind .le. df_pn ) then

          i = df_ind
          do it = 1,nite
             i = i - fite(it)
             if ( i .le. 0 ) exit
          end do

          it1 = 3 * (it - 1) + 1
          it2 = i + fite(it)

          call getvert(df_x(it1),df_x(it1 + 2),lite(it),fite(it),it2,v)
          call cregion(v,1,df_c)

       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .ge. 1 .and. df_ind .le. df_pn ) then

          i = df_ind
          do it = 1,nite
             i = i - fite(it)
             if ( i .le. 0 ) exit
          end do

          it1 = 3 * (it - 1) + 1
          it2 = i + fite(it)

          call getvert(df_x(it1),df_x(it1 + 2),lite(it),fite(it),it2,v)
          
          call jregion(v,1,df_jcnnz,df_jcvar,df_jcval)
          
          df_jcnnz = df_jcnnz + 1
          do i = 1,2
             df_jcvar(i) = it1 - 1 + i
          end do
          df_jcvar(3) = it1 + 2
          
          call getjvert(df_x(it1),df_x(it1 + 2),lite(it),fite(it),it2,v)
          df_jcval(3) = df_jcval(1) * v(1) + df_jcval(2) * v(2)
          
       else
          flag = - 1
       end if

    case(5) ! Finalization

       deallocate(df_l,df_u,df_x,df_jcval,df_jcvar)

    case default ! Incorrect mode
       flag = - 1
    end select

  end subroutine P7

  ! --------------------------------------------------------------- !
  ! Packing polygons in arbitrary nonconvex regions                 !
  !                                                                 !
  !                                                                 !
  ! --------------------------------------------------------------- !

  subroutine P8(wtd,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,wtd

    intent(in ) :: wtd
    intent(out) :: flag

    ! INTERFACES
    interface
       subroutine mct(ct1,theta1,l1,nf1,ct2,theta2,l2,nf2,over)
         integer :: nf1,nf2
         real(8) :: ct1(2),ct2(2),l1,l2,over,theta1,theta2

         intent(in ) :: ct1,ct2,l1,l2,nf1,nf2,theta1,theta2
         intent(out) :: over
       end subroutine mct

       subroutine mctf(ct,theta,l,nf,over)
         integer :: nf
         real(8) :: ct(2),l,over,theta

         intent(in ) :: ct,l,nf,theta
         intent(out) :: over
       end subroutine mctf

       subroutine getvert(ct,theta,l,nf,vi,v)
         real(8) :: ct(2),v(2)
         real(8) :: l,theta
         integer :: nf,vi

         intent(in ) :: ct,l,nf,vi,theta
         intent(out) :: v
       end subroutine getvert

       subroutine getjvert(ct,theta,l,nf,vi,jv)
         real(8) :: ct(2),jv(2)
         real(8) :: l,theta
         integer :: nf,vi

         intent(in ) :: ct,l,nf,vi,theta
         intent(out) :: jv
       end subroutine getjvert

    end interface


    ! LOCAL SCALARS
    integer :: error,i,j,nite,it,it1,it2
    real(8) :: tmp

    ! LOCAL ARRAYS
    integer :: fite(10)
    real(8) :: lite(10),v(2)

    flag = 0

    nite = 3
    lite(1:nite) = 2.0D0
    fite(1:nite) = (/ 3, 4, 5 /)

    select case(wtd)

    case(1) ! Initialization

       df_n  = 3 * nite
       df_pl = 0
       df_pn = nite
       df_ml = 0
       df_mn = 0
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       call RANDOM_NUMBER(df_x(1:df_n))

       df_x(1:df_n) = - 1.0D0 + 2.0D0 * df_x(1:df_n)

!       df_x(1:df_n) = (/ 2.0D0, 2.0D0, 3.1415D0 / 4.0D0, &
!            1.0D0, 2.0D0, 3.0D0 * 3.1415D0 / 4.0D0, &
!            2.0D0, 1.5D0, 3.1415D0, -1.0D0, -1.0D0, 3.1415D0 /)
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = 0.0D0

       do i = 1,nite - 1
          it1 = 3 * (i - 1) + 1
          do j = i + 1,nite

             it2 = 3 * (j - 1) + 1
             
             call mct(df_x(it1),df_x(it1 + 2),lite(i),fite(i), &
                  df_x(it2),df_x(it2 + 2),lite(j),fite(j),tmp)
             
             df_f = df_f + tmp
          end do
       end do

    case(3) ! Constraints

       if ( df_ind .ge. 1 .and. df_ind .le. df_pn ) then

          it1 = 3 * (df_ind - 1) + 1

          call mctf(df_x(it1),df_x(it1 + 2),lite(df_ind),fite(df_ind),df_c)

       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       flag = - 1

    case(5) ! Finalization

       deallocate(df_l,df_u,df_x,df_jcval,df_jcvar)

    case default ! Incorrect mode
       flag = - 1
    end select

  end subroutine P8





  subroutine P9(wtd,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,wtd

    intent(in ) :: wtd
    intent(out) :: flag


    ! LOCAL SCALARS
    integer :: error,i
    real(8) :: pi


    pi = acos(- 1.0D0)

    flag = 0

    select case(wtd)

    case(1) ! Initialization

       df_n  = 1
       df_pl = 0
       df_pn = 0
       df_ml = 0
       df_mn = 0
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       df_x(1:df_n) = 0.0D0
       
       df_l(1:df_n) = - 30.0D0
       df_u(1:df_n) =   30.0D0

    case(2) ! Objective function

          df_f = df_x(1) ** 2 * cos(df_x(1) / 3.0D0)

    case(3) ! Constraints

       flag = - 1

    case(4) ! df_ind-th Jacobian of constraints

       flag = - 1

    case(5) ! Finalization

       deallocate(df_l,df_u,df_x,df_jcval,df_jcvar)

    case default ! Incorrect mode
       flag = - 1
    end select

  end subroutine P9

end module dftests
