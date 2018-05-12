! --------------------------------------------------------------- !
! GREEDINESS.F90                                                  !
!                                                                 !
! This is a collection of problems in which simple Augmented      !
! Lagrangian methods suffer the greediness effect.                !
!                                                                 !
! min f(x)                                                        !
! s.t. a_i^T x - b_i <= 0  (i = 1,...,df_pl)                      !
!      g_i(x)        <= 0  (i = df_pl + 1,...,df_pl + df_pn)      !
!      a_i^T x - b_i  = 0 (i = df_pl + df_pn + 1,...,             !
!                              df_pl + df_pn + df_ml)             !
!      h_i(x)         = 0 (i = df_pl + df_pn + df_ml + 1,...,df_m)!
!                                                                 !
! --------------------------------------------------------------- !


module dftests

  implicit none

  ! COMMON SCALARS
  integer :: df_ind,df_flag,df_ml,df_mn,df_n,df_NP,df_pl,df_pn
  real(8) :: df_c,df_f,df_jcnnz

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

    select case(df_NP)
    case(1)
       call P1(wtd,callp)
    case(2)
       call P2(wtd,callp)
    case(3)
       call P3(wtd,callp)
    case(4)
       call P4(wtd,callp)
    case(5)
       call P5(wtd,callp)
    case(6)
       call P6(wtd,callp)
    case(7)
       call P7(wtd,callp)
    case default
       callp = - 1
    end select

  end function callp

  ! -------- !
  ! PROBLEMS !
  ! -------- !


  ! --------------------------------------------------------------- !
  ! Problem 1 from                                                  !
  !                                                                 !
  ! Castelani, Martinez, Martinez and Svaiter, "Addressing the      !
  ! greediness phenomenon in Nonlinear Programming by means of      !
  ! Proximal Augmented Lagrangians", Computational Optimization and !
  ! Applications, Volume 46, Number 2, 229-245, 2010                !
  !                                                                 !
  ! --------------------------------------------------------------- !

  subroutine P1(wtd,flag)

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

       df_n  = 100
       df_pl = 100
       df_pn =   0
       df_ml =   0
       df_mn =   0
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),     STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       df_x(1:df_n) = - 7.0D0
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = 0.0D0
       do i = 1,df_n
          df_f = df_f + df_x(i) ** 3.0D0
       end do

    case(3) ! Constraints

       if ( df_ind .ge. 1 .and. df_ind .le. df_n ) then
          df_c = - df_x(df_ind)
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .ge. 1 .and. df_ind .le. df_n ) then
          df_jcnnz = 1
          df_jcvar(1) =  df_ind
          df_jcval(1) = - 1.0D0          
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
  ! Problem 56 from Hock and Schittkowski, 1981.                    !
  ! Initial point taken from Castelani, Martinez, Martinez and      !
  ! Svaiter, 2010. (Problem 2)                                      !
  ! --------------------------------------------------------------- !

  subroutine P2(wtd,flag)

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

       df_n  = 7
       df_pl = 0
       df_pn = 0
       df_ml = 0
       df_mn = 4
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(4),df_jcvar(4), STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       do i = 1,df_n
          df_x(i) = DBLE(i)
       end do
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = - df_x(1) * df_x(2) * df_x(3)

    case(3) ! Constraints

       if ( df_ind .eq. 1 ) then
          df_c = df_x(1) - 4.2D0 * sin(df_x(4)) ** 2
       else if ( df_ind .eq. 2 ) then
          df_c = df_x(2) - 4.2D0 * sin(df_x(5)) ** 2
       else if ( df_ind .eq. 3 ) then
          df_c = df_x(3) - 4.2D0 * sin(df_x(6)) ** 2
       else if ( df_ind .eq. 4 ) then
          df_c = df_x(1) + 2.0D0 * (df_x(2) + df_x(3)) - &
                 7.2D0 * sin(df_x(7)) ** 2
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .eq. 1 ) then
          df_jcnnz = 2

          df_jcvar(1) = 1
          df_jcval(1) = 1.0D0
          df_jcvar(2) = 4
          df_jcval(2) = - 8.4D0 * sin(df_x(4)) * cos(df_x(4))
       else if ( df_ind .eq. 2 ) then
          df_jcnnz = 2

          df_jcvar(1) = 2
          df_jcval(1) = 1.0D0
          df_jcvar(2) = 5
          df_jcval(2) = - 8.4D0 * sin(df_x(5)) * cos(df_x(5))
       else if ( df_ind .eq. 3 ) then
          df_jcnnz = 2

          df_jcvar(1) = 3
          df_jcval(1) = 1.0D0
          df_jcvar(2) = 6
          df_jcval(2) = - 8.4D0 * sin(df_x(6)) * cos(df_x(6))
       else if ( df_ind .eq. 4 ) then
          df_jcnnz = 4

          df_jcvar(1) = 1
          df_jcval(1) = 1.0D0
          df_jcvar(2) = 2
          df_jcval(2) = 2.0D0
          df_jcvar(3) = 3
          df_jcval(3) = 2.0D0
          df_jcvar(4) = 7
          df_jcval(4) = - 14.4D0 * sin(df_x(7)) * cos(df_x(7))
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
  ! Problem 3 from Castelani, Martinez, Martinez and Svaiter, 2010  !
  ! --------------------------------------------------------------- !

  subroutine P3(wtd,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,wtd

    intent(in ) :: wtd
    intent(out) :: flag


    ! LOCAL SCALARS
    integer :: error

    flag = 0

    select case(wtd)

    case(1) ! Initialization

       df_n  = 2
       df_pl = 0
       df_pn = 0
       df_ml = 0
       df_mn = 1
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),     STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       df_x(1:df_n) = 1.0D0
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = - df_x(1) * df_x(2) ** 3.0D0

    case(3) ! Constraints

       if ( df_ind .eq. 1 ) then
          df_c = df_x(1) * df_x(2) - 4.0D0 * sin(df_x(1)) ** 2
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .eq. 1 ) then
          df_jcnnz = 2
          df_jcvar(1) = 1
          df_jcval(1) = df_x(2) - 8.0D0 * sin(df_x(1)) * cos(df_x(1))
          df_jcvar(2) = 2
          df_jcval(2) = df_x(1)
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
  ! Problem 4 from Castelani, Martinez, Martinez and Svaiter, 2010  !
  ! --------------------------------------------------------------- !

  subroutine P4(wtd,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,wtd

    intent(in ) :: wtd
    intent(out) :: flag


    ! LOCAL SCALARS
    integer :: error

    flag = 0

    select case(wtd)

    case(1) ! Initialization

       df_n  = 2
       df_pl = 0
       df_pn = 0
       df_ml = 0
       df_mn = 1
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),     STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       df_x(1) =   1.0D0
       df_x(2) = - 1.5D0
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = - df_x(1) * exp(- df_x(1) * df_x(2))

    case(3) ! Constraints

       if ( df_ind .eq. 1 ) then
          df_c = - (df_x(1) + 1.0D0) ** 3.0D0 + &
               3.0D0 * (df_x(1) + 1.0D0) ** 2.0D0 - 1.5D0 + df_x(2)
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .eq. 1 ) then
          df_jcnnz = 2
          df_jcvar(1) = 1
          df_jcval(1) = - 3.0D0 * (df_x(1) + 1.0D0) ** 2 + &
               6.0D0 * (df_x(1) + 1.0D0)
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

  end subroutine P4


  ! --------------------------------------------------------------- !
  ! Problem 5 from Castelani, Martinez, Martinez and Svaiter, 2010  !
  ! --------------------------------------------------------------- !

  subroutine P5(wtd,flag)

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

       df_n  = 50
       df_pl =  0
       df_pn =  1
       df_ml =  0
       df_mn =  0
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),     STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       df_x(1:df_n) = 1.0D-01
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = 0.0D0
       do i = 1,df_n
          df_f = df_f + df_x(i) ** 8.0D0 + df_x(i)
       end do
       df_f = - df_f

    case(3) ! Constraints

       if ( df_ind .eq. 1 ) then
          df_c = - 1.0D0
          do i = 1,df_n
             df_c = df_c + df_x(i) ** 2.0D0
          end do
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .eq. 1 ) then
          df_jcnnz = df_n
          do i = 1,df_n
             df_jcvar(i) =               i
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

  end subroutine P5


  ! --------------------------------------------------------------- !
  ! Problem 6 from Castelani, Martinez, Martinez and Svaiter, 2010  !
  ! --------------------------------------------------------------- !

  subroutine P6(wtd,flag)

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

       df_n  = 100
       df_pl =   0
       df_pn =   1
       df_ml =   0
       df_mn =   0
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(df_n),df_jcvar(df_n),     STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       df_x(1:df_n) = 1.0D0 / DBLE(df_n)
       
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = 0.0D0
       do i = 1,df_n
          df_f = df_f + phi6(df_x(i))
       end do

    case(3) ! Constraints

       if ( df_ind .eq. 1 ) then
          df_c = - 1.0D0
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
             df_jcvar(i) =               i
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

  end subroutine P6

  ! Auxiliar function

  function phi6(t)
    
    ! RETURN VALUE
    real(8) :: phi6

    ! SCALAR ARGUMENTS
    real(8) :: t

    phi6 = cos(t)
    if ( phi6 .gt. 0.0D0 ) then
       phi6 = log(phi6)
    else
       phi6 = - 1.0D+30
    end if

  end function phi6


  ! --------------------------------------------------------------- !
  ! Problem 40 from Hock and Schittkowski, 1981.                    !
  ! --------------------------------------------------------------- !

  subroutine P7(wtd,flag)

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

       df_n  = 4
       df_pl = 0
       df_pn = 0
       df_ml = 0
       df_mn = 3
       
       allocate(df_l(df_n),df_u(df_n),df_x(df_n), &
                df_jcval(4),df_jcvar(4), STAT=error)
       if ( error .ne. 0 ) then
          flag = - 1
          return
       end if

       df_x(1:df_n) =   8.0D-01
       df_l(1:df_n) = - 1.0D+20
       df_u(1:df_n) =   1.0D+20

    case(2) ! Objective function

       df_f = - df_x(1) * df_x(2) * df_x(3) * df_x(4)

    case(3) ! Constraints

       if ( df_ind .eq. 1 ) then
          df_c = df_x(1) ** 3.0D0 + df_x(2) ** 2.0D0 - 1.0D0
       else if ( df_ind .eq. 2 ) then
          df_c = df_x(4) * df_x(1) ** 2.0D0 - df_x(3)
       else if ( df_ind .eq. 3 ) then
          df_c = df_x(4) ** 2.0D0 - df_x(2)
       else
          flag = - 1
       end if

    case(4) ! df_ind-th Jacobian of constraints

       if ( df_ind .eq. 1 ) then
          df_jcnnz = 2

          df_jcvar(1) = 1
          df_jcval(1) = 3.0D0 * df_x(1) ** 2.0D0
          df_jcvar(2) = 2
          df_jcval(2) = 2.0D0 * df_x(2)
       else if ( df_ind .eq. 2 ) then
          df_jcnnz = 3

          df_jcvar(1) = 1
          df_jcval(1) = 2.0D0 * df_x(1) * df_x(4)
          df_jcvar(2) = 3
          df_jcval(2) = - 1.0D0
          df_jcvar(3) = 4
          df_jcval(3) = df_x(1) ** 2.0D0
       else if ( df_ind .eq. 3 ) then
          df_jcnnz = 2

          df_jcvar(1) = 2
          df_jcval(1) = - 1.0D0
          df_jcvar(2) = 4
          df_jcval(2) = 2.0D0 * df_x(4)
       else
          flag = - 1
       end if

    case(5) ! Finalization

       deallocate(df_l,df_u,df_x,df_jcval,df_jcvar)

    case default ! Incorrect mode
       flag = - 1
    end select

  end subroutine P7

end module dftests
