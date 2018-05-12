module ellipsedata

  ! PARAMETERS
  real(8),parameter :: PI   = 3.141592653589793D0
  real(8),parameter :: ZERO =             1.0D-16

!!$    real(8),dimension(2 * 5),parameter :: P = (/ 1.0D0, 1.0D0, 3.0D0, 1.0D0,&
!!$         2.0D0, 5.0D0, 2.0D0, - 1.0D0, 3.0D0, 3.0D0 /)
!!$  !real(8),dimension(2 * 5),parameter :: P = (/ 1.0D0, 1.0D0, 3.0D0, 1.0D0,&
!!$  !     2.0D0, 0.0D0, 2.0D0, 2.0D0, (2.0D0 + sqrt(2.0D0) / 2.0D0),&
!!$  !     (1.0D0 + sqrt(2.0D0) / 2.0D0) /)

  real(8),dimension(2 * 4),parameter :: P = (/ 1.0D0, 1.0D0, 3.0D0, 1.0D0,&
       2.0D0, 5.0D0, 2.0D0, - 1.0D0 /)

contains

  !------------------------------------------------------------!
  ! FUNCTION EVALELLIPSE                                       !
  !                                                            !
  ! This function is used to verify the relation between a     !
  ! point 'p' and the (probably) ellipse defined by 'x'.       !
  !                                                            !
  !------------------------------------------------------------!

  function evalellipse(x,p)

    implicit none

    ! ARRAY ARGUMENTS
    real(8),intent(in) :: p(:),x(:)

    ! RETURN VALUE
    real(8) :: evalellipse

!!$    evalellipse = x(1) * p(1) ** 2 + 2.0D0 * x(2) * p(1) * p(2) + &
!!$         x(3) * p(2) ** 2 + 2.0D0 * x(4) * p(1) + &
!!$         2.0D0 * x(5) * p(2) + x(6)

    evalellipse = - 2.0D0 * x(5)
    evalellipse = evalellipse + SQRT((x(1) - p(1)) ** 2 + (x(2) - p(2)) ** 2)
    evalellipse = evalellipse + SQRT((x(3) - p(1)) ** 2 + (x(4) - p(2)) ** 2)

  end function evalellipse

  !------------------------------------------------------------!
  ! SUBROUTINE EVALGELLIPSE                                    !
  !                                                            !
  ! This subroutine evaluates the gradient of the above        !
  ! function.                                                  !
  !                                                            !
  !------------------------------------------------------------!

  subroutine evalgellipse(x,p,jcvar,jcval,jcnnz)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: jcnnz

    ! ARRAY ARGUMENTS
    real(8) :: jcval(:),p(:),x(:)
    integer :: jcvar(:)

    intent(in ) :: p,x
    intent(out) :: jcnnz,jcval,jcvar

    ! LOCAL SCALARS
    real(8) :: tmp

    jcnnz = 5

    tmp = SQRT((x(1) - p(1)) ** 2 + (x(2) - p(2)) ** 2)

    jcvar(1) = 1
    jcval(1) = (x(1) - p(1)) / tmp
    
    jcvar(2) = 2
    jcval(2) = (x(2) - p(2)) / tmp

    tmp = SQRT((x(3) - p(1)) ** 2 + (x(4) - p(2)) ** 2)

    jcvar(3) = 3
    jcval(3) = (x(3) - p(1)) / tmp

    jcvar(4) = 4
    jcval(4) = (x(4) - p(2)) / tmp

    jcvar(5) = 5
    jcval(5) = - 2.0D0

!!$    jcnnz = 6
!!$
!!$    jcvar(1) = 1
!!$    jcval(1) = p(1) ** 2
!!$
!!$    jcvar(2) = 2
!!$    jcval(2) = 2.0D0 * p(1) * p(2)
!!$
!!$    jcvar(3) = 3
!!$    jcval(3) = p(2) ** 2
!!$
!!$    jcvar(4) = 4
!!$    jcval(4) = 2.0D0 * p(1)
!!$
!!$    jcvar(5) = 5
!!$    jcval(5) = 2.0D0 * p(2)
!!$
!!$    jcvar(6) = 6
!!$    jcval(6) = 1.0D0
  end subroutine evalgellipse

  !------------------------------------------------------------!
  ! SUBROUTINE CELLIPSE                                        !
  !                                                            !
  !------------------------------------------------------------!

  subroutine cellipse(x,ox,oy,a,b,angle)

    implicit none

    ! SCALAR ARGUMENTS
    real(8),intent(out) :: a,angle,b,ox,oy

    ! ARRAY ARGUMENTS
    real(8),intent(in ) :: x(:)

    ! LOCAL SCALARS
    real(8) :: dist

    dist = SQRT((x(1) - x(3)) ** 2 + (x(2) - x(4)) ** 2)

    ox = (x(1) + x(3)) / 2.0D0
    oy = (x(2) + x(4)) / 2.0D0

    a  = x(5)
    b  = SQRT(x(5) ** 2 - (dist / 2.0D0) ** 2)

    angle = ASIN(ABS(x(2) - x(4)) / dist)

!!$    ox = (x(3) * x(4) - x(2) * x(5)) / (x(2) ** 2 - x(1) * x(3))
!!$    oy = (x(1) * x(5) - x(2) * x(4)) / (x(2) ** 2 - x(1) * x(3))
!!$
!!$    a = sqrt(2.0D0 * (x(1) * x(5) ** 2 + x(3) * x(4) ** 2 + &
!!$         x(6) * x(2) ** 2 - 2.0D0 * x(2) * x(4) * x(5) - &
!!$         x(1) * x(3) * x(6)) / ((x(2) ** 2 - x(1) * x(3)) * &
!!$         (sqrt((x(1) - x(3)) ** 2 + 4.0D0 * x(2) ** 2) - &
!!$         (x(1) + x(3)))) )
!!$
!!$    b = sqrt(2.0D0 * (x(1) * x(5) ** 2 + x(3) * x(4) ** 2 + &
!!$         x(6) * x(2) ** 2 - 2.0D0 * x(2) * x(4) * x(5) - &
!!$         x(1) * x(3) * x(6)) / ((x(2) ** 2 - x(1) * x(3)) * &
!!$         (- sqrt((x(1) - x(3)) ** 2 + 4.0D0 * x(2) ** 2) - &
!!$         (x(1) + x(3)))) )
!!$
!!$    if ( abs(x(2)) .le. ZERO ) then
!!$       if ( x(1) .lt. x(3) ) then
!!$          angle = 0.0D0
!!$       else
!!$          angle = pi / 2.0D0
!!$       end if
!!$    else
!!$       if ( x(1) .lt. x(3) ) then
!!$          angle = ATAN(2.0D0 * x(2) / (x(1) - x(3))) / 2.0D0
!!$       else
!!$          angle = pi / 2.0D0 + ATAN(2.0D0 * x(2) / (x(1) - x(3))) / 2.0D0
!!$       end if
!!$    end if

    ! From rad to deg
    angle = angle * 180.0D0 / PI 

  end subroutine cellipse

end module ellipsedata
