module engdata

  implicit none

  integer         :: numfeval, & ! Numb. of func. evaluations
                     numceval    ! Numb. of constr. evaluations
  integer,private :: m_original,n_original,return_flag
  logical,private :: jccheck = .false.
  real(8),private :: MYINFINITY = 1.0D+99

  real(8),private :: fmin = 1.0D+99

  real(8),allocatable         :: engXPrev(:)
  real(8),allocatable,target,private :: gamma_prob(:)
  real(8),pointer,private     :: l(:),u(:)

  ! true  = calling from restoration. 
  ! false = calling from BOBYQA.
  logical :: ENG_FORRES 
                        
  ! Fixed penalty parameter to the extreme barrier functions
  real(8) :: EXBARRHO = 1.0D+00

contains

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  subroutine engUpdateFmin(fnew)

    ! SCALAR ARGUMENT
    real(8),intent(in) :: fnew

    fmin = min(fmin,fnew)

  end subroutine engUpdateFmin

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  function engGetFmin()

    ! RETURN VALUE
    real(8) :: engGetFmin

    engGetFmin = fmin
  end function engGetFmin

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  function engGetM()
    
    ! RETURN VALUE
    integer :: engGetM

    engGetM = m_original
  end function engGetM

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  subroutine engSetM(m)

    ! SCALAR ARGUMENTS
    integer,intent(in) :: m

    m_original = m

  end subroutine engSetM

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  function engGetGamma()

    ! RETURN VALUE
    real(8),pointer :: engGetGamma(:)

    !write(*,*) 'engdata',gamma_prob

    engGetGamma => gamma_prob
  end function engGetGamma

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  subroutine engSetGamma(m,gamma)

    ! SCALAR ARGUMENTS
    integer,intent(in) :: m
    real(8),intent(in) :: gamma(m)

    if ( allocated(gamma_prob) .and. m .ne. m_original ) then
       deallocate(gamma_prob)
    end if

    if ( .not. allocated(gamma_prob) ) allocate(gamma_prob(m))

    gamma_prob = gamma

    !write(*,*) 'Entrou engSetGamma!',m

  end subroutine engSetGamma

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  function engGetN()

    ! RETURN VALUE
    integer :: engGetN

    engGetN = n_original
  end function engGetN

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  subroutine engSetN(n)

    ! SCALAR ARGUMENTS
    integer,intent(in) :: n

    n_original = n
  end subroutine engSetN

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  function engGetFlag()

    ! RETURN VALUE
    integer :: engGetFlag

    engGetFlag = return_flag
  end function engGetFlag

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  subroutine engSetFlag(flag)

    ! SCALAR ARGUMENTS
    integer,intent(in) :: flag

    return_flag = flag
  end subroutine engSetFlag

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  function engGetInfinity()

    ! RETURN VALUE
    real(8) :: engGetInfinity

    engGetInfinity = MYINFINITY
  end function engGetInfinity

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  subroutine engSetInfinity(infty)

    ! SCALAR ARGUMENTS
    real(8),intent(in) :: infty

    MYINFINITY = infty
  end subroutine engSetInfinity

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  function engJccheck()

    ! RETURN VALUE
    logical :: engJccheck

    engJccheck = jccheck
  end function engJccheck

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  subroutine engSetBounds(xl,xu)

    ! ARRAY ARGUMENTS
    real(8),target,intent(in) :: xl(:),xu(:)

    l => xl
    u => xu

  end subroutine engSetBounds

  !------------------------------------------------------------!
  !------------------------------------------------------------!

  function engVerifyBounds(n,x)

    ! SCALAR ARGUMENTS
    integer :: n

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    ! RETURN VALUE
    real(8) :: engVerifyBounds

    ! LOCAL SCALARS
    integer :: i

    engVerifyBounds = 0.0D0

    do i = 1,n
       engVerifyBounds = max(engVerifyBounds,l(i) - x(i))
       engVerifyBounds = max(engVerifyBounds,x(i) - u(i))
    end do

  end function engVerifyBounds

end module engdata
