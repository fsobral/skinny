module svmdata

  implicit none

  integer :: dim,nhyper,npbo,npin,npot
  real(8) :: scale,xmax,xmin,ymax,ymin
  real(8),dimension(:,:),allocatable :: pntbo,pntin,pntot

  ! PARAMETERS
  real(8),parameter :: SVMEPSIN   = 5.0D-01
  integer,parameter :: SVMMAXDOTS =    1000

contains

  subroutine load()

    ! Format of input file
    !
    ! DIM
    ! N_HYPER
    ! N_BORDER N_IN N_OUT N_OTHER
    !
    ! BORDER POINTS
    ! :
    !
    ! IN POINTS
    ! :
    !
    ! OUT POINTS
    ! :
    !
    ! OTHER POINTS
    ! :

    implicit none

    ! LOCAL SCALARS
    integer :: i,j

    ! LOCAL ARRAYS
    character(LEN=100)               :: fin
    real(8),dimension(:),allocatable :: maxp,minp
      

!    write(*,FMT=0010)
    read( *,FMT=0020)fin

!    fin = 'SVM3.IN'
!    fin = '../../personal/doutorado/engordamento/SVM2.IN'

    OPEN(60,FILE=fin)

    read(60,*) dim
    read(60,*) nhyper
    read(60,*) npbo,npin,npot

    nhyper = max(1,nhyper)

    allocate(pntbo(dim,npbo),pntin(dim,npin),pntot(dim,npot))
    allocate(maxp(dim),minp(dim))

    maxp = - 1.0D+20
    minp = + 1.0D+20

    do j = 1,npbo
       read(60,*) pntbo(:,j)
       maxp = max(maxp,pntbo(:,j))
       minp = min(minp,pntbo(:,j))
    end do
    do j = 1,npin
       read(60,*) pntin(:,j)
       maxp = max(maxp,pntin(:,j))
       minp = min(minp,pntin(:,j))
    end do
    do j = 1,npot
       read(60,*) pntot(:,j)
       maxp = max(maxp,pntot(:,j))
       minp = min(minp,pntot(:,j))
    end do

    CLOSE(60)
    
    ! X-scale (A4)
    scale = 10.0D0 / (maxp(1) - minp(1))
    ! Y-scale (A4)
    scale = min(scale,20.0D0 / (maxp(2) - minp(2)))

    xmax = maxp(1) + 0.5D0 / scale
    ymax = maxp(2) + 0.5D0 / scale
    xmin = minp(1) - 0.5D0 / scale
    ymin = minp(2) - 0.5D0 / scale

    deallocate(maxp,minp)

    ! NON-EXECUTABLE STATEMENTS
    
0010 FORMAT('Please, enter the name of the file containing the data.')
0020 FORMAT(A100)

  end subroutine load
  
  !------------------------------------------------------------!
  !                                                            !
  !------------------------------------------------------------!

  subroutine unload()

    deallocate(pntbo,pntin,pntot)

  end subroutine unload


end module svmdata
