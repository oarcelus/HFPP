module utilities

  implicit none 

  public :: utility_recip_lattice 
  public :: utility_cart_to_frac
  public :: utility_frac_to_cart
  public :: utility_spin_rotation
  public :: utility_gr_rotation
  public :: cutility_diag
  public :: rutility_diag
  public :: call_invert
  public :: heavyside
  public :: trace
  public :: gauss

  external  :: zheev
  external  :: dsyev 
  external  :: zgetrf
  external  :: zgetri

  contains

  subroutine utility_gr_rotation(M,U,N,MFLAT,mode)

    use constants,     only : nspin, cmplx_0

    implicit none

    complex*16, intent(in)       :: U(N,N,nspin,nspin), M(N,N,nspin)
    integer, intent(in)          :: N
    character(len=1), intent(in) :: mode ! f = forward rotation, b = backward rotation
    complex*16, intent(out)      :: MFLAT(N*nspin,N*nspin)
    complex*16                   :: UFLAT(N*nspin,N*nspin)
    integer                      :: iorb,iorb1,ispin,ispin1,isite,coun,coun1

    ! Flaten U
    UFLAT(1:N,1:N) = U(:,:,1,1)
    UFLAT(1:N,N+1:N*nspin) = U(:,:,1,2)
    UFLAT(N+1:N*nspin,1:N) = U(:,:,2,1)
    UFLAT(N+1:N*nspin,N+1:N*nspin) = U(:,:,2,2)
   
    ! Flatten M
    MFLAT = cmplx_0
    MFLAT(1:N,1:N) = M(:,:,1)
    MFLAT(N+1:N*nspin,N+1:N*nspin) = M(:,:,2)

    if (mode == 'f') then
       MFLAT(:,:) = matmul(UFLAT,matmul(MFLAT(:,:),conjg(transpose(UFLAT))))
    else if (mode == 'b') then
       MFLAT(:,:) = matmul(conjg(transpose(UFLAT)),matmul(MFLAT(:,:),UFLAT))
    else
       print*, 'Please insert rotation mode'
       stop
    end if

  end subroutine utility_gr_rotation

  subroutine utility_spin_rotation(M,U,N,NAT,mode,typ)

    use constants,     only : nspin

    implicit none

    complex*16, intent(inout)    :: M(N,N,NAT,nspin,nspin)
    complex*16, intent(in)       :: U(N,N,nspin,nspin)
    integer, intent(in)          :: N, NAT
    character(len=1), intent(in) :: mode ! f = U O U_dagg      , b = U_dagg O U
    character(len=1), intent(in) :: typ  ! o = operators       , n = occupation matrices
    complex*16                   :: MFLAT(N*nspin,N*nspin,NAT), UFLAT(N*nspin,N*nspin)
    integer                      :: iorb,iorb1,ispin,ispin1,isite,coun,coun1

    ! Flaten U
    UFLAT(1:N,1:N) = U(:,:,1,1)
    UFLAT(1:N,N+1:N*nspin) = U(:,:,1,2)
    UFLAT(N+1:N*nspin,1:N) = U(:,:,2,1)
    UFLAT(N+1:N*nspin,N+1:N*nspin) = U(:,:,2,2)

    ! Flatten M
    do isite = 1, NAT
       MFLAT(1:N,1:N,isite) = M(:,:,isite,1,1)
       MFLAT(1:N,N+1:N*nspin,isite) = M(:,:,isite,1,2)
       MFLAT(N+1:N*nspin,1:N,isite) = M(:,:,isite,2,1)
       MFLAT(N+1:N*nspin,N+1:N*nspin,isite) = M(:,:,isite,2,2)
    end do

    if (mode == 'f'.and.typ == 'o') then
       do isite = 1, NAT
          MFLAT(:,:,isite) = matmul(UFLAT,matmul(MFLAT(:,:,isite),conjg(transpose(UFLAT))))
       end do
    else if (mode == 'b'.and.typ == 'o') then
       do isite = 1, NAT
          MFLAT(:,:,isite) = matmul(conjg(transpose(UFLAT)),matmul(MFLAT(:,:,isite),UFLAT))
       end do
    else if (mode == 'b'.and.typ == 'n') then
       UFLAT = transpose(UFLAT)
       do isite = 1, NAT
          MFLAT(:,:,isite) = matmul(conjg(transpose(UFLAT)),matmul(MFLAT(:,:,isite),UFLAT))
       end do
    else if (mode == 'f'.and.typ == 'n') then
       UFLAT = transpose(UFLAT)
       do isite = 1, NAT
          MFLAT(:,:,isite) = matmul(UFLAT,matmul(MFLAT(:,:,isite),conjg(transpose(UFLAT))))
       end do
    else
       print*, 'Please insert rotation mode'
       stop
    end if

    ! Re order indeces M
    do isite = 1, NAT
       M(:,:,isite,1,1) = MFLAT(1:N,1:N,isite)
       M(:,:,isite,1,2) = MFLAT(1:N,N+1:N*nspin,isite)
       M(:,:,isite,2,1) = MFLAT(N+1:N*nspin,1:N,isite)
       M(:,:,isite,2,2) = MFLAT(N+1:N*nspin,N+1:N*nspin,isite)
    end do

  end subroutine utility_spin_rotation

  function heavyside(x) result (hside)
    
    real*8, intent(in)  :: x
    real*8              :: hside

    if (x < 0) then
       hside = 0.0
    else if (x >= 1) then
       hside = 1.0
    end if
 
  end function heavyside

  function gauss(x,x0,sigma) result (gaussfunc)

    use constants,    only : twopi

    real*8, intent(in)  :: x, x0, sigma
    real*8              :: gaussfunc

    if ((x-x0) > -6*sigma .and. (x-x0) < 6*sigma) then
       gaussfunc = exp(-((x-x0)**2)/(2*sigma**2))/sqrt(twopi*sigma**2)
    else
       gaussfunc = 0.0
    end if

  end function gauss

  function trace(A,N) result (tr)

    integer    , intent(in) :: N
    complex*16 , intent(in) :: A(N,N)
    complex*16              :: tr
    integer                 :: i  

    tr = 0.0
    do i = 1, N
       tr = tr + A(i,i)
    end do
  
  end function trace

  subroutine utility_recip_lattice (real_lat,recip_lat,volume)  !
    !==================================================================!
    !                                                                  !
    !!  Calculates the reciprical lattice vectors and the cell volume
    !                                                                  !
    !===================================================================

    implicit none
    real*8, intent(in)  :: real_lat (3, 3)
    real*8, intent(out) :: recip_lat (3, 3)
    real*8, intent(out) :: volume

    recip_lat(1,1)=real_lat(2,2)*real_lat(3,3)-real_lat(3,2)*real_lat(2,3)
    recip_lat(2,1)=real_lat(3,2)*real_lat(1,3)-real_lat(3,3)*real_lat(1,2)
    recip_lat(3,1)=real_lat(1,2)*real_lat(2,3)-real_lat(2,2)*real_lat(1,3)
    recip_lat(1,2)=real_lat(2,3)*real_lat(3,1)-real_lat(3,3)*real_lat(2,1)
    recip_lat(2,2)=real_lat(3,3)*real_lat(1,1)-real_lat(1,3)*real_lat(3,1)
    recip_lat(3,2)=real_lat(2,1)*real_lat(1,3)-real_lat(2,3)*real_lat(1,1)
    recip_lat(1,3)=real_lat(2,1)*real_lat(3,2)-real_lat(3,1)*real_lat(2,2)
    recip_lat(2,3)=real_lat(3,1)*real_lat(1,2)-real_lat(1,1)*real_lat(3,2)
    recip_lat(3,3)=real_lat(2,2)*real_lat(1,1)-real_lat(2,1)*real_lat(1,2)

    volume=real_lat(1,1)*recip_lat(1,1) + &
         real_lat(2,1)*recip_lat(2,1) + &
         real_lat(3,1)*recip_lat(3,1)

    recip_lat=2.*acos(-1.0)*recip_lat/volume
    volume=abs(volume)

  end subroutine utility_recip_lattice

  subroutine utility_cart_to_frac(cart,frac,recip_lat)
    !==================================================================!
    !                                                                  !
    !!  Convert from Cartesian to fractional coordinates
    !                                                                  !
    !===================================================================
    implicit none

    real*8, intent(in)  :: recip_lat(3,3)
    real*8, intent(out) :: frac(3)
    real*8, intent(in)  :: cart(3)

    integer :: i

    do i=1,3
       frac(i)=recip_lat(i,1)*cart(1) + recip_lat(i,2)*cart(2) + recip_lat(i,3)*cart(3)
    end do

    frac=frac/(2.0*acos(-1.0))


  end subroutine utility_cart_to_frac

  subroutine utility_frac_to_cart(frac,cart,real_lat)
    !==================================================================!
    !                                                                  !
    !!  Convert from fractional to Cartesian coordinates
    !                                                                  !
    !===================================================================
    implicit none

    real*8, intent(in)  :: real_lat(3,3)
    real*8, intent(in)  :: frac(3)
    real*8, intent(out) :: cart(3)

    integer :: i

    do i=1,3
       cart(i)=real_lat(i,1)*frac(1) + real_lat(i,2)*frac(2) + real_lat(i,3)*frac(3)
    end do

    return

  end subroutine utility_frac_to_cart   

  subroutine cutility_diag(MAT,EIG,N)

    implicit none

    integer, intent(in) :: N
    integer :: INF, LWORK
    integer, parameter :: LWMAX = 1000
    real*8, dimension(N), intent(out) :: EIG
    complex*16, dimension(N,N), intent(inout) :: MAT
    complex*16 :: W(LWMAX)
    real*8, dimension(3*N-2) :: RW

    call zheev('V','U',N,MAT,N,EIG,W,-1,RW,INF)
    LWORK = min(LWMAX, int(W(1)))
    call zheev('V','U',N,MAT,N,EIG,W,LWORK,RW,INF)

  end subroutine cutility_diag

  subroutine rutility_diag(MAT,EIG,N)

    implicit none

    integer, intent(in) :: N
    integer :: INF, LWORK
    integer, parameter :: LWMAX = 1000
    real*8, dimension(N), intent(out) :: EIG
    real*8, dimension(N,N), intent(inout) :: MAT
    real*8 :: W(LWMAX)

    call dsyev('V','U',N,MAT,N,EIG,W,-1,INF)
    LWORK = min(LWMAX, int(W(1)))
    call dsyev('V','U',N,MAT,N,EIG,W,LWORK,INF)

  end subroutine rutility_diag

  subroutine call_invert(MAT,N)

    integer, intent(in)                     :: N
    integer                                 :: LWORK, INFO, IPIV(N)
    integer, parameter                      :: LMAX = 1000
    complex*16, intent(inout)               :: MAT(N,N)
    complex*16                              :: WORK(LMAX)

    ! Factorize
    call zgetrf(N,N,MAT,N,IPIV,INFO)

    if (INFO == 0) then
       !Invert matrix
       call zgetri(N,MAT,N,IPIV,WORK,-1,INFO)
       LWORK = min(LMAX, int(WORK(1)))
       call zgetri(N,MAT,N,IPIV,WORK,LWORK,INFO)
    else
       print*, 'zgetrf failed'
    end if

  end subroutine call_invert

end module utilities
