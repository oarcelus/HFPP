module fractions

  real*8, allocatable, save, public   :: polres(:,:)  

! Using the continuous fraction aproximation to the fermi dirac function

  public :: fermi_dirac
  public :: cfermi_dirac
  public :: matsubara  
  public :: expfrac
  public :: get_residues
  public :: frac_fd
 
  contains

  function fermi_dirac(x) result (fd)

    use constants,    only : kboltz
    use param,        only : fdtemp

    implicit none

    real*8, intent(in)  :: x
    real*8              :: fd, y

    y = x/(kboltz*fdtemp)
    fd = 1.0/(exp(y) + 1.0)

  end function fermi_dirac

  function cfermi_dirac(x) result (fd)

    use constants,    only : kboltz
    use param,        only : fdtemp
 
    implicit none

    complex*16, intent(in)  :: x
    complex*16              :: fd, y

    y = x/(kboltz*fdtemp)
    fd = 1.0/(exp(y) + 1.0)

  end function cfermi_dirac

  function matsubara(x) result (mats)

    use constants,    only : kboltz, twopi
    use param,        only : fdtemp, npoles

    implicit none

    integer                 :: ipol
    complex*16, intent(in)  :: x
    complex*16              :: mats, y
    real*8                  :: pi
 
    y = x/(kboltz*fdtemp)
    pi = twopi/2.0
    mats = 0.5
    do ipol = 1, npoles
       mats = mats - 2.0*y/(y**2.0+((2.0*ipol-1.0)**2)*(pi**2.0))
    end do

  end function matsubara

  function expfrac(x) result (efr)

    use constants,    only : kboltz
    use param,        only : fdtemp, npoles

    implicit none
 
    complex*16, intent(in)  :: x
    complex*16              :: efr, y

    y = x/(kboltz*fdtemp)
    efr = 1.0/(1.0 + (1.0+(y/npoles))**npoles)

  end function expfrac

  subroutine get_residues()
    
    use constants,          only : kboltz
    use utilities,          only : rutility_diag 
    use param,              only : npoles, fdtemp

    implicit none

    real*8          :: A(npoles*2,npoles*2), B(npoles*2,npoles*2), &
                       D(npoles*2,npoles*2), eig(npoles*2), tmp(2)
    integer         :: ipol, ipol1

    A = 0.0
    B = 0.0
    D = 0.0

    do ipol = 1, npoles*2 
       A(ipol,ipol) = real(2*ipol - 1)
    end do

    do ipol = 1, npoles*2 - 1
       B(ipol+1,ipol) = 0.5
       B(ipol,ipol+1) = 0.5
    end do

    do ipol = 1, npoles*2
       A(ipol,ipol) = 1.0/sqrt(A(ipol,ipol))
    end do

    D = matmul(A,matmul(B,A))

    call rutility_diag(D,eig,npoles*2)

    allocate(polres(2,npoles*2))

    do ipol = 1, npoles*2
       polres(1,ipol) = 1/eig(ipol)
       polres(2,ipol) = -(1.0/4.0)*(D(1,ipol)/eig(ipol))**2
    end do

    do ipol = 1, 2*npoles - 1
       do ipol1 = ipol + 1, 2*npoles
          if (polres(1,ipol) > polres(1,ipol1)) then
             tmp(:) = polres(:,ipol)
             polres(:,ipol) = polres(:,ipol1)
             polres(:,ipol1) = tmp(:)
          end if
       end do
    end do
 
  end subroutine get_residues

  function frac_fd(x) result (frac)

    use constants,             only : kboltz, cmplx_i
    use param,                 only : npoles, fdtemp

    complex*16, intent(in) :: x
    complex*16             :: frac, y
    integer                :: ipol

    y = x/(kboltz*fdtemp)

    frac = 0.5 
    do ipol = 1, npoles
       frac = frac + polres(2,ipol+npoles)/(y-cmplx_i*polres(1,ipol+npoles)) + &
                     polres(2,ipol+npoles)/(y+cmplx_i*polres(1,ipol+npoles))
    end do

  end function frac_fd

end module fractions
