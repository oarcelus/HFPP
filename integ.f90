module integral

  complex*16, allocatable, save, public   :: Cpath(:)

  public :: define_path

  contains

  subroutine define_path()

    use param,              only : efermi, npath, npoles, fdtemp
    use constants,          only : twopi, cmplx_i, kboltz
    use fractions,          only : polres

    implicit none

    integer    :: ipath
    real*8     :: rad, theeta, dang
 
    rad = polres(1,2*npoles)*kboltz*fdtemp + 1d+8

    theeta = twopi/2.0
 
    allocate(Cpath(0:npath))

    Cpath(0) = efermi + rad
    do ipath = 1, npath
       dang = ipath*theeta/npath
       Cpath(ipath) = efermi + rad*exp(cmplx_i*dang)
    end do

  end subroutine define_path

end module integral
