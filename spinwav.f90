module spindispersion

  real*8, allocatable, public, save    :: wq(:,:)
  real*8,              public, save    :: cw

  public :: spin_wav

  contains

  subroutine get_cw()

    use constants,      only : kboltz
    use param,          only : sspin, atom_label, nrpt, nsite, r_bond

    implicit none

    logical                 :: io_find
    integer                 :: irpt
    complex*16              :: isoj(nrpt)
    real*8                  :: dum, suma

    inquire(file='lf.dat', exist=io_find)
    if (io_find) then
       open(2, file='lf.dat', status='old')
          do irpt = 1, nrpt
             read(2,*) dum, dum, dum, dum, dum, isoj(irpt)
          end do
       close(2)
    else
       print*, 'lf.dat file not found'
       stop
    end if

    suma = 0.0
    do irpt = 1, nrpt
       if (norm2(r_bond(:,irpt))<=1d-6) cycle
       if (atom_label(1,irpt) == 1) then
          suma = suma + real(isoj(irpt))/1000.0
       end if
    end do

    cw = (1.0+1.0/sspin)/(3.0*kboltz)*suma

  end subroutine get_cw

  subroutine spin_wav()

    use param,          only : kpt_pathr, nrpt, totpth, nsite, npth, npnts, atom_label, &
                               r_bond, norb, sspin, maglatt, atom_cart
    use constants,      only : cmplx_i, cmplx_0
    use utilities,      only : cutility_diag

    logical                 :: io_find
    integer                 :: ipth, irpt, isite, isite1, num, iline
    real*8                  :: rdotk, rkpt, r(3), dum
    complex*16              :: fac
    complex*16, allocatable :: jllk(:,:,:), jl(:,:), omg(:,:,:), temp(:,:), isoj(:)

    allocate(jllk(nsite,nsite,totpth))
    allocate(jl(nsite,nsite))
    allocate(isoj(nrpt))

    inquire(file='lf.dat', exist=io_find)
    if (io_find) then
       open(2, file='lf.dat', status='old')
          do irpt = 1, nrpt
             read(2,*) dum, dum, dum, dum, dum, isoj(irpt)
          end do
       close(2)
    else
       print*, 'lf.dat file not found'
       stop
    end if

    do irpt = 1, nrpt
       if (maglatt(atom_label(1,irpt)) + maglatt(atom_label(2,irpt)) == 0) then
          isoj(irpt) = -1.0*isoj(irpt)
       end if
    end do

    jllk = cmplx_0
    do ipth = 1, totpth
       do irpt = 1, nrpt
          rdotk = dot_product(kpt_pathr(:,ipth),r_bond(:,irpt))!(atom_cart(:,atom_label(2,irpt))-atom_cart(:,atom_label(1,irpt))))!r_bond(:,irpt))
          fac = exp(cmplx_i*rdotk)
          jllk(atom_label(1,irpt),atom_label(2,irpt),ipth) = jllk(atom_label(1,irpt),atom_label(2,irpt),ipth) + &
                                                             isoj(irpt)*fac
       end do
    end do

    jl = cmplx_0
    do isite = 1, nsite
       do isite1 = 1, nsite
          jl(isite,isite) = jl(isite,isite) + jllk(isite,isite1,1)
       end do
    end do

    allocate(omg(nsite,nsite,totpth))
    allocate(temp(nsite,nsite))
    omg = cmplx_0
    do ipth = 1, totpth
       omg(:,:,ipth) = (jl(:,:) - jllk(:,:,ipth))/sspin
    end do

    allocate(wq(nsite,totpth))
    do ipth = 1, totpth
       temp(:,:) = omg(:,:,ipth)
       call cutility_diag(temp(:,:),wq(:,ipth),nsite)
    end do

    open(unit = 23, file = 'spinwav.dat', status = 'replace')
    do isite = 1, nsite
       rkpt = 0.0
       write(23,*) rkpt, wq(isite,1)
       do ipth = 1, totpth - 1
          r(:) = kpt_pathr(:,ipth+1) - kpt_pathr(:,ipth)
          rkpt = rkpt + sqrt(r(1)**2 + r(2)**2 + r(3)**2)
          write(23,*) rkpt, wq(isite,ipth+1)
       end do
       write(23,*) ' '
    end do
    close(23)

  end subroutine spin_wav

end module spindispersion 
