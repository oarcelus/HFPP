module dminter

  real*8, save, public, allocatable     :: dmi(:,:)

  contains
  
  subroutine calc_dmi()

    use param,               only : nrpt, norb, nsite, atom_label, r_bond
    use constants,           only : nspin, cmplx_i, cmplx_0

    implicit none

    integer                 :: irpt, isite
    real*8                  :: dum, devec(3,nsite), isoj(nrpt), forces(3,nrpt), &
                               dmix(3,nrpt), dmiy(3,nrpt), dmiz(3,nrpt)
    logical                 :: io_find, io_findz, io_findx, io_findy

    allocate(dmi(3,nrpt))

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

    dmix = 0.0
    dmiy = 0.0
    dmiz = 0.0
    dmi = 0.0

    inquire(file='forces_z.dat', exist=io_findz)
    inquire(file='forces_x.dat', exist=io_findx)
    inquire(file='forces_y.dat', exist=io_findy)

    if (.not.io_findx.or..not.io_findy.or..not.io_findz) then 
       print*, 'forces_*.dat file is missing please run SCLR in all three directions x, y, z'   
       stop
    end if
    open(2, file='forces_z.dat', status='old')   
       do isite = 1, nsite
          read(2,*) devec(:,isite)
       end do
       read(2,*)
       do irpt = 1, nrpt
          read(2,*) dum, dum, dum, dum, dum, forces(:,irpt)
          dmiz(2,irpt) = forces(1,irpt) - 0.5*isoj(irpt)*(devec(1,atom_label(2,irpt))-devec(1,atom_label(1,irpt)))
          dmiz(1,irpt) = -1.0*forces(2,irpt) + 0.5*isoj(irpt)*(devec(2,atom_label(2,irpt))-devec(2,atom_label(1,irpt)))
       end do
    close(2)
    open(2, file='forces_x.dat', status='old')
       do isite = 1, nsite
          read(2,*) devec(:,isite)
       end do
       read(2,*) 
       do irpt = 1, nrpt
          read(2,*) dum, dum, dum, dum, dum, forces(:,irpt)
          dmix(3,irpt) = forces(2,irpt) - 0.5*isoj(irpt)*(devec(2,atom_label(2,irpt))-devec(2,atom_label(1,irpt)))
          dmix(2,irpt) = -1.0*forces(3,irpt) + 0.5*isoj(irpt)*(devec(3,atom_label(2,irpt))-devec(3,atom_label(1,irpt)))
       end do
    close(2)
    open(2, file='forces_y.dat', status='old')
       do isite = 1, nsite
          read(2,*) devec(:,isite)
       end do
       read(2,*)
       do irpt = 1, nrpt
          read(2,*) dum, dum, dum, dum, dum, forces(:,irpt)
          dmiy(1,irpt) = forces(3,irpt) - 0.5*isoj(irpt)*(devec(3,atom_label(2,irpt))-devec(3,atom_label(1,irpt)))
          dmiy(3,irpt) = -1.0*forces(1,irpt) + 0.5*isoj(irpt)*(devec(1,atom_label(2,irpt))-devec(1,atom_label(1,irpt)))
       end do
    close(2)

    ! Check dmi values between diiferent direction
    do irpt = 1, nrpt
       dmi(1,irpt) = dmiz(1,irpt)
       dmi(2,irpt) = dmiz(2,irpt)
       dmi(3,irpt) = dmix(3,irpt)
    end do

  end subroutine calc_dmi

end module dminter
