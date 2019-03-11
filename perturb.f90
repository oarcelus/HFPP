module perturbations

  complex*16, public, allocatable, save           :: dn(:,:,:,:,:), dvp(:,:,:,:,:), dv(:,:,:,:,:)

  public :: calc_dvp
  public :: calc_dn
  public :: calc_dv

  contains

  subroutine calc_dvp(pt)

    use param,           only : norb, nsite
    use constants,       only : nspin, cmplx_0
    use responsefunc,    only : invur 

    implicit none

    integer                 :: iorb, iorb1, iorb2, iorb3, ispin, ispin1, isite, isite1, counter, i4, j4
    complex*16              :: pt4(norb,norb,nsite,nspin*2), dvp4(norb,norb,nsite,nspin*2)
    complex*16, intent(in)  :: pt(norb,norb,nsite,nspin,nspin)

    allocate(dvp(norb,norb,nsite,nspin,nspin))
  
    counter = 0
    do ispin = 1, nspin
       do ispin1 = 1, nspin
          counter = counter + 1
          pt4(:,:,:,counter) = pt(:,:,:,ispin1,ispin)
       end do
    end do

    dvp4 = cmplx_0
    do i4 = 1, 4
       do j4 = 1, 4
          do isite = 1, nsite
             do iorb = 1, norb
                do iorb1 = 1, norb
                   do isite1 = 1, nsite
                      do iorb2 = 1, norb
                         do iorb3 = 1, norb
                            dvp4(iorb,iorb1,isite,i4) = dvp4(iorb,iorb1,isite,i4) + &
                                                        invur(iorb,iorb1,isite,iorb2,iorb3,isite1,i4,j4)*pt4(iorb2,iorb3,isite1,j4)
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    dvp = cmplx_0
    counter = 0
    do ispin = 1, nspin
       do ispin1 = 1, nspin
          counter = counter + 1
          dvp(:,:,:,ispin1,ispin) = dvp4(:,:,:,counter)
       end do
    end do

  end subroutine calc_dvp

  subroutine calc_dn()

    use param,           only : norb, nsite
    use constants,       only : nspin, cmplx_0
    use responsefunc,    only : r

    implicit none

    integer                 :: iorb, iorb1, iorb2, iorb3, ispin, ispin1, isite, isite1

    allocate(dn(norb,norb,nsite,nspin,nspin))

    dn = cmplx_0
    do isite = 1, nsite
       do iorb = 1, norb
          do iorb1 = 1, norb
             do isite1 = 1, nsite
                do iorb2 = 1, norb
                   do iorb3 = 1, norb
                      dn(iorb,iorb1,isite,1,1) = dn(iorb,iorb1,isite,1,1) + &
                                                 r(iorb,iorb1,isite,iorb2,iorb3,isite1,1)*dvp(iorb2,iorb3,isite1,1,1)
                      dn(iorb,iorb1,isite,2,1) = dn(iorb,iorb1,isite,2,1) + &
                                                 r(iorb,iorb1,isite,iorb2,iorb3,isite1,2)*dvp(iorb2,iorb3,isite1,1,2)
                      dn(iorb,iorb1,isite,1,2) = dn(iorb,iorb1,isite,1,2) + &
                                                 r(iorb,iorb1,isite,iorb2,iorb3,isite1,3)*dvp(iorb2,iorb3,isite1,2,1)
                      dn(iorb,iorb1,isite,2,2) = dn(iorb,iorb1,isite,2,2) + &
                                                 r(iorb,iorb1,isite,iorb2,iorb3,isite1,4)*dvp(iorb2,iorb3,isite1,2,2)
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine calc_dn

  subroutine calc_dv()

    use param,           only : norb, nsite
    use constants,       only : nspin, cmplx_0
    use responsefunc,    only : ur

    implicit none

    integer                 :: iorb, iorb1, iorb2, iorb3, ispin, ispin1, isite, isite1, counter, i4, j4

    allocate(dv(norb,norb,nsite,nspin,nspin))

    do isite = 1, nsite
       do iorb = 1, norb
          do iorb1 = 1, norb
             do isite1 = 1, nsite
                do iorb2 = 1, norb
                   do iorb3 = 1, norb
                      dv(iorb,iorb1,isite,1,1) = dv(iorb,iorb1,isite,1,1) + &
                                                 ur(iorb,iorb1,isite,iorb2,iorb3,isite1,1,1)*dvp(iorb2,iorb3,isite1,1,1) + &
                                                 ur(iorb,iorb1,isite,iorb2,iorb3,isite1,1,4)*dvp(iorb2,iorb3,isite1,2,2)
                      dv(iorb,iorb1,isite,2,1) = dv(iorb,iorb1,isite,2,1) + &
                                                 ur(iorb,iorb1,isite,iorb2,iorb3,isite1,2,2)*dvp(iorb2,iorb3,isite1,2,1)
                      dv(iorb,iorb1,isite,1,2) = dv(iorb,iorb1,isite,1,2) + &
                                                 ur(iorb,iorb1,isite,iorb2,iorb3,isite1,3,3)*dvp(iorb2,iorb3,isite1,1,2)
                      dv(iorb,iorb1,isite,2,2) = dv(iorb,iorb1,isite,2,2) + &
                                                 ur(iorb,iorb1,isite,iorb2,iorb3,isite1,4,1)*dvp(iorb2,iorb3,isite1,1,1) + &
                                                 ur(iorb,iorb1,isite,iorb2,iorb3,isite1,4,4)*dvp(iorb2,iorb3,isite1,1,1)
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine calc_dv

end module perturbations
