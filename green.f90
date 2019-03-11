module greenfunction

  complex*16, allocatable, save, public  :: greenR(:,:,:), gr1(:,:,:), gr2(:,:,:), gr1b(:,:), gr2b(:,:)
  complex*16,              save, public  :: lf, fifj(3)


  public :: gfR
  public :: lich
  public :: lfdmi

  contains

  subroutine gfR(eps,bond) 

    use param,          only : hf_vec, hf_val, r_bond, nkpt, ntot, norb,    &
                               atom_label, kpt_cart, nsite, fdtemp, efermi, &
                               nrpt, nsub
    use constants,      only : nspin, cmplx_0, cmplx_i, delta

    implicit none

    integer, intent(in)      :: bond
    integer                  :: ikpt, itot, ispin, iorb, iorb1, irpt
    complex*16, intent(in)   :: eps
    real*8                   :: rdotk
    complex*16               :: fac, div, prod

    greenR = cmplx_0
    do ispin = 1, nspin
       do itot = 1, nsub
          do ikpt = 1, nkpt
             rdotk = dot_product(kpt_cart(:,ikpt),r_bond(:,bond))
             fac = exp(cmplx_i*rdotk)
             div = 1.0/(eps - hf_val(itot,ispin,ikpt) + cmplx_i*delta)
             do iorb = 1, norb
                do iorb1 = 1, norb
                   prod =       hf_vec((atom_label(1,bond)-1)*norb+iorb,itot,ispin,ikpt) *&
                          conjg(hf_vec((atom_label(2,bond)-1)*norb+iorb1,itot,ispin,ikpt))
                   greenR(iorb,iorb1,ispin) = greenR(iorb,iorb1,ispin) + prod*fac*div/nkpt 
                end do
             end do
          end do
       end do
    end do

  end subroutine gfR

  subroutine lich(bond)

    use param,          only : hf_pot, nrpt, norb, nsite, r_bond, atom_label, dpot
    use utilities,      only : trace
    use constants,      only : twopi

    implicit none 

    integer, intent(in) :: bond
    integer             :: isite, irpt
    complex*16          :: matr(norb,norb)

    matr(:,:) = matmul(matmul(dpot(:,:,atom_label(1,bond)),gr1(:,:,1)), &
                       matmul(dpot(:,:,atom_label(2,bond)),gr2(:,:,2)))

    lf = trace(matr,norb)

  end subroutine lich

  subroutine lfdmi(bond,dvper,exfield)

    use param,          only : norb, nsite, atom_label
    use utilities,      only : trace
    use constants,      only : twopi, nspin, cmplx_0

    implicit none

    integer, intent(in)    :: bond
    integer                :: isite, irpt,idir, iorb, ispin, iorb1, ispin1, coun, coun1
    complex*16, intent(in) :: dvper(norb,norb,nsite,nspin,nspin), exfield(norb,norb,nsite,nspin,nspin,3)
    complex*16             :: dvperx(norb*nspin,norb*nspin,nsite), term1(norb*nspin,norb*nspin,3), &
                              term2(norb*nspin,norb*nspin,3), exfieldx(norb*nspin,norb*nspin,nsite,3)
   
    ! Flatten dvp
    do isite = 1, nsite
       coun = 0
       do ispin = 1, nspin
          do iorb = 1, norb
             coun1 = 0
             coun = coun + 1
             do ispin1 = 1, nspin
                do iorb1 = 1, norb
                   coun1 = coun1 + 1
                   dvperx(coun,coun1,isite) = dvper(iorb,iorb1,isite,ispin,ispin1)
                end do
             end do
          end do
       end do
    end do

    ! Flatten exchange field
    do idir = 1, 3
    do isite = 1, nsite
       coun = 0
       do ispin = 1, nspin
          do iorb = 1, norb
             coun1 = 0
             coun = coun + 1
             do ispin1 = 1, nspin
                do iorb1 = 1, norb
                   coun1 = coun1 + 1
                   exfieldx(coun,coun1,isite,idir) = exfield(iorb,iorb1,isite,ispin,ispin1,idir)
                end do
             end do
          end do
       end do
    end do
    end do

    do idir = 1, 3
       term1(:,:,idir) = matmul(matmul(gr1b,dvperx(:,:,atom_label(2,bond))), &
                                matmul(gr2b,exfieldx(:,:,atom_label(1,bond),idir)))
       term2(:,:,idir) = matmul(matmul(gr1b,exfieldx(:,:,atom_label(2,bond),idir)), &
                                matmul(gr2b,dvperx(:,:,atom_label(1,bond))))
       fifj(idir) = trace((term1(:,:,idir)-term2(:,:,idir)),norb*nspin)
    end do

  end subroutine lfdmi


end module greenfunction
