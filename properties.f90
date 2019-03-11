module prop

  complex*16, save, public                 ::  de, e3, e0
  complex*16, save, public, allocatable    ::  mut0(:,:), muts(:,:), mutl(:,:)

  public :: calc_mut0
  public :: calc_muts
  public :: calc_mutl
  public :: calc_de
  public :: calc_3de
  public :: calc_e0

  contains

  subroutine calc_mut0(dnst)

    use param,         only : nsite, norb, sigma_x, sigma_y, sigma_z
    use constants,     only : nspin, cmplx_0, cmplx_i
    use utilities,     only : trace

    implicit none

    integer                :: isite
    complex*16, intent(in) :: dnst(norb,norb,nsite,nspin,nspin)

    allocate(mut0(3,nsite))

    mut0 = cmplx_0
    do isite = 1, nsite
       mut0(1,isite) = mut0(1,isite) +         trace(dnst(:,:,isite,1,2)+dnst(:,:,isite,2,1),norb)
       mut0(2,isite) = mut0(2,isite) + cmplx_i*trace(dnst(:,:,isite,2,1)-dnst(:,:,isite,1,2),norb)
       mut0(3,isite) = mut0(3,isite) +         trace(dnst(:,:,isite,1,1)-dnst(:,:,isite,2,2),norb)
    end do

  end subroutine calc_mut0

  subroutine calc_muts(ddens)

    use param,          only : nsite, norb, sigma_x, sigma_y, sigma_z
    use constants,      only : nspin, cmplx_0, cmplx_i
    use utilities,      only : trace

    implicit none

    integer                :: isite    
    complex*16, intent(in) :: ddens(norb,norb,nsite,nspin,nspin)
    complex*16             :: dmuts(3,nsite)   
 
    allocate(muts(3,nsite))

    muts = cmplx_0
    dmuts = cmplx_0
    do isite = 1, nsite
       dmuts(1,isite) = dmuts(1,isite) +         trace(ddens(:,:,isite,1,2)+ddens(:,:,isite,2,1),norb)
       dmuts(2,isite) = dmuts(2,isite) + cmplx_i*trace(ddens(:,:,isite,2,1)-ddens(:,:,isite,1,2),norb)
       dmuts(3,isite) = dmuts(3,isite) +         trace(ddens(:,:,isite,1,1)-ddens(:,:,isite,2,2),norb)
    end do

    muts(:,:) = mut0(:,:) + dmuts(:,:)

  end subroutine calc_muts 

  subroutine calc_mutl(ddens)

    use param,          only : nsite, norb, lx, ly, lz
    use constants,      only : nspin, cmplx_0
    use utilities,      only : trace

    implicit none

    integer                :: isite, iorb, iorb1, ispin, ispin1
    complex*16, intent(in) :: ddens(norb,norb,nsite,nspin,nspin)

    allocate(mutl(3,nsite))

    mutl = cmplx_0
    do isite = 1, nsite
       do iorb = 1, norb
          do iorb1 = 1, norb
                mutl(1,isite) = mutl(1,isite) + lx(iorb,iorb1)*(ddens(iorb,iorb1,isite,1,1) + ddens(iorb,iorb1,isite,2,2))
                mutl(2,isite) = mutl(2,isite) + ly(iorb,iorb1)*(ddens(iorb,iorb1,isite,1,1) + ddens(iorb,iorb1,isite,2,2))
                mutl(3,isite) = mutl(3,isite) + lz(iorb,iorb1)*(ddens(iorb,iorb1,isite,1,1) + ddens(iorb,iorb1,isite,2,2))
          end do
       end do
    end do

  end subroutine calc_mutl

  subroutine calc_e0(dens)

    use param,          only : nsite, norb, hf_val, hf_pot, nkpt, efermi
    use constants,      only : nspin, cmplx_0
    use utilities,      only : trace

    implicit none

    integer                :: itot, ispin, ispin1, ikpt, iorb, iorb1, isite
    complex*16, intent(in) :: dens(norb,norb,nsite,nspin,nspin)

    e0 = 0.0
    do itot = 1, norb*nsite
       do ispin = 1, nspin
          do ikpt = 1, nkpt
             if (hf_val(itot,ispin,ikpt) > efermi) then
                cycle
             end if

             e0 = e0 + hf_val(itot,ispin,ikpt)/real(nkpt)
          end do
       end do
    end do

    do isite = 1, nsite
       do iorb = 1, norb
          do iorb1 = 1, norb
             e0 = e0 - 0.5*hf_pot(iorb,iorb1,isite,1)*dens(iorb,iorb1,isite,1,1)
             e0 = e0 - 0.5*hf_pot(iorb,iorb1,isite,2)*dens(iorb,iorb1,isite,2,2)
          end do
       end do
    end do

  end subroutine calc_e0

  subroutine calc_de(dvpot,ddens)

    use param,          only : nsite, norb
    use constants,      only : nspin, cmplx_0
    use utilities,      only : trace

    implicit none

    integer                :: isite, ispin, ispin1, iorb, iorb1, icomb, icomb1
  
    complex*16, intent(in) :: dvpot(norb,norb,nsite,nspin,nspin),&
                              ddens(norb,norb,nsite,nspin,nspin)
    de = 0.0
    do ispin = 1, nspin
       do ispin1 = 1, nspin
          do isite = 1, nsite
             do iorb = 1, norb
                do iorb1 = 1, norb
                   de = de + 0.5*dvpot(iorb,iorb1,isite,ispin,ispin1)*ddens(iorb,iorb1,isite,ispin,ispin1)
                end do
             end do
          end do
       end do
    end do


  end subroutine calc_de

  subroutine calc_3de(dvpot,pot)
 
    use param,                only : nrpt, ntot, norb, nsite, r_bond, nkpt, kpt_cart, &
                                     fdtemp, efermi, nelec
    use constants,            only : nspin, hamname, cmplx_i, kboltz, cmplx_0
    use utilities,            only : cutility_diag
    use fractions,            only : fermi_dirac

    implicit none

    integer                :: i, ikpt, irpt, itot, iat, iat1, dum, ispin, ispin1, isite, isite1, iorb, iorb1
    logical                :: io_find
    real*8                 :: eigval3(ntot,nkpt), emax, emin, numax, numin, num, e
    complex*16             :: ham(norb,norb,nspin,nspin), ham_k(norb,norb,nsite,nsite,nkpt,nspin,nspin), &
                              vnew(norb,norb,nsite,nspin,nspin), ham_full(ntot,ntot,nkpt), temp(ntot,ntot), &
                              eigvec3(ntot,ntot,nkpt), newdens(norb,norb,nsite,nspin,nspin)                                
    complex*16             :: rdotk, fac
    complex*16, intent(in) :: dvpot(norb,norb,nsite,nspin,nspin), pot(norb,norb,nsite,nspin,nspin)

    inquire(file=hamname, exist=io_find)
    if (io_find) then
       ham_k = cmplx_0
       do ikpt = 1, nkpt
          open(2, file=hamname, status='old')
          ham = cmplx_0
          do irpt = 1, nrpt
             read(2,*) iat, iat1, dum, dum
             read(2,*) 
             do iorb = 1, norb
                read(2,*) ham(iorb,:,1,1)
             end do

             !Expand to spindown subspace
             ham(:,:,2,2) = ham(:,:,1,1)

             rdotk = dot_product(kpt_cart(:,ikpt),r_bond(:,irpt))
             fac = exp(cmplx_i*rdotk)
             ham_k(:,:,iat,iat1,ikpt,1,1) = ham_k(:,:,iat,iat1,ikpt,1,1) + fac*ham(:,:,1,1)
             ham_k(:,:,iat,iat1,ikpt,2,2) = ham_k(:,:,iat,iat1,ikpt,2,2) + fac*ham(:,:,2,2)
             ham_k(:,:,iat,iat1,ikpt,1,2) = ham_k(:,:,iat,iat1,ikpt,1,2) + fac*ham(:,:,1,2)
             ham_k(:,:,iat,iat1,ikpt,2,1) = ham_k(:,:,iat,iat1,ikpt,2,1) + fac*ham(:,:,2,1)
          end do
          close(2)
       end do
    else
       print*, 'No sclr_ham.in file found'
       stop
    end if

    ! New hf potential as V_hf + dvp

    vnew = pot + dvpot

    !Sum kinetic part to the hamiltonian
    ham_full = cmplx_0
    do ikpt = 1, nkpt
       do ispin = 1, nspin
          do ispin1 = 1, nspin
             do isite = 1, nsite
                do isite1 = 1, nsite
                   ham_full((isite-1)*norb+1+(ispin-1)*norb*nsite:isite*norb+(ispin-1)*norb*nsite, &
                            (isite1-1)*norb+1+(ispin1-1)*norb*nsite:isite1*norb+(ispin1-1)*norb*nsite,ikpt) = &
                            ham_k(:,:,isite,isite1,ikpt,ispin,ispin1)
                end do
             end do
          end do
       end do
    end do

    !Sum hf_potential
    do ikpt = 1, nkpt
       do ispin = 1, nspin
          do ispin1 = 1, nspin
             do isite = 1, nsite
                ham_full((isite-1)*norb+1+(ispin-1)*norb*nsite:isite*norb+(ispin-1)*norb*nsite, &
                         (isite-1)*norb+1+(ispin1-1)*norb*nsite:isite*norb+(ispin1-1)*norb*nsite,ikpt) =   &
                ham_full((isite-1)*norb+1+(ispin-1)*norb*nsite:isite*norb+(ispin-1)*norb*nsite, &
                         (isite-1)*norb+1+(ispin1-1)*norb*nsite:isite*norb+(ispin1-1)*norb*nsite,ikpt) +   &
                vnew(:,:,isite,ispin,ispin1)
             end do
          end do
       end do
    end do

    ! Diagonalize
    do ikpt = 1, nkpt
       temp(:,:) = ham_full(:,:,ikpt)
       call cutility_diag(temp(:,:),eigval3(:,ikpt),ntot)
       eigvec3(:,:,ikpt) = temp(:,:)
    end do

    ! Calculate new efermi
    !emax = maxval(eigval3(ntot,:))
    !emin = minval(eigval3(1,:))
    !numax = 0.0
    !numin = 0.0
    !do itot = 1, ntot
    !   do ikpt = 1, nkpt
    !      numax = numax + fermi_dirac((eigval3(itot,ikpt)-emax))/real(nkpt)
    !      numin = numin + fermi_dirac((eigval3(itot,ikpt)-emin))/real(nkpt)
    !   end do
    !end do
!
!    if ((numax > nelec) .and. (numin < nelec)) then
!       do i = 1, 30
!          num = 0
!          e = (emax-emin)/2.0 + emin
!          do itot = 1, ntot
!             do ikpt = 1, nkpt
!                num = num + fermi_dirac((eigval3(itot,ikpt)-e))/real(nkpt)
!             end do
!          end do
!          if (abs(num-nelec) < 1e-6) then
!             efermi = e
!             exit
!          else if (num < nelec) then
!             emin = e
!          else if (num > nelec) then
!             emax = e
!          end if
!       end do
!    end if

    ! Calculate new density matrix
    newdens = cmplx_0
    do itot = 1, ntot
       do ikpt = 1, nkpt
          !if (eigenval(itot,ikpt) > efermi) then
          !   cycle
          !end if
          do isite = 1, nsite
             do ispin = 1, nspin
                do ispin1 = 1, nspin
                   do iorb = 1, norb
                      do iorb1 = 1, norb
                         newdens(iorb,iorb1,isite,ispin,ispin1) = &
                         newdens(iorb,iorb1,isite,ispin,ispin1) + fermi_dirac(eigval3(itot,ikpt)-efermi)*  &
                                                                   conjg(eigvec3((isite-1)*norb+iorb+(ispin-1)*norb*nsite,itot,ikpt))* &
                                                                   eigvec3((isite-1)*norb+iorb1+(ispin1-1)*norb*nsite,itot,ikpt)/nkpt
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    ! Calculate energy
    e3 = 0.0
    do itot = 1, ntot
       do ikpt = 1, nkpt
          if (eigval3(itot,ikpt) > efermi) then
             cycle
          end if

          e3 = e3 + eigval3(itot,ikpt)/nkpt

       end do
    end do 

    do ispin = 1, nspin
       do ispin1 = 1, nspin
          do isite = 1, nsite
             do iorb = 1, norb
                do iorb1 = 1, norb
                   e3 = e3 - 0.5*vnew(iorb,iorb1,isite,ispin,ispin1)*newdens(iorb,iorb1,isite,ispin,ispin1)
                end do
             end do
          end do
       end do
    end do

  end subroutine calc_3de

end module prop
