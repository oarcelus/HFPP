module param 

  implicit none

  character(len=50), save, public                :: atomname

  ! Calculation parameters
  complex*16     , save, public, allocatable     :: hf_pot(:,:,:,:), hf_vec(:,:,:,:), dpot(:,:,:), dvext(:,:,:,:,:), u_mat(:,:,:,:,:), &
                                                    sigma_x(:,:,:,:), sigma_y(:,:,:,:), sigma_z(:,:,:,:), lx(:,:), ly(:,:), lz(:,:), &
                                                    urot(:,:,:,:), sigpot(:,:,:,:,:,:), hf_potex(:,:,:,:,:), hf_den(:,:,:,:)
  integer        , save, public, allocatable     :: atom_label(:,:), maglatt(:)!, ml(:)
  integer        , save, public                  :: nsite, npth, nsub, norb, nkpt, ntot, nrpt, npoles, npath, totpth, npnts, nelec
  real*8         , save, public                  :: lattvec(3,3), recipvec(3,3)
  real*8         , save, public, allocatable     :: kpt_cart(:,:), atom_dir(:,:), atom_cart(:,:), r_bond(:,:), hf_val(:,:,:), &
                                                    kpt_pathr(:,:)
  real*8         , save, public                  :: volume, efermi, fdtemp, sspin, th, phi
  logical        , save, public                  :: lorb, lexch, ldmi, lmgnon, lsph, lsclr

  character(len=5), save, public                 :: soctype

  ! Subroutines
  public    :: get_names
  public    :: get_pauli_matrix
  public    :: get_spin_rot
  public    :: get_l_operator
  public    :: get_latt_parameters
  public    :: get_maglatt_parameters
  public    :: get_wav_parameters
  public    :: get_path_parameters
  public    :: get_eig_parameters
  public    :: get_bond_parameters
  public    :: get_pot_parameters
  public    :: get_so_potential
  public    :: get_coulomb_matrix
  public    :: spherical_param

  contains

  subroutine get_names()

    open(unit=1, file= 'gflr.in', status='old') 
       read(1,*) ! Space for separating file names and parameters
       read(1,*) npoles ! Number of poles for the FD aproximation
       read(1,*) npath  ! Number of points in the integration path
       read(1,*) sspin  ! Spin number of hubbard sites
       read(1,*) fdtemp ! Temperature of the FD approximation
       read(1,*) nelec  ! number of electrons
       read(1,*) lexch  ! .true. if magnetic coupling calculated
       read(1,*) lmgnon ! .true. if spin wave dispersion calculated
       read(1,*) ! Below here SCLR module parameters
       read(1,*) lsclr  ! .true. if SCLR module activated
          if (lsclr) then
             read(1,*) lsph ! .true. if HF calculation were done using spherical parameters (carefull with the spherical values)
             read(1,*) lorb ! .true. if orbital operators are readed and orbital magnetization claculated
             read(1,*) ldmi ! .true. if DMI are calculated. Make sure to calculate lexch = .true. first
             read(1,*) soctype ! 'r': As H_so - H in same wannier basis (local part). 'f': From file sclr_vso.in
             read(1,*) th, phi ! read direction of magnetization
!             read(1,*) nmom ! nmom = 1 for p; nmom = 2 for d electrons
!             allocate(ml(2*nmom+1))
!             read(1,*) ml(:)! ml values, norb of them, give them in order from small to high,e.g., -2 -1 0 1 2
          end if
    close(1)

  end subroutine get_names

  subroutine get_latt_parameters()

    use utilities,                only : utility_recip_lattice, utility_frac_to_cart
    use constants,                only : lattname

    implicit none 
    logical                       :: io_find
    integer                       :: ilatt, isite

    inquire(file=lattname, exist=io_find)
    if (io_find) then
       open(2, file=lattname, status='old')
          read(2,*)
          do ilatt = 1, 3
             read(2,*) lattvec(:,ilatt)
          end do
          read(2,*) atomname, nsite
          allocate(atom_dir(3,nsite))
          allocate(atom_cart(3,nsite))
          do isite = 1, nsite
             read(2,*) atom_dir(:,isite)
          end do
       close(2)
    else
       print*, 'No lattvec file found' 
       stop
    end if

    call utility_recip_lattice(lattvec, recipvec, volume)
    do isite = 1, nsite
       call utility_frac_to_cart(atom_dir(:,isite),atom_cart(:,isite),lattvec)
    end do

  end subroutine get_latt_parameters

  subroutine get_maglatt_parameters()

    use constants,                only : maglattname

    implicit none

    logical                       :: io_find
    integer                       :: isite

    allocate(maglatt(nsite))

    inquire(file=maglattname, exist=io_find)
    if (io_find) then
       open(2, file=maglattname, status='old')
          do isite = 1, nsite
             read(2,*) maglatt(isite)
          end do
       close(2)
    end if

  end subroutine get_maglatt_parameters

  subroutine get_wav_parameters()

    use constants,                only : wavname, nspin

    implicit none

    integer                       :: ikpt, ispin, itot, idum, idum1
    logical                       :: io_find
    real*8                        :: A, B

    inquire(file=wavname, exist=io_find)
    if (io_find) then
       open(2, file=wavname, status='old')
          read(2,*) nkpt
          read(2,*) nsub
          allocate(kpt_cart(3,nkpt))
          allocate(hf_vec(nsub,nsub,nspin,nkpt))
          do ikpt = 1, nkpt
             read(2,*) kpt_cart(:,ikpt)
             do ispin = 1, nspin
                read(2,*) 
                do itot = 1, nsub*nsub
                   read(2,*) idum, idum1, A, B
                   hf_vec(idum,idum1,ispin,ikpt) = cmplx(A,B)
                end do
                read(2,*)
             end do
          end do
       close(2)
    else
       print*, 'No wav.out file found'
       stop
    end if

    ntot = nsub*nspin  
    norb = nsub/nsite
    allocate(hf_pot(norb,norb,nsite,nspin))
    allocate(hf_den(norb,norb,nsite,nspin))

  end subroutine get_wav_parameters


  subroutine get_pauli_matrix()

    use constants,                only : cmplx_0, cmplx_i, nspin

    implicit none

    integer              :: ispin, ispin1, iorb, iorb1

    allocate(sigma_x(norb,norb,nspin,nspin))
    allocate(sigma_y(norb,norb,nspin,nspin))
    allocate(sigma_z(norb,norb,nspin,nspin))

    sigma_x = cmplx_0
    sigma_y = cmplx_0
    sigma_z = cmplx_0

    allocate(urot(norb,norb,nspin,nspin))

    do iorb = 1, norb
       do ispin = 1, nspin
          do ispin1 = 1, nspin
             if (ispin == 1 .and. ispin1 == 1) then
                sigma_x(iorb,iorb,ispin,ispin1) = cmplx_0
                sigma_y(iorb,iorb,ispin,ispin1) = cmplx_0
                sigma_z(iorb,iorb,ispin,ispin1) = 1.0
             else if (ispin == 2 .and. ispin1 == 2) then
                sigma_x(iorb,iorb,ispin,ispin1) = cmplx_0
                sigma_y(iorb,iorb,ispin,ispin1) = cmplx_0
                sigma_z(iorb,iorb,ispin,ispin1) = -1.0
             else if (ispin == 1 .and. ispin1 == 2) then
                sigma_x(iorb,iorb,ispin,ispin1) = 1.0
                sigma_y(iorb,iorb,ispin,ispin1) = -1.0*cmplx_i
                sigma_z(iorb,iorb,ispin,ispin1) = cmplx_0
             else if (ispin == 2 .and. ispin1 == 1) then
                sigma_x(iorb,iorb,ispin,ispin1) = 1.0 
                sigma_y(iorb,iorb,ispin,ispin1) = cmplx_i
                sigma_z(iorb,iorb,ispin,ispin1) = cmplx_0
             end if
          end do
       end do
    end do

  end subroutine get_pauli_matrix

  subroutine get_angles()

    use constants,        only : twopi

    implicit none
 

    th = th*twopi/360.0
    phi = phi*twopi/360.0

  end subroutine get_angles

  subroutine get_spin_rot(theeta,phii)

    use constants,               only : cmplx_i, nspin, cmplx_0

    implicit none

    integer              :: iorb
    real*8, intent(in)   :: theeta, phii

    do iorb = 1, norb
       urot(iorb,iorb,1,1) = cos(theeta/2.0)*exp(-1.0*cmplx_i*phii/2.0) 
       urot(iorb,iorb,2,2) = cos(theeta/2.0)*exp(cmplx_i*phii/2.0)
       urot(iorb,iorb,1,2) = -1.0*sin(theeta/2.0)*exp(-1.0*cmplx_i*phii/2.0)
       urot(iorb,iorb,2,1) = sin(theeta/2.0)*exp(cmplx_i*phii/2.0)
    end do

  end subroutine get_spin_rot

  subroutine get_l_operator()

  ! THE FORM OF THE L OPERATORS ARE TAKEN FROM GAVS PAPER OF SERGEY ONLY FOR
  ! TESTING, ORDER OF ORBITALS IS Z2, XY, X2-Y2
  
    use constants,               only : cmplx_i, cmplx_0, orbitname

    implicit none

    integer              :: iorb, iorb1
    logical              :: io_find
    real*8               :: i(norb,norb), j(norb,norb)

    allocate(lx(norb,norb))
    allocate(ly(norb,norb))
    allocate(lz(norb,norb))

    lx = cmplx_0
    ly = cmplx_0
    lz = cmplx_0
    inquire(file=orbitname, exist=io_find)
    if (io_find) then
       open(40, file=orbitname, status='old')
          ! x-direction L
          read(40,*)
          do iorb = 1, norb
             read(40,*) i(iorb,:)
          end do
          read(40,*)
          do iorb = 1, norb
             read(40,*) j(iorb,:)
          end do
          lx = cmplx(i,j)
          read(40,*)
          ! y-direction L
          read(40,*)
          do iorb = 1, norb
             read(40,*) i(iorb,:)
          end do
          read(40,*)
          do iorb = 1, norb
             read(40,*) j(iorb,:)
          end do
          ly = cmplx(i,j)
          read(40,*)
          ! z-direction L
          read(40,*)
          do iorb = 1, norb
             read(40,*) i(iorb,:)
          end do
          read(40,*)
          do iorb = 1, norb
             read(40,*) j(iorb,:)
          end do
          lz = cmplx(i,j)
       else
          print*, 'No sclr_l.out found'
          stop
       end if
 
  end subroutine get_l_operator

  ! Get path kpoitnsgg

  subroutine get_path_parameters()

    use constants,                only : pathname
    use utilities,                only : utility_frac_to_cart

    implicit none

    logical               :: io_find
    integer               :: ipth, iline, num
    real*8, allocatable   :: kpt_pathl(:,:), pathpoints(:,:)

    ! Read paths
    inquire(file=pathname, exist=io_find)
    if (io_find) then
       open(40, file=pathname, status='old')
          read(40,*) npth, npnts
          allocate(pathpoints(3,npth))
          totpth = 0
          do ipth = 1, npth
             read(40,*) pathpoints(:,ipth)
          end do
       close(40)
    else
       print*, 'No path file found'
       stop
    end if

    ! Linearly interpolate the end path points to obtain kpoints along high sym
    ! lines
    totpth = (npth-1)*npnts + 1
    allocate(kpt_pathl(3,totpth))
    allocate(kpt_pathr(3,totpth))

    num = 0
    do ipth = 1, npth - 1
       do iline = 1, npnts
          num = num + 1
          kpt_pathl(:,num) = (iline-1)*(pathpoints(:,ipth+1) - pathpoints(:,ipth))/npnts + pathpoints(:,ipth)
       end do
    end do

    kpt_pathl(:,totpth) = pathpoints(:,npth)

    do ipth = 1, totpth
       call utility_frac_to_cart(kpt_pathl(:,ipth),kpt_pathr(:,ipth),recipvec)
    end do

  end subroutine get_path_parameters

  ! Get parameters

  subroutine get_eig_parameters()

    use constants,                only : eigname, nspin

    implicit none

    integer                       :: ikpt, isub, idum, idum1, ispin
    logical                       :: io_find
    real*8                        :: kpt(3,nkpt)

    inquire(file=eigname, exist=io_find)
    if (io_find) then
       open(2, file=eigname, status='old')
          read(2,*) efermi
          allocate(hf_val(nsub,nspin,nkpt))
          do ikpt = 1, nkpt
             read(2,*) kpt(:,ikpt)
             if (kpt(1,ikpt) /= kpt_cart(1,ikpt) .and. &
                 kpt(2,ikpt) /= kpt_cart(2,ikpt) .and. &
                 kpt(3,ikpt) /= kpt_cart(3,ikpt)) then
                 print*, 'k-point mismatch in wav.out and eig.out'
                 stop
             end if           
             do ispin = 1, nspin       
                read(2,*) 
                do isub = 1, nsub
                   read(2,*) idum1, hf_val(isub,ispin,ikpt)
                end do
                read(2,*)
             end do
          end do
       close(2)
    else
       print*, 'No eig.out file found'
       stop
    end if

  end subroutine get_eig_parameters

  subroutine get_bond_parameters()

    use constants,                only : bondname
    use utilities,                only : utility_frac_to_cart

    implicit none

    integer                       :: irpt
    logical                       :: io_find

    inquire(file=bondname, exist=io_find)
    if (io_find) then
       open(2, file=bondname, status='old')
          read(2,*) nrpt
          allocate(atom_label(2,nrpt))
          allocate(r_bond(3,nrpt))
          do irpt = 1, nrpt
             read(2,*) atom_label(:,irpt)
             read(2,*) r_bond(:,irpt)
          end do
       close(2)
    else
       print*, 'No bond.out file found'
       stop
    end if

  end subroutine get_bond_parameters

  subroutine get_dens_parameters()

    use constants,                only : densname, nspin, cmplx_0

    implicit none

    integer                       :: ind, ind1, iorb, iorb1, ispin, ispin1, isite
    logical                       :: io_find
    real*8                        :: A, B

    inquire(file=densname, exist=io_find)
    if (io_find) then
       open(2, file=densname, status='old')
          do ind = 1, nspin*nsite
             read(2,*) ispin, isite
             do ind1 = 1, norb*norb
                read(2,*) iorb, iorb1, A, B
                hf_den(iorb,iorb1,isite,ispin) = cmplx(A,B)
             end do
             read(2,*)
          end do
       close(2)
    else
       print*, 'No dens.out file found'
       stop
    end if

  end subroutine get_dens_parameters

  subroutine get_pot_parameters()

    use constants,                only : potname, nspin, cmplx_0

    implicit none

    integer                       :: ind, ind1, iorb, iorb1, ispin, ispin1, isite
    logical                       :: io_find
    real*8                        :: A, B

    allocate(u_mat(norb,norb,norb,norb,nsite))
    allocate(hf_potex(norb,norb,nsite,nspin,nspin))
    inquire(file=potname, exist=io_find)
    if (io_find) then
       open(2, file=potname, status='old')
          do ind = 1, nspin*nsite
             read(2,*) ispin, isite
             do ind1 = 1, norb*norb
                read(2,*) iorb, iorb1, A, B
                hf_pot(iorb,iorb1,isite,ispin) = cmplx(A,B)
             end do
             read(2,*)
          end do
       close(2)
    else
       print*, 'No pot.out file found'
       stop
    end if

    hf_potex = cmplx_0
    do isite = 1, nsite
       hf_potex(:,:,isite,1,1) = hf_pot(:,:,isite,1)
       hf_potex(:,:,isite,2,2) = hf_pot(:,:,isite,2)
    end do

    allocate(dpot(norb,norb,nsite))
    do isite = 1, nsite
       dpot(:,:,isite) = hf_pot(:,:,isite,1) - hf_pot(:,:,isite,2)
    end do

    allocate(sigpot(norb,norb,nsite,nspin,nspin,3))
    ! Build vector of \sigma*dpot for DMI interactions
    do isite = 1, nsite
       do ispin = 1, nspin
          do ispin1 = 1, nspin
             ! x-direction
             sigpot(:,:,isite,ispin,ispin1,1) = matmul(dpot(:,:,isite),sigma_x(:,:,ispin,ispin1))
             ! y-direction
             sigpot(:,:,isite,ispin,ispin1,2) = matmul(dpot(:,:,isite),sigma_y(:,:,ispin,ispin1))
             ! z-direction
             sigpot(:,:,isite,ispin,ispin1,3) = matmul(dpot(:,:,isite),sigma_z(:,:,ispin,ispin1))
          end do
       end do
    end do

  end subroutine get_pot_parameters

  subroutine get_so_potential()

    use constants,             only : socname, hamname, hamsoname, cmplx_0, nspin
    use utilities,             only : utility_spin_rotation

    implicit none

    integer                       :: isite, iorb, iorb1
    logical                       :: io_find
    real*8         , allocatable  :: rtempnn(:,:), itempnn(:,:)
    complex*16     , allocatable  :: hamR0(:,:,:,:,:), hamR0_so(:,:,:,:,:), tempnn(:,:)


    allocate(rtempnn(norb*nspin,norb*nspin))
    allocate(itempnn(norb*nspin,norb*nspin))
    allocate(tempnn(norb*nspin,norb*nspin))

    if (soctype == 'r') then
 
       allocate(hamR0(norb,norb,nsite,nspin,nspin))
       allocate(hamR0_so(norb,norb,nsite,nspin,nspin))

       hamR0 = cmplx_0
       hamR0_so = cmplx_0
       inquire(file=hamname, exist=io_find)
       if (io_find) then
          open(2, file=hamname, status='old')
             do isite = 1, nsite
                read(2,*)
                read(2,*)
                do iorb = 1, norb
                   read(2,*) hamR0(iorb,:,isite,1,1)
                end do
                hamR0(:,:,isite,2,2) = hamR0(:,:,isite,1,1)
             end do
          close(2)
       else
          print*, 'No sclr_ham.in file found'
          stop
       end if

       inquire(file=hamsoname, exist=io_find)
       if (io_find) then
          open(2, file=hamsoname, status='old')
             do isite = 1, nsite
                read(2,*)
                read(2,*)
                do iorb = 1, norb*nspin
                   read(2,*) rtempnn(iorb,:)
                end do
                read(2,*)
                do iorb = 1, norb*nspin
                   read(2,*) itempnn(iorb,:)
                end do
                ! reduce from 2n-n to n-n
                do iorb = 1, norb*nspin
                   do iorb1 = 1, norb*nspin
                      tempnn(iorb,iorb1) = cmplx(rtempnn(iorb,iorb1),itempnn(iorb,iorb1))
                   end do
                end do
                ! Format the hamiltonian to the ham_r variable
                hamR0_so(:,:,isite,1,1) = tempnn(1:norb,1:norb)
                hamR0_so(:,:,isite,2,2) = tempnn(norb+1:nspin*norb,norb+1:nspin*norb)
                hamR0_so(:,:,isite,1,2) = tempnn(1:norb,norb+1:nspin*norb)
                hamR0_so(:,:,isite,2,1) = tempnn(norb+1:nspin*norb,1:norb)
             end do
          close(2)
       else
          print*, 'No sclr_hamso.in file found'
          stop
       end if

       allocate(dvext(norb,norb,nsite,nspin,nspin))

       ! Calculate diference between Hamiltonian with LDA and LDA+SO.
       do isite = 1, nsite
          dvext(:,:,isite,1,1) = hamR0_so(:,:,isite,1,1) - hamR0(:,:,isite,1,1)
          dvext(:,:,isite,2,2) = hamR0_so(:,:,isite,2,2) - hamR0(:,:,isite,2,2)
          dvext(:,:,isite,1,2) = hamR0_so(:,:,isite,1,2) - hamR0(:,:,isite,1,2)
          dvext(:,:,isite,2,1) = hamR0_so(:,:,isite,2,1) - hamR0(:,:,isite,2,1)
       end do

    else if (soctype == 'f') then
       inquire(file=socname, exist=io_find)
       if (io_find) then
          open(2, file=socname, status='old')
             do isite = 1, nsite
                read(2,*)
                read(2,*)
                do iorb = 1, norb*nspin
                   read(2,*) rtempnn(iorb,:)
                end do
                read(2,*)
                do iorb = 1, norb*nspin
                   read(2,*) itempnn(iorb,:)
                end do
                ! reduce from 2n-n to n-n
                do iorb = 1, norb*nspin
                   do iorb1 = 1, norb*nspin
                      tempnn(iorb,iorb1) = cmplx(rtempnn(iorb,iorb1),itempnn(iorb,iorb1))
                   end do
                end do
                allocate(dvext(norb,norb,nsite,nspin,nspin))
                ! Format the hamiltonian to the ham_r variable
                dvext(:,:,isite,1,1) = tempnn(1:norb,1:norb)
                dvext(:,:,isite,2,2) = tempnn(norb+1:nspin*norb,norb+1:nspin*norb)
                dvext(:,:,isite,1,2) = tempnn(1:norb,norb+1:nspin*norb)
                dvext(:,:,isite,2,1) = tempnn(norb+1:nspin*norb,1:norb)
             end do
          close(2)
       else
          print*, 'No sclr_vso.in file found'
          stop
       end if
    else
      
       print*, 'No soctype tag in sclr.in, please write how you want to get dvext' 
       stop

    end if
 
  end subroutine get_so_potential

  subroutine get_coulomb_matrix()

    use constants,            only : nspin, coulombname, cmplx_0

    integer                   :: isite, ispin, i, j, k, l, itot
    complex*16  , allocatable :: u_tmp(:,:,:,:)
    real*8                    :: u_real, u_imag
    logical                   :: io_find

    allocate(u_tmp(ntot/2,ntot/2,ntot/2,ntot/2))

    ! Read screened.dat provisionally
    u_mat = cmplx_0
    u_tmp = cmplx_0

    inquire(file=coulombname, exist=io_find)
    if (io_find) then
       open(unit=1, file=coulombname, status='old')
          do itot = 1, ntot*ntot*ntot*ntot/(2*2*2*2)
             read(1,*) i, j, k, l, u_real, u_imag
             u_tmp(i,j,k,l) = cmplx(u_real,u_imag)
          end do
       close(1)
    end if

    do isite = 1, nsite
       u_mat(:,:,:,:,isite) =  u_tmp((isite-1)*norb+1:isite*norb, &
                                     (isite-1)*norb+1:isite*norb, &
                                     (isite-1)*norb+1:isite*norb, &
                                     (isite-1)*norb+1:isite*norb)
    end do

    deallocate(u_tmp)

  end subroutine get_coulomb_matrix

  subroutine spherical_param()

    use constants,            only : cmplx_0

    complex*16, allocatable   :: usph(:), jsph(:)
    integer                   :: isite, iorb, iorb1

    allocate(usph(nsite))
    allocate(jsph(nsite))

    u_mat = cmplx_0
    do isite = 1, nsite
       do iorb = 1, norb
          do iorb1 = 1, norb
             u_mat(iorb,iorb,iorb1,iorb1,isite) = 2.9 !usph(isite)
          end do
      end do
    end do

    do isite = 1, nsite
       do iorb = 1, norb
          do iorb1 = 1, norb
             if (iorb /= iorb1) then
                u_mat(iorb,iorb1,iorb1,iorb,isite) = 0.9 !jsph(isite)
             end if
          end do
       end do
    end do

!  u_mat = cmplx_0
!  u_mat(1,1,1,1,1) = 0.707
!  u_mat(1,1,2,2,1) = 0.509
!  u_mat(1,1,3,3,1) = 0.509
!  u_mat(2,2,1,1,1) = 0.509
!  u_mat(2,2,2,2,1) = 0.676
!  u_mat(2,2,3,3,1) = 0.505
!  u_mat(3,3,1,1,1) = 0.509
!  u_mat(3,3,2,2,1) = 0.505
!  u_mat(3,3,3,3,1) = 0.676

!  u_mat(1,1,1,1,1) = 0.707
!  u_mat(1,2,2,1,1) = 0.086
!  u_mat(1,3,3,1,1) = 0.086
!  u_mat(2,1,1,2,1) = 0.086
!  u_mat(2,2,2,2,1) = 0.676
!  u_mat(2,3,3,2,1) = 0.086
!  u_mat(3,1,1,3,1) = 0.086
!  u_mat(3,2,2,3,1) = 0.086
!  u_mat(3,3,3,3,1) = 0.676

!  u_mat(1,1,1,1,1) = 0.707
!  u_mat(1,2,1,2,1) = 0.086
!  u_mat(1,3,1,3,1) = 0.086
!  u_mat(2,1,2,1,1) = 0.086
!  u_mat(2,2,2,2,1) = 0.676
!  u_mat(2,3,2,3,1) = 0.086
!  u_mat(3,1,3,1,1) = 0.086
!  u_mat(3,2,3,2,1) = 0.086
!  u_mat(3,3,3,3,1) = 0.676


  end subroutine spherical_param

end module param
