program localforce 
  use :: param
  use :: constants
  use :: utilities
  use :: fractions
  use :: integral
  use :: greenfunction 
  use :: spindispersion
  use :: responsefunc
  use :: perturbations
  use :: prop
  use :: dminter

  implicit none

  integer                 :: ipol, ipath, irpt, irpt1, isite, ispin, ispin1, iorb, iorb1, idos
  complex*16              :: e, jbond, lf1, lf2, force(3), fifj1(3), fifj2(3)
  complex*16, allocatable :: n(:,:,:,:,:), jll(:), resforc(:,:)
  real*8                  :: emin, emax, pdos1, pdos2, dum, ez(3), ey(3), ex(3)
  real*8, allocatable     :: occ(:), mag(:), ev0(:,:), ev(:,:), dev(:,:)
  character(len=50)       :: filename
  logical                 :: io_find


  open(unit=12, file='gf.out', status='replace')
    write(12,*) '!-------------------------------------------------!'
    write(12,*) '!       Green-Function extraction and local       !'
    write(12,*) '!       force theorem for magnetic coupling       !'
    write(12,*) '!----------------- Oier Arcelus ------------------!'
    write(12,*) ' '
    write(12,*) ' '
    write(12,*) ' Reading parameters from Hartree-Fock outputs...'
    write(12,*) ' ------------------------------------------------'
    write(12,*) ' ------------------------------------------------'
    write(12,*) ' '
    call get_names()
    call get_latt_parameters()
    if (lmgnon) then
    call get_maglatt_parameters() 
    end if
    call get_wav_parameters()
    call get_pauli_matrix()
    if (lmgnon) then
    call get_path_parameters()
    end if
    call get_eig_parameters()
    call get_bond_parameters()
    call get_pot_parameters()
    call get_dens_parameters()
    write(12,*) ' Calculation parameters:'
    write(12,*) ' '
    write(12,*) ' norb = ', norb
    write(12,*) ' nsite = ', nsite
    write(12,*) ' nkpt = ', nkpt
    write(12,*) ' nspin = ', nspin
    write(12,*) ' ntot = ', ntot
    write(12,*) ' npoles = ', npoles
    write(12,*) ' npath = ', npath
    write(12,*) ' nrpt = ', nrpt
    write(12,*) ' fdtemp = ', fdtemp
    write(12,*) ' lexch = ', lexch
    write(12,*) ' lmgnon = ', lmgnon
    write(12,*) ' '
    write(12,*) ' Calculation of SCLR: '
    write(12,*) ' lsclr = ', lsclr
    if (lsclr) write(12,*) ' Find the results from SCLR in file sclr.out' 
    write(12,*) ' '    

    call get_residues()

    write(12,*) ' Residues and poles calculated, check poles.dat, separation.dat and frac.dat '
    write(12,*) ' to verify that separation between poles is correct and FD function looks    '
    write(12,*) ' correct when compared with other aproximations'
    write(12,*) ' '

    open(unit=14, file='poles.dat', status='replace')
    do ipol = 1, 2*npoles
       write(14,*) 0.0, polres(1,ipol), polres(2,ipol)
    end do
    close(14)
    open(unit=14, file='separation.dat', status='replace')
    do ipol = 1, npoles-1
       write(14,*) ipol, polres(1,ipol+npoles+1)-polres(1,ipol+npoles), twopi
    end do
    close(14)

    emax = maxval(hf_val(:,:,:))
    emin = minval(hf_val(:,:,:))

    open(unit=14, file='frac.dat', status='replace')
    do idos = 1, 10000
       e = emin + idos*(emax-emin)/10000.0
       write(14,'(5F10.4)') real(e), real(cfermi_dirac(e-efermi)), real(matsubara(e-efermi)), &
                   & real(frac_fd(e-efermi)), real(expfrac(e-efermi))
    end do
    close(14) 

    write(12,*) ' Check error.dat to see absolute error of the FD aproximation'
    write(12,*) ' is low enough. '
    open(unit=14, file='error.dat', status='replace')
    do idos = 1, 10000
       e = emin + idos*(emax-emin)/10000.0
       write(14,*) real(e), abs(real(cfermi_dirac(e-efermi))-real(frac_fd(e-efermi)))
    end do
    close(14)
    write(12,*) ' '
    
    call define_path()

    write(12,*) ' Semicircle for line integral defined in path.dat '
    open(unit=14, file='path.dat', status='replace')
    do ipath = 0, npath
       write(14,*) real(Cpath(ipath)), aimag(Cpath(ipath)) 
    end do
    close(14)

    allocate(gr1(norb,norb,nspin))
    allocate(gr2(norb,norb,nspin))
    allocate(gr1b(norb*nspin,norb*nspin))
    allocate(gr2b(norb*nspin,norb*nspin))
    allocate(greenR(norb,norb,nspin))

    allocate(n(norb,norb,nsite,nspin,nspin))

    allocate(jll(nrpt))

    jll = cmplx_0

    if (lexch) then
       write(12,*) ' Calculating exchange parameters and printing them to lf.dat '
       write(12,'(a)') ' '
       open(unit=14, file='lf.dat', status='replace')
       do irpt = 1, nrpt
          jbond = cmplx_0
          do irpt1 = 1, nrpt
             if ((atom_label(1,irpt).eq.atom_label(2,irpt1)).and.(atom_label(2,irpt).eq.atom_label(1,irpt1)) .and. &
                 (norm2(r_bond(:,irpt)+r_bond(:,irpt1))<1d-6)) then

                do ipath = 0, npath-1
                   call gfR(Cpath(ipath),irpt)
                   gr1 = greenR
                   call gfR(Cpath(ipath),irpt1)
                   gr2 = greenR
                   call lich(irpt)
                   lf1 = lf*frac_fd(Cpath(ipath)-efermi)
                   call gfR(Cpath(ipath+1),irpt)
                   gr1 = greenR
                   call gfR(Cpath(ipath+1),irpt1)
                   gr2 = greenR
                   call lich(irpt)
                   lf2 = lf*frac_fd(Cpath(ipath+1)-efermi)
                   jbond = jbond - (1.0/(2.0*twopi))*(Cpath(ipath+1)-Cpath(ipath))*(lf2+lf1)
                end do
                do ipol = 1, npoles
                   call gfR(efermi+cmplx_i*kboltz*fdtemp*polres(1,ipol),irpt)
                   gr1 = greenR
                   call gfR(efermi+cmplx_i*kboltz*fdtemp*polres(1,ipol),irpt1)
                   gr2 = greenR
                   call lich(irpt)
                   jbond = jbond + cmplx_i*kboltz*fdtemp*polres(2,ipol)*lf
                end do
                jll(irpt) = aimag(jbond)
                write(14,'(i1,1x,i1,2x,3F10.4,2x,F12.6)') atom_label(1,irpt),atom_label(2,irpt),r_bond(:,irpt),1000*aimag(jbond)
             end if
          end do
       end do
       close(14)
    end if

    call get_cw()
    write(12,*) ' '
    write(12,*) ' CW temp:', cw

    if (lmgnon) then
       write(12,*) ' Spin wave dispersion calculation (MAKE SURE THE LF.DAT FILE IS COMPLETE!!!)'
       write(12,*) ' '
        
       call spin_wav() 
    
       write(12,*) ' Spin wave dispersion DONE '
       write(12,*) '  '
    end if

  close(12)

  if (lsclr) then
  open(unit=12, file='sclr.out', status='replace')
    write(12,*) '!-------------------------------------------------!'
    write(12,*) '!          SSS      CCC      LL      RRRRR        !'
    write(12,*) '!         SS       CC        LL      RR  RR       !'
    write(12,*) '!          SS     CC         LL      RRRRRR       !'
    write(12,*) '!           SS     CC        LL      RR RR        !'
    write(12,*) '!         SSS       CCC      LLLL    RR  RR       !'
    write(12,*) '!----------------- Oier Arcelus ------------------!'
    write(12,*) ' '
    write(12,*) ' '
    write(12,*) ' Get inputs for SCLR'
    write(12,*) ' '
    write(12,*) ' lsph = ', lsph
    write(12,*) ' '
    ! get dvext

    write(12,*) ' Getting the SOC potential...'
    write(12,*) ' '
  
    call get_so_potential()
    if (lorb) call get_l_operator()

    do ispin = 1, nspin
       do ispin1 = 1, nspin
          write(12,*) 'spin1 = ', ispin, ' spin2 = ', ispin1
          do isite = 1, nsite
             write(12,*) 'site = ', isite
             write(12,'(a)') ' '
             do iorb = 1, norb
                write(12,'(<2*norb>F9.4)')(dvext(iorb,iorb1,isite,ispin,ispin1), iorb1 = 1, norb)
             end do
             write(12,'(a)') ' '
          end do
       end do
    end do

    n = cmplx_0

    n(:,:,:,1,1) = hf_den(:,:,:,1)
    n(:,:,:,2,2) = hf_den(:,:,:,2)

    allocate(mag(nsite))
    do isite = 1, nsite
       mag(isite) = real(trace(hf_den(:,:,isite,1)-hf_den(:,:,isite,2),norb))
    end do

    write(12,*) ' SOC potential generated'
    write(12,*) ' '
    call calc_e0(n)
    ! get umat
    call get_coulomb_matrix()
    ! transform if spherical parameters are used
    if (lsph) call spherical_param()
    ! See if file for r exists, if not create it
    ! CAREFULL IF R DOES NOT EXIST: FIRST GENERATE IT
    ! WITH MAGNETIZATION COLLINEAR TO Z
    inquire(file=rname, exist=io_find)
    if (.not.io_find) then    
       ! Calculate r and solve rotate instabilities
       call calc_response_tensor()
       call calc_response_instable(mag)
       ! Write response tensor to sclr_r.out
       call io_response('w')
    else
       ! Read response tensor from sclr_r.out
       call io_response('r')
       ! R IS DEFINED IN THE Z-AXIS, BUT MAGNETIZATION IS DEFINED IN WHATEVER
       ! DIRECTIONS FROM THE GREENS FUNCTIONS ABOVE
       ! Calculate UR 4x4 matrix
       call calc_uporr()
       ! Calculate [1-UR]^-1
       call calc_invur()
       ! get theeta and phi angles
       call get_angles()
       ! get spin rotation matrix
       call get_spin_rot(th,phi)
       ! rotate SO potential
       call utility_spin_rotation(dvext,urot,norb,nsite,'b','o')
       ! Calculate dvp, dn and dv from previous matrices
       call calc_dvp(dvext)
       call calc_dn()
       call calc_dv()
       ! rotate n, dvp, dn and dv back
       call utility_spin_rotation(dvext,urot,norb,nsite,'f','o')
       call utility_spin_rotation(dvp,urot,norb,nsite,'f','o')
       call utility_spin_rotation(dn,urot,norb,nsite,'b','n')  
       call utility_spin_rotation(dv,urot,norb,nsite,'f','o')
       call utility_spin_rotation(n,urot,norb,nsite,'b','n')
       call utility_spin_rotation(hf_potex,urot,norb,nsite,'f','o') ! For the 3rd order correction

       write(12,*) ' -----------------------------------'
       write(12,*) ' '
       write(12,*) ' '
       write(12,*) ' -----------------------------------'
       write(12,*) ' Energy Values'
       write(12,'(a,F12.6,a)') ' 0th order of SO: ', real(e0), ' eV'
       call calc_de(dv,n)
       write(12,'(a,F12.6,a)') ' 1st order of SO: ',real(e0)+ real(de), ' eV'
       write(12,'(a,F12.6,a)') '      Correction: ', real(de), ' eV'
       call calc_de(dvext,dn)
       write(12,'(a,F12.6,a)') ' 2nd order of SO: ', real(e0) +real(de), ' eV'
       write(12,'(a,F12.6,a)') '      Correction: ', real(de), ' eV'
       call calc_3de(dvp,hf_potex)
       write(12,'(a,F12.6,a)') ' 3rd order of SO: ', real(e3), ' eV <--- Still experimental'
       write(12,'(a,F12.6,a)') '      Correction: ', real(e3)-real(e0), ' eV <--- Still experimental'
       write(12,*) ' -----------------------------------'


       ! Calculate properties
       call calc_mut0(n)
       write(12,*) ' '
       write(12,*) ' '
       write(12,*) ' -----------------------------------'
       write(12,*) ' Magnetization values'
       write(12,*) ' 0. Order of SO: '
       do isite = 1, nsite
          write(12,'(a,i2,a,F10.4,a,F10.4,a,F10.4)') '   Atom: ', isite, ' Sx: ', real(mut0(1,isite)), &
                                                        ' Sy: ', real(mut0(2,isite)), &
                                                        ' Sz: ', real(mut0(3,isite))
       end do
       call calc_muts(dn)
       if (lorb) call calc_mutl(dn)
       write(12,*) ' 1st Order of SO: '
       do isite = 1, nsite
          write(12,'(a,i2,a,F10.4,a,F10.4,a,F10.4)') '   Atom: ', isite, ' Sx: ', real(muts(1,isite)), &
                                                        ' Sy: ', real(muts(2,isite)), &
                                                        ' Sz: ', real(muts(3,isite))
          if (lorb) then
          write(12,'(a,F10.4,a,F10.4,a,F10.4)')      '            Lx: ', real(mutl(1,isite)), &
                                                        ' Ly: ', real(mutl(2,isite)), &
                                                        ' Lz: ', real(mutl(3,isite))
          end if
       end do
       write(12,*) ' -----------------------------------'
       write(12,*) ' '
       write(12,*) ' Calculating force components due to dvp '
       if (ldmi) then
          allocate(resforc(3,nrpt))
          allocate(ev0(3,nsite))
          allocate(ev(3,nsite))
          allocate(dev(3,nsite))
          ! Normalize direction of magnetization 
          do isite = 1, nsite
             ev0(:,isite) = real(mut0(:,isite))/norm2(real(mut0(:,isite)))
             ev(:,isite)  = real(muts(:,isite))/norm2(real(muts(:,isite)))
             dev(:,isite) = ev(:,isite) - ev0(:,isite)
          end do
          ez = 0.0
          ex = 0.0
          ey = 0.0
          ez(3) = 1.0
          ey(2) = 1.0
          ex(1) = 1.0
          if (norm2(ev0(:,1) - ez(:)) < 1d-7) then
             open(unit=14, file='forces_z.dat', status='replace')
          else if (norm2(ev0(:,1) - ex(:)) < 1d-7) then
             open(unit=14, file='forces_x.dat', status='replace')
          else if (norm2(ev0(:,1) - ey(:)) < 1d-7) then
             open(unit=14, file='forces_y.dat', status='replace')
          else
             print*, 'You activated ldmi = .true., but the direction of the magnetization'
             print*, 'has to be along z, y, or x axis, please set the appropiate th and phi'
             stop
          end if
          do isite = 1, nsite
             write(14,'(3F12.6)') dev(:,isite)
          end do
          write(14,*) ' '

          do irpt = 1, nrpt
             force = cmplx_0
             do irpt1 = 1, nrpt
             if ((atom_label(1,irpt).eq.atom_label(2,irpt1)).and.(atom_label(2,irpt).eq.atom_label(1,irpt1)).and. &
                 (norm2(r_bond(:,irpt)+r_bond(:,irpt1))<1d-6)) then

                do ipath = 0, npath-1
                   ! Value of fifj in Cpath
                   call gfR(Cpath(ipath),irpt)
                   call utility_gr_rotation(greenR,urot,norb,gr1b,'f')
                   call gfR(Cpath(ipath),irpt1)
                   call utility_gr_rotation(greenR,urot,norb,gr2b,'f')
                   call lfdmi(irpt,dvp,sigpot)
                   fifj1 = fifj*frac_fd(Cpath(ipath)-efermi)
                   ! Value of fifj in Cpath+1
                   call gfR(Cpath(ipath+1),irpt)
                   call utility_gr_rotation(greenR,urot,norb,gr1b,'f')
                   call gfR(Cpath(ipath+1),irpt1)
                   call utility_gr_rotation(greenR,urot,norb,gr2b,'f')
                   call lfdmi(irpt,dvp,sigpot)
                   fifj2 = fifj*frac_fd(Cpath(ipath+1)-efermi)
                   force = force + (1.0/(2.0*twopi))*(Cpath(ipath+1)-Cpath(ipath))*(fifj2+fifj1)
                end do
                do ipol = 1, npoles
                   call gfR(efermi+cmplx_i*kboltz*fdtemp*polres(1,ipol),irpt)
                   call utility_gr_rotation(greenR,urot,norb,gr1b,'f')
                   call gfR(efermi+cmplx_i*kboltz*fdtemp*polres(1,ipol),irpt1)
                   call utility_gr_rotation(greenR,urot,norb,gr2b,'f')
                   call lfdmi(irpt,dvp,sigpot)
                   force = force - cmplx_i*kboltz*fdtemp*polres(2,ipol)*fifj
                end do
                resforc(:,irpt) = aimag(force(:))
                write(14,'(i1,1x,i1,2x,3F10.4,2x,3F12.6)') atom_label(1,irpt),atom_label(2,irpt),r_bond(:,irpt),1000*aimag(force(:))
             end if
             end do
          end do
          close(14) 
       end if
    end if
    if (ldmi) then
       call calc_dmi()
       open(unit=14, file='dmi.dat', status='replace')
       do irpt = 1, nrpt
          write(14,'(i1,1x,i1,2x,3F10.4,2x,3F12.6)') atom_label(1,irpt),atom_label(2,irpt),r_bond(:,irpt), dmi(:,irpt)
       end do
       close(14)
    end if
  close(12)
  end if

end program localforce
