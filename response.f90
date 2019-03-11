module responsefunc

  complex*16, save, public, allocatable        ::  r(:,:,:,:,:,:,:), ur(:,:,:,:,:,:,:,:), invur(:,:,:,:,:,:,:,:)

  public  :: calc_response_tensor
  public  :: calc_response_instable
  public  :: io_response
  public  :: calc_uporr
  public  :: calc_invur

  contains

  subroutine calc_response_tensor()

    use param,           only : norb, nsub, nsite, ntot, hf_vec, hf_val, nkpt, efermi
    use constants,       only : nspin, cmplx_0, cmplx_i, delta 

    implicit none

    integer        :: isub, isub1, ikpt, iorb, iorb1, iorb2, iorb3, ispin, ispin1, isite, isite1, is, is1
    complex*16     :: term1(norb,norb,nsite,norb,norb,nsite), term2(norb,norb,nsite,norb,norb,nsite)
    real*8         :: diff, fac 

    allocate(r(norb,norb,nsite,norb,norb,nsite,nspin*nspin))

    ! up-up ispin = 1
    ! dw-up ispin = 2
    ! up-dw ispin = 3
    ! dw-dw ispin = 4

    r = cmplx_0

    do ispin = 1, nspin
       do ispin1 = 1, nspin

          do ikpt = 1, nkpt
          do isub = 1, nsub
             if ((hf_val(isub,ispin,ikpt)-efermi) .ge. 0.0) cycle
          do isub1 = 1, nsub
             if ((hf_val(isub1,ispin1,ikpt)-efermi) .le. 0.0) cycle

             do isite = 1, nsite
             do isite1 = 1, nsite
                do iorb = 1, norb
                do iorb1 = 1, norb
                do iorb2 = 1, norb
                do iorb3 = 1, norb

                   term1(iorb,iorb1,isite,iorb2,iorb3,isite1) = &
                         dconjg(hf_vec(iorb+(isite-1)*norb,isub,ispin,ikpt))    *&
                                hf_vec(iorb1+(isite-1)*norb,isub1,ispin1,ikpt)  *&
                         dconjg(hf_vec(iorb2+(isite1-1)*norb,isub1,ispin1,ikpt))*&
                                hf_vec(iorb3+(isite1-1)*norb,isub,ispin,ikpt)

                   term2(iorb,iorb1,isite,iorb2,iorb3,isite1) = &
                         dconjg(hf_vec(iorb+(isite-1)*norb,isub1,ispin1,ikpt))   *&
                                hf_vec(iorb1+(isite-1)*norb,isub,ispin,ikpt)   *&
                         dconjg(hf_vec(iorb2+(isite1-1)*norb,isub,ispin,ikpt)) *&
                                hf_vec(iorb3+(isite1-1)*norb,isub1,ispin1,ikpt)
  
                end do
                end do
                end do
                end do
             end do
             end do

             diff = hf_val(isub,ispin,ikpt) - hf_val(isub1,ispin1,ikpt)
             fac = (1.0/real(nkpt))*(1.0/(diff + cmplx_i*delta))

             if (ispin .eq. ispin1) then
                if (ispin .eq. 1) is = 1
                if (ispin .eq. 2) is = 4
                r(:,:,:,:,:,:,is) = r(:,:,:,:,:,:,is) + term1(:,:,:,:,:,:)*fac
                r(:,:,:,:,:,:,is) = r(:,:,:,:,:,:,is) + term2(:,:,:,:,:,:)*fac
             else
                if (ispin .eq. 1) then
                   is = 3
                   is1 = 2
                end if
                if (ispin .eq. 2) then
                   is = 2
                   is1 = 3
                end if
                r(:,:,:,:,:,:,is) = r(:,:,:,:,:,:,is) + term1(:,:,:,:,:,:)*fac
                r(:,:,:,:,:,:,is1) = r(:,:,:,:,:,:,is1) + term2(:,:,:,:,:,:)*fac
             end if

          end do     
          end do
          end do

       end do
    end do

  end subroutine calc_response_tensor

  subroutine calc_response_instable(mag0)

    use param,           only : norb, nsite, ntot, nkpt, sigma_x, sigma_y, sigma_z
    use constants,       only : nspin, cmplx_0, cmplx_i
    use utilities,       only : trace

    implicit none

    integer                  :: itot, itot1, ikpt, iorb, iorb1, iorb2, iorb3, ispin, ispin1, isite, isite1
    complex*16               :: cx,cy,c2,c3,trx, try, trxmu, trymu, mu2, lmbx, lmby, fact(norb,norb,nspin,nspin)
    real*8, intent(in)       :: mag0(nsite)
   
    mu2 = cmplx_0
    do isite = 1, nsite
       mu2 = mu2 + mag0(isite)*mag0(isite)
    end do

    mu2 = mu2*real(norb)*2.0

    do isite1 = 1, nsite
       do iorb2 = 1, norb
          do iorb3 = 1, norb
             lmbx = cmplx_0
             lmby = cmplx_0
             trxmu = cmplx_0
             trymu = cmplx_0
             do isite = 1, nsite
                trx = cmplx_0
                try = cmplx_0
                do iorb = 1, norb
                trx = trx + r(iorb,iorb,isite,iorb2,iorb3,isite1,2) + &
                            r(iorb,iorb,isite,iorb2,iorb3,isite1,3)
                try = try - cmplx_i*(r(iorb,iorb,isite,iorb2,iorb3,isite1,2) - &
                                     r(iorb,iorb,isite,iorb2,iorb3,isite1,3))  
                end do
                trymu = trymu + trx*mag0(isite)
                trxmu = trxmu - try*mag0(isite)
             end do
             lmby = trxmu/mu2
             lmbx = trymu/mu2
             do isite = 1, nsite
                cx = lmby*mag0(isite)
                cy = -lmbx*mag0(isite)
                c2 = cmplx(real(cx),real(-cy))
                c3 = cmplx(real(cx),real(cy))
                do iorb = 1, norb
                   r(iorb,iorb,isite,iorb2,iorb3,isite1,2) = r(iorb,iorb,isite,iorb2,iorb3,isite1,2) - c2
                   r(iorb,iorb,isite,iorb2,iorb3,isite1,3) = r(iorb,iorb,isite,iorb2,iorb3,isite1,3) - c3
                end do
             end do
          end do
       end do
    end do                  

  end subroutine calc_response_instable

  subroutine io_response(mode)

    use param,           only : norb, nsite
    use constants,       only : nspin, rname, cmplx_0

    implicit none

    character(len=1), intent(in)      :: mode
    integer                           :: iorb,iorb1,iorb2,iorb3,isite,isite1,ispin,ispin1, i
    real*8                            :: A,B

    ! Read and write, the FASTEST INDEX IS THE FIRST ONE 
    ! THE SLOWEST IS THE LAST ONE
    if (mode=='w') then
       open(unit=11, file=rname, status='replace')
           do ispin = 1, nspin*nspin
              do isite = 1, nsite
                 do isite1 = 1, nsite
                    do iorb = 1, norb
                       do iorb1 = 1, norb
                          do iorb2 = 1, norb
                             do iorb3 = 1, norb
                                write(11,'(7i,2ES18.10)') iorb, iorb1, isite, iorb2, iorb3, isite1, ispin, &
                                   & real(r(iorb,iorb1,isite,iorb2,iorb3,isite1,ispin)), &
                                   & aimag(r(iorb,iorb1,isite,iorb2,iorb3,isite1,ispin))
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        close(11)
     else if (mode=='r') then
       allocate(r(norb,norb,nsite,norb,norb,nsite,nspin*nspin))
       r = cmplx_0
       open(unit=11, file=rname, status='old')
          do i = 1, (norb**4)*(nsite**2)*(nspin**2)
             read(11,*) iorb, iorb1, isite, iorb2, iorb3, isite1, ispin, A, B
             r(iorb,iorb1,isite,iorb2,iorb3,isite1,ispin) = cmplx(A,B)  
          end do
       close(11)
     end if
  
  end subroutine io_response

  subroutine calc_uporr()

    use param,           only : norb, nsite, u_mat
    use constants,       only : nspin, cmplx_0

    implicit none

    integer              :: iorb, iorb1, iorb2, iorb3, iorb4, iorb5, ispin, ispin1, isite, isite1
    complex*16           :: u(norb,norb,norb,norb,nsite), j(norb,norb,norb,norb,nsite)

    allocate(ur(norb,norb,nsite,norb,norb,nsite,nspin*2,nspin*2))

    do isite = 1, nsite
       do iorb = 1, norb
          do iorb1 = 1, norb
             do iorb2 = 1, norb
                do iorb3 = 1, norb
                   u(iorb,iorb1,iorb2,iorb3,isite) = u_mat(iorb,iorb1,iorb2,iorb3,isite)
                   j(iorb,iorb1,iorb2,iorb3,isite) = u_mat(iorb,iorb3,iorb2,iorb1,isite)
                end do
             end do
          end do
       end do
    end do

    ur = cmplx_0
    do isite = 1, nsite
       ! \sum_cd U_abcd_i * R_cd_i_ef_j
       do iorb = 1, norb
          do iorb1 = 1, norb
             do isite1 = 1, nsite
                do iorb4 = 1, norb
                   do iorb5 = 1, norb  
                   ! Begin summation over cd  
                   do iorb2 = 1, norb
                   do iorb3 = 1, norb
                      ur(iorb,iorb1,isite,iorb4,iorb5,isite1,1,1) = ur(iorb,iorb1,isite,iorb4,iorb5,isite1,1,1) + &
                                                                   (u(iorb,iorb1,iorb2,iorb3,isite)-j(iorb,iorb1,iorb2,iorb3,isite))*&
                                                                    r(iorb2,iorb3,isite,iorb4,iorb5,isite1,1)
                      ur(iorb,iorb1,isite,iorb4,iorb5,isite1,1,4) = ur(iorb,iorb1,isite,iorb4,iorb5,isite1,1,4) + &
                                                                    u(iorb,iorb1,iorb2,iorb3,isite)*r(iorb2,iorb3,isite,iorb4,iorb5,isite1,4)
                      ur(iorb,iorb1,isite,iorb4,iorb5,isite1,2,2) = ur(iorb,iorb1,isite,iorb4,iorb5,isite1,2,2) - &
                                                                    j(iorb,iorb1,iorb2,iorb3,isite)*r(iorb2,iorb3,isite,iorb4,iorb5,isite1,3)
                      ur(iorb,iorb1,isite,iorb4,iorb5,isite1,3,3) = ur(iorb,iorb1,isite,iorb4,iorb5,isite1,3,3) - &
                                                                    j(iorb,iorb1,iorb2,iorb3,isite)*r(iorb2,iorb3,isite,iorb4,iorb5,isite1,2)
                      ur(iorb,iorb1,isite,iorb4,iorb5,isite1,4,1) = ur(iorb,iorb1,isite,iorb4,iorb5,isite1,4,1) + &
                                                                    u(iorb,iorb1,iorb2,iorb3,isite)*r(iorb2,iorb3,isite,iorb4,iorb5,isite1,1)
                      ur(iorb,iorb1,isite,iorb4,iorb5,isite1,4,4) = ur(iorb,iorb1,isite,iorb4,iorb5,isite1,4,4) + &
                                                                   (u(iorb,iorb1,iorb2,iorb3,isite)-j(iorb,iorb1,iorb2,iorb3,isite))*&
                                                                    r(iorb2,iorb3,isite,iorb4,iorb5,isite1,4)
                   end do
                   end do
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine calc_uporr

  subroutine calc_invur()

    use param,           only : norb, nsite
    use constants,       only : nspin, cmplx_0 
    use utilities,       only : call_invert

    implicit none

    integer             :: iorb, iorb1, iorb2, iorb3, isite, isite1, ispin, ispin1, flat1, flat2
    complex*16          :: urflat(norb*norb*nsite*nspin*2,norb*norb*nsite*nspin*2), imat(norb*norb*nsite*nspin*2,norb*norb*nsite*nspin*2)

    allocate(invur(norb,norb,nsite,norb,norb,nsite,nspin*2,nspin*2))

    urflat = cmplx_0
    imat = cmplx_0
    invur = cmplx_0     

    ! Reshape to flatten rank6 tensors into rank2 matrices
    flat1 = 0
    do ispin = 1, nspin*2
       do isite = 1, nsite
          do iorb = 1, norb
             do iorb1 = 1, norb
                flat1 = flat1 + 1
                flat2 = 0
                do ispin1 = 1, nspin*2
                   do isite1 = 1, nsite
                      do iorb2 = 1, norb
                         do iorb3 = 1, norb
                            flat2 = flat2 + 1
                            urflat(flat1,flat2) = ur(iorb,iorb1,isite,iorb2,iorb3,isite1,ispin,ispin1)
                         end do
                      end do    
                   end do
                end do
             end do
          end do
       end do
    end do

    flat1 = 0
    do flat1 = 1, norb*norb*nsite*nspin*2
       imat(flat1,flat1) = 1.0
    end do

    urflat(:,:) = imat(:,:) - urflat(:,:)  

    call call_invert(urflat,norb*norb*nsite*nspin*2)

    ! Reconvert rank2 matrix to rank6 tensor

    flat1 = 0
    do ispin = 1, nspin*2
       do isite = 1, nsite
          do iorb = 1, norb
             do iorb1 = 1, norb
                flat1 = flat1 + 1
                flat2 = 0
                do ispin1 = 1, nspin*2
                   do isite1 = 1, nsite
                      do iorb2 = 1, norb
                         do iorb3 = 1, norb
                            flat2 = flat2 + 1
                            invur(iorb,iorb1,isite,iorb2,iorb3,isite1,ispin,ispin1) = urflat(flat1,flat2) 
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    
  end subroutine calc_invur

end module responsefunc
