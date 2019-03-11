module constants

implicit none

complex*16       , public, parameter   :: cmplx_i = (0.0,1.0), cmplx_0 = (0.0,0.0)
real*8           , public, parameter   :: twopi = 2*acos(-1.0), kboltz = 8.6173324d-5, delta = 1d-8
integer          , public, parameter   :: nspin = 2
character(len=50), public, parameter   :: wavname     = 'wav.out',      &
                                          eigname     = 'eig.out',      &
                                          lattname    = 'struct.in',    &
                                          bondname    = 'bond.out',     &
                                          maglattname = 'maglatt.in',   &
                                          potname     = 'pot.out',       &
                                          densname    = 'dens.out',      &
                                          pathname    = 'bandpath.in',   &
                                          hamname     = 'sclr_ham.in',   &
                                          hamsoname   = 'sclr_hamso.in', &
                                          coulombname = 'sclr_u.in',     &
                                          rotname     = 'sclr_rot.in',   &
                                          rname       = 'sclr_r.out',    &
                                          socname     = 'sclr_vso.in',   &
                                          orbitname   = 'sclr_l.in'
end module constants
