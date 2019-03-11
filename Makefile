FC     = ifort

default: hf.x 

hf.x: 
	$(FC) -mkl const.f90 util.f90 parameter.f90 frac.f90 integ.f90 green.f90 spinwav.f90 response.f90 perturb.f90 properties.f90 dmi.f90 main.f90

clean: 
	rm -f  *.mod 

