#----------------------------------------------------------------------------
#          Makefile for the Efit package
#----------------------------------------------------------------------------

ifdef CPUNAME
  # do nothing -- CPUNAME defined in the environment vars
else
  # Try to regognize the machine automatically
  CPUNAME = $(shell uname -n)
endif


#ifneq ($(findstring sin,$(CPUNAME)),)  #singe.llnl.gov
#  FF = gfortran
#  LINKER = gfortran
#endif

FFILES = readg.f90
EXECFIL = efit


$(EXECFIL) : $(FFILES)
	python -m numpy.f2py -c -m $(EXECFIL) $(FFILES)



clean: 
	rm -f *.so *.pyc
	rm -r *.so.dSYM
