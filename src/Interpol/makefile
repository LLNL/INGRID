#----------------------------------------------------------------------------              
#          Makefile for the Interpol package                                                   
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

FFILES = bicubic2.f90
EXECFIL = bicubicf


$(EXECFIL) : $(FFILES)
	python -m numpy.f2py -c -m $(EXECFIL) $(FFILES)



clean:
	rm -f *.so *.pyc
	rm -r *.so.dSYM