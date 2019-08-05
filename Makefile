#=============================================#
#       Upper level makefile for Ingrid       #
#=============================================#

all:
	cd src; make; cd -
	cd Docs; make html; cd -

clean:
	cd src; make clean; cd -
	cd Docs; make clean; cd -
