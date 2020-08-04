all:
	cd Efit; make; cd -
	cd Uegrid; make; cd -

clean:
	cd Efit; make clean; cd -
	cd GUI; make clean; cd -; cd Topologies; make clean; cd -
	cd Interpol; make clean; cd -

	rm -f *.so *.pyc *.rst~ *.nml~
	rm -rf __pycache__
