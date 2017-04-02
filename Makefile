V = 0.2.6
EXEDIR =.
#
.PHONY: lib setup discrete
#
all: setup discrete
	
lib: 
	cd lib; make

discrete: lib
	cd discrete; make

install: all
	cp discrete/discrete $(EXEDIR)/discrete.exe

clean:
	cd lib;make clean
	cd discrete; make clean
	rm $(EXEDIR)/discrete.exe

        	
