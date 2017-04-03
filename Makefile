EXEDIR = exe
#
.PHONY: lib setup discrete
#
all: install 
	
lib: 
	cd lib; make

discrete: lib
	cd discrete; make

install: discrete
	cp discrete/discrete $(EXEDIR)/godmd

clean:
	cd lib;make clean
	cd discrete; make clean
	rm $(EXEDIR)/godmd

        	
