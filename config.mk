DEBUG =
#DEBUG = -g -fbacktrace -fbounds-check 
#DEBUG = -g -fbacktrace -fbounds-check -DDEBUG
DEBUGNOCHECK = 
#DEBUGNOCHECK = -g -fbacktrace -DDEBUG
#PROFILE =
PROFILE = -pg
#
OPTIONS = -O3 -fimplicit-none -funroll-loops -mtune=generic -Wall -ftree-vectorizer-verbose=1 -ffast-math  
LOPTIONS= -static-libgfortran  
F77= /usr/local/bin/gfortran
CPP= cpp
LIBDIR = ../lib
CFLAGS = -c $(DEBUG) $(PROFILE) -I$(LIBDIR)
CFLAGSNODEBUG = -c $(PROFILE) -I$(LIBDIR)
CFLAGSNOCHECK = -c $(DEBUGNOCHECK) $(PROFILE) -I$(LIBDIR)
LFLAGS = $(PROFILE) -L$(LIBDIR) 
ENDFLAG = -fconvert=little-endian
