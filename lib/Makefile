include ../config.mk
#
OBJS = constants.o commLine.o geometry.o geometryDP.o stepPotentials.o intList.o txtUtils.o
LIBFILE = libDiscrete.a

libDiscrete: $(OBJS)
	ar r $(LIBFILE) $(OBJS) 
	
txtUtils.o: $<
	$(F77) -ffree-form $(OPTIONS) $(CFLAGSNOCHECK) txtUtils.f
.f.o: $(INCS) $<
	$(F77) -ffree-form $(OPTIONS) $(CFLAGS) $<
.F.o: $(INCS) $<
	$(F77) -ffree-form $(OPTIONS) $(CFLAGS) $<
clean: 
	rm *.o *.mod $(LIBFILE)
