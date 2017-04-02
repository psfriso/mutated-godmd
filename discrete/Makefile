include ../config.mk
#
OBJS = IO.o random.o energy.o colisio.o DMD.o dims_utils.o least_sq_fit.o ANM.o
MODS = paramSet.o  Structure.o doubleDMD.o ls_rmsd.o
EXEFILE = discrete

discrete: ${MODS} $(OBJS) $(LIBDIR) main.o
	$(F77) $(LFLAGS) $(ENDFLAG) -o $(EXEFILE) $(LOPTIONS) $(OBJS) $(MODS) main.o -l Discrete 

main.o: vars.h

Structure.o: Structure.f 
	$(F77) -ffree-form $(OPTIONS) $(CFLAGSNOCHECK) Structure.f

.f.o: $<
	$(F77) -ffree-form $(OPTIONS) $(CFLAGS) $<
.F.o: $<
	$(F77) -ffree-form $(OPTIONS) $(CFLAGS) $<
clean: 
	rm *.o *.mod $(EXEFILE)
