COMPIL = gfortran
CFLAGS=-O0 
CFLAGS=-O2
LFLAGS=-o
#CFLAGS=-O2 
#LFLAGS=-O2 -L$(HOME)
COPTS = $(CFLAGS)
LOPTS = $(LFLAGS)
OBJS =	deftype.o\
	solveur_openmp.o\
	maillage.o\
	operateur.o\
  maillage_sphere.o\
	main.o
PROG = solve
#
.SUFFIXES: .f90 .o
SUFF = .f90.o
$(SUFF):
	$(COMPIL) -c $(CFLAGS) $*.f90 -I.
#
$(PROG): $(OBJS)
	@echo "Compilation terminee"
	@echo "Creation de l'executable $(PROG)"
	$(COMPIL) $(LFLAGS) $(OBJS) -o $(PROG)
#
clean:
	@echo "Nettoyage"
	@rm -f core *.o *.mod ./solve




