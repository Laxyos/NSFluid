COMPIL = gfortran
CFLAGS=-O0 
#CFLAGS=-O2
LFLAGS=-o
#CFLAGS=-O2 
#LFLAGS=-O2 -L$(HOME)
COPTS = $(CFLAGS)
LOPTS = $(LFLAGS)
OBJS =	deftype.o\
	solveur_openmp.o\
	operateur.o\
	sortie.o\
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
	rm *.plt
	@echo "Nettoyage"
	@rm -f core *.o *.mod ./solve




