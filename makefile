#--Top level trees makefile--
OBJS    = variables.o main.o define_variables.o \
				  power_parser.o io.o metryc.o shapes.o \
					distribution.o fuels_create.o \
			 		fuel_read.o trees_read.o \
					In_DUET.o VarRoutines_DUET.o \
  				MainLoop_DUET.o support_DUET.o Out_DUET.o
VPATH = .:Modular/DUET

DUET_DIR  = ModularDUET

FC = gfortran   ## compiler
FFLAGS = -O2 -ffixed-line-length-none

Inputs/trees.exe: ${OBJS}
	${FC} -o $@ ${FFLAGS} $^

clean: 
	rm Inputs/trees.exe *.o *.mod Inputs/trees*.dat

# optimized suffix rules
.SUFFIXES: .f .f90

%.o: ${DUET_DIR}/%.f90
	${FC} ${FFLAGS} -c $<
%.o: %.f
	${FC} ${FFLAGS} -c $<
%.o: %.f90
	${FC} ${FFLAGS} -c $<
