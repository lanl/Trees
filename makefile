#--Top level trees makefile--
OBJS    = variables.o main.o define_variables.o \
				  power_parser.o io.o metryc.o shapes.o \
					distribution.o baseline.o duet_inputs.o \
			 		fuel_read.o trees_read.o \

FC = gfortran   ## compiler
FFLAGS = -O2 -ffixed-line-length-none

Inputs/trees.exe: ${OBJS}
	${FC} -o $@ ${FFLAGS} ${OBJS}

clean: 
	rm Inputs/trees.exe *.o *.mod Inputs/trees*.dat

# optimized suffix rules
.SUFFIXES: .f .f90

.f.o:
	${FC} ${FFLAGS} -c $<
.f90.o:
	${FC} ${FFLAGS} -c $<
