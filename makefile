#--Top level trees makefile--
OBJS    = variables.o main.o define_variables.o \
				  io.o metryc.o shapes.o distribution.o \
			 		baseline.o duet.o treatments.o \
			 		fuel_read.o trees_read.o \
#OBJS    = $(patsubst %,$(ODIR)/%,$(OBJS_))

FC = gfortran   ## compiler
FFLAGS = -O2 -ffixed-line-length-none

trees: ${OBJS}
	${FC} -o $@ ${FFLAGS} ${OBJS}

clean: 
	rm trees *.o *.mod *.dat

# optimized suffix rules
.SUFFIXES: .f .f90

.f.o:
	${FC} ${FFLAGS} -c $<
.f90.o:
	${FC} ${FFLAGS} -c $<
