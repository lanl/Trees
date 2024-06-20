#--Top level trees makefile--
OBJS    = MakeFFInputs_DUET.o variables_DUET.o main_DUET.o define_variables_DUET.o \
				  io_DUET.o full_DUET.o \
#OBJS    = $(patsubst %,$(ODIR)/%,$(OBJS_))

FC = gfortran   ## compiler
FFLAGS = -O2 -ffixed-line-length-none

duet: ${OBJS}
	${FC} -o $@ ${FFLAGS} ${OBJS}

clean: 
	rm duet *.o *.mod duet*.dat

# optimized suffix rules
.SUFFIXES: .f .f90

.f.o:
	${FC} ${FFLAGS} -c $<
.f90.o:
	${FC} ${FFLAGS} -c $<
