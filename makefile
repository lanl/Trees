#--Top level trees makefile--
JSONDIR:= ./jsonfortran
VPATH   = ${JSONDIR}
ODIR    = obj
OBJS    = variables.o main.o define_variables.o \
					json_kinds.o json_parameters.o \
			 		json_string_utilities.o json_value_module.o \
			 		json_file_module.o json_module.o \
				  io.o metryc.o shapes.o distribution.o \
			 		baseline.o treatments.o \
			 		fuel_read.o trees_read.o \
#OBJS    = $(patsubst %,$(ODIR)/%,$(OBJS_))

FC = gfortran   ## compiler
FFLAGS = -O2 -ffixed-line-length-none

trees: ${OBJS}
	${FC} -o $@ ${FFLAGS} ${OBJS}

clean: 
	rm trees *.o *.mod *.dat

# optimized suffix rules
.SUFFIXES: .f .F90

.f.o:
	${FC} ${FFLAGS} -c $<
.F90.o:
	${FC} ${FFLAGS} -c $<
