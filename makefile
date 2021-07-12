## Commands below for pink - must have loaded modules with 
#              module load lampi intel

## Linux cluster (ie, pink) options - select only one FFLAGS definition
#----------------------------------------
  FC = gfortran   ## compiler
  #MPILIB =  -L${MPI_ROOT}/lib64 -lmpi -I${MPI_ROOT}/include/
  #MPILIB =  -I${MPI_ROOT}/include/
  FFLAGS = -O2 -ffixed-line-length-none 
#  FFLAGS = -g -extend_source 
#   FFLAGS = -g -extend_source -warn all -check bounds -check overflow 

JSONDIR = json-fortran

OBJS = variables.o main.o define_variables.o \
			 io.o metryc.o shapes.o distribution.o \
			 baseline.o treatments.o \
			 json_kinds.o json_parameters.o \
			 json_string_utilities.o json_value_module.o \
			 json_file_module.o json_module.o fuel_read.o \

trees: ${OBJS}
	${FC} -o $@  ${FFLAGS} ${OBJS} ${MPILIB}

clean: 
	rm trees *.o *.mod *.dat

# optimized suffix rules
.SUFFIXES: .f .F90

.f.o:
	${FC} ${FFLAGS} -c $<
.F90.o:
	${FC} ${FFLAGS} -c $<
