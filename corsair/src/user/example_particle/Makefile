include ../../../Makefile.${ARCH}

# Names of files that should be compiled and 
# inserted to libgridbuilder.a:

OBJS = fields_B_constant.o\
	particle_species.o\
	simcontrol.o user.o

# Make targets

default:
	${MAKE} lib${SIM}.a

clean:
	rm -rf *.o *.a *~
	rm -f ../lib/lib${SIM}.a

lib${SIM}.a: ${OBJS}
	${AR} r lib${SIM}.a ${OBJS}
	ln -f -s ${CURDIR}/lib${SIM}.a -t ../../lib

# Dependencies

override INCS+=${INC} ${INC_BOOST} ${INC_PARGRID} ${INC_VLSV} ${INC_ZOLTAN}
override INCS+=-I../../include -I../../dataoperator -I../../particleinjector -I../../particlepropagator -I../../gridbuilder

INCS_REG=${INCS} -I../../particleinjector -I../../dataoperator
INCS_REG+=-I../../gridbuilder

DEPS_ACC=particle_accumulator.h
DEPS_INJ=default_injector.h test_injector.h
DEPS_FIELDS=fields_B_constant.h fields_B_constant.cpp
DEPS_PARTICLE_SPECIES=particle_species.h particle_species.cpp
DEPS_PROP=particle_propagator_rk2_gc.h
DEPS_SIMCONTROL=simcontrol.h simcontrol.cpp
DEPS_USER=../../include/user.h user.cpp ../../include/particle_list_skeleton.h ${DEPS_ACC} ${DEPS_INJ} ${DEPS_PROP} ${DEPS_FIELD}

# Compilation rules

.cpp.o:
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< ${INCS}

fields_B_constant.o: ${DEPS_FIELDS}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c fields_B_constant.cpp ${INCS}

particle_accumulator.o: ${DEPS_ACC}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c particle_accumulator.cpp ${INCS}

particle_species.o: ${DEPS_PARTICLE_SPECIES}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c particle_species.cpp ${INCS}

simcontrol.o: ${DEPS_SIMCONTROL}
	${CMP} ${CXXFLAGS} ${FLAGS} -c simcontrol.cpp ${INCS}

user.o: ${DEPS_USER}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c user.cpp ${INCS}
