include ../../../Makefile.${ARCH}

# Names of files that should be compiled and 
# inserted to lib${SIM}.a:

OBJS = user.o advection_propagator.o example_advection.o

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
override INCS+=-I../../include

INCS_REG=${INCS} -I../../dataoperator -I../../gridbuilder

DEPS_USER=../../include/user.h user.cpp
DEPS_ADV_PROP=example_advection.h advection_propagator.h advection_propagator.cpp
DEPS_EX_ADV=example_advection.h example_advection.cpp

# Compilation rules

.cpp.o:
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< ${INCS}

user.o: ${DEPS_USER}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c user.cpp ${INCS} ${INCS_REG}

advection_propagator.o: ${DEPS_ADV_PROP}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c advection_propagator.cpp ${INCS}

example_advection.o: ${DEPS_EX_ADV}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c example_advection.cpp ${INCS}

