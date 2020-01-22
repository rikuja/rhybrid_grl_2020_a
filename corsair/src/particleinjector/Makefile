include ../../Makefile.${ARCH}

# Names of files that should be compiled and 
# inserted to libparticleinjector.a:

OBJS=particleinjector.o particleinjectorfactory.o

# Make targets
default:
	${MAKE} libparticleinjector.a

clean:
	rm -rf *.o *.a *~
	rm -f ../lib/libparticleinjector.a

# Adjust include paths

override INCS+=${INC} ${INC_BOOST} ${INC_PARGRID} ${INC_VLSV} ${INC_ZOLTAN}
override INCS+=-I../include

# Dependencies

DEPS = ../include/simulation.h ../include/constants.h ../include/randomnumber.h ../include/mpilogger.h\
     ../include/definitions.h ../include/cell.h ../include/profiler.h ../include/boundary_particle.h
DEPS_FACTORY = ${DEPS} ../include/particleinjector.h ../include/particleinjectorfactory.h particleinjectorfactory.cpp
DEPS_INJECTOR = ${DEPS} ../include/particleinjector.h particleinjector.cpp

# Compilation rules
.cpp.o:
	${CMP} ${CXXFLAGS} ${FLAGS} -c $< ${INCS}

particleinjector.o: ${DEPS_INJECTOR}
	${CMP} ${CXXFLAGS} ${FLAGS} -c particleinjector.cpp ${INCS}

particleinjectorfactory.o: ${DEPS_FACTORY}
	${CMP} ${CXXFLAGS} ${FLAGS} -c particleinjectorfactory.cpp ${INCS}

libparticleinjector.a: ${OBJS}
	${AR} r libparticleinjector.a ${OBJS}
	ln -f -s ${CURDIR}/libparticleinjector.a -t ../lib
