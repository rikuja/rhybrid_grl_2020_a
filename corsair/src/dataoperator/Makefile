include ../../Makefile.${ARCH}

# Names of files that should be compiled and 
# inserted to libdataoperator.a:

OBJS = dataoperator.o dataoperatorcontainer.o\
     operator_cellid.o operator_load.o operator_macroparticles.o operator_mpirank.o\
	operator_pargrid_array.o\
	operator_particle.o operator_timeseries_load.o

# Make targets

default:
	${MAKE} libdataoperator.a -f Makefile

clean:
	rm -rf *.o *.a *~
	rm -f ../lib/libdataoperator.a

libdataoperator.a: ${OBJS}
	${AR} r libdataoperator.a ${OBJS}
	ln -f -s ${CURDIR}/libdataoperator.a -t ../lib

# Adjust include paths

override INCS+=${INC} ${INC_BOOST} ${INC_PARGRID} ${INC_VLSV} ${INC_ZOLTAN}
override INCS+=-I../include

# Dependencies

DEPS = ../include/simulation.h ../include/constants.h ../include/randomnumber.h ../include/mpilogger.h\
     ../include/definitions.h ../include/cell.h ../include/profiler.h

DEPS_DATACONT = ${DEPS} ../include/dataoperator.h\
	      ../include/dataoperatorcontainer.h dataoperatorcontainer.cpp
DEPS_DATAOPERATOR = ${DEPS} ../include/dataoperator.h dataoperator.cpp
DEPS_OP_CELLID = ${DEPS} ../include/dataoperator.h operator_cellid.h operator_cellid.cpp
DEPS_OP_LOAD = ${DEPS} ../include/dataoperator.h operator_load.h operator_load.cpp
DEPS_OP_TSERIES_LOAD = ${DEPS} ../include/dataoperator.h operator_timeseries_load.h operator_timeseries_load.cpp
DEPS_OP_MACROP = ${DEPS} ../include/dataoperator.h operator_macroparticles.h operator_macroparticles.cpp 
DEPS_OP_MPIRANK = ${DEPS} ../include/dataoperator.h operator_mpirank.h operator_mpirank.cpp
DEPS_OP_PARGRID_ARRAY=${DEPS} ../include/dataoperator.h operator_pargrid_array.h operator_pargrid_array.cpp
DEPS_OP_PARTICLE=${DEPS} ../include/dataoperator.h operator_particle.h operator_particle.cpp

# Compilation rules

.cpp.o:
	${CMP} ${CXXFLAGS} ${FLAGS} -c $< ${INCS}

dataoperator.o: ${DEPS_DATAOPERATOR}
	${CMP} ${CXXFLAGS} ${FLAGS} -c dataoperator.cpp ${INCS}

dataoperatorcontainer.o: ${DEPS_DATACONT}
	${CMP} ${CXXFLAGS} ${FLAGS} -c dataoperatorcontainer.cpp ${INCS}

operator_cellid.o: ${DEPS_CELLID}
	${CMP} ${CXXFLAGS} ${FLAGS} -c operator_cellid.cpp ${INCS}

operator_load.o: ${DEPS_LOAD}
	${CMP} ${CXXFLAGS} ${FLAGS} -c operator_load.cpp ${INCS}

operator_macroparticles.o: ${DEPS_OP_MACROP}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c operator_macroparticles.cpp ${INCS}

operator_mpirank.o: ${DEPS_OP_MPIRANK}
	${CMP} ${CXXFLAGS} ${FLAGS} -c operator_mpirank.cpp ${INCS}

operator_pargrid_array.o: ${DEPS_OP_PARGRID_ARRAY}
	${CMP} ${CXXFLAGS} ${FLAGS} -c operator_pargrid_array.cpp ${INCS}

operator_particle.o: ${DEPS_OP_PARTICLE}
	${CMP} ${CXXFLAGS} ${FLAGS} -c operator_particle.cpp ${INCS}

operator_timeseries_load.o: ${DEPS_OP_TSERIES_LOAD}
	${CMP} ${CXXFLAGS} ${FLAGS} -c operator_timeseries_load.cpp ${INCS}
