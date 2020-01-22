include ../../Makefile.${ARCH}

# Names of files that should be compiled and 
# inserted to libgridbuilder.a:

OBJS = gridbuilder.o rect_cuboid_builder.o restart_builder.o\
	logically_cartesian_builder.o

# Make targets

default:
	${MAKE} libgridbuilder.a

clean:
	rm -rf *.o *.a *~
	rm -f ../lib/libgridbuilder.a

libgridbuilder.a: ${OBJS}
	${AR} r libgridbuilder.a ${OBJS}
	ln -f -s ${CURDIR}/libgridbuilder.a -t ../lib

# Adjust include paths, override forces these values to be added 
# to the value(s) passed from earlier Makefile(s)

override INCS+=${INC} ${INC_BOOST} ${INC_PARGRID} ${INC_VLSV} ${INC_ZOLTAN}
override INCS+=-I../include

# Dependencies

DEPS = ../include/simulation.h ../include/constants.h ../include/randomnumber.h ../include/mpilogger.h\
     ../include/definitions.h ../include/cell.h ../include/profiler.h ../include/object_factory_generic.h
DEPS_GRIDBUILDER = ${DEPS} ../include/gridbuilder.h gridbuilder.cpp
DEPS_LCB = ${DEPS} ../include/gridbuilder.h logically_cartesian_builder.h logically_cartesian_builder.cpp
DEPS_RCB = ${DEPS} ../include/gridbuilder.h rect_cuboid_builder.h rect_cuboid_builder.cpp
DEPS_RESTART = ${DEPS} ../include/gridbuilder.h restart_builder.h restart_builder.cpp

# Compilation rules

.cpp.o:
	${CMP} ${CXXFLAGS} ${FLAGS} -c $< ${INCS}

gridbuilder.o: ${DEPS_GRIDBUILDER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c gridbuilder.cpp ${INCS}

logically_cartesian_builder.o: ${DEPS_LCB}
	${CMP} ${CXXFLAGS} ${FLAGS} -c logically_cartesian_builder.cpp ${INCS}

rect_cuboid_builder.o: ${DEPS_RCB}
	${CMP} ${CXXFLAGS} ${FLAGS} -c rect_cuboid_builder.cpp ${INCS}

restart_builder.o: ${DEPS_RESTART}
	${CMP} ${CXXFLAGS} ${FLAGS} -c restart_builder.cpp ${INCS}
