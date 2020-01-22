include ../Makefile.${ARCH}

override LIBS += -lsimkernel -l${SIM} -lplistbase -ldataoperator -lgridbuilder
override LIBS += ${LIB_OPENMP} ${LIB_BOOST} ${LIB_PARGRID} ${LIB_VLSV} ${LIB_ZOLTAN}

default:
	${MAKE} "ARCH=${ARCH}" "INCS=${INCS}" -C dataoperator
	${MAKE} "ARCH=${ARCH}" "INCS=${INCS}" -C gridbuilder
	${MAKE} "ARCH=${ARCH}" "INCS=${INCS}" -C kernel
	${MAKE} "ARCH=${ARCH}" "INCS=${INCS}" -C user/${SIM}

	${LNK} -o corsair_${SIM} -Llib ${LIBS}
	mv corsair_${SIM} ../

clean:
	rm -f *~
	rm -f lib/*
	${MAKE} clean -C include
	${MAKE} clean "ARCH=${ARCH}" -C dataoperator
	${MAKE} clean "ARCH=${ARCH}" -C gridbuilder
	${MAKE} clean "ARCH=${ARCH}" -C kernel
	${MAKE} clean "ARCH=${ARCH}" -C user/${SIM}

simclean:
	rm -f lib/lib${SIM}.a
	${MAKE} clean "ARCH=${ARCH}" -C user/${SIM}
