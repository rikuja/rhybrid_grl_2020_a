### MAKEFILE FOR CREATING PARGRID DISTRIBUTION ###
###               DO NOT EDIT                  ###

# Distribution package name
DIR=pargrid
DIST=pargrid_v01_100.tar

# Build targets

clean:
	rm -rf *~ *.tar *.tar.gz

dist:

	ln -s ${CURDIR} ${DIR}
	tar -rf ${DIST} ${DIR}/Makefile
	tar -rf ${DIST} ${DIR}/COPYING*
	tar -rf ${DIST} ${DIR}/*.h
	tar -rf ${DIST} ${DIR}/*.pdf
	gzip -9 ${DIST}
	rm ${DIR}
