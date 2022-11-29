main = searchdirac.f90
OBJmain = searchdirac.o
FFLAGS = -O -ffree-line-length-0

 
FC = gfortran
LIB=-llapack -lblas

searchdirac : ${OBJmain} 
	${FC} ${FFLAGS} -o searchdirac searchdirac.o ${LIB}

${OBJmain} : ${@:.o=.f90} ${main}
	${FC} ${FFLAGS} -c ${@:.o=.f90}

clean:
	rm *.o *.mod searchdirac *.c
