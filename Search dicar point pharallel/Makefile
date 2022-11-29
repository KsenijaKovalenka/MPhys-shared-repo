main = main.f90
OBJmain = main.o
FFLAGS = -O -ffree-line-length-0

 
FC = gfortran
LIB=-llapack -lblas

main : ${OBJmain} 
	${FC} ${FFLAGS} -o main main.o ${LIB}

${OBJmain} : ${@:.o=.f90} ${main}
	${FC} ${FFLAGS} -c ${@:.o=.f90}

clean:
	rm *.o *.mod main *.c
