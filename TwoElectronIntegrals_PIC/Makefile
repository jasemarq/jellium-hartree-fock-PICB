LIB        = -L. -lm
INCLUDE    = -I.
CFLAGS     = -O2
EXEC       = pq_int.exe
CXX        = g++

${EXEC}: pq_int.c gauss_legendre.o
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} pq_int.c gauss_legendre.o -o ${EXEC}

Radial.o: gauss_legendre.c gauss_legendre.h
	${CXX} ${LIB} -c gauss_legendre.c ${CFLAGS}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<
