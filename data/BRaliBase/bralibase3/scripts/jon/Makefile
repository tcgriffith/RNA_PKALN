LIBS     =   -lm

CFLAGS   =   -O3

CC       =   gcc

OBJECTS_1  =   main.o ancstates.o errorexit.o gamma.o linalg.o complex.o matrices.o structures.o tokenize.o
OBJECTS_2  =   stomb.o errorexit.o

PROGS    =   rnanc stomb

all:		$(PROGS)

rnanc:		$(OBJECTS_1)
		$(CC) $(CFLAGS) $(OBJECTS_1) $(LIBS) -o rnanc

stomb:		$(OBJECTS_2)
		$(CC) $(CFLAGS) $(OBJECTS_2) $(LIBS) -o stomb

main.o:	main.c
		$(CC) $(CFLAGS) -c main.c
		
ancstates.o:	ancstates.c
		$(CC) $(CFLAGS) -c ancstates.c

errorexit.o:	errorexit.c
		$(CC) $(CFLAGS) -c errorexit.c

gamma.o:	gamma.c
		$(CC) $(CFLAGS) -c gamma.c

linalg.o:	linalg.c
		$(CC) $(CFLAGS) -c linalg.c

complex.o:	complex.c
		$(CC) $(CFLAGS) -c complex.c

matrices.o:	matrices.c
		$(CC) $(CFLAGS) -c matrices.c

structures.o:	structures.c
		$(CC) $(CFLAGS) -c structures.c

tokenize.o:	tokenize.c
		$(CC) $(CFLAGS) -c tokenize.c

stomb.o:	stomb.c
		$(CC) $(CFLAGS) -c stomb.c

clean:		
		rm -f *.o
