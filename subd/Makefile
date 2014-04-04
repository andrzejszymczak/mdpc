
CC = gcc
OPT =  -O3 -I. 
LIBOPT = -lm 

all : vfsubd

%.o: %.cpp *.h Makefile
	$(CC) $(OPT) -c -o $@ $< 

vfsubd : mesh.o vfsubd.o Makefile
	$(CC) $(OPT) -o vfsubd mesh.o vfsubd.o $(LIBOPT)

clean :
	rm *.o vfsubd *~
