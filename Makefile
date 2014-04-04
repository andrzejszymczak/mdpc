
CC = g++
OPT =  -O3 -I. 
LIBOPT = -lm -lGL -lglut -lGLEW

all : mdpc msvis
	make -C subd

%.o: %.cpp *.h Makefile
	$(CC) $(OPT) -c -o $@ $< 

mdpc : pcenv.o pcvf.o vfield_base.o mdpc.o mesh_base.o mesh.o tgraph.o primitive.o pcstable.o mstype.o tskel.o pchull.o Makefile
	$(CC) $(OPT) -o mdpc pcenv.o pcvf.o vfield_base.o mdpc.o mesh_base.o mesh.o tgraph.o primitive.o pcstable.o mstype.o tskel.o pchull.o $(LIBOPT)

msvis : msvis.o program.o trackball.o pcvfdisplay.o primset.o primitive.o pcvf.o mesh.o mesh_base.o vfield_base.o Makefile
	$(CC) $(OPT) -o msvis msvis.o program.o trackball.o pcvfdisplay.o primset.o primitive.o pcvf.o mesh.o mesh_base.o vfield_base.o $(LIBOPT)


clean :
	rm *.o msvis mdpc *~
	cd subd
	make -C subd clean
