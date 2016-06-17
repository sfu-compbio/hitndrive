CPLEXDIR=
CPLEX_BUILD=

CPLEXINC=$(CPLEXDIR)/cplex/include
CPLEXLIB=$(CPLEXDIR)/cplex/lib/$(CPLEX_BUILD)/static_pic
CONCERTINC=$(CPLEXDIR)/concert/include
CONCERTLIB=$(CPLEXDIR)/concert/lib/$(CPLEX_BUILD)/static_pic
CPLEXFLAGS=$(CPLEXLIB)/libilocplex.a $(CPLEXLIB)/libcplex.a $(CONCERTLIB)/libconcert.a 

CC:=g++

FLAGS=-O3 -std=gnu++0x -Iinclude -pthread
FLAGSWITHCPLEX=-O3 -std=gnu++0x -DIL_STD -I$(CPLEXINC) -I$(CONCERTINC)

EXE1=buildGraph
SRC1=buildGraph.cpp
OBJ1=buildGraph.o

EXE2=getHTMatrixInversion
SRC2=getHTMatrixInversion.cpp
OBJ2=getHTMatrixInversion.o

EXE3=hitndrive
SRC3=hitndrive.cpp
OBJ3=hitndrive.o

all: $(EXE1) $(EXE2) $(EXE3)
	

$(EXE1): $(OBJ1)
	$(CC) $(OBJ1) $(FLAGS) -o $@

$(EXE2): $(OBJ2)
	$(CC) $(OBJ2) $(FLAGS) -o $@

$(EXE3): $(OBJ3)
	$(CC) $(OBJ3) $(FLAGS) $(CPLEXFLAGS) -o $@

.cpp.o:
	$(CC) -c $(FLAGSWITHCPLEX) $^ -o $@

clean:
	rm -f *.o $(EXE1) $(EXE2) $(EXE3)

