CPLEXDIR=/path_to_cplex_rppt/cplex128/cplex
SYSTEM=x86-64_linux
LIBFORMAT=static_pic

CONCERTDIR= $(CPLEXDIR)/../concert
CPLEXBINDIR=$(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIBDIR=$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR=$(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread -ldl

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
CCC = g++ -std=gnu++0x -O4
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

###################################

all: hitndrive getHTMatrixInversion buildGraph

hitndrive: hitndrive.o
        $(CCC) $(CCFLAGS) $(CCLNDIRS) -o $@ hitndrive.o $(CCLNFLAGS)
hitndrive.o: hitndrive.cpp
        $(CCC) -c $(CCFLAGS) hitndrive.cpp -o $@

getHTMatrixInversion: getHTMatrixInversion.o
        $(CCC) $(CCOPT) -o $@ getHTMatrixInversion.o
getHTMatrixInversion.o: getHTMatrixInversion.cpp
        $(CCC) -c $(CCOPT) getHTMatrixInversion.cpp -o $@

buildGraph: buildGraph.o
        $(CCC) $(CCOPT) -o $@ buildGraph.o
buildGraph.o: buildGraph.cpp
        $(CCC) -c $(CCOPT) buildGraph.cpp -o $@

clean:
        rm -f *.o hitndrive getHTMatrixInversion buildGraph
