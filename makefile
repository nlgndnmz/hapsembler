#--------------------------------------------------------+++++++++++++++++++++
#  makefile for Hapsembler-2.0
#
#  "make"                 -> build all the binaries and move them to ./bin
#  "make clean"           -> remove all *.o files and binaries
#
#--------------------------------------------------------+++++++++++++++++++++


# change the following line to use a different C++ compiler
CC = g++

ifdef VDEBUG
TOPT = -Wall -g -O1
else ifdef NO_OMP
TOPT = -Wall -O3
else
TOPT = -Wall -O3 -fopenmp -DUSEOPENMP
endif

ifndef NO_LP
COPT = $(TOPT) -DLPSOLVER
else
COPT = $(TOPT)
endif

default: all
	mv overlappr encore consensr preprocr hapsemblr scarpa ./bin
	rm -f *.o

DNAOBJECTIBLES = Dna DnaRead HapUtils KmerHash KmerIter KmerList Naive ReadAln ReadMatch SmithWaterman
DNAOBJECTS=$(DNAOBJECTIBLES:=.o)
DNAHEADERS=$(DNAOBJECTIBLES:=.h) KmerNode.h GraphDef.h HapSuite.h

HAPOBJECTIBLES = BiEdge Containee EdgeTuple HapUtils MateEdge OvlEdge OvlGraph PathEdge
HAPOBJECTS=$(HAPOBJECTIBLES:=.o)
HAPHEADERS=$(HAPOBJECTIBLES:=.h) BiNode.h DnaTempl.h GraphDef.h HapSuite.h LinkedIter.h LinkedList.h LinkedNode.h MatePair.h OvlGraphTempl.h PathNode.h ReadNode.h

SCAOBJECTIBLES = BiEdge ContigEdge DirNode DnaRead HapUtils Naive ReadAln Scaffer SmithWaterman UEdge UGraph UNode
SCAOBJECTS=$(SCAOBJECTIBLES:=.o)
SCAHEADERS=$(SCAOBJECTIBLES:=.h) BiNode.h ContigNode.h GraphDef.h HapSuite.h LinkedIter.h LinkedList.h LinkedNode.h lp_lib.h

all: overlappr encore consensr preprocr hapsemblr scarpa

overlappr: overlappr.o $(DNAOBJECTS)
	$(CC) $(COPT) overlappr.o $(DNAOBJECTS) -o overlappr
overlappr.o: DnaAlignment.cpp $(DNAHEADERS)
	$(CC) -c $(COPT) -DOVERLAPPR DnaAlignment.cpp -o overlappr.o

encore: encore.o $(DNAOBJECTS)
	$(CC) $(COPT) encore.o $(DNAOBJECTS) -o encore
encore.o: DnaAlignment.cpp $(DNAHEADERS)
	$(CC) -c $(COPT) -DCORRECTR DnaAlignment.cpp -o encore.o

consensr: consensr.o $(DNAOBJECTS)
	$(CC) $(COPT) consensr.o $(DNAOBJECTS) -o consensr
consensr.o: DnaAlignment.cpp $(DNAHEADERS)
	$(CC) -c $(COPT) -DCONSENSR DnaAlignment.cpp -o consensr.o

preprocr: preprocr.o $(DNAOBJECTS)
	$(CC) $(COPT) preprocr.o $(DNAOBJECTS) -o preprocr
preprocr.o: DnaAlignment.cpp $(DNAHEADERS)
	$(CC) -c $(COPT) -DPREPROCR DnaAlignment.cpp -o preprocr.o

hapsemblr: hapsemblr.o $(HAPOBJECTS)
	$(CC) $(COPT) hapsemblr.o $(HAPOBJECTS) -o hapsemblr
hapsemblr.o: HapsemblerCore.cpp $(HAPHEADERS)
	$(CC) -c $(COPT) HapsemblerCore.cpp -o hapsemblr.o

scarpa: scarpa.o $(SCAOBJECTS)
	$(CC) $(COPT) scarpa.o $(SCAOBJECTS) liblpsolve55.a -o scarpa -lm -ldl
scarpa.o: ScafferMain.cpp $(SCAHEADERS)
	$(CC) -c $(COPT) ScafferMain.cpp -o scarpa.o

%.o: %.cpp $(DNAHEADERS) $(HAPHEADERS) $(SCAHEADERS)
	$(CC) -c $(COPT) $< -o $@

clean:
	rm -f *.o ./bin/overlappr ./bin/encore ./bin/consensr ./bin/preprocr ./bin/hapsemblr ./bin/scarpa

