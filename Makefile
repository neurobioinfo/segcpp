CC=g++
#CC=icc
CDEBUG_FLAGS=-g
CFLAGS=-c -O2 -Wall $(CDEBUG_FLAGS) -I /usr/include/boost 
#-I /usr/include/boost
#LDFLAGS=-lboost_regex -lboost_program_options -lboost_thread
LDFLAGS=-lboost_program_options -lboost_system
#SOURCES=segregation.cpp parsers.cpp main.cpp Variant.cpp  VariantCallerParser.cpp Chromosome.cpp Sample.cpp Family.cpp VariantInSample.cpp
#SOURCES=CPedFile.cpp
SOURCES=CChromosome.cpp  CGene.cpp  CFamily.cpp  CSample.cpp  CVariant.cpp  CVariantInSample.cpp  UPedFileParser.cpp UStringUtils.cpp  CRunSettings.cpp CVCFFile.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=segregation
STATIC=segregation.static
GFF2IBS=gff2ibs

all: $(SOURCES) $(EXECUTABLE)
#all: $(SOURCES) $(STATIC)
	
$(EXECUTABLE): $(OBJECTS) PAnalysis.o
	$(CC) $(OBJECTS) PAnalysis.o -o $@ $(LDFLAGS)

$(GFF2IBS): $(OBJECTS) GFF2IBS.o
	$(CC) $(LDFLAGS) $(OBJECTS) GFF2IBS.o -o $@

vcfparser: $(OBJECTS) PVCFParser.o
	$(CC) $(LDFLAGS) $(OBJECTS) PVCFParser.o -o $@

static:	$(SOURCES) $(STATIC)

$(STATIC): $(OBJECTS) 
	$(CC) $(LDFLAGS) -static $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS) PAnalysis.o
