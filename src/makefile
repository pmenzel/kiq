CXX=g++
CXXFLAGS=-O3 -pthread -std=c++11 -DNDEBUG
LDLIBS=-lpthread -lz
INCLUDES=-isystem ../include/
INCLUDES_SRA=-I$(NCBI_DIR)/include

LDLIBS_SRA=-L$(NCBI_DIR)/lib64 \
					 -lncbi-ngs-c++ \
					 -lngs-c++ \
					 -lncbi-vdb \
					 -lpthread -lz -ldl \
					 -Wl,-rpath,$(NCBI_DIR)/lib64

all: makefile kiq testread testwrite convert1to2
	mkdir -p ../bin && cp kiq ../bin/

sra: CXXFLAGS:=$(CXXFLAGS) -D KIQ_SRA
sra: makefile main.o kdump.o kindex.o kquery.o kdb.o kmodify.o util.o ReadItem.o CountThread.o ksra.o
	$(CXX) $(LDFLAGS) -o kiq kdump.o main.o kindex.o kdb.o kquery.o util.o kmodify.o ReadItem.o CountThread.o ksra.o $(LDLIBS_SRA)
	mkdir -p ../bin && cp kiq ../bin/

kiq: makefile main.o kdump.o kindex.o kquery.o kdb.o kmodify.o util.o ReadItem.o CountThread.o
	$(CXX) $(LDFLAGS) -o kiq kdump.o main.o kindex.o kdb.o kquery.o kmodify.o util.o ReadItem.o CountThread.o $(LDLIBS)

testread: util.o testread.o
	$(CXX) $(LDFLAGS) -o testread util.o testread.o $(LDLIBS)
testwrite: util.o testwrite.o
	$(CXX) $(LDFLAGS) -o testwrite util.o testwrite.o $(LDLIBS)
convert1to2: util.o convert1to2.o
	$(CXX) $(LDFLAGS) -o convert1to2 util.o convert1to2.o $(LDLIBS)

ksra.o: ksra.cpp version.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(INCLUDES_SRA) -c -o ksra.o ksra.cpp

%.o : %.cpp version.hpp ../include/ProducerConsumerQueue/ProducerConsumerQueue.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -f -v kiq ../bin/*
	find . -name "*.o" -delete

static: LDFLAGS = -static
static: LDLIBS = -Wl,--whole-archive -lpthread -lz -Wl,--no-whole-archive
static: all

debug: CXXFLAGS = -O3 -pthread -std=c++11 -ggdb3 -Wall -Wpedantic -Wextra -Wconversion -Wno-unused-function -fno-omit-frame-pointer
debug: all

.PHONY: clean debug static

