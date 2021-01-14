CC = g++
LD = $(CC)
CFLAGS = -c -Wall -Os
#IFLAGS =  -I./ -I$(BDXRECO_ROOT)/src/libraries -I$(BDXRECO_ROOT)/src/external/jana_0.7.7p1/src
ROOTINC := $(shell root-config --cflags)
ROOTLIB := $(shell root-config --glibs)
#LIBS = -L$(BDXRECO_ROOT)/lib -lbdxReco -lbdxRecoExt -lJANA -lProof
#LFLAGS = -Wl


SRC = $(wildcard *.C)

TARGET = main

OBJECTS = $(patsubst %.C, %.o, $(wildcard *.C))

all:   dict $(TARGET)

dict: Cluster_Ecal_Selector.h
	@echo "Generating dictionary $@..."
	@rootcint -v -f Cluster_Ecal_Selector_Dict.C -c -p Cluster_Ecal_Selector.h Cluster_Ecal_Selector_LinkDef.h

$(TARGET): $(OBJECTS)
	$(LD) -shared -fPIC -o libCluster_Ecal_Selector.so $(ROOTLIB) $^
	$(LD) -o $@ $^ $(ROOTLIB)

%.o: %.C
	$(CC) -fPIC -g $(CFLAGS) $(ROOTINC) $^ -o $@

clean:
	rm $(TARGET) $(OBJECTS)
	rm libCluster_Ecal_Selector.so
