#
# Simple Makefile for compiling C++ code 
#

CXX := g++
CXXFLAGS := -O3 -m64 -Wall -std=c++11
CXXLINKFLAGS := -lm 
CXXINC := -Iinclude

OBJECTDIR := obj
BINDIR := bin
SRCDIR := src

OBJECTS := \
	${OBJECTDIR}/DataSource.o \
	${OBJECTDIR}/EM.o \
	${OBJECTDIR}/FileIO.o \
	${OBJECTDIR}/GraphTools.o \
	${OBJECTDIR}/Layer.o \
	${OBJECTDIR}/Network.o \
	${OBJECTDIR}/NetworkFactory.o \
	${OBJECTDIR}/Node.o \
	${OBJECTDIR}/PDStruct.o \
	${OBJECTDIR}/Utils.o \
	


all: directories ${BINDIR}/trainEM ${BINDIR}/MPEsynth

clean: 
	rm -rf $(OBJECTS) ${BINDIR}/trainEM ${BINDIR}/MPEsynth

directories: 
	mkdir -p ${OBJECTDIR}
	mkdir -p ${BINDIR}


#
# binaries
#
${BINDIR}/trainEM: ${SRCDIR}/trainEM.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CXXINC) -Isrc/ ${SRCDIR}/trainEM.cpp -o ${BINDIR}/trainEM $(CXXLINKFLAGS) $(OBJECTS)

${BINDIR}/MPEsynth: ${SRCDIR}/MPEsynth.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CXXINC) -Isrc/ ${SRCDIR}/MPEsynth.cpp -o ${BINDIR}/MPEsynth $(CXXLINKFLAGS) $(OBJECTS)


#
# object files
#
${OBJECTDIR}/DataSource.o:${SRCDIR}/DataSource.hpp ${SRCDIR}/DataSource.cpp
	$(CXX) $(CXXFLAGS) $(CXXINC) -c -Isrc/ ${SRCDIR}/DataSource.cpp -o ${OBJECTDIR}/DataSource.o

${OBJECTDIR}/EM.o:${SRCDIR}/EM.hpp ${SRCDIR}/EM.cpp
	$(CXX) $(CXXFLAGS) $(CXXINC) -c -Isrc/ ${SRCDIR}/EM.cpp -o ${OBJECTDIR}/EM.o

${OBJECTDIR}/FileIO.o:${SRCDIR}/FileIO.hpp ${SRCDIR}/FileIO.cpp
	$(CXX) $(CXXFLAGS) $(CXXINC) -c -Isrc/ ${SRCDIR}/FileIO.cpp -o ${OBJECTDIR}/FileIO.o

${OBJECTDIR}/GraphTools.o:${SRCDIR}/GraphTools.cpp ${SRCDIR}/GraphTools.hpp
	$(CXX) $(CXXFLAGS) $(CXXINC) -c -Isrc/ ${SRCDIR}/GraphTools.cpp -o ${OBJECTDIR}/GraphTools.o

${OBJECTDIR}/Layer.o: ${SRCDIR}/Layer.hpp ${SRCDIR}/Layer.cpp ${SRCDIR}/Node.hpp
	$(CXX) $(CXXFLAGS) $(CXXINC) -c -Isrc/ ${SRCDIR}/Layer.cpp -o ${OBJECTDIR}/Layer.o

${OBJECTDIR}/Network.o:${SRCDIR}/Network.hpp ${SRCDIR}/Network.cpp ${SRCDIR}/Node.hpp ${SRCDIR}/Layer.hpp
	$(CXX) $(CXXFLAGS) $(CXXINC) -c -Isrc/ ${SRCDIR}/Network.cpp -o ${OBJECTDIR}/Network.o

${OBJECTDIR}/NetworkFactory.o:${SRCDIR}/NetworkFactory.hpp ${SRCDIR}/NetworkFactory.cpp 
	$(CXX) $(CXXFLAGS) $(CXXINC) -c -Isrc/ ${SRCDIR}/NetworkFactory.cpp -o ${OBJECTDIR}/NetworkFactory.o

${OBJECTDIR}/Node.o:${SRCDIR}/Node.hpp ${SRCDIR}/Node.cpp ${SRCDIR}/Utils.hpp 
	$(CXX) $(CXXFLAGS) $(CXXINC) -c -Isrc/ ${SRCDIR}/Node.cpp -o ${OBJECTDIR}/Node.o

${OBJECTDIR}/PDStruct.o:${SRCDIR}/PDStruct.hpp ${SRCDIR}/PDStruct.cpp  
	$(CXX) $(CXXFLAGS) $(CXXINC) -c -Isrc/ ${SRCDIR}/PDStruct.cpp -o ${OBJECTDIR}/PDStruct.o

${OBJECTDIR}/Utils.o: ${SRCDIR}/Utils.hpp ${SRCDIR}/Utils.cpp
	$(CXX) $(CXXFLAGS) $(CXXINC) -c -Isrc/ ${SRCDIR}/Utils.cpp -o ${OBJECTDIR}/Utils.o
