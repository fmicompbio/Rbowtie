CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

all: SpliceMap

clean: 
	rm *.o SpliceMap runSpliceMap_QuasR sortsam amalgamateSAM getSpliceMapUnmapped fuseReorder

SpliceMap: runSpliceMap_QuasR sortsam amalgamateSAM getSpliceMapUnmapped fuseReorder params.o cfgfile.o SpliceMap_utils_QuasR.o main.cpp main.h
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o SpliceMap main.cpp SpliceMap_utils_QuasR.o params.o cfgfile.o

runSpliceMap_QuasR: SpliceMap_utils_QuasR.o params.o cfgfile.o runSpliceMap_QuasR.cpp runSpliceMap_QuasR.h 
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o runSpliceMap_QuasR SpliceMap_utils_QuasR.o cfgfile.o runSpliceMap_QuasR.cpp params.o

sortsam: sortsam.cpp sortsam.h SpliceMap_utils_QuasR.o params.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o sortsam sortsam.cpp SpliceMap_utils_QuasR.o params.o

amalgamateSAM: amalgamateSAM.cpp amalgamateSAM.h SpliceMap_utils_QuasR.o params.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o amalgamateSAM amalgamateSAM.cpp SpliceMap_utils_QuasR.o params.o

getSpliceMapUnmapped: getSpliceMapUnmapped.cpp
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o getSpliceMapUnmapped getSpliceMapUnmapped.cpp

fuseReorder: fuseReorder.cpp
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o fuseReorder fuseReorder.cpp

SpliceMap_utils_QuasR.o: SpliceMap_utils_QuasR.cpp SpliceMap_utils_QuasR.h
	$(CXX) $(CXXFLAGS) -c SpliceMap_utils_QuasR.cpp

params.o: params.cpp params.h
	$(CXX) $(CXXFLAGS) -c params.cpp

cfgfile.o: cfgfile.cpp cfgfile.h
	$(CXX) $(CXXFLAGS) -c cfgfile.cpp

