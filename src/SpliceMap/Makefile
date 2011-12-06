#R_HOME=`R RHOME`
#if test -z "${R_HOME}"
#then
#  echo "could not determine R_HOME"
#  exit 1
#fi
#CXX=`"${R_HOME}/bin/R" CMD config CXX`
#CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CXX=`R CMD config CXX`
CXXFLAGS=`R CMD config CXXFLAGS`

all: SpliceMap

clean: 
	rm *.o SpliceMap runSpliceMap_QuasR sortsam amalgamateSAM

SpliceMap: runSpliceMap_QuasR sortsam amalgamateSAM params.o cfgfile.o SpliceMap_utils_QuasR.o main.cpp main.h
	$(CXX) $(CXXFLAGS) -o SpliceMap main.cpp SpliceMap_utils_QuasR.o params.o cfgfile.o

runSpliceMap_QuasR: SpliceMap_utils_QuasR.o params.o cfgfile.o runSpliceMap_QuasR.cpp runSpliceMap_QuasR.h 
	$(CXX) $(CXXFLAGS) -o runSpliceMap_QuasR SpliceMap_utils_QuasR.o cfgfile.o runSpliceMap_QuasR.cpp params.o

sortsam: sortsam.cpp sortsam.h SpliceMap_utils_QuasR.o params.o
	$(CXX) $(CXXFLAGS) -o sortsam sortsam.cpp SpliceMap_utils_QuasR.o params.o

amalgamateSAM: amalgamateSAM.cpp amalgamateSAM.h SpliceMap_utils_QuasR.o params.o
	$(CXX) $(CXXFLAGS) -o amalgamateSAM amalgamateSAM.cpp SpliceMap_utils_QuasR.o params.o

SpliceMap_utils_QuasR.o: SpliceMap_utils_QuasR.cpp SpliceMap_utils_QuasR.h
	$(CXX) $(CXXFLAGS) -c SpliceMap_utils_QuasR.cpp

params.o: params.cpp params.h
	$(CXX) $(CXXFLAGS) -c params.cpp

cfgfile.o: cfgfile.cpp cfgfile.h
	$(CXX) $(CXXFLAGS) -c cfgfile.cpp

