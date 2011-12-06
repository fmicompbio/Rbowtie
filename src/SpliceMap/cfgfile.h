

#include <map>
#include <string>
#include <set>
#include <vector>
#include <fstream>
#include <sstream>
#include "params.h"
#include "SpliceMap_utils_QuasR.h"
using namespace std;



class cfgfile{
public:
	cfgfile(string file_name);
	string getVal(string key);
	vector<string> getList(string key);
	
private:
	map<string,vector<string> > data; 
};

