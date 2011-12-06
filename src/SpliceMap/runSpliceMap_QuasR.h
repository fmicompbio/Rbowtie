#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <sstream>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <limits>
#include <climits>
#include <cstdlib>
#include <unistd.h>
#include "SpliceMap_utils_QuasR.h"
#include "cfgfile.h"

using namespace std;

typedef pair<int, int> coord_t;

struct sam_mismatch_t{
	int num_mismatch;
	std::string data;
};


struct dot_t_t {
	int8_t mismatch_dir; 
	string chr_name;
	uint_fast32_t location;
};



inline void print_usage_and_exit();
pair<string,string> get_path_and_filename(const string& path);

inline void processbowtie(string &input_filename, string &output_filename);
inline int sam_char2int(char c);
inline bool sam_is_mapped(unsigned int flag);
inline string sam_get_direction(unsigned int flag);

void output_index(string reads_filename, ofstream &map_out);





