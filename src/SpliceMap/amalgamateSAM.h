

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <time.h>
#include <sys/time.h>
#include "SpliceMap_utils_QuasR.h"

struct posline_t{
	int pos;
	int file_index;
	std::string data;
};

struct comparison_t {
	bool operator() (posline_t i,posline_t j) { return (i.pos<j.pos);}
} compare_posline;

struct comparison_reverse_t {
	bool operator() (posline_t i,posline_t j) { return (i.pos>j.pos);}
} compare_posline_reverse;

void print_usage_and_exit();
void update_jun_coverage(map<string,jundict_ij_t> &jundict_chr, string line, bool multi_mapped);
void addnNR(nNR_t &store, vector<uint_fast32_t> exon_start_list, vector<uint_fast32_t> exon_end_list, int num_clip);
bool compareIntVector(vector<uint_fast32_t> &a, vector<uint_fast32_t> &b);


