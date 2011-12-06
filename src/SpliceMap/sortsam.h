


#include "params.h"
#include "SpliceMap_utils_QuasR.h"
#include <algorithm>
#include <queue>

struct posline_t{
	uint_fast32_t pos;
	std::string data;
};

struct comparison_t {
	bool operator() (posline_t i,posline_t j) { return (i.pos<j.pos);}
} compare_posline;

struct comparison_reverse_t {
	bool operator() (posline_t i,posline_t j) { return (i.pos>j.pos);}
} compare_posline_reverse;



void print_usage_and_exit();


