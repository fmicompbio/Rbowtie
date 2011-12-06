

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include <map>
#include <set>
#include <exception>
#include <stdexcept>
#include <algorithm>

#include "SpliceMap_utils_QuasR.h"
#include "cfgfile.h"

using namespace std;











const uint_fast32_t MAX_HASH = 1048576;

typedef vector<uint_fast32_t>* seq_hash;





struct line_t {
	string chr_name;
	bool unique; 
	int_fast32_t errors;   
	uint_fast32_t coor;
	bool direction;  
};


struct mix_t {
	string full_line;  
	string seed_line;
	bool unique;
	int_fast32_t errors;
	string chr_name;
	uint_fast32_t coor;
	bool direction; 
};



void print_usage_and_exit();
vector<string> make_onesegment_namelist( string &reads_filename,int fullread_length);
string get_filename(const string& path);
inline void DelLocalHalfMappableHits(list<line_t> &line_list, uint_fast32_t distance);
inline void FindExonicRead(good_vec_t& good, list<line_t>  line_list[2]);
inline uint_fast32_t cal_half_distance(uint_fast32_t length, uint_fast32_t intron_curve[3]);
inline pair<vector<uint_fast32_t>::iterator, vector<uint_fast32_t>::iterator> checkindexlistFwd(uint_fast32_t ref,vector<uint_fast32_t> &indexlist,uint_fast32_t range1, uint_fast32_t range2);
inline pair<vector<uint_fast32_t>::iterator, vector<uint_fast32_t>::iterator> checkindexlistRes(uint_fast32_t ref,vector<uint_fast32_t> &indexlist,uint_fast32_t range1, uint_fast32_t range2);
inline bool isgoodtequal(good_t &first, good_t &second);

inline bool checkPairInfo(list<good_t> &left_jun_list, list<good_t> &right_jun_list, int distance);
vector<good_t> combine_result_list(vector<good_t> &l1, vector<good_t> &l2);
good_vec_t checkgoodlists(good_vec_t &good1, good_vec_t &good2);
inline int get_base_val(char c);
inline int base2int(string &s);


vector< list<good_t> > check_singleread_group_result(vector<vector<good_t> > &singleread_group_result, pair<vector<coord_t>,vector<coord_t> > &range_list);


bool check_segment_overlap(good_t &junctiona, good_t &junctionb,coord_t &lena, coord_t &lenb);
good_t good_seq(good_t &junction, int local_a, int local_b);
void recursive_record(vector<vector<good_t> > &singleread_group_result, map<int,map<int,list<int> > > &result
					  , int ref_group_index, int good_index, list<good_t> &one_output,vector< list<good_t> > &all_output );
void InsertGoodExtend(good_vec_t &good_list,good_vec_t &good_exonic_list,good_t contents);
bool compare_line_t(line_t &first, line_t &second);


string count_mismatch(start_end_nano_t &alignment, string &read, string &genome,int_fast32_t num_read_mismatch,int_fast32_t max_clip_allowed);
	

