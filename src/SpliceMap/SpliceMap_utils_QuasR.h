


#ifndef _SPLICEMAP_UTILS_H
#define  _SPLICEMAP_UTILS_H

#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <limits>
#include <climits>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "params.h"
//#include "ostools.h"


using namespace std;

bool make_DNA_upper(string &str);

string LongintToStr( uint_fast64_t x );
string IntToStr(int x);  
string DoubleToStr( double x );
void trim2(string& str);
vector<string> split(const string &s, char delim);
vector<string> split(const string &s);
list<string> split_list(const string &s, char delim);
list<string> split_list(const string &s);
double diffclock(struct timeval &start_tv,struct timeval &tv);
vector<coord_t> designsuffix(int read_length);
void rtrim(string& str);
void ltrim(string& str);


void print_jun_dict_ij(jundict_ij_t &jundict, string bedname, string chr_name,string track_name);
string getDirection(int direc_type);
string directranslate(bool direc_type);
int rangedict(nNR_t &dict);
void load_bed_ref_file(ifstream &primary_file, chrdict_t &primary_jundict, int primary_file_type);
void compleseq(string &seq);
bool sameSign(int x, int y);
void print_good(ostream& out, good_t& good);



coord_t junction2boundary(good_t junction);


int read_full_reads_file(ifstream &in_file, vector<string> &full_reads);
void print_int_vector(ofstream &out, vector<int> &vec,char delim);

void phred642phred33(string &quals);
void solexa2phred33(string &quals);
vector<pair<uint_fast32_t,uint_fast32_t> > SAM2updown(string line);
vector<pair<uint_fast32_t,uint_fast32_t> > cigar2updown(uint_fast32_t start_loc, string cigar);
bool read_reference_map(string chromosome_path,map<string,reference_t> &ref_map, vector<string> chr_file_list);
void add_good_t(start_end_nano_t &start_end, list<good_t> good_seg_list, vector<coord_t> suffix);
pair<pair<int,int>,string> roll_cigar(start_end_nano_t &start_end, int full_read_len, bool cufflinks);

#endif

