

#ifndef _PARAMS_H
#define  _PARAMS_H

#define __STDC_LIMIT_MACROS

#include <iostream>
#include <string>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <climits>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdint.h>






extern const int Iblank;
extern const int Iextend;
extern const int Iexonic;
extern const int Ijunpositive;   
extern const int Ijunnegative;
extern const int Mjunpositive;   
extern const int Mjunnegative;

extern const int_fast32_t read_length;  
extern const std::string read_len_s[2];
extern const int_fast32_t read_halflen;







extern const std::string debug_path; 
extern const std::string reference_filename; 

extern const std::string good_dict_list_filename; 


extern const std::string map_filename; 

typedef std::pair<int, int> coord_t;

extern const coord_t read_len[2];

extern const std::string VERSION;
extern const std::string WEBSITE;

extern const bool ORIGINAL; 

struct nNR_inside_t {
	std::vector<uint_fast32_t> chr_start;  
	std::vector<uint_fast32_t> chr_end;
	int num_clip; 
};

typedef std::vector<nNR_inside_t> nNR_t;  



 
 

struct jun_store_ij { 
	
	bool direction; 
	int i;  
	int j;  
	
	
	
	nNR_t nNR; 
};

typedef std::map<int,jun_store_ij > jundict_b_ij_t;  
typedef std::map<int,jundict_b_ij_t> jundict_ij_t;  


typedef std::map<int,int> endpos_t; 
typedef std::map<int,endpos_t> beginpos_t;
typedef std::map<std::string,beginpos_t> chrdict_t;





struct good_t {
	int_fast32_t a; 
	int_fast32_t b; 
	int_fast32_t c; 
	int_fast32_t d; 
};


typedef std::vector<good_t> good_vec_t;


typedef std::pair<std::vector<std::list<good_t> >,std::pair<std::vector<coord_t>,std::vector<coord_t> >  >  good_seg_vec_t;


typedef std::map< int,std::vector<good_t> > good_list_t;

typedef std::map< int,std::pair<std::vector<std::list<good_t> >,int > > good_seg_list_t;





typedef std::vector<std::pair<std::string,std::list<good_t> > > good_chr_vec_t;
typedef std::map< std::pair<int,int>,good_chr_vec_t > good_chr_list_t;


struct start_end_nano_t { 
	float coverage; 
	
	bool strand; 
	bool direction; 
	
	std::vector<int> chr_start;  
	std::vector<int> chr_end;
	std::vector<short> read_start; 
	std::vector<short> read_end;
};

typedef std::list< start_end_nano_t > start_end_vec_t;
typedef std::map< std::pair<int,int>, start_end_vec_t> start_end_list_t;

typedef std::map<int,float> coverage_t;

struct coverage_entry_t {
	
	float coverage; 
	std::vector<int> chr_start;  
	std::vector<int> chr_end;
};

struct reference_t {
	std::string file_name;
	std::string file_path;
	uint_fast64_t file_index_start; 
	uint_fast64_t file_index_end; 
};

#endif


