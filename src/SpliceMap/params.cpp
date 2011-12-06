#include "params.h"



const std::string VERSION = "3.3.5.2 (55)";
const std::string WEBSITE = "http://www.stanford.edu/group/wonglab/SpliceMap/";


const coord_t read_len[2] = {coord_t(1,25),coord_t(26,50)};
const std::string read_len_s[2] = {"1-25","26-50"};




const std::string debug_path = "debug_logs/";
const std::string reference_filename = "ref_list";

const std::string good_dict_list_filename = "gooddict_list";


const std::string map_filename = "25mers.map";




const int_fast32_t read_halflen = 25;
const int_fast32_t read_length = 50;  




const int Iblank = 7;
const int Iextend = 6;
const int Iexonic = 5;
const int Ijunpositive = 1;   
const int Ijunnegative = 2;
const int Mjunpositive = 3;   
const int Mjunnegative = 4;

const bool ORIGINAL = true; 


