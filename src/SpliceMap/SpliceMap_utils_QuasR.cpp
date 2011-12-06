

#include "SpliceMap_utils_QuasR.h"



char solexa2phred32_table[68] = 
{1,1,2,2,3,3,4,4,5,5,6,7,
8,9,10,10,11,12,13,14,15,16,17,18,
19,20,21,22,23,24,25,26,27,28,29,30,
31,32,33,34,35,36,37,38,39,40,41,42,
43,44,45,46,47,48,49,50,51,52,53,54,
55,56,57,58,59,60,61,62};

bool make_DNA_upper(string &str)
{
	
	string::iterator it = str.begin();
	while (it != str.end()) {
		switch (*it) {
			case 'a':
				*it = 'A';
				break;
			case 'c':
				*it = 'C';
				break;
			case 't':
				*it = 'T';
				break;
			case 'g':
				*it = 'G';
				break;
			case 'n':
				*it = 'N';
				break;
			case 'A':
				break;
			case 'C':
				break;
			case 'T':
				break;
			case 'G':
				break;
			case 'N':
				break;
			default:
				*it = 'N';
				break;
		}
		it++;
	}
	
	return true;
}


string LongintToStr( uint_fast64_t x ) 
{
	stringstream o; 
	o << x;
	return o.str();
}

string IntToStr( int x )
{
	stringstream o; 
	o << x;
	return o.str();
}

string DoubleToStr( double x )
{
	stringstream o; 
	o << x;
	return o.str();
}


void ltrim(string& str) 
{
	string::size_type pos = 0;
	while (pos < str.size() && (isspace(str[pos]))) pos++;
	str.erase(0, pos);
}
void rtrim(string& str) 
{
	string::size_type pos = str.size();
	while (pos > 0 && (isspace(str[pos - 1]))) pos--;
	str.erase(pos);
}
void trim2(string& str) 
{
	ltrim(str);
	rtrim(str);
}


vector<string> split(const string &s, char delim) 
{
    vector<string> elems;
	stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


vector<string> split(const string &s) 
{
    vector<string> elems;
	stringstream ss(s);
    string item;
    while(ss >> item) {
        elems.push_back(item);
    }
    return elems;
}

list<string> split_list(const string &s, char delim) 
{
    list<string> elems;
	stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


list<string> split_list(const string &s) 
{
    list<string> elems;
	stringstream ss(s);
    string item;
    while(ss >> item) {
        elems.push_back(item);
    }
    return elems;
}

double diffclock(struct timeval &start_tv,struct timeval &tv)
{
	double elapsed;
	
	elapsed = (tv.tv_sec - start_tv.tv_sec) + (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
	
	return elapsed;
}




vector<coord_t> designsuffix(int full_read_length)
{
	vector<coord_t> result;
	int Nsegment = 0;
	
	if (full_read_length  == 50) {
		Nsegment = 1;
	}else if (full_read_length <= 59) {
		Nsegment = 2;
	}else {
		Nsegment = (full_read_length - 60)/30 + 2;  
	}

	result.push_back(coord_t(1,50));
	int i = 1;
	int start_pt = 1;
	int end_pt = 50;
	
	while (i < Nsegment) {
		int distance = (full_read_length - 50)/(Nsegment - i);
		start_pt = distance + start_pt;
		end_pt = distance + end_pt;
		result.push_back(coord_t(start_pt,end_pt));
		full_read_length = full_read_length - distance;
		i++;
	}

	return result;
}


void compleseq(string &seq)
{
	string temp = seq;
	string::iterator it = seq.begin();
	string::reverse_iterator rit = temp.rbegin();
	while (it != seq.end()) {
		switch (*rit) {
			case 'A':
				*it = 'T';
				break;
			case 'C':
				*it = 'G';
				break;
			case 'G':
				*it = 'C';
				break;
			case 'T':
				*it = 'A';
				break;
			case 'a':
				*it = 't';
				break;
			case 'c':
				*it = 'g';
				break;
			case 'g':
				*it = 'c';
				break;
			case 't':
				*it = 'a';
				break;
			default:
				*it = *rit;
				break;
		}
		it++;
		rit++;
	}
}

void print_good(ostream& out, good_t& good)
{
	out << good.a << "," << good.b << "," << good.c << "," << good.d;
}

bool sameSign(int x, int y)
{
    return (x >= 0) ^ (y < 0);
}

int read_full_reads_file(ifstream &in_file, vector<string> &full_reads)
{
	int len = 0;
	string line = "";
	
	
	
	while (!in_file.eof()) {
		
		
		
		getline(in_file, line);
		trim2(line);
		if(line.length() == 0)
			continue;
		full_reads.push_back(line);
		
		len++;
	}
	
	return len;
}

void print_int_vector(ofstream &out, vector<int> &vec,char delim)
{
	
	
	
	vector<int>::iterator vec_it = vec.begin();
	out << *vec_it;
	
	vec_it++;
	while (vec_it != vec.end()) {
		out << delim << *vec_it;
		
		vec_it++;
	}
}












void print_jun_dict_ij(jundict_ij_t &jundict, string bedname, string chr_name,string track_name)
{
	
	ofstream junction;
	
	if (track_name.length() == 0) {
		junction.open(bedname.c_str(), ios::app);
	}else{
		junction.open(bedname.c_str(), ios::out);
		junction << "track\tname=junctions\tdescription=\"" << track_name << "\"" << endl;
	}
	
	int number = 1;
	
	jundict_ij_t::iterator leftpos = jundict.begin();
	
	
	
	while (leftpos != jundict.end()) {
		
		jundict_b_ij_t::iterator rightpos = (leftpos->second).begin();
		while (rightpos != (leftpos->second).end()) {
			int count = 1; 
			string str_leftpos = IntToStr(leftpos->first - count);
			string str_rightpos = IntToStr(rightpos->first + count);
			string str_count = IntToStr(count);
			string str_length = IntToStr(rightpos->first - leftpos->first + count);
			string str_right_length = IntToStr(rangedict((rightpos->second).nNR)) 
			+ "_" + IntToStr((int)(rightpos->second).nNR.size());
			string str_direc = directranslate((rightpos->second).direction); 
			
			int i_val = (rightpos->second).i;
			int j_val = (rightpos->second).j;
			int nR = i_val + j_val;
			
			string ij = IntToStr(i_val) + "/" + IntToStr(j_val);
			
			junction << chr_name << '\t';
			junction << str_leftpos << '\t';
			junction << str_rightpos << '\t';
			
			junction << "(" << nR << ")";
			junction << "[" << str_right_length << "]";
			junction << "(" << ij << ")\t";
			junction << str_count << '\t';
			junction << str_direc << '\t';
			junction << str_leftpos << '\t';
			junction << str_rightpos << '\t';
			junction << "255,0,0\t2\t";  
			junction << str_count << "," << str_count << '\t';
			junction << "0," << str_length << "\n";
			
			number++;
			rightpos++;
		}
		
		leftpos++;
	}
	junction << flush;
	
	junction.close();
}


string getDirection(int direc_type)
{
	if (direc_type > 0) {
		return "F";
	}if (direc_type < 0) {
		return "R";
	}
	return "error"; 
}


string directranslate(bool direc_type)
{
	direc_type = abs(direc_type);
	if (direc_type) { 
		return "+";
	}else { 
		return "-";
	}
	return "#"; 
}


int rangedict(nNR_t &dict)
{
	int smallest_left = 0;
	int biggest_left = 0;
	
	nNR_t::iterator it = dict.begin();
	
	smallest_left = it->chr_end.back();
	biggest_left = it->chr_end.back();
	
	it++;
	
	while (it != dict.end()) {
		int temp = it->chr_end.back();
		if (temp < smallest_left) {
			smallest_left = temp;
		}
		
		if (temp > biggest_left) {
			biggest_left = temp;
		}
		
		it++;
	}
	
	
	return biggest_left - smallest_left + 1;
}


void load_bed_ref_file(ifstream &primary_file, chrdict_t &primary_jundict, int primary_file_type)
{
	string line;
	string temp;
	
	
	while (!primary_file.eof()) {
		getline(primary_file, line);
		trim2(line);
		
		if (line.length() == 0){
			continue;
		}
		
		if(primary_file_type == 0){ 
			
			vector<string> refline_list = split(line);
			
			if (atoi(refline_list[8].c_str()) == 1) {
				continue; 
			}
			
			list<string> temp_list = split_list(refline_list[9], ',');
			list<int> exon_start_list;
			for (list<string>::iterator it = temp_list.begin(); it!=temp_list.end(); it++) {
				exon_start_list.push_back(atoi(it->c_str()));
			}
			exon_start_list.pop_front(); 
			
			temp_list = split_list(refline_list[10], ',');
			list<int> exon_end_list;
			for (list<string>::iterator it = temp_list.begin(); it!=temp_list.end(); it++) {
				exon_end_list.push_back(atoi(it->c_str()));
			}
			exon_end_list.pop_back(); 
			
			string line_chr = refline_list[2]; 
			chrdict_t::iterator primarychr_it = primary_jundict.find(line_chr);
			
			
			list<int>::iterator start_it = exon_end_list.begin();
			list<int>::iterator end_it = exon_start_list.begin(); 
			
			if(primarychr_it == primary_jundict.end()){
				
				endpos_t end_contents;
				end_contents.insert(pair<int,int>(*end_it,1));
				beginpos_t begin_contents;
				begin_contents.insert(pair<int,endpos_t>(*start_it,end_contents));
				pair<chrdict_t::iterator,bool> out = primary_jundict.insert(pair<string,beginpos_t>(line_chr,begin_contents));
				primarychr_it = out.first;
				
				
				
				start_it++;
				end_it++;
			}
			
			while (start_it != exon_end_list.end()) {
				
				
				beginpos_t::iterator begin_it = (primarychr_it->second).find(*start_it);
				
				if (begin_it != (primarychr_it->second).end()) {
					endpos_t::iterator finish_it = (begin_it->second).find(*end_it);
					
					if (finish_it != (begin_it->second).end()) {
						finish_it->second++;
					}else {
						(begin_it->second).insert(pair<int,int>(*end_it,1));
						
					}
					
				}else {
					endpos_t contents;
					contents.insert(pair<int,int>(*end_it,1));
					
					(primarychr_it->second).insert(pair<int,endpos_t>(*start_it,contents));
					
					
				}
				
				start_it++;
				end_it++;
			}
			
		}else if(primary_file_type == 1){
			
			
			if (line.substr(0, 5).compare("track") == 0) {
				continue;
			}
			
			vector<string> bedline_list = split(line);
			
			vector<string> thickness = split(bedline_list[10], ',');
			int leftpos = atoi(bedline_list[1].c_str()) + atoi(thickness[0].c_str());
			int rightpos = atoi(bedline_list[2].c_str()) - atoi(thickness[1].c_str());
			
			string line_chr = bedline_list[0]; 
			chrdict_t::iterator primarychr_it = primary_jundict.find(line_chr);
			
			
			if(primarychr_it == primary_jundict.end()){
				
				endpos_t end_contents;
				end_contents.insert(pair<int,int>(rightpos,1));
				beginpos_t begin_contents;
				begin_contents.insert(pair<int,endpos_t>(leftpos,end_contents));
				pair<chrdict_t::iterator,bool> out = primary_jundict.insert(pair<string,beginpos_t>(line_chr,begin_contents));
				primarychr_it = out.first;
				
				
				
			}else {
				beginpos_t::iterator begin_it = (primarychr_it->second).find(leftpos);
				
				if (begin_it != (primarychr_it->second).end()) {
					endpos_t::iterator finish_it = (begin_it->second).find(rightpos);
					
					if (finish_it != (begin_it->second).end()) {
						finish_it->second++;
					}else {
						(begin_it->second).insert(pair<int,int>(rightpos,1));
						
					}
					
				}else {
					endpos_t contents;
					contents.insert(pair<int,int>(rightpos,1));
					
					(primarychr_it->second).insert(pair<int,endpos_t>(leftpos,contents));
				}
			}
			
		}
		
		
	}
	
}







pair<pair<int,int>,string> roll_cigar(start_end_nano_t &start_end, int full_read_len, bool cufflinks)
{
	bool DNA_mode = false;
	
	string cigar;
	pair<int,int> clipped;
	clipped.first = 0;
	clipped.second = 0;
	
	vector<short>::iterator start_it     = start_end.read_start.begin();
	vector<short>::iterator end_it       = start_end.read_end.begin();
	vector<int>::iterator chr_start_it = start_end.chr_start.begin();
	vector<int>::iterator chr_end_it   = start_end.chr_end.begin();
	
	int num_sec = (int)start_end.read_start.size();
	
	
	if (*start_it != 1) {
		if (!cufflinks){
			cigar += (IntToStr(*start_it - 1) + "S");
		}
		clipped.first = *start_it - 1;
	}
	
	if (num_sec > 1) {
		while (start_it != --(start_end.read_start.end())) {
			vector<short>::iterator next_start_it = start_it + 1;
			vector<int>::iterator next_chr_start_it = chr_start_it + 1;
			short M_end = 0;
			if (*next_chr_start_it - *chr_end_it - 1 < 0) {
				M_end = *end_it + (*next_chr_start_it - *chr_end_it - 1);
			}else {
				M_end = *end_it;
			}
			
			cigar += (IntToStr(M_end - *start_it + 1) + "M");
			start_it = next_start_it;
			chr_start_it = next_chr_start_it;
			if (!DNA_mode) {
				cigar += (IntToStr(*chr_start_it - *chr_end_it - 1) + "N");
			}else{
				if(*chr_start_it - *chr_end_it - 1 > 0)
					cigar += (IntToStr(*chr_start_it - *chr_end_it - 1) + "D");
				else {
					cigar += (IntToStr(*chr_end_it + 1 - *chr_start_it) + "I");
				}
			}
			
			end_it++;
			chr_end_it++;			
		}
	}
	
	
	cigar += (IntToStr(*end_it - *start_it + 1) + "M");
	if (*end_it != full_read_len) {
		if (!cufflinks){
			cigar += (IntToStr(full_read_len - *end_it) + "S");
		}
		clipped.second = full_read_len - *end_it;
	}
	
	return pair<pair<int,int>,string>(clipped,cigar);
}




void add_good_t(start_end_nano_t &start_end, list<good_t> good_seg_list, vector<coord_t> suffix)
{
	bool strand = false; 
	int type = 0;
	int suffix_idx = 0;
	
	int back_pos = 0;  
	int start_pos = INT_MAX; 
	
	int read_back_pos = 0;  
	int read_start_pos = INT_MAX; 
	
	
	if (good_seg_list.front().c > 0) { 
		good_seg_list.reverse();
		start_end.direction = true;
	}else { 
		reverse(suffix.begin(),suffix.end());
		start_end.direction = false;
	}
	
	
	list<good_t>::iterator good_it = good_seg_list.begin();       
	list<good_t>::iterator good_it_back = good_seg_list.end();  
	good_it_back--;
	
	
	
	
	
	type = abs(good_it->c);
	start_end.direction = (good_it->c)>0; 
	
	
	
	
	
	if (type == Iextend ) {
		start_end.chr_start.push_back(good_it->a);
		if (good_it->d > 0) {
			start_end.read_start.push_back(suffix[suffix_idx].first);
		}else {
			start_end.read_start.push_back(suffix[suffix_idx].second + good_it->d + 1); 
		}
	}else if(type == Iexonic){
		
		int temp_back_pos = good_it->b;
		if (temp_back_pos > back_pos) {
			back_pos = temp_back_pos;
			read_back_pos = suffix[suffix_idx].second;
		}
		
		start_end.chr_start.push_back(good_it->a);
		start_end.read_start.push_back(suffix[suffix_idx].first);
	}else { 
		
		int temp_back_pos = good_it->b + good_it->d;
		if (temp_back_pos > back_pos) {
			back_pos = temp_back_pos;
			read_back_pos = suffix[suffix_idx].second;
		}
		
		int jun_start = good_it->a - (read_length - good_it->d) + 1;
		
		start_end.chr_start.push_back(jun_start);
		start_end.read_start.push_back(suffix[suffix_idx].first);
		
		start_end.chr_end.push_back(good_it->a);
		start_end.read_end.push_back(suffix[suffix_idx].second - good_it->d);  
		
		start_end.chr_start.push_back(good_it->b+1);
		start_end.read_start.push_back(suffix[suffix_idx].second - good_it->d + 1);
		
		
		strand = (type==1); 
	}
	
	
	if (good_it != good_it_back) {
		suffix_idx++;   
		good_it++;
	}
	
	
	
	
	
	
	
	while (good_it != good_it_back) {
		
		
		
		type = abs(good_it->c);
		
		if (type == Iblank) {
			
		}else if (type == Iextend ) {
			
		}else if(type == Iexonic){
			
		}else { 
			
			int temp_start_pos = good_it->a - (read_length-good_it->d) + 1;
			if(temp_start_pos < start_pos){
				start_pos = temp_start_pos;
				read_start_pos = suffix[suffix_idx].first;
			}
			
			int temp_back_pos = good_it->b + good_it->d;
			if (temp_back_pos > back_pos) {
				back_pos = temp_back_pos;
				read_back_pos = suffix[suffix_idx].second;
			}
			
			
			if (start_end.chr_end.size() > 0 && good_it->a == start_end.chr_end.back()) {
				
			}else {
				start_end.chr_end.push_back(good_it->a);
				start_end.read_end.push_back(suffix[suffix_idx].second - good_it->d);  
				
				start_end.chr_start.push_back(good_it->b + 1);
				start_end.read_start.push_back(suffix[suffix_idx].second - good_it->d + 1);
			}
			
			strand = (type==1); 
		}
		
		good_it++;
		suffix_idx++;
	}
	
	
	
	
	
	type = abs(good_it->c);
	
	
	
	
	
	
	if (type == Iextend ) {
		start_end.chr_end.push_back(good_it->b);
		if (good_it->d > 0) {
			start_end.read_end.push_back(suffix[suffix_idx].first + good_it->d - 1);
		}else {
			start_end.read_end.push_back(suffix[suffix_idx].second); 
		}
	}else if(type == Iexonic){
		int temp_start_pos = good_it->a;
		if(temp_start_pos < start_pos){
			start_pos = temp_start_pos;
			read_start_pos = suffix[suffix_idx].first;
		}
		
		start_end.chr_end.push_back(good_it->b);
		start_end.read_end.push_back(suffix[suffix_idx].second);
	}else { 
		
		int temp_start_pos = good_it->a - (read_length-good_it->d) + 1;
		if(temp_start_pos < start_pos){
			start_pos = temp_start_pos;
			read_start_pos = suffix[suffix_idx].first;
		}		
		
		
		int jun_end = good_it->b + good_it->d;
		
		if (start_end.chr_end.size() > 0 && good_it->a == start_end.chr_end.back()) {
			
		}else {
			start_end.chr_end.push_back(good_it->a);
			start_end.read_end.push_back(suffix[suffix_idx].second - good_it->d);  
			
			start_end.chr_start.push_back(good_it->b + 1);
			start_end.read_start.push_back(suffix[suffix_idx].second - good_it->d + 1);
		}
		
		start_end.chr_end.push_back(jun_end);
		start_end.read_end.push_back(suffix[suffix_idx].second);  
		
		
		strand = (type==1); 
	}
	
	if (start_end.chr_start.front() > start_pos) {
		start_end.chr_start.front() = start_pos;
		start_end.read_start.front() = read_start_pos;
	}
	
	
	if (start_end.chr_end.back() < back_pos) {
		start_end.chr_end.back() = back_pos;
		start_end.read_end.back() = read_back_pos;
	}
	
	
	start_end.strand = strand;
	
}


coord_t junction2boundary(good_t junction)
{
	coord_t result = coord_t(-1,-1);
	int readtype = abs(junction.c);
	if (readtype == Ijunpositive || readtype == Ijunnegative) {
		result.first = junction.a - (read_length-junction.d) + 1;  
		
		result.second = junction.b + junction.d;
	}else if (readtype == Iexonic) {
		result.first = junction.a ;
		result.second = junction.b;
	}else if (readtype == Iextend) {
		result.first = junction.a;
		result.second = junction.b;
	}
	
	return result;
}


void phred642phred33(string &quals)
{
	string::iterator sit = quals.begin();
	while (sit != quals.end()) {
		*sit = *sit - 31;
		sit++;
	}
}

void solexa2phred33(string &quals)
{
	string::iterator sit = quals.begin();
	while (sit != quals.end()) {
		int val = *sit - 64; 
		
		if (val < -5 || val > 62) {
			cout << "ERROR: Solexa quality score out of range, should be between -5 and 62 (inclusive)" << endl;
			cout << quals << endl;
			exit(2);
		}
		
		*sit = solexa2phred32_table[(val)+5]+33; 
		
		sit++;
	}
}

vector<pair<uint_fast32_t,uint_fast32_t> > cigar2updown(uint_fast32_t start_loc, string cigar)
{

	vector<pair<uint_fast32_t,uint_fast32_t> > output;
	
	
	size_t read_ptr = 0;
	
	while (read_ptr != string::npos) {
		size_t temp_read_ptr = cigar.find_first_of("SMN", read_ptr);
		if(temp_read_ptr != string::npos){
			if(cigar[temp_read_ptr] == 'M'){
				uint_fast32_t range = (uint_fast32_t)atoi(cigar.substr(read_ptr, temp_read_ptr-read_ptr).c_str());
				
				output.push_back(pair<uint_fast32_t,uint_fast32_t>(start_loc,start_loc+range-1));
				start_loc = start_loc + range;
			}else if(cigar[temp_read_ptr] == 'N'){ 
				uint_fast32_t range = (uint_fast32_t)atoi(cigar.substr(read_ptr, temp_read_ptr-read_ptr).c_str());
				
				start_loc = start_loc + range;
			}
			
			read_ptr = temp_read_ptr + 1;
			
			if (read_ptr == cigar.length()) {
				break;
			}
		}else {
			read_ptr = temp_read_ptr;
		}
		
		
	}
	
	
	
	
	return output;
}

vector<pair<uint_fast32_t,uint_fast32_t> > SAM2updown(string line) 
{
	vector<string> line_list = split(line, '\t');
	
	vector<pair<uint_fast32_t,uint_fast32_t> > output;
	
	
	string cigar = line_list[5];
	uint_fast32_t start_loc = (uint_fast32_t)atoi(line_list[3].c_str());
	size_t read_ptr = 0;
	
	while (read_ptr != string::npos) {
		size_t temp_read_ptr = cigar.find_first_of("SMN", read_ptr);
		if(temp_read_ptr != string::npos){
			if(cigar[temp_read_ptr] == 'M'){
				uint_fast32_t range = (uint_fast32_t)atoi(cigar.substr(read_ptr, temp_read_ptr-read_ptr).c_str());
				
				output.push_back(pair<uint_fast32_t,uint_fast32_t>(start_loc,start_loc+range-1));
				start_loc = start_loc + range;
			}else if(cigar[temp_read_ptr] == 'N'){ 
				uint_fast32_t range = (uint_fast32_t)atoi(cigar.substr(read_ptr, temp_read_ptr-read_ptr).c_str());
				
				start_loc = start_loc + range;
			}
			
			read_ptr = temp_read_ptr + 1;
			
			if (read_ptr == cigar.length()) {
				break;
			}
		}else {
			read_ptr = temp_read_ptr;
		}
		
		
	}
	
	
	return output;
}


bool read_reference_map(string chromosome_path,map<string,reference_t> &ref_map, vector<string> chr_file_list)
{
        // ostools is not portable -> get chr_file_list directly from caller in runSpliceMap.cpp
        //ostools os(ostools::NIX);
	//vector<string> chr_file_list;
	//chr_file_list = os.list_dir(chromosome_path);
	
	vector<string>::iterator chr_file_list_it = chr_file_list.begin();
	while (chr_file_list_it != chr_file_list.end()) {
		ifstream reference_file;
		string line;
		
		
		size_t split_loc = chr_file_list_it->rfind('/');  
		string name = chr_file_list_it->substr(split_loc+1);
		string path = chr_file_list_it->substr(0,split_loc+1);
		
		reference_file.open(chr_file_list_it->c_str(), ios::in);
		
		if (!reference_file.is_open()) {
			cout << "ERROR: failed to open " << *chr_file_list_it << endl;
			return false;
		}
		
		string chr_name = "";
		uint_fast64_t file_index_start = 0;
		uint_fast64_t file_index_end = 0;
		
		uint_fast64_t curr_pos = 0;
		
		while (!reference_file.eof()) {
			getline(reference_file, line);
			curr_pos = curr_pos + line.length() + 1; 
			trim2(line);
			if (line.length() == 0) {
				continue;
			}
			
			if (line[0] == '>') { 

				
				if (chr_name.length() != 0) {
					reference_t contents;
					contents.file_name = name;
					contents.file_path = path;
					contents.file_index_start = file_index_start;
					contents.file_index_end = file_index_end;
					
					pair<map<string,reference_t>::iterator,bool> out = ref_map.insert(pair<string,reference_t>(chr_name,contents));
					if (!out.second) {
						cout << "ERROR: duplicate chromosome (" << chr_name << ") in file: " << path << name << endl;
						cout << "Please ensure that all genome files are unique." << endl;
						exit(2);
					}
					
					
				}
				

				file_index_start = curr_pos;
				
				size_t first_blank = line.length();
				
				for (size_t i = 1; i<line.length(); i++) {
					if(isspace(line[i])){
						first_blank = i;
						break;
					}
				}
				
				chr_name = line.substr(1,first_blank-1);  
				
				
			}
			
			
			file_index_end = curr_pos;
		}
		
		
		if (chr_name.length() != 0) {
			reference_t contents;
			contents.file_name = name;
			contents.file_path = path;
			contents.file_index_start = file_index_start;
			contents.file_index_end = file_index_end;
			
			pair<map<string,reference_t>::iterator,bool> out = ref_map.insert(pair<string,reference_t>(chr_name,contents));
			if (!out.second) {
				cout << "ERROR: duplicate chromosome (" << chr_name << ") in file: " << path << name << endl;
				cout << "Please ensure that all genome files are unique." << endl;
				exit(2);
			}
			
			
		}
		
		reference_file.close();
		reference_file.clear();
		
		chr_file_list_it++;
	}
	
	return true;
}
