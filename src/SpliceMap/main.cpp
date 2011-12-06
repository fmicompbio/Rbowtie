

#include "main.h"

using namespace std;


int main (int argc, char * const argv[]) 
{	
	string genome_path; 
	string genome_filename; 
	string chr_name;  
	string genome_start_loc_str;
	uint_fast64_t genome_file_index_start = 0;
	uint_fast64_t genome_file_index_end = 0;
	
	string temp_path;
	
	ofstream debug_out;
	
	
	vector<string> reads_filename[2];  
	int num_pair = 0; 
	string temp("temp"); 
	string line("");
	
	bool DNA_mode = false;
	int min_intron_search_distance = 20;
	int extra_length = 0; 
	
	
	bool cufflinks = false;
	
	seq_hash* chrdict_10 = new seq_hash[MAX_HASH]; 
	for (uint_fast32_t i=0; i<MAX_HASH; i++) {
		(chrdict_10[i]) = NULL; 
	}
	
	uint_fast32_t chrsize = 0;   
	string chrline; 
	
	ifstream genome_filename_file;
	
	ofstream good_SAM_file;
	
	
	
	
	
	
	
	
	
	
	
	
	int_fast32_t num_seed_mismatch = 1; 
	int_fast32_t num_read_mismatch = 2;  
	int_fast32_t max_clip_allowed = 40;  
	uint_fast32_t intron_curve[3] = {20000,100000,400000}; 
	uint_fast32_t min_extend_req = 0; 
	
	
	
	
	
	
	struct timeval tv;
	struct timeval start_tv;
	struct timeval temp_start_tv;
	
	gettimeofday(&start_tv, NULL);
	 
	
	
	
	
	
	if (argc == 3) {
		
		
		temp = argv[1];
		chr_name = argv[2];
		
		
		cfgfile *run_cfg = new cfgfile(temp);
		
		
		temp_path = run_cfg->getVal("temp_path");
		if (temp_path.length() == 0) {
			
			temp_path = "temp/";
		}else {
			if (temp_path[temp_path.length()-1] != '/') {
				temp_path = temp_path + '/';
			}
		}
	
		
		ifstream ref_map_file;
		temp = temp_path+reference_filename;
		ref_map_file.open(temp.c_str(), ios::in);
		
		if (!ref_map_file.is_open()) {
			cout << "ERROR: cannot read reference list " << temp << endl;
			exit(1);
		}
		
		while (!ref_map_file.eof()) {
			getline(ref_map_file, line);
			trim2(line);
			if (line.length() == 0) {
				continue;
			}
			
			vector<string> line_list = split(line, '\t');
			
			if (line_list[0] == chr_name) {
				genome_path = line_list[2];
				genome_filename = line_list[1];
				stringstream ss;
				uint_fast64_t pos;
				genome_start_loc_str = line_list[3]; 
				ss.str(line_list[3]);
				ss >> pos;
				genome_file_index_start = pos;
				ss.clear();
				ss.str(line_list[4]);
				ss >> pos;
				genome_file_index_end = pos;
				
				
			}
		}
		
		ref_map_file.close();
		ref_map_file.clear();
		
		if (genome_filename.length() == 0) {
			cout << "ERROR: Failed to read chromosome " << chr_name << endl;
			exit(1);
		}
		
		
		temp = temp_path+debug_path+genome_filename+"_"+genome_start_loc_str+".log";
		
		debug_out.open(temp.c_str(), ios::out);
		
		
		reads_filename[0] = run_cfg->getList("reads_list1");
		reads_filename[1] = run_cfg->getList("reads_list2");
		
		
		
		if(reads_filename[0].size() == 0){
			cout << "ERROR: There should be at least one reads file in the first list. " << endl;
			debug_out << "ERROR: There should be at least one reads file in the first list. " << endl;
			
		
			
			exit(2);
		}
		
		if(reads_filename[1].size() == 0){
			num_pair = 1;
		}else {
			num_pair = 2;
		}
		
		reads_filename[0].clear();
		reads_filename[1].clear();
		
		
		
		
		
		
		temp = run_cfg->getVal("sam_file");
		if (temp.compare("sam") == 0) {
		}else if(temp.compare("cuff") == 0){
			cufflinks = true;
		}else if(temp.length() == 0){
			cufflinks = false;
		}else {
			cout << "ERROR: The choices for SAM output are \"sam\", \"cuff\" ... not  " << temp << endl;
			debug_out << "ERROR: The choices for SAM output are \"sam\", \"cuff\" ... not  " << temp << endl;
		}
		
		
		
		temp = run_cfg->getVal("seed_mismatch");
		if (temp.length() > 0) {
			num_seed_mismatch = atoi(temp.c_str());
		}
		
		
		temp = run_cfg->getVal("read_mismatch");
		if (temp.length() > 0) {
			num_read_mismatch = atoi(temp.c_str());
		}
		
		
		temp = run_cfg->getVal("max_clip_allowed");
		if (temp.length() > 0) {
			max_clip_allowed = atoi(temp.c_str());
		}
		
		temp = run_cfg->getVal("max_intron");
		int max_intron = atoi(temp.c_str());
		
		temp = run_cfg->getVal("min_intron");
		int min_intron = atoi(temp.c_str());
		
		
		
		
		if (max_intron > 0 && min_intron > 0) {
			
			if (min_intron > max_intron) {
				max_intron = min_intron;
			}
			
			intron_curve[0] = min_intron;
			intron_curve[1] = (uint_fast32_t)floor(max_intron-(max_intron-min_intron)*0.75);
			intron_curve[2] = max_intron;
		}else if (max_intron > 0){
			intron_curve[0] = (uint_fast32_t)floor(max_intron*0.10);
			intron_curve[1] = (uint_fast32_t)floor(max_intron*0.75);
			intron_curve[2] = max_intron;
		}else if (min_intron > 0) {
			intron_curve[0] = min_intron;
			intron_curve[1] = max_intron*5;
			intron_curve[2] = min_intron*20;
		}else {
			
		}
		
		
		temp = run_cfg->getVal("min_extend_req");
		min_extend_req = atoi(temp.c_str()); 

		
		delete run_cfg;
		
		
		
	} else {
		print_usage_and_exit();
	}
	
	
	
	
	debug_out << "Reading read-file-list files!... " << endl;
	string reads_filename_lists[2] = {temp_path+"reads_list1",temp_path+"reads_list2"}; 

	for (int i = 0; i<num_pair; i++) {
		ifstream reads_filename_lists_file;
		
		reads_filename_lists_file.open(reads_filename_lists[i].c_str(), ios::in);
		if (reads_filename_lists_file.is_open()) {
			while (!reads_filename_lists_file.eof()) {
				getline(reads_filename_lists_file, temp);
				trim2(temp);
				if (temp.length()>0) {
					reads_filename[i].push_back(temp);
				}
			}
			reads_filename_lists_file.close();
			reads_filename_lists_file.clear();
		}else {
			cout << "ERROR: I'm sorry, cannot open the reads filename list file, " << reads_filename_lists[i] << endl;
			debug_out << "ERROR: I'm sorry, cannot open the reads filename list file, " << reads_filename_lists[i] << endl;
			
			
			
			exit(1);
		}
	}
	
	 
	 
	
	if (num_seed_mismatch > 2){
		cout << "I'm sorry, SpliceMap " << VERSION << " supports a maximum of 2 mismatchs in the seeds" << endl;
		debug_out << "I'm sorry, SpliceMap " << VERSION << " supports a maximum of 2 mismatchs in the seeds" << endl;
		
		
		
		exit(1);
	}
	
	
	debug_out << "Creating output SAM file" << endl;
	temp = temp_path+genome_filename + "_" + genome_start_loc_str + ".sam";
	good_SAM_file.open(temp.c_str(), ios::out);
	
	if (!good_SAM_file.is_open()) {
		cout << "ERROR: something went wrong creating SAM file :( -- " << temp << endl;
		debug_out << "ERROR: something went wrong creating SAM file :( -- " << temp << endl;
	}
	
	

	
	
	
	
	debug_out << "Reading and indexing reference genome!..." << endl;

	temp = genome_path+genome_filename;
	genome_filename_file.open(temp.c_str(), ios::in);
	
	
	stringstream chr_buf (stringstream::in | stringstream::out);  
	
	if (genome_filename_file.is_open()) {
		
		genome_filename_file.seekg(genome_file_index_start); 
		
		
		
		debug_out << "Chromosome name: " << chr_name << endl;
		
		uint_fast64_t curr_pos = genome_file_index_start;
		while (curr_pos < genome_file_index_end) {
			getline(genome_filename_file, temp);
			curr_pos =curr_pos + temp.length() + 1;
			trim2(temp);
			
			string prev = temp;
			if(!make_DNA_upper(temp)){
				cout << "Warning: Unexpected character in genome (changed to N):" << endl;
				cout << prev << endl;
				debug_out << "Warning: Unexpected character in genome (changed to N):" << endl;
				debug_out << prev << endl;
				
			}			
			chr_buf << temp;
		}
		genome_filename_file.close();
		genome_filename_file.clear();

	}else {
		cout << "ERROR: I'm sorry, cannot open the genome file, " << genome_path+genome_filename << endl;
		debug_out << "ERROR: I'm sorry, cannot open the genome file, " << genome_path+genome_filename << endl;
		
		
		
		exit(1);
	}
	
	
	chrline = chr_buf.str();
	
	
	chrsize = (uint_fast32_t)chrline.length();
	
	
	uint_fast32_t key_num = 0; 
	int pos = -9;  
	uint_fast32_t num_read = 0; 
	
	while(num_read < 10){
		char c;
		chr_buf.get(c);
		int val = get_base_val(c);
		if (val < 0) {
			num_read = 0;
			key_num = 0;
		}else {
			
			key_num = key_num + (val<<(2*num_read));
			num_read++;
		}
		pos++;
	}
	
	
	
	seq_hash v = new vector<uint_fast32_t>();
	v->push_back(pos);
	chrdict_10[key_num] = v;
	
	
	while (!chr_buf.eof()) {
		char c;

		chr_buf.get(c);
		int val = get_base_val(c);
		
		if (val >= 0) {
			key_num = key_num>>2;  
			
			key_num = key_num + (val<<18);  
			pos++;
		}else {
			num_read = 0; 
			key_num = 0;
			pos++; 
			
			while(num_read < 10 && !chr_buf.eof()){
				chr_buf.get(c);
				val = get_base_val(c);
				if (val < 0) {
					num_read = 0;
					key_num = 0;
				}else {
					
					
					key_num = key_num + (val<<(2*num_read));
										 
					num_read++;
				}
				pos++;
			}
			
		}
		
		if(chr_buf.eof()){
			break; 
		}
		
		
		
		
		
		if(chrdict_10[key_num]!=NULL){
			chrdict_10[key_num]->push_back(pos);
		}else {
			seq_hash vec = new vector<uint_fast32_t>();
			vec->push_back(pos);
			chrdict_10[key_num] = vec;
		}
		
		
	}
	
	
	debug_out << "Chromosome size: " << chrsize << endl;

	
	chr_buf.str(string());
	
	
	

	gettimeofday(&tv, NULL);

	
	debug_out<<"Index creation time: "<< diffclock(start_tv,tv) <<" s."<<endl;
	
	gettimeofday(&temp_start_tv, NULL);
	
	for (uint_fast32_t i = 0; i<MAX_HASH; i++) {
		
		if (chrdict_10[i] != NULL) {  
			
			sort(chrdict_10[i]->begin(),chrdict_10[i]->end());
		}
	}
	
	
	gettimeofday(&tv, NULL);
	
	
	debug_out<<"Index sorting time: "<< diffclock(temp_start_tv,tv) <<" s."<<endl;
	
	
	
	
	
	
	
	
	good_SAM_file << "@HD\tVN:1.0\tSO:coordinate" << "\n";
	good_SAM_file << "@SQ\tSN:" << chr_name << "\tLN:" << chrsize << "\n";
	good_SAM_file << "@PG\tID:SpliceMap\tVN:" << VERSION << "\n";
	

	
	
	debug_out << "SpliceMapping on "<< chr_name << "!... please wait" << endl;

	gettimeofday(&temp_start_tv, NULL);
	
	
	
	
	
	
	uint_fast32_t len = (uint_fast32_t)reads_filename[0].size(); 
	
	for(uint_fast32_t file_idx = 0;file_idx< len;file_idx++){
		ifstream fullread[2]; 
		ifstream fullqual[2]; 
		ifstream fullname[2]; 
		
		
		
		uint_fast32_t fullread_length[2];
		
		
		ifstream read_index_tfile[2];
		
		for (uint_fast8_t pair_idx = 0;pair_idx<num_pair;pair_idx++){
			temp = reads_filename[pair_idx][file_idx] + ".index.t";
			read_index_tfile[pair_idx].open(temp.c_str());
			
			
			
			fullread[pair_idx].open(reads_filename[pair_idx][file_idx].c_str(), ios::in);
			if(!fullread[pair_idx].is_open()){
				cout << "ERROR: I'm sorry, cannot open the full-reads file, "<< reads_filename[pair_idx][file_idx] << endl;
				debug_out << "ERROR: I'm sorry, cannot open the full-reads file, "<< reads_filename[pair_idx][file_idx] << endl;
				
				
				
				exit(1);
			}
			
			
			temp = reads_filename[pair_idx][file_idx] + ".quals";
			fullqual[pair_idx].open(temp.c_str(), ios::in);
			if(!fullqual[pair_idx].is_open()){
				cout << "ERROR: I'm sorry, cannot open the full-quals file, "<< temp << endl;
				debug_out << "ERROR: I'm sorry, cannot open the full-quals file, "<< temp << endl;
				
				
				
				exit(1);
			}
			
			temp = reads_filename[pair_idx][file_idx] + ".names";
			fullname[pair_idx].open(temp.c_str(), ios::in);
			if(!fullname[pair_idx].is_open()){
				cout << "ERROR: I'm sorry, cannot open the full-name file, "<< temp << endl;
				debug_out << "ERROR: I'm sorry, cannot open the full-name file, "<< temp << endl;
				
				
				
				exit(1);
			}
			
			
			
		}
		
		
		
		uint_fast32_t curr_read_idx = 1;
		
		
		
		
		while (!read_index_tfile[0].eof()) { 
			
			
			string full_line_read[2] = ""; 
			string full_line_qual[2] = ""; 
			string full_line_name[2] = ""; 
			
			
			
			getline(read_index_tfile[0], line);
			trim2(line);
			
			if (line.length() == 0) {
				continue;
			}
			
			fullread_length[0] = atoi(line.c_str());
			
			if (num_pair == 2) {
				getline(read_index_tfile[1], line);
				trim2(line);
				
				if (line.length() == 0) {
					continue;
				}
				
				fullread_length[1] = atoi(line.c_str());
			}
			

			vector<coord_t> suffix_list[2];
			
			
			for (int pair_idx = 0; pair_idx<num_pair; pair_idx++) {
				suffix_list[pair_idx] = designsuffix(fullread_length[pair_idx]);
			}
			
			good_seg_vec_t goodpair_list[2];  

			
			
			for (int pair_idx = 0; pair_idx<num_pair; pair_idx++) {
				
				
				
				getline(fullread[pair_idx], full_line_read[pair_idx]); 
				
				if(full_line_read[pair_idx].length() == 0){
					cout << "ERROR: blank line? -- " << reads_filename[pair_idx][file_idx] << "," << curr_read_idx << endl;
					exit(2);
				}
				
				
				if(!make_DNA_upper(full_line_read[pair_idx])){
					cout << "ERROR: Unexpected character in read:" << endl;
					cout << full_line_read[pair_idx] << endl;
					exit(2);
				}
				
				
				getline(fullqual[pair_idx], full_line_qual[pair_idx]);
				if (full_line_qual[pair_idx][0] == 0) {
					full_line_qual[pair_idx] = "";
				}
				getline(fullname[pair_idx], full_line_name[pair_idx]);
				
				
				list<good_vec_t> group_result; 
				
				uint_fast8_t group_len = (uint_fast8_t)suffix_list[pair_idx].size();  
				
				
				for (uint_fast8_t seg = 0;seg<group_len;seg++){
					good_vec_t good_list_temp[2];
					good_vec_t good_list;
					good_vec_t good_exonic_list;
					
					uint_fast32_t local_start = (suffix_list[pair_idx][seg]).first;
					uint_fast32_t local_end = (suffix_list[pair_idx][seg]).second;

					vector<mix_t> mix_list[2];
					
					
					

					list<line_t> temp_line_list[2]; 
					
					for (int k = 0;k<2;k++){
						
						
						
						getline(read_index_tfile[pair_idx], line);
						uint_fast32_t curr_num_matches = (uint_fast32_t)atoi(line.c_str());
						
						for (uint_fast32_t match_idx = 0; match_idx < curr_num_matches; match_idx++) {
							getline(read_index_tfile[pair_idx], line);
							
							vector<string> line_list = split(line, '\t');
							
							
							
							uint_fast32_t chr_pos = (uint_fast32_t)atoi(line_list[1].c_str());
							int_fast32_t mismatch_direction = (int_fast32_t)atoi(line_list[2].c_str()); 
							
							
							
							
							if(chr_name == line_list[0]){
								line_t temp_line;
								temp_line.chr_name = line_list[0];
								temp_line.unique = (curr_num_matches == 1);
								temp_line.errors = abs(mismatch_direction)-1;
								
								
								
								
								temp_line.coor = chr_pos;
								temp_line.direction = mismatch_direction>0;
	
								
								temp_line_list[k].push_back(temp_line);
								
							}
						}
						
						
						
						
						
						
						
						
						
						if (temp_line_list[k].size() > 1) {
							
							DelLocalHalfMappableHits(temp_line_list[k], intron_curve[2]);
						}
					}
					
					string linefullreads = "";
					
					if(temp_line_list[0].size() == 0 && temp_line_list[1].size() == 0 ){ 
						
					}else {
						
						
						
						
						linefullreads = full_line_read[pair_idx].substr(local_start-1, local_end-local_start+1); 
						
						if (temp_line_list[0].size() > 0 && temp_line_list[1].size() > 0 ) {
							FindExonicRead(good_list,temp_line_list);
						}
					}
					
					
					
					
					for (int k = 0;k<2;k++){
						list<line_t>::iterator it = temp_line_list[k].begin();
						while (it != temp_line_list[k].end()) {
							
							
							
							if (it->errors <= num_seed_mismatch) { 
								
								mix_t m;
								m.full_line = linefullreads;
								m.seed_line = linefullreads.substr(k * read_halflen, read_halflen);
								m.unique = it->unique;
								m.errors = it->errors;
								m.chr_name = it->chr_name;
								m.coor = it->coor;
								m.direction = it->direction;
								
								
								
								
								
								
								
								
								
								if (it->direction) { 
									
									mix_list[k].push_back(m);
									
								}else if (!it->direction) { 
									
									mix_list[(k+1)%2].push_back(m);
									
								}
							}
							it++;
						}
					}
					
					
					
					
					good_exonic_list = good_list;
					

					
					

					
					
					for (int k = 0; k<2; k++) {
						
						good_list_temp[k].clear();
						
						vector<mix_t>::iterator mix_it = mix_list[k].begin();
						while (mix_it != mix_list[k].end()) {
							
							mix_t mix_temp_line = *mix_it;
							
							int IRF = 1; 
							
							if (!(mix_temp_line.direction)) {
								
								compleseq(mix_temp_line.full_line);
								compleseq(mix_temp_line.seed_line);
								
								IRF = -1;
							}
							
							
							uint_fast32_t chr_pos = mix_temp_line.coor - 1 + ((k+1)%2)*(read_halflen-1); 
							
							
							int_fast32_t read_pos = read_halflen - ((k+1)%2);
							
							
							list< coord_t > coord_list[2]; 
							
							
							
							
							if(k==0){  
								
								while (read_pos < 50 && chr_pos < (uint_fast32_t)(chrsize-read_halflen) 
									   && mix_temp_line.full_line[read_pos] == chrline[chr_pos]) {  
									
									read_pos++;
									chr_pos++;
									
									if(!DNA_mode){
										
										if (chrline[chr_pos] == 'G' && chrline[chr_pos+1] == 'T') {
											coord_list[0].push_back(coord_t(read_pos,chr_pos));
										}
										if (chrline[chr_pos] == 'C' && chrline[chr_pos+1] == 'T') {
											coord_list[1].push_back(coord_t(read_pos,chr_pos));
										}
									}else {
										
										coord_list[0].push_back(coord_t(read_pos,chr_pos));
									}
									
									
									
								}
							}else {  
								
								
								
								
								
								
								
								
								
								
								while (read_pos >= 0 && chr_pos > (uint_fast32_t)(read_halflen-1) && mix_temp_line.full_line[read_pos] == chrline[chr_pos]) {
									
									
									

									read_pos--;
									chr_pos--;
									
									
									if(!DNA_mode){
										
										if (chrline[chr_pos-1] == 'A' && chrline[chr_pos] == 'G') {
											coord_list[0].push_back(coord_t(read_pos,chr_pos));
										}
										if (chrline[chr_pos-1] == 'A' && chrline[chr_pos] == 'C') {
											coord_list[1].push_back(coord_t(read_pos,chr_pos));
										}
									}else {
										coord_list[1].push_back(coord_t(read_pos,chr_pos));
									}
									
									
								}
								
							}
							
							int longest_read_extend_pos = read_pos;
							int furthest_read_extend_chr_pos = chr_pos;
							
							
							
							
							
							
							vector<good_t> result_list;
							
							
							int direc_type = 0;

							
							for (int n=0; n<2; n++) { 
								direc_type++;  
								
								list< coord_t >::iterator it = coord_list[n].begin();
								while (it != coord_list[n].end()) {
									read_pos = it->first;  
									chr_pos = it->second; 
									
									if (k == 0){  
										
										string line_left = mix_temp_line.full_line.substr(read_pos); 
										
										if (read_pos < (read_length-(9+extra_length))){ 
											int i = 10;
											int line_left_len = (int)line_left.length();
											string line_left_10 = line_left.substr(0, 10);
											
											int key_val = base2int(line_left_10);

											
											if(key_val>=0 && chrdict_10[key_val] != NULL){
												uint_fast32_t half_distance = cal_half_distance(line_left_len,intron_curve);
												
												pair<vector<uint_fast32_t>::iterator, vector<uint_fast32_t>::iterator> left_indexlist = checkindexlistFwd(chr_pos, (*chrdict_10[key_val]), half_distance, min_intron_search_distance);

												vector<uint_fast32_t>::iterator index_it = left_indexlist.first;
												while (index_it != left_indexlist.second   && index_it != chrdict_10[key_val]->end()) {
													
													
													if (DNA_mode && *index_it == chr_pos - (read_pos-8)) {
														index_it++;
														
														continue;
													}
													
													string tempread;
													
													try{
														tempread = chrline.substr((*index_it)-1, read_halflen);  
													}catch (std::out_of_range& e) {
														
														
														index_it++;
														continue; 
													}
													
													i = 10;
													
													
													
													if (DNA_mode||(chrline[(*index_it)-3]=='A' && chrline[(*index_it)-2]=='G' && direc_type == 1)
														||(chrline[(*index_it)-3]=='A' && chrline[(*index_it)-2]=='C' && direc_type == 2)) {
														

														while (i < line_left_len && tempread[i] == line_left[i]) { 
															i++;
														}
														
														
														
														if (i == line_left_len) {
															
															
															good_t contents;
															contents.a = chr_pos;
															contents.b = (*index_it);
															contents.c = IRF*direc_type;
															contents.d = line_left_len;
															
															result_list.push_back(contents);
															
														}
													}
													index_it++;
												}

											}
										}
									}else {  
										
										string line_left = mix_temp_line.full_line.substr(0, read_pos+1);
										
										if (read_pos > (8+extra_length)) { 
											int i = 10;
											int line_left_len = (int)line_left.length();
											string line_left_10 = line_left.substr(line_left_len-10, 10);

											int key_val = base2int(line_left_10);

											
											if(key_val>=0 && chrdict_10[key_val] != NULL){
												
												uint_fast32_t half_distance = cal_half_distance(line_left_len,intron_curve);
												
												pair<vector<uint_fast32_t>::iterator, vector<uint_fast32_t>::iterator> left_indexlist = checkindexlistRes(chr_pos, (*chrdict_10[key_val]), half_distance, min_intron_search_distance);
												
												vector<uint_fast32_t>::iterator index_it = left_indexlist.first;
												while (index_it != left_indexlist.second  && index_it != chrdict_10[key_val]->end()) {
													
													if (DNA_mode && *index_it == chr_pos - (read_pos-10)) {
														index_it++;
														
														continue;
													}
													
													
													string tempread;
													
													try{
														tempread = chrline.substr((*index_it)-read_halflen+9, read_halflen); 
													} catch (std::out_of_range& e) {
														
														
														index_it++;
														continue; 
													}
													
													i = 10;
													
													
													
													if (DNA_mode||(chrline[(*index_it)+9]=='G' && chrline[(*index_it)+10]=='T' && direc_type == 1)
														||(chrline[(*index_it)+9]=='C' && chrline[(*index_it)+10]=='T' && direc_type == 2)) {
														
														
														while (i < line_left_len && tempread[read_halflen-i-1] == line_left[line_left_len-i-1]) { 
															i++;
														}
														if (i == line_left_len) {
															
															good_t contents;
															contents.a = chr_pos;
															contents.b = (*index_it);
															contents.c = IRF*direc_type;
															contents.d = line_left_len;
															
															result_list.push_back(contents);
														}
													}
													index_it++;
												}
											}
										}
									}
									
									it++;
								}

							}
							
							
							
							
							
							if(k==0){  
								
								
								
								if (longest_read_extend_pos != 50 && result_list.size() > 0) {
									
									
									if(result_list.size() == 1){  
										good_t contents = result_list.front();
										
										contents.b = contents.b - 1;
										
										if(contents.a != contents.b){
											
											vector<good_t>::iterator good_it = good_list.begin();
											int Indicator_exist = 0;
											while (good_it != good_list.end()) {
												if(isgoodtequal(contents,(*good_it))){  
													Indicator_exist++;
													break;
												}
												good_it++;
											}
											if (Indicator_exist == 0) {
												good_list.push_back(contents);
											}
										}
									}else { 
										vector<good_t>::iterator result_it = result_list.begin();
										while (result_it != result_list.end()) {
											
											result_it->b = result_it->b - 1;
											
											if(result_it->c > 0){
												if (abs(result_it->c) == Ijunpositive) {
													result_it->c = Mjunpositive;
												}else {
													result_it->c = Mjunnegative;
												}

											}else {
												if (abs(result_it->c) == Ijunpositive) {
													result_it->c = -Mjunpositive;
												}else {
													result_it->c = -Mjunnegative;
												}
											}

											
											if(result_it->a != result_it->b){
												
												vector<good_t>::iterator good_it = good_list.begin();
												int Indicator_exist = 0;
												while (good_it != good_list.end()) {
													if(isgoodtequal(*result_it,(*good_it))){  
														Indicator_exist++;
														break;
													}
													good_it++;
												}
												if (Indicator_exist == 0) {
													good_list.push_back(*result_it);
												}
											}
											
											result_it++;
										}
									}

									
								}else if (longest_read_extend_pos!=50) { 
									good_t contents = {furthest_read_extend_chr_pos-longest_read_extend_pos+1,
										furthest_read_extend_chr_pos,IRF*Iextend,longest_read_extend_pos};
									
									
									InsertGoodExtend(good_list,good_exonic_list,contents);
									
								}
								
								
								
							}else {  
								
								
								
								
								
								if (longest_read_extend_pos != -1 && result_list.size() > 0) {
									
									if(result_list.size() == 1){  
										
										good_t contents = result_list.front();
										
										contents.a = contents.a + 1;
										contents.b = contents.b + 9;
										int temp_int = contents.a;
										contents.a = contents.b;
										contents.b = temp_int;
										contents.d = read_length-contents.d;
										
										if(contents.a != contents.b){
											vector<good_t>::iterator good_it = good_list.begin();
											int Indicator_exist = 0;
											while (good_it != good_list.end()) {
												if(isgoodtequal(contents,(*good_it))){
													
													Indicator_exist++;
												}
												good_it++;
											}
											if (Indicator_exist == 0) {
												good_list.push_back(contents);
											}
										}
									}else {
										vector<good_t>::iterator result_it = result_list.begin();
										while (result_it != result_list.end()) {
											
											result_it->a = result_it->a + 1;
											result_it->b = result_it->b + 9;
											int temp_int = result_it->a;
											result_it->a = result_it->b;
											result_it->b = temp_int;
											result_it->d = read_length-result_it->d;
											
											if(result_it->c > 0){
												if (abs(result_it->c) == Ijunpositive) {
													result_it->c = Mjunpositive;
												}else {
													result_it->c = Mjunnegative;
												}
												
											}else {
												if (abs(result_it->c) == Ijunpositive) {
													result_it->c = -Mjunpositive;
												}else {
													result_it->c = -Mjunnegative;
												}
											}
											
											
											if(result_it->a != result_it->b){
												vector<good_t>::iterator good_it = good_list.begin();
												int Indicator_exist = 0;
												while (good_it != good_list.end()) {
													if(isgoodtequal(*result_it,(*good_it))){
														
														Indicator_exist++;
													}
													good_it++;
												}
												if (Indicator_exist == 0) {
													good_list.push_back(*result_it);
												}
											}
											
											result_it++;
										}
									}
									
									
								}else if (longest_read_extend_pos!=-1) {  
									good_t contents = {furthest_read_extend_chr_pos+2,
										50 - longest_read_extend_pos+furthest_read_extend_chr_pos,IRF*Iextend,
										longest_read_extend_pos-49}; 
									
									
									InsertGoodExtend(good_list,good_exonic_list, contents);
									
								}
								
								

							}
							

							mix_it++;
						}
						
						

						
						
						good_list_temp[k] = good_list;
						good_list.clear();
						
					}
					
					
					good_exonic_list.clear();
					
					
					
					
					good_list = checkgoodlists(good_list_temp[0], good_list_temp[1]);
					
					group_result.push_back(good_list);
					
					
					
				}

				
				
				
				
				
				good_seg_vec_t group_filtered_results_list;
				
				bool empty_list = true;
				vector< vector<good_t> > singleread_group_result;
				pair<vector<coord_t>,vector<coord_t> > seg_suffix_list; 
				int max_contig = 0;
				int curr_contig = 0;
				
				
				
				
				
				
				vector<coord_t>::iterator suffix_it = suffix_list[pair_idx].begin();
				vector<coord_t>::reverse_iterator suffix_rev_it = suffix_list[pair_idx].rbegin();
				
				
				for (list<good_vec_t>::iterator good_it = group_result.begin(); good_it != group_result.end(); good_it++) {

					if (good_it->size() > 0) {
						
						singleread_group_result.push_back(*good_it);
						
						
						
						
						
						seg_suffix_list.first.push_back(*suffix_it);
						seg_suffix_list.second.push_back(*suffix_rev_it);
						
						curr_contig++;
						if (curr_contig > max_contig) {
							max_contig = curr_contig;
						}
						
						empty_list = false;
						
						
					}else {
						
						
						vector<good_t> contents;
						good_t good;
						good.a = 0;
						good.b = 0;
						good.c = Iblank;
						good.d = 0;
						
						contents.push_back(good);
						
						singleread_group_result.push_back(contents);
						seg_suffix_list.first.push_back(*suffix_it);
						seg_suffix_list.second.push_back(*suffix_rev_it);
						
						curr_contig = 0;
						
					}
					
					suffix_it++;
					suffix_rev_it++;
				}
				
				if (!empty_list){
					
					
					
					
					
					vector< vector<good_t> >::iterator temp_it = singleread_group_result.begin();
					pair<vector<coord_t>::iterator,vector<coord_t>::iterator> temp2_it(seg_suffix_list.first.begin(),seg_suffix_list.second.begin());
					while (temp_it->front().c == Iblank) {
						temp_it++;
						(temp2_it.first)++;
						(temp2_it.second)++;
					}
					if (temp_it != singleread_group_result.begin()) {
						singleread_group_result.erase(singleread_group_result.begin(),temp_it); 
						seg_suffix_list.first.erase(seg_suffix_list.first.begin(),temp2_it.first);
						seg_suffix_list.second.erase(seg_suffix_list.second.begin(),temp2_it.second);
					}
					
					int first_contig_length = 0;
					temp_it = singleread_group_result.begin();
					temp2_it = pair<vector<coord_t>::iterator,vector<coord_t>::iterator>(seg_suffix_list.first.begin(),seg_suffix_list.second.begin());
					vector< vector<good_t> >::iterator begin_temp_it = singleread_group_result.begin();
					pair<vector<coord_t>::iterator,vector<coord_t>::iterator> begin_temp2_it(seg_suffix_list.first.begin(),seg_suffix_list.second.begin());
					while (temp_it != singleread_group_result.end()) {
						if (temp_it->front().c != Iblank) {
							first_contig_length++;
						}else {
							if (first_contig_length>=2) {
								break;
							}else {
								first_contig_length = 0;
								begin_temp_it = temp_it+1;
								begin_temp2_it.first = temp2_it.first+1;
								begin_temp2_it.second = temp2_it.second+1;
							}
							
						}
						
						
						temp_it++;
						temp2_it.first++;
						temp2_it.second++;
					}
					
					singleread_group_result.erase(temp_it,singleread_group_result.end()); 
					seg_suffix_list.first.erase(temp2_it.first,seg_suffix_list.first.end());
					seg_suffix_list.second.erase(temp2_it.second,seg_suffix_list.second.end());
					
					singleread_group_result.erase(singleread_group_result.begin(),begin_temp_it); 
					seg_suffix_list.first.erase(seg_suffix_list.first.begin(),begin_temp2_it.first);
					seg_suffix_list.second.erase(seg_suffix_list.second.begin(),begin_temp2_it.second);
					
					
					
					
					bool read_is_ok = false;
					if (suffix_list[pair_idx].size() == 1) {
						if(singleread_group_result.size() > 0){
							read_is_ok = true;
						}
					}else if(suffix_list[pair_idx].size() == 2){
						if(singleread_group_result.size() > 1){  
							read_is_ok = true;
						}
					}else if(suffix_list[pair_idx].size() > 2){
						
						
						if(singleread_group_result.size() > 1){
							read_is_ok = true;
						}
					}
					
					
					
					
					
					if (read_is_ok){
						
						
						
						
						
						
						
						
						
						vector< list<good_t> > temp_group = check_singleread_group_result(singleread_group_result, seg_suffix_list);
						if(temp_group.size() == 0 && singleread_group_result.size() > 2){
							
							int num_pop = 1;
							vector< vector<good_t> > temp_result;
							pair<vector<coord_t>,vector<coord_t> > temp_suffix;
							
							while (temp_group.size() == 0 && (singleread_group_result.size()-num_pop) >= 2) {
								temp_result = vector< vector<good_t> >(singleread_group_result.begin()+num_pop,singleread_group_result.end());
								
								temp_suffix.first = vector<coord_t>(seg_suffix_list.first.begin()+num_pop,seg_suffix_list.first.end());
								temp_suffix.second = vector<coord_t>(seg_suffix_list.second.begin(),seg_suffix_list.second.end()-num_pop);
						
									
								temp_group = check_singleread_group_result(temp_result, temp_suffix);
								
								if (temp_group.size() == 0) {
									temp_result = vector< vector<good_t> >(singleread_group_result.begin(),singleread_group_result.end()-num_pop);
									temp_suffix.first = vector<coord_t>(seg_suffix_list.first.begin(),seg_suffix_list.first.end()-num_pop);
									temp_suffix.second = vector<coord_t>(seg_suffix_list.second.begin()+num_pop,seg_suffix_list.second.end());
									
									
									temp_group = check_singleread_group_result(temp_result, temp_suffix);
									
								}
								
								
									

								num_pop++;
							}
							
							if(temp_group.size() > 0){
								singleread_group_result = temp_result;
								seg_suffix_list = temp_suffix;	
							}
							 
						}
						
						
						
						
						
						
						
						
						if (temp_group.size() > 0) {
							
							vector<coord_t>::iterator suffix_end_it = seg_suffix_list.first.begin() + temp_group.front().size();
							seg_suffix_list.first.erase(suffix_end_it,seg_suffix_list.first.end());
							
							
							
							suffix_end_it = seg_suffix_list.second.begin() + temp_group.front().size();
							seg_suffix_list.second.erase(suffix_end_it,seg_suffix_list.second.end());

							
							
							
						}
						
						
						
						
						
						
						if (temp_group.size() > 1) { 
							vector<vector< list<good_t> >::iterator> erase_list;
							list< pair<start_end_nano_t,vector< list<good_t> >::iterator> > start_end_vec;
							for (vector< list<good_t> >::iterator it = temp_group.begin(); it != temp_group.end(); it++) {
								start_end_nano_t contents;
								
								if(it->front().c > 0){
									add_good_t(contents, *it, seg_suffix_list.first); 
								} else {
									add_good_t(contents, *it, seg_suffix_list.second);
								}

								start_end_vec.push_back(pair<start_end_nano_t,vector< list<good_t> >::iterator>(contents,it));
							}
							
							
							
							list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator start_end_it  = start_end_vec.begin();
							
							
							
							map<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator> start_map;
							
							start_map.insert(pair<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator>((start_end_it->first).chr_start.front(),start_end_it));
							start_end_it++;
							
							
							while (start_end_it != start_end_vec.end()) {
								pair<map<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator>::iterator,bool> out 
								= start_map.insert(pair<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator>((start_end_it->first).chr_start.front(),start_end_it));
								
								
								
								
								
								if (!out.second && ((start_end_it->first).chr_end.back() - (start_end_it->first).chr_start.front() + 1) != (int)fullread_length[pair_idx]) {
									
									int first_len = (int)(((out.first)->second)->first).chr_start.size();
									int second_len = (int)(start_end_it->first).chr_start.size();
									
									
									
									
									
									if (first_len > second_len && second_len == 1) {
										erase_list.push_back(start_end_it->second);
									}else if (second_len > first_len && first_len == 1) {
										erase_list.push_back(((out.first)->second)->second);
										start_map.erase(out.first);
										start_map.insert(pair<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator>((start_end_it->first).chr_start.front(),start_end_it));
									}else {
										
										
										
										
										
										
										
										int first_range = ((out.first)->second)->first.read_end.back() - ((out.first)->second)->first.read_start.front();
										int second_range = start_end_it->first.read_end.back() - start_end_it->first.read_start.front();
										if (first_range > second_range) { 
											erase_list.push_back(start_end_it->second);
										}else if(first_range < second_range){
											erase_list.push_back((out.first)->second->second);
											start_map.erase(out.first);
											start_map.insert(pair<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator>((start_end_it->first).chr_start.front(),start_end_it));
										}
										
										
										
									}
									
								}
								
								start_end_it++;
							}
							
							
							
							vector<vector< list<good_t> >::iterator>::reverse_iterator erase_it = erase_list.rbegin();
							
							while (erase_it != erase_list.rend()) {
								temp_group.erase(*erase_it);
								erase_it++;
							}
							
							erase_list.clear();
							start_end_vec.clear();
							
							
							if (temp_group.size()>1) {
								
								for (vector< list<good_t> >::iterator it = temp_group.begin(); it != temp_group.end(); it++) {
									start_end_nano_t contents;
									
									if(it->front().c > 0){
										add_good_t(contents, *it, seg_suffix_list.first); 
									} else {
										add_good_t(contents, *it, seg_suffix_list.second);
									}
									
									start_end_vec.push_back(pair<start_end_nano_t,vector< list<good_t> >::iterator>(contents,it));
								}
								
								
								start_end_it  = start_end_vec.begin();
								map<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator> end_map;
								end_map.insert(pair<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator>((start_end_it->first).chr_end.back(),start_end_it));
								start_end_it++;
								
								while (start_end_it != start_end_vec.end()) {
									pair<map<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator>::iterator,bool> out 
									= end_map.insert(pair<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator>((start_end_it->first).chr_end.back(),start_end_it));
									
									
									
									
									
									if (!out.second && ((start_end_it->first).chr_end.back() - (start_end_it->first).chr_start.front() + 1) != (int)fullread_length[pair_idx]) {
										
										
										int first_len = (int)(((out.first)->second)->first).chr_end.size();
										int second_len = (int)(start_end_it->first).chr_end.size();
										
										
										
										
										
										if (first_len > second_len && second_len == 1) {
											erase_list.push_back(start_end_it->second);
										}else if (second_len > first_len && first_len == 1) {
											erase_list.push_back(((out.first)->second)->second);
											end_map.erase(out.first);
											out = end_map.insert(pair<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator>((start_end_it->first).chr_end.back(),start_end_it));
										}else {
											
											
											
											
											
											
											int first_range = ((out.first)->second)->first.read_end.back() - ((out.first)->second)->first.read_start.front();
											int second_range = start_end_it->first.read_end.back() - start_end_it->first.read_start.front();
											if (first_range > second_range) { 
												erase_list.push_back(start_end_it->second);
											}else if(first_range < second_range){
												erase_list.push_back((out.first)->second->second);
												end_map.erase(out.first);
												out = end_map.insert(pair<int,list< pair<start_end_nano_t,vector< list<good_t> >::iterator> >::iterator>((start_end_it->first).chr_end.back(),start_end_it));
											}
											
											
											
										}
										
									}
									
									start_end_it++;
								}
							}
							
							
							
							erase_it = erase_list.rbegin();
							
							while (erase_it != erase_list.rend()) {
								temp_group.erase(*erase_it);
								erase_it++;
							}
							
							
						} 
						
						
						
						if ((temp_group.size() > 0) && ((suffix_list[pair_idx].size() == 1 && temp_group.front().size() > 0) || temp_group.front().size() > 1)) {
							
							
							group_filtered_results_list = pair<vector<list<good_t> >,pair<vector<coord_t>,vector<coord_t> > >(temp_group,seg_suffix_list);
							
							
						}else {
							
						}
						
					}
					
				}

				goodpair_list[pair_idx] = group_filtered_results_list;
				
				
			}
			
			
			
			
			
			
			
			
			
			if (num_pair == 1) {  
				
				
				vector<list<good_t> >::iterator side_seg_jun_list;
				
				
				vector<vector<list<good_t> >::iterator>  matched;
				
				
				
				side_seg_jun_list = (goodpair_list[0].first).begin();
				
				
				
				while (side_seg_jun_list != (goodpair_list[0].first).end()) {
					
					int curr_num_non_extend = 0;
					int max_extend_distance = -1;
					int max_non_extend_contig = 0;
					for (list<good_t>::iterator it = side_seg_jun_list->begin(); it != side_seg_jun_list->end(); it++) {
						if (abs(it->c) < Iextend) {
							
							curr_num_non_extend++;
						}else {
							if (abs(it->d) > max_extend_distance) {
								max_extend_distance = abs(it->d);
							}
							if (curr_num_non_extend > max_non_extend_contig) {
								max_non_extend_contig = curr_num_non_extend;
							} 
							curr_num_non_extend = 0;
						}
						
					}
					if (curr_num_non_extend > max_non_extend_contig) {
						max_non_extend_contig = curr_num_non_extend;
					} 
					if (max_extend_distance == -1) {
						max_extend_distance = 50;
					}
					
					
					bool all_extends = !((max_non_extend_contig>=2)||(((uint_fast32_t)max_extend_distance)>(25+min_extend_req) && max_non_extend_contig >= 1)||(side_seg_jun_list->size()<=2 && max_non_extend_contig >= 1));
					
					
					
					
					if(!all_extends){
						
						matched.push_back( side_seg_jun_list );
						
						
						
						
						
						
						
						
						
					}

					(side_seg_jun_list)++;
				}
				
				
				
				
				
				
				int num_multi_jun = 0;
				
				vector<vector<list<good_t> >::iterator>::iterator match_it = matched.begin();
				while (match_it != matched.end()) {
					list<good_t>::iterator list_it = (*match_it)->begin();
					while (list_it != (*match_it)->end()) {
						if (abs(list_it->c) == Mjunnegative || abs(list_it->c) == Mjunpositive) {
							num_multi_jun++;
						}
						list_it++;
					}
					match_it++;
				}
				
				match_it = matched.begin();
				while (match_it != matched.end()) {
					bool add_jun = true;
					
					if (!DNA_mode && num_multi_jun > 1) {

						list<good_t>::iterator list_it = (*match_it)->begin();
						while (list_it != (*match_it)->end() && add_jun) {
							if (abs(list_it->c) == Mjunnegative || abs(list_it->c) == Mjunpositive) {
								add_jun = false;
							}
							list_it++;
						}
						
					}
					
					if(add_jun){
						
						
						
						
						vector<coord_t> seg_suffix_list;
						
						
						
						if ((*match_it)->front().c > 0) {
							seg_suffix_list = goodpair_list[0].second.first;
						}else {
							seg_suffix_list = goodpair_list[0].second.second;
						}
						
						
						
						
						
						
						
						
						
						start_end_nano_t curr_start_end;
						add_good_t(curr_start_end,*(*match_it),seg_suffix_list);
						
						string mismatch_str;
						string read_out;
						pair<pair<int,int>,string> cigar = roll_cigar(curr_start_end,fullread_length[0],cufflinks);
						
						if (curr_start_end.direction) {
							mismatch_str = count_mismatch(curr_start_end,full_line_read[0],chrline,num_read_mismatch,max_clip_allowed);
							
							if(mismatch_str.length() > 0){
								if (!cufflinks){
									read_out = full_line_read[0] + "\t";   
								}else {
									read_out = full_line_read[0].substr(cigar.first.first, fullread_length[0] - cigar.first.first - cigar.first.second) + "\t";   
								}
							}
							
						}else {
							temp = full_line_read[0];
							compleseq(temp);
							mismatch_str = count_mismatch(curr_start_end,temp,chrline,num_read_mismatch,max_clip_allowed);
							
							if(mismatch_str.length() > 0){
								if (!cufflinks){
									read_out = temp + "\t";  
								}else {
									read_out = temp.substr(cigar.first.first, fullread_length[0] - cigar.first.first - cigar.first.second) + "\t"; 
								}
							}
						} 
						
						
						
						
						if(mismatch_str.length() > 0){
							
							good_SAM_file << full_line_name[0] << '\t';
							
							
							int flag = 0;
							if (!curr_start_end.direction) {
								flag += 16;
							}
							
							good_SAM_file << flag << '\t';
							
							
							
							
							good_SAM_file << chr_name << "\t";
							
							
							good_SAM_file << (curr_start_end.chr_start).front() << "\t";
							
							
							good_SAM_file << "255" << "\t";
							
							
							
							good_SAM_file << cigar.second << "\t";
							
							
							good_SAM_file << "*" << "\t";  
							
							
							good_SAM_file << 0 << "\t";  
							
							
							good_SAM_file << "0" << "\t";  
							
							
							good_SAM_file << read_out;
							

							
							if(full_line_qual[0].length() == 0){ 
								if (!cufflinks){
									good_SAM_file << string(fullread_length[0],'I') ;
								}else {
									good_SAM_file << string(fullread_length[0],'I').substr(cigar.first.first, fullread_length[0] - cigar.first.first - cigar.first.second) ;
								}
							}else {
								
								if (!cufflinks){
									good_SAM_file << full_line_qual[0];   
								}else {
									good_SAM_file << full_line_qual[0].substr(cigar.first.first, fullread_length[0] - cigar.first.first - cigar.first.second);   
								}
								
							}
							
							
							if (num_multi_jun > 0){
								good_SAM_file << '\t' << "XM:i:" << num_multi_jun;
							}
							
							
							if (curr_start_end.chr_start.size() > 1){
								if (curr_start_end.strand) {
									good_SAM_file << '\t' << "XS:A:+";
								}else{
									good_SAM_file << '\t' << "XS:A:-";
								}
								
							}
							
							good_SAM_file << mismatch_str;
							
							good_SAM_file << '\n';
							
						}
						
					}
					
					match_it++;
				}
				

				
			}else {
				

				
				
				vector<list<good_t> >::iterator side_seg_jun_list[2];
				
				
				
				
				if (goodpair_list[0].first.size() != 0 && goodpair_list[1].first.size() != 0) {
					
					
					
					vector<pair<vector<list<good_t> >::iterator,vector<list<good_t> >::iterator> > matched_pairs;
					
					
					
					
					
					
					
					side_seg_jun_list[0] = (goodpair_list[0].first).begin();
					
					
					
					while (side_seg_jun_list[0] != (goodpair_list[0].first).end()) {
						

						
						side_seg_jun_list[1] = (goodpair_list[1].first).begin();
						
						
						while (side_seg_jun_list[1] != (goodpair_list[1].first).end()) {
							
							
							
							
							
							if (checkPairInfo(*(side_seg_jun_list[0]), *(side_seg_jun_list[1]), intron_curve[2])) {
								
								
								
								
								
								
								
								int left_curr_num_non_extend = 0;
								int left_max_extend_distance = -1;
								int left_max_non_extend_contig = 0;
								for (list<good_t>::iterator it = side_seg_jun_list[0]->begin(); it != side_seg_jun_list[0]->end(); it++) {
									if (abs(it->c) < Iextend) {
										
										left_curr_num_non_extend++;
									}else {
										if (abs(it->d) > left_max_extend_distance) {
											left_max_extend_distance = abs(it->d);
										}
										if (left_curr_num_non_extend > left_max_non_extend_contig) {
											left_max_non_extend_contig = left_curr_num_non_extend;
										} 
										left_curr_num_non_extend = 0;
									}
									
								}
								if (left_curr_num_non_extend > left_max_non_extend_contig) {
									left_max_non_extend_contig = left_curr_num_non_extend;
								} 
								if (left_max_extend_distance == -1) {
									left_max_extend_distance = 50;
								}
								
								
								bool left_all_extends = !((left_max_non_extend_contig>=2)||(((uint_fast32_t)left_max_extend_distance)>(25+min_extend_req) && left_max_non_extend_contig >= 1)||(side_seg_jun_list[0]->size()<=2 && left_max_non_extend_contig >= 1));
								
								
								
								int right_curr_num_non_extend = 0;
								int right_max_extend_distance = -1;
								int right_max_non_extend_contig = 0;
								for (list<good_t>::iterator it = side_seg_jun_list[1]->begin(); it != side_seg_jun_list[1]->end(); it++) {
									if (abs(it->c) < Iextend) {
										
										right_curr_num_non_extend++;
									}else {
										if (abs(it->d) > right_max_extend_distance) {
											right_max_extend_distance = abs(it->d);
										}
										if (right_curr_num_non_extend > right_max_non_extend_contig) {
											right_max_non_extend_contig = right_curr_num_non_extend;
										} 
										right_curr_num_non_extend = 0;
									}
								}
								if (right_curr_num_non_extend > right_max_non_extend_contig) {
									right_max_non_extend_contig = right_curr_num_non_extend;
								} 
								if (right_max_extend_distance == -1) {
									right_max_extend_distance = 50;
								}
								
								
								bool right_all_extends = !((right_max_non_extend_contig>=2) || (((uint_fast32_t)right_max_extend_distance)>(25+min_extend_req) && right_max_non_extend_contig>=1)||(side_seg_jun_list[1]->size()<=2 && right_max_non_extend_contig >= 1));
								
								
								
								
								
								
								if((!left_all_extends && !right_all_extends)||(side_seg_jun_list[0]->size() >= 3 && !right_all_extends)||(!left_all_extends && side_seg_jun_list[1]->size()>=3)){

									matched_pairs.push_back(pair<vector<list<good_t> >::iterator,vector<list<good_t> >::iterator >( side_seg_jun_list[0],side_seg_jun_list[1]));

								}else {
									 if(!left_all_extends){ 

										 matched_pairs.push_back(pair<vector<list<good_t> >::iterator,vector<list<good_t> >::iterator >( side_seg_jun_list[0],goodpair_list[1].first.end()));

									 }else if(!right_all_extends){ 

										 matched_pairs.push_back(pair<vector<list<good_t> >::iterator,vector<list<good_t> >::iterator >( (goodpair_list[0].first).end(),side_seg_jun_list[1]));
										 

									 }
									
								}

							}
							
							
							
							(side_seg_jun_list[1])++;
							
						}
						
						(side_seg_jun_list[0])++;
						
					}
					
					
					
					
					
					
					
					int num_multi_jun[2] = {0,0};
					
					vector<pair<vector<list<good_t> >::iterator,vector<list<good_t> >::iterator> >::iterator match_it = matched_pairs.begin();
					while (match_it != matched_pairs.end()) {
						
						if(match_it->first != (goodpair_list[0].first).end()){
							list<good_t>::iterator list_it = (match_it->first)->begin();
							while (list_it != (match_it->first)->end()) {
								if (abs(list_it->c) == Mjunnegative || abs(list_it->c) == Mjunpositive) {
									num_multi_jun[0]++;
								}
								list_it++;
							}
						}
						
						
						if(match_it->second != (goodpair_list[1].first).end()){
							list<good_t>::iterator list_it = (match_it->second)->begin();
							while (list_it != (match_it->second)->end()) {
								if (abs(list_it->c) == Mjunnegative || abs(list_it->c) == Mjunpositive) {
									num_multi_jun[1]++;
								}
								list_it++;
							}
						}
						
						match_it++;
					}
					
					
					
					match_it = matched_pairs.begin();
					while (match_it != matched_pairs.end()) {
						
						bool add_jun[2] = {true,true};
						
						if (!DNA_mode && num_multi_jun[0] > 1 && match_it->first != (goodpair_list[0].first).end()) {
							
							list<good_t>::iterator list_it = (match_it->first)->begin();
							while (list_it != (match_it->first)->end() && add_jun[0]) {
								if (abs(list_it->c) == Mjunnegative || abs(list_it->c) == Mjunpositive) {
									add_jun[0] = false;
								}
								list_it++;
							}
						}
						if (!DNA_mode && num_multi_jun[1] > 1 && match_it->second != (goodpair_list[1].first).end()){
							list<good_t>::iterator list_it = (match_it->second)->begin();
							while (list_it != (match_it->second)->end() && add_jun[1]) {
								if (abs(list_it->c) == Mjunnegative || abs(list_it->c) == Mjunpositive) {
									add_jun[1] = false;
								}
								list_it++;
							}
						}
						
						
						vector<coord_t> seg_suffix_list[2];
						start_end_nano_t curr_start_end[2];
						string mismatch_str[2] = {"",""};
						string read_out[2] = {"",""};
						pair<pair<int,int>,string> cigar[2];
						
						if(match_it->first != (goodpair_list[0].first).end() && match_it->second != (goodpair_list[1].first).end() && add_jun[0] && add_jun[1]){
							
							
							
							
							
							if (match_it->first->front().c > 0) {
								seg_suffix_list[0] = goodpair_list[0].second.first;
							}else {
								seg_suffix_list[0] = goodpair_list[0].second.second;
							}
							
							if (match_it->second->front().c > 0) {
								seg_suffix_list[1] = goodpair_list[1].second.first;
							}else {
								seg_suffix_list[1] = goodpair_list[1].second.second;
							}
							
							
							
						
							add_good_t(curr_start_end[0],*(match_it->first),seg_suffix_list[0]);
							add_good_t(curr_start_end[1],*(match_it->second),seg_suffix_list[1]);
							
							cigar[0] = roll_cigar(curr_start_end[0],fullread_length[0],cufflinks);
							cigar[1] = roll_cigar(curr_start_end[1],fullread_length[1],cufflinks);
							
							if (curr_start_end[0].direction) {
								mismatch_str[0] = count_mismatch(curr_start_end[0],full_line_read[0],chrline,num_read_mismatch,max_clip_allowed);
								
								temp = full_line_read[1];
								compleseq(temp);
								mismatch_str[1] = count_mismatch(curr_start_end[1],temp,chrline,num_read_mismatch,max_clip_allowed);
								
								if (!cufflinks){
									read_out[0] = full_line_read[0] + "\t";   
									read_out[1] = temp + "\t";   
								}else {
									read_out[0] = full_line_read[0].substr(cigar[0].first.first, fullread_length[0] - cigar[0].first.first - cigar[0].first.second) + "\t";   
									read_out[1] = temp.substr(cigar[1].first.first, fullread_length[1] - cigar[1].first.first - cigar[1].first.second) + "\t";   
								}
								
							}else {
								temp = full_line_read[0];
								compleseq(temp);
								mismatch_str[0] = count_mismatch(curr_start_end[0],temp,chrline,num_read_mismatch,max_clip_allowed);
								
								mismatch_str[1] = count_mismatch(curr_start_end[1],full_line_read[1],chrline,num_read_mismatch,max_clip_allowed);
								
								if (!cufflinks){
									read_out[0] = temp + "\t";  
									read_out[1] = full_line_read[1] + "\t";  
								}else {
									read_out[0] = temp.substr(cigar[0].first.first, fullread_length[0] - cigar[0].first.first - cigar[0].first.second) + "\t"; 
									read_out[1] = full_line_read[1].substr(cigar[1].first.first, fullread_length[1] - cigar[1].first.first - cigar[1].first.second) + "\t"; 
								}
							} 
							

						}else if (match_it->first != (goodpair_list[0].first).end() && add_jun[0]){
							
							
							
							
							if (match_it->first->front().c > 0) {
								seg_suffix_list[0] = goodpair_list[0].second.first;
							}else {
								seg_suffix_list[0] = goodpair_list[0].second.second;
							}
							
							
							
							
							
							
							add_good_t(curr_start_end[0],*(match_it->first),seg_suffix_list[0]);
							
							cigar[0] = roll_cigar(curr_start_end[0],fullread_length[0],cufflinks);
							
							if (curr_start_end[0].direction) {
								mismatch_str[0] = count_mismatch(curr_start_end[0],full_line_read[0],chrline,num_read_mismatch,max_clip_allowed);
								
								if (!cufflinks){
									read_out[0] = full_line_read[0] + "\t";   
								}else {
									read_out[0] = full_line_read[0].substr(cigar[0].first.first, fullread_length[0] - cigar[0].first.first - cigar[0].first.second) + "\t";   
								}
								
							}else {
								temp = full_line_read[0];
								compleseq(temp);
								mismatch_str[0] = count_mismatch(curr_start_end[0],temp,chrline,num_read_mismatch,max_clip_allowed);
								
								if (!cufflinks){
									read_out[0] = temp + "\t";  
								}else {
									read_out[0] = temp.substr(cigar[0].first.first, fullread_length[0] - cigar[0].first.first - cigar[0].first.second) + "\t"; 
								}
							} 
							
							
						}else if (match_it->second != (goodpair_list[1].first).end() && add_jun[1]){
							
							
							
							
							if (match_it->second->front().c > 0) {
								seg_suffix_list[1] = goodpair_list[1].second.first;
							}else {
								seg_suffix_list[1] = goodpair_list[1].second.second;
							}
							
							
							
							
							
							
							add_good_t(curr_start_end[1],*(match_it->second),seg_suffix_list[1]);

							cigar[1] = roll_cigar(curr_start_end[1],fullread_length[1],cufflinks);
							
							if (curr_start_end[1].direction) {
								mismatch_str[1] = count_mismatch(curr_start_end[1],full_line_read[1],chrline,num_read_mismatch,max_clip_allowed);
								
								
								if (!cufflinks){
									read_out[1] = full_line_read[1] + "\t";   
								}else {
									read_out[1] = full_line_read[1].substr(cigar[1].first.first, fullread_length[1] - cigar[1].first.first - cigar[1].first.second) + "\t";   
								}
								
							}else {
								
								temp = full_line_read[1];
								compleseq(temp);
								mismatch_str[1] = count_mismatch(curr_start_end[1],temp,chrline,num_read_mismatch,max_clip_allowed);
								if (!cufflinks){
									read_out[1] = temp + "\t";  
								}else {
									read_out[1] = temp.substr(cigar[1].first.first, fullread_length[1]  - cigar[1].first.first - cigar[1].first.second) + "\t"; 
								}
							} 
							
						}
						
						
						
						
						
						if(mismatch_str[0].length() > 0 && mismatch_str[1].length() > 0){
							
							int insert_size = 0;
							if (curr_start_end[0].direction) { 
								insert_size = curr_start_end[1].chr_start.front() - curr_start_end[0].chr_end.back();
							}else {
								insert_size = curr_start_end[0].chr_start.front() - curr_start_end[1].chr_end.back();
							}

							
							
							
							
							good_SAM_file << full_line_name[0] << '\t';
							
							
							int flag = 1; 
							if (!curr_start_end[0].direction) { 
								flag += 16;
								flag += 32; 
							}
							flag += 64; 
							
							
							good_SAM_file << flag << '\t';
							
							
							
							
							good_SAM_file << chr_name << "\t";
							
							
							good_SAM_file << (curr_start_end[0].chr_start).front() << "\t";
							
							
							good_SAM_file << "255" << "\t";
							
							
							
							good_SAM_file << cigar[0].second << "\t";
							
							
							good_SAM_file  << "=\t"; 
							
							
							good_SAM_file << (curr_start_end[1].chr_start).front() << "\t";
							
							
							good_SAM_file << IntToStr(insert_size) << "\t";  
							
							

							good_SAM_file << read_out[0];

							
							if(full_line_qual[0].length() == 0){ 
								if (!cufflinks){
									good_SAM_file << string(fullread_length[0],'I') ;
								}else {
									good_SAM_file << string(fullread_length[0],'I').substr(cigar[0].first.first, fullread_length[0] - cigar[0].first.first - cigar[0].first.second) ;
								}
							}else {
								
								if (!cufflinks){
									good_SAM_file << full_line_qual[0];   
								}else {
									good_SAM_file << full_line_qual[0].substr(cigar[0].first.first, fullread_length[0] - cigar[0].first.first - cigar[0].first.second);   
								}
								
							}
							
							
							if (num_multi_jun[0] > 0){
								good_SAM_file << '\t' << "XM:i:" << num_multi_jun[0];
							}
							
							
							if (curr_start_end[0].chr_start.size() > 1){
								if (curr_start_end[0].strand) {
									good_SAM_file << '\t' << "XS:A:+";
								}else{
									good_SAM_file << '\t' << "XS:A:-";
								}
								
							}
							
							good_SAM_file << mismatch_str[0];
							
							good_SAM_file << '\n';
							
							
							
							
							
							
							
							
							
							
							
							
							
							good_SAM_file << full_line_name[1] << '\t';
							
							
							flag = 1; 
							if (!curr_start_end[0].direction) { 
								flag += 16;
								flag += 32; 
							}
							flag += 128; 
							
							
							good_SAM_file << flag << '\t';
							
							
							
							
							good_SAM_file << chr_name << "\t";
							
							
							good_SAM_file << (curr_start_end[1].chr_start).front() << "\t";
							
							
							good_SAM_file << "255" << "\t";
							
							
							
							good_SAM_file << cigar[1].second << "\t";
							
							
							good_SAM_file  << "=\t"; 
							
							
							good_SAM_file << (curr_start_end[0].chr_start).front() << "\t";
							
							
							good_SAM_file <<  IntToStr(insert_size)  << "\t"; 
							
							
	
							good_SAM_file << read_out[1];
							
							
							
							if(full_line_qual[1].length() == 0){ 
								if (!cufflinks){
									good_SAM_file << string(fullread_length[1],'I') ;
								}else {
									good_SAM_file << string(fullread_length[1],'I').substr(cigar[1].first.first, fullread_length[1] - cigar[1].first.first - cigar[1].first.second) ;
								}
							}else {
								
								if (!cufflinks){
									good_SAM_file << full_line_qual[1];   
								}else {
									good_SAM_file << full_line_qual[1].substr(cigar[1].first.first, fullread_length[1] - cigar[1].first.first - cigar[1].first.second);   
								}
								
							}
							
							
							if (num_multi_jun[1] > 0){
								good_SAM_file << '\t' << "XM:i:" << num_multi_jun[1];
							}
							
							
							if (curr_start_end[1].chr_start.size() > 1){
								if (curr_start_end[1].strand) {
									good_SAM_file << '\t' << "XS:A:+";
								}else{
									good_SAM_file << '\t' << "XS:A:-";
								}
								
							}
							
							good_SAM_file << mismatch_str[1];
							
							good_SAM_file << '\n';
							
							
							
							
							
							
						}else if (mismatch_str[0].length() > 0){
							
							
							
							good_SAM_file << full_line_name[0] << '\t';
							
							
							int flag = 1; 
							flag += 8;    
							if (!curr_start_end[0].direction) { 
								flag += 16;
								flag += 32; 
							}
							flag += 64; 
							
							
							good_SAM_file << flag << '\t';
							
							
							
							
							good_SAM_file << chr_name << "\t";
							
							
							good_SAM_file << (curr_start_end[0].chr_start).front() << "\t";
							
							
							good_SAM_file << "255" << "\t";
							
							
							
							good_SAM_file << cigar[0].second << "\t";
							
							
							good_SAM_file << "*" << "\t";  
							
							
							good_SAM_file << 0 << "\t";  
							
							
							good_SAM_file << "0" << "\t";  
							
							
							good_SAM_file << read_out[0];
							
							
							
							if(full_line_qual[0].length() == 0){ 
								if (!cufflinks){
									good_SAM_file << string(fullread_length[0],'I') ;
								}else {
									good_SAM_file << string(fullread_length[0],'I').substr(cigar[0].first.first, fullread_length[0] - cigar[0].first.first - cigar[0].first.second) ;
								}
							}else {
								
								if (!cufflinks){
									good_SAM_file << full_line_qual[0];   
								}else {
									good_SAM_file << full_line_qual[0].substr(cigar[0].first.first, fullread_length[0] - cigar[0].first.first - cigar[0].first.second);   
								}
								
							}
							
							
							if (num_multi_jun[0] > 0){
								good_SAM_file << '\t' << "XM:i:" << num_multi_jun[0];
							}
							
							
							if (curr_start_end[0].chr_start.size() > 1){
								if (curr_start_end[0].strand) {
									good_SAM_file << '\t' << "XS:A:+";
								}else{
									good_SAM_file << '\t' << "XS:A:-";
								}
								
							}
							
							good_SAM_file << mismatch_str[0];
							
							good_SAM_file << '\n';
						
						
						}else if (mismatch_str[1].length() > 0){
							
							
							
							good_SAM_file << full_line_name[1] << '\t';
							
							
							int flag = 1; 
							flag += 8;    
							if (curr_start_end[1].direction) { 
								flag += 16;
								flag += 32; 
							}
							flag += 128; 
							
							
							good_SAM_file << flag << '\t';
							
							
							
							
							good_SAM_file << chr_name << "\t";
							
							
							good_SAM_file << (curr_start_end[1].chr_start).front() << "\t";
							
							
							good_SAM_file << "255" << "\t";
							
							
							
							good_SAM_file << cigar[1].second << "\t";
							
							
							good_SAM_file << "*" << "\t";  
							
							
							good_SAM_file << 0 << "\t";  
							
							
							good_SAM_file << "0" << "\t";  
							
							
							
							good_SAM_file << read_out[1];
							
							
							
							if(full_line_qual[1].length() == 0){ 
								if (!cufflinks){
									good_SAM_file << string(fullread_length[1],'I') ;
								}else {
									good_SAM_file << string(fullread_length[1],'I').substr(cigar[1].first.first, fullread_length[1]  - cigar[1].first.first - cigar[1].first.second) ;
								}
							}else {
								
								if (!cufflinks){
									good_SAM_file << full_line_qual[1];   
								}else {
									good_SAM_file << full_line_qual[1].substr(cigar[1].first.first, fullread_length[1]  - cigar[1].first.first - cigar[1].first.second);   
								}
								
							}
							
							
							if (num_multi_jun[1] > 0){
								good_SAM_file << '\t' << "XM:i:" << num_multi_jun[1];
							}
							
							
							if (curr_start_end[1].chr_start.size() > 1){
								if (curr_start_end[1].strand) {
									good_SAM_file << '\t' << "XS:A:+";
								}else{
									good_SAM_file << '\t' << "XS:A:-";
								}
								
							}
							
							good_SAM_file << mismatch_str[1];
							
							
							good_SAM_file << '\n';
							
						}
						
						
						
						match_it++;
						
					}

				}
				
			}
			
			
			
			
			
			
			
			curr_read_idx++;
		}
		
		
		for (int k = 0;k<num_pair;k++){
			read_index_tfile[k].close();
			read_index_tfile[k].clear();
			
			fullread[k].close();
			fullread[k].clear();
			
			fullqual[k].close();
			fullqual[k].clear();
			
			fullname[k].close();
			fullname[k].clear();
		}
		

	}
	
	
	gettimeofday(&tv, NULL);
	
	debug_out<<"SpliceMapping time: "<< diffclock(temp_start_tv,tv) <<" s."<<endl;
	 
	

	
	
	for (uint_fast32_t i = 0; i<MAX_HASH; i++) {
		if (chrdict_10[i] != NULL) {
			delete chrdict_10[i];
		}
	}
	
	
	delete [] chrdict_10;
	
	 
	
	
	gettimeofday(&tv, NULL);
	
	
	cout<<"Total " << chr_name << " execution time: "<<diffclock(start_tv,tv)<<" s."<<endl;
	debug_out << "______________________" << endl;
	debug_out<<"Total " << chr_name << " execution time: "<<diffclock(start_tv,tv)<<" s."<<endl;

	
	debug_out.close();
	debug_out.clear();
	
    return 0;
}








vector<good_t> combine_result_list(vector<good_t> &l1, vector<good_t> &l2)
{
	vector<good_t> result;
	for (vector<good_t>::iterator it = l1.begin(); it != l1.end(); it++) {
		bool exists = false;
		for (vector<good_t>:: iterator it2 = l2.begin(); it2 != l2.end(); it2++) {
			if (isgoodtequal(*it, *it2)) {
				exists = true;
			}
		}
		if (!exists) {
			result.push_back(*it);
		}
	}
	for (vector<good_t>:: iterator it2 = l2.begin(); it2 != l2.end(); it2++) {
		result.push_back(*it2);
	}
	
	return result;
}

good_vec_t checkgoodlists(good_vec_t &good1, good_vec_t &good2)
{
	if (good1.size() > 0 && good2.size() > 0) {
		return combine_result_list(good1, good2);
	}else if(good1.size() > 0){
		return good1;
	}else {
		return good2;
	}
}


inline bool checkPairInfo(list<good_t> &left_seg_jun_list, list<good_t> &right_seg_jun_list, int distance)
{
		
	if (left_seg_jun_list.front().c > 0 && right_seg_jun_list.front().c < 0) {
		coord_t left_range = junction2boundary(left_seg_jun_list.back());
		coord_t right_range = junction2boundary(right_seg_jun_list.front());
		
		if (right_range.first < (distance + left_range.second) && right_range.first > left_range.first) {
			return true;
		}
		
	}else if (left_seg_jun_list.front().c < 0 && right_seg_jun_list.front().c > 0) {
		coord_t left_range = junction2boundary(left_seg_jun_list.front());
		coord_t right_range = junction2boundary(right_seg_jun_list.back());
		
		if (left_range.first < (distance + right_range.second) && left_range.first > right_range.first) {
			return true;
		}
	}
	
	return false;
}




inline bool isgoodtequal(good_t &first, good_t &second)
{
	if (first.a == second.a && first.b == second.b && first.c == second.c && first.d == second.d) {
		return true;
	}else {
		return false;
	}

}


inline pair<vector<uint_fast32_t>::iterator, vector<uint_fast32_t>::iterator> checkindexlistFwd(uint_fast32_t ref,vector<uint_fast32_t> &indexlist,uint_fast32_t range1, uint_fast32_t range2)
{
	pair<vector<uint_fast32_t>::iterator, vector<uint_fast32_t>::iterator> result;
	uint_fast32_t uplim = ref+range1;
	uint_fast32_t lowlim = ref+range2+25;
	
	

	result.first =  lower_bound(indexlist.begin(),indexlist.end(),lowlim+1);
	result.second = lower_bound(indexlist.begin(),indexlist.end(),uplim);
	
	return result;
}



inline pair<vector<uint_fast32_t>::iterator, vector<uint_fast32_t>::iterator> checkindexlistRes(uint_fast32_t ref,vector<uint_fast32_t> &indexlist,uint_fast32_t range1, uint_fast32_t range2)
{
	pair<vector<uint_fast32_t>::iterator, vector<uint_fast32_t>::iterator> result;
	uint_fast32_t uplim;
	if (ref>range2) {
		uplim = ref-range2;
	}else {
		uplim = 0;
	}
	
	
	uint_fast32_t lowlim;
	if (ref>range1) {
		lowlim = ref-range1;
	}else {
		lowlim = 0;
	}

	
	
	

	result.first =  lower_bound(indexlist.begin(),indexlist.end(),lowlim+1);
	
	
	
	
	
	result.second = lower_bound(indexlist.begin(),indexlist.end(),uplim);

	return result;
}


inline uint_fast32_t cal_half_distance(uint_fast32_t length, uint_fast32_t intron_curve[3])
{
	int result = intron_curve[0];
	if (length > 9 && length < 16) {
		result = intron_curve[0];
	}else if (length >= 16 && length < 21) {
		result = intron_curve[1];
	}else if (length >= 21) {
		result = intron_curve[2];
	}
	
	return result;
}







inline void FindExonicRead(good_vec_t& good, list<line_t> line_list[2])
{
	list<line_t>::iterator it1 = line_list[0].begin();
	
	while (it1 != line_list[0].end()) {
		list<line_t>::iterator it2 = line_list[1].begin();
		int pos1 = it1->coor;
		int dir1 = it1->direction;
		while (it2 != line_list[1].end()) {
			int pos2 = it2->coor;
			int dir2 = it2->direction;
			
			if (dir1 && dir2 ) {  
				if (pos2 - pos1 == 25) {
					
					good_t exon = {pos1,pos1+49,Iexonic,0};
					good.push_back(exon);
				}
			}else if ((!dir1)  && (!dir2) ) { 
				if (pos1 - pos2 == 25) {
					
					good_t exon = {pos2,pos2+49,-Iexonic,0};
					good.push_back(exon);
				}
			}
			it2++;
		}
		it1++;
	}
}


void InsertGoodExtend(good_vec_t &good_list,good_vec_t &good_exonic_list,good_t contents)
{	
	bool good_extend = true;
	
	vector<good_t>::iterator cur_list = good_exonic_list.begin();
	while (cur_list != good_exonic_list.end()) {
		
		if ((contents.a == cur_list->a || contents.b == cur_list->b) && sameSign(contents.c, cur_list->c)) {
			good_extend = false;
		}
		cur_list++;
	}
	
	
	if (good_extend) {
		good_list.push_back(contents);
	}
}


bool compare_line_t(line_t &first, line_t &second)
{
	return second.coor > first.coor;
}

inline void DelLocalHalfMappableHits(list<line_t> &line_list, uint_fast32_t distance)
{
	
	line_list.sort(compare_line_t);
	
	list<line_t>::iterator it = line_list.begin();
	list<list<line_t>::iterator> erase_list;
	list<list<line_t>::iterator>::iterator erase_it;
	set<uint_fast32_t> erased_coors;
	
	uint_fast32_t prev_val = it->coor;
	list<line_t>::iterator prev_it = it;
	it++;
	
	
	
	while (it != line_list.end()) {
		
		if ((it->coor - prev_val) < distance) {
			erase_list.push_back(it);
			erased_coors.insert(it->coor);
			
			pair<set<uint_fast32_t>::iterator,bool> out = erased_coors.insert(prev_val);
			
			if (out.second){
				erase_list.push_back(prev_it);
			}
		}
		prev_val = it->coor;
		prev_it = it;
		it++;
	}
	
	
	erase_it = erase_list.begin();
	
	while (erase_it != erase_list.end()) {
	    
		line_list.erase(*erase_it);
		erase_it++;
		
	}
	
}

 

inline int get_base_val(char c)
{
	int val = -1;
	switch (c) {
		case 'A':
			val = 0;
			break;
		case 'C':
			val = 1;
			break;
		case 'G':
			val = 2;
			break;
		case 'T':
			val = 3;
			break;
		default:
			break;
	}
	
	return val;
}

string get_filename(const string& path)
{
	string::size_type lastslash = path.rfind('/');
	if (lastslash == string::npos) {
		lastslash = path.rfind('\\');
	}
	
	if (lastslash == string::npos) {
		
		return path;
	}
	
	return path.substr(lastslash+1,path.length()-lastslash);
}






inline int base2int(string &s)
{
	int out = 0;
	for (int i=0;i<10;i++){
		int val = get_base_val(s[i]);
		if (val < 0){
			return -1;  
		}

		
		out = out + (val<<(2*i));
	}
	return out;
}




vector< list<good_t> > check_singleread_group_result(vector<vector<good_t> > &singleread_group_result, pair<vector<coord_t>,vector<coord_t> > &range_list)
{
	
	
	vector< list<good_t> > all_output;
	map<int,map<int,list<int> > > result;
	int ref_group_index = 0; 
	map<int,int> ref_index_list;
	
	int length_segment = (int)singleread_group_result.size();
	
	if (length_segment == 1) {
		for (vector<good_t>::iterator it = singleread_group_result.front().begin(); it != singleread_group_result.front().end(); it++) {
			list<good_t> contents;
			contents.push_back(*it);
			all_output.push_back(contents);
		}
		
		return all_output;  
	}
	
	while (ref_group_index < length_segment) {

		if (ref_group_index == 0) {
			int len_list = (int)singleread_group_result[ref_group_index].size();
			for (int i = 0; i<len_list; i++) {
				ref_index_list.insert(pair<int,int>(i,1));
			}
		}else {
			pair<map<int,map<int,list<int> > >::iterator,bool> result_ref = result.insert(pair<int,map<int,list<int> > >(ref_group_index,map<int,list<int> >()));
			map<int,int> temp_ref_index_list;
			int index = 0;
			
			
			
			for (vector<good_t>::iterator junction = singleread_group_result[ref_group_index].begin()
				 ; junction != singleread_group_result[ref_group_index].end(); junction++) {
				
				pair<map<int,list<int> >::iterator,bool> result_ref_index = ((result_ref.first)->second).insert(pair<int,list<int> >(index,list<int>()));
				
				
				
				for (map<int,int>::iterator ref_index = ref_index_list.begin(); ref_index!=ref_index_list.end(); ref_index++) {
					
					bool check_overlap;
					
					if(junction->c > 0){
						check_overlap = check_segment_overlap(singleread_group_result[ref_group_index-1][ref_index->first],
															   *junction, (range_list.first)[ref_group_index-1], (range_list.first)[ref_group_index]);
					}else {
						check_overlap = check_segment_overlap(singleread_group_result[ref_group_index-1][ref_index->first],
															  *junction, (range_list.second)[ref_group_index], (range_list.second)[ref_group_index-1]);
					}

					if (check_overlap) {
						
						((result_ref_index.first)->second).push_back(ref_index->first);
						
						if(temp_ref_index_list.find(index) == temp_ref_index_list.end()){
							temp_ref_index_list.insert(pair<int,int>(index,1));
						}
					}
					
				}
				
				index++;
			}
			
			
			if (temp_ref_index_list.size() == 0) {
				break;
			}
			
			ref_index_list = temp_ref_index_list;
		}
		ref_group_index++;
	}
	
	
	
	
	ref_group_index = ref_group_index - 1;
	
	
	
	
	for (map<int,list<int> >::iterator index = result[ref_group_index].begin(); 
		 index != result[ref_group_index].end(); index++) {
		
		for (list<int>::iterator good_index = (index->second).begin(); good_index != (index->second).end(); good_index++) {
			list<good_t> one_output;
			one_output.push_back(singleread_group_result[ref_group_index][index->first]);
			int new_ref_group_index = ref_group_index - 1;
			
			
			
			
			
			recursive_record(singleread_group_result, result, new_ref_group_index, *good_index, one_output, all_output);
			
			
			 
			
		}
		
	}
	
	

	return all_output;
}


void recursive_record(vector<vector<good_t> > &singleread_group_result, map<int,map<int,list<int> > > &result
					  , int ref_group_index, int good_index, list<good_t> &one_output,vector< list<good_t> > &all_output )
{
	if (ref_group_index == 0) {
		one_output.push_back(singleread_group_result[ref_group_index][good_index]);
		all_output.push_back(one_output);
		
		one_output.pop_back();
	}else {
		one_output.push_back(singleread_group_result[ref_group_index][good_index]);
		int new_ref_group_index = ref_group_index - 1;
		int index = good_index;
		
		
		list<int>::iterator list_start = result[ref_group_index][index].begin();
		list<int>::iterator list_end = result[ref_group_index][index].end();
		
		
		
		for (list<int>::iterator list_it = list_start; list_it != list_end; list_it++) {
			
			
			recursive_record(singleread_group_result,result,new_ref_group_index,*list_it,one_output,all_output);
		}
		one_output.pop_back();
	}

}



bool check_segment_overlap(good_t &junctiona, good_t &junctionb,coord_t &lena, coord_t &lenb)
{
	
	
	
	
	
	
	
	const int EXON_BUFFER = 5;
	
	int local_a = 0;
	int local_b = 0;
	
	if (junctiona.c < 0) {
		
		
		local_a = 1;
		local_b = lena.second - lenb.first + 1;
	}else {
		local_a = lenb.first - lena.first + 1;
		local_b = lena.second - lena.first + 1;
	}
	
	

	
	good_t juna_good_seq = good_seq(junctiona, local_a, local_b);  
	
	
	
	if (junctionb.c < 0) {
		
		
		
		local_a = lenb.first - lena.first + 1;
		local_b = lena.second - lena.first + 1;
	}else {
		
		local_a = 1;
		local_b = lena.second - lenb.first + 1;
	}
	
	
	
	good_t junb_good_seq = good_seq(junctionb, local_a, local_b);  
	
	
	

	int typea = abs(junctiona.c);
	int typeb = abs(junctionb.c);
	
	
	if (typea == Iblank || typeb == Iblank){
		return true;
	}
	
	if (typea < Iextend && typeb < Iextend) {
		
		if((typea == Iexonic && typeb == Iexonic) || (typea < Iexonic && typeb < Iexonic)){
			if (isgoodtequal(juna_good_seq, junb_good_seq)) {
				return true;
			}else {
				return false;
			}
		}else if(typea == Iexonic && typeb < Iexonic){
			if (juna_good_seq.a >= junb_good_seq.a && juna_good_seq.b <= junb_good_seq.b + EXON_BUFFER){
				return true;
			}else if (junb_good_seq.c!=0 && junb_good_seq.d!=0 && juna_good_seq.a >= junb_good_seq.c && juna_good_seq.b <= junb_good_seq.d + EXON_BUFFER){  
				return true;
			}else{
				return false;
			}
		}else if(typea < Iexonic && typeb == Iexonic){
			if (junb_good_seq.a >= juna_good_seq.a && junb_good_seq.b <= juna_good_seq.b + EXON_BUFFER){
				return true;
			}else if (juna_good_seq.c!=0 && juna_good_seq.d!=0 && junb_good_seq.a >= juna_good_seq.c && junb_good_seq.b <= juna_good_seq.d + EXON_BUFFER){
				return true;
			}else{
				return false;
			}
		}

	}else if(typea == Iextend && typeb < Iextend){
		if (juna_good_seq.a >= junb_good_seq.a && juna_good_seq.b <= junb_good_seq.b){
            return true;
		}else if (junb_good_seq.c!=0 && junb_good_seq.d!=0 && juna_good_seq.a >= junb_good_seq.c && juna_good_seq.b <= junb_good_seq.d){  
            return true;
		}else{
			return false;
		}
	}else if(typea < Iextend && typeb == Iextend){
		if (junb_good_seq.a >= juna_good_seq.a && junb_good_seq.b <= juna_good_seq.b){
            return true;
		}else if (juna_good_seq.c!=0 && juna_good_seq.d!=0 && junb_good_seq.a >= juna_good_seq.c && junb_good_seq.b <= juna_good_seq.d){
            return true;
		}else{
			return false;
		}
	}else if(typea == Iextend && typeb == Iextend){
		if ((junb_good_seq.b - juna_good_seq.a) == (local_b - local_a)){
            return true;
		}else{
			return false;
		}
	}
	
	return false; 
}


good_t good_seq(good_t &junction, int local_a, int local_b)
{
	int splice_local_pt = read_length - junction.d;
	good_t result = {0,0,0,0};
	
	
	if (abs(junction.c) < Iexonic) {
		int jun_start = junction.a - (read_length - junction.d) + 1;
		int jun_end = junction.b + junction.d;
		
		if (splice_local_pt >= local_b) {
			result.a = local_a + jun_start - 1;
			result.b = local_b + jun_start - 1;
			return result;
		}else if (splice_local_pt >= local_a) {
			result.a = local_a + jun_start - 1;
            result.b = junction.a;
            result.c = junction.b + 1;
            result.d = local_b + jun_end - read_length;
            return result;			
		}else if (splice_local_pt <= local_a - 1){  
			result.a = jun_end - (read_length - local_a);
            result.b = jun_end - (read_length - local_b);
            return result;			
		}
	}else if (abs(junction.c) == Iexonic) {  
		result.a = local_a + junction.a - 1;
        result.b = local_b + junction.a - 1;
        return result;
	}else if (abs(junction.c) == Iextend) {
		if (junction.d > 0){
			if (junction.d >= local_b) {
				result.a = local_a + junction.a - 1;
                result.b = local_b + junction.a - 1;
			}else if (junction.d >= local_a) {
				result.a = local_a + junction.a - 1;
                result.b = junction.b;
			}else if (junction.d < local_a) {
				result.a = 0;  
                result.b = 0;
			}
		}else if (junction.d < 0) { 
			int left_arm = read_length + junction.d;
            int local_left_start = left_arm + 1;
			
            if (local_left_start <= local_a){
                result.a = local_a + junction.a - local_left_start;
                result.b = local_b + junction.a - local_left_start;
			}else if (local_left_start <= local_b){
                result.a = junction.a;
                result.b = local_b + junction.a - local_left_start;
			}else if (local_left_start > local_b){
                result.a = 0;
				result.b = 0;   
			}
		}
		return result;
	}
	
	return result; 
}





string count_mismatch(start_end_nano_t &alignment, string &read, string &genome,int_fast32_t num_read_mismatch,int_fast32_t max_clip_allowed){
	int num_mismatch = 0;
	int num_clip = 0;
	
	string output;
	
	num_clip = (alignment.read_start.front()-1) + ((int)read.length() - alignment.read_end.back());
	
	if (num_clip > max_clip_allowed) {
		return ""; 
	}
	
	vector<short>::iterator read_start_it = alignment.read_start.begin();
	vector<short>::iterator read_end_it = alignment.read_end.begin();
	vector<int>::iterator chr_start_it = alignment.chr_start.begin();
	
	int match_count = 0;
	string MD_str = "MD:Z:";
	
	while (read_start_it != alignment.read_start.end()) {
		int r = *read_start_it;
		int c = *chr_start_it;
		
		while (r <= *read_end_it) {

			if (read[r-1] != genome[c-1]) {
				if(match_count!=0){
					MD_str += IntToStr(match_count) + genome[c-1];
				}else {
					MD_str += genome[c-1];
				}

					
				num_mismatch++;
				match_count = 0;
			}else {
				match_count++;
			}

			
			r++;
			c++;
		}
		
		read_start_it++;
		read_end_it++;
		chr_start_it++;
	}
	
	if (num_mismatch > num_read_mismatch) {
		return ""; 
	}
	
	if (match_count != ((int)read.length() - num_clip)) {
		if(match_count!=0){
			output += "\t" + MD_str + IntToStr(match_count);
		}else {
			output += "\t" + MD_str;
		}

	}
	
	output += "\tNM:i:" + IntToStr(num_mismatch); 
	
	output += "\tXC:i:" + IntToStr(num_clip); 
	
	return output;
}


inline void print_usage_and_exit()
{
	
	
	
	cout << "WARNING: This program is for internal use only, I suggest you use \"runSpliceMap\" instead. " << endl;
	
	
	
	
	exit(2);
}







