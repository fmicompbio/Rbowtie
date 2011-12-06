

#include "amalgamateSAM.h"



using namespace std;


int main (int argc, char * const argv[]) 
{
	ifstream curr_SAM_file;
	string temp_path;
	map<string,string> chr_list; 
	
	map<string,jundict_ij_t> jundict_chr;   
	
	string SAM_header;
	
	string outfile_name;
	
	string temp;
	string line;

	struct timeval tv;
	struct timeval start_tv;
	
	gettimeofday(&start_tv, NULL);
	
	
	if (argc == 3) {
		temp_path = argv[1];
		outfile_name = argv[2];
	}else {
		print_usage_and_exit();
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
		
		temp = temp_path+line_list[1] + "_" + line_list[3] + ".sam";
		chr_list.insert(pair<string,string>(line_list[0],temp));

	}
	
	ref_map_file.close();
	ref_map_file.clear();
	
	cout << "Reading SAM headers" << endl;
	
	SAM_header = "@HD\tVN:1.0\tSO:coordinate\n";
	
	
	map<string,string>::iterator chr_list_it = chr_list.begin();
	while (chr_list_it != chr_list.end()) { 
		curr_SAM_file.open(chr_list_it->second.c_str(), ios::in);
		if (!curr_SAM_file.is_open()) {
			cout << "ERROR: could not open " << chr_list_it->second << endl;
			exit(2);
		}
		
		bool found = false;
		while (!found && !curr_SAM_file.eof()) {
			getline(curr_SAM_file, line);
			trim2(line);
			if (line.length() == 0) {
				continue;
			}
			
			vector<string> line_list = split(line, '\t');
			
			if (line_list[0] == "@SQ") {
				found = true;
				
				SAM_header += line + '\n';
			}
		}
		
		curr_SAM_file.close();
		curr_SAM_file.clear();
		
		chr_list_it++;
	}
	
	SAM_header += "@PG\tID:SpliceMap\tVN:"+VERSION+"\n";
	
	
	
	
	
	cout << "Processing SAM files..." << endl;

	
	
	
	vector<ifstream*> insam_list;
	vector<ifstream*>::iterator insam_list_it;
	string *curr_line = new string[chr_list.size()];
	int file_idx = 0;
	chr_list_it = chr_list.begin();
	
	while (chr_list_it != chr_list.end()) { 
		ifstream *temp_infile = new ifstream(chr_list_it->second.c_str(), ios::in);
		if (!temp_infile->is_open()) {
			cout << "ERROR: could not open " << chr_list_it->second << "   :("<< endl;
			exit(2);
		}
		
		insam_list.push_back(temp_infile);

		getline(*temp_infile,curr_line[file_idx]);
		trim2(curr_line[file_idx]);
		
		file_idx++;
		chr_list_it++;
	}
	
	
	vector<ofstream*> outsam_list;
	vector<ofstream*>::iterator outsam_list_it;
	
	chr_list_it = chr_list.begin();
	while (chr_list_it != chr_list.end()) { 
		temp = chr_list_it->second + "_uniq";
		ofstream *temp_outfile = new ofstream(temp.c_str(), ios::out);

		if (!temp_outfile->is_open()) {
			cout << "ERROR: could not write " << temp << " :("<< endl;
			exit(2);
		}
		
		outsam_list.push_back(temp_outfile);
		
		chr_list_it++;
	}
	
	
	
	priority_queue<posline_t,vector<posline_t>, comparison_reverse_t> data_queue;
	
	
	bool done = false;
	
	
	while(!done){
		done = true;
		
		
		int min_file_index = INT_MAX;
		file_idx = 0;
		insam_list_it = insam_list.begin();
		while (insam_list_it != insam_list.end()) {

			int first_idx = 0;
			bool keep_reading = true;
			while (keep_reading && !(*insam_list_it)->eof()) {
				
				if (curr_line[file_idx].length() == 0) {
					getline(**insam_list_it,curr_line[file_idx]);
					trim2(curr_line[file_idx]);
					continue;
				}
				
				
				
				done = false; 
				
				size_t tabloc = curr_line[file_idx].find('\t');
				string name_str = curr_line[file_idx].substr(0, tabloc);
				size_t start_pos = name_str.rfind('[')+1;
				size_t end_pos = name_str.rfind(']');
				int curr_idx = atoi(name_str.substr(start_pos, end_pos-start_pos).c_str());
				
				if (first_idx == 0) {
					first_idx = curr_idx;
					posline_t contents;
					contents.pos = curr_idx;
					contents.file_index = file_idx;
					contents.data = curr_line[file_idx];
					
					data_queue.push(contents);
					
					
					if (curr_idx < min_file_index) {
						min_file_index = curr_idx;
					}
					
					getline(**insam_list_it,curr_line[file_idx]);
					trim2(curr_line[file_idx]);
				}else {
					if (curr_idx == first_idx) {
						posline_t contents;
						contents.pos = curr_idx;
						contents.file_index = file_idx;
						contents.data = curr_line[file_idx];
						
						data_queue.push(contents);
						
						
						if (curr_idx < min_file_index) {
							min_file_index = curr_idx;
						}
						
						getline(**insam_list_it,curr_line[file_idx]);
						trim2(curr_line[file_idx]);
					}else {
						keep_reading = false;
					}

				}
				
				
			}
			
			file_idx++;
			insam_list_it++;
		}
		
		
		
		while (data_queue.size() > 0 && data_queue.top().pos <= min_file_index) {
			
			done = false;
			
			int first_idx = data_queue.top().pos;
			vector<posline_t> curr_read_group;
			
			
			while (data_queue.size() > 0 && data_queue.top().pos == first_idx) {
				curr_read_group.push_back(data_queue.top());
				data_queue.pop();
				
			}
			
			
			int map_count[2] = {0,0}; 
			if (curr_read_group.size() > 1) {
				
				vector<posline_t>::iterator curr_read_group_it = curr_read_group.begin();
				while(curr_read_group_it != curr_read_group.end()){
					
					size_t tabloc  = curr_read_group_it->data.find('\t'); 
					size_t tabloc2 = curr_read_group_it->data.find('\t',tabloc+1); 
					
					int flag = atoi(curr_read_group_it->data.substr(tabloc+1,tabloc2 - tabloc - 1).c_str());
					
					if (flag & 1){ 
						if (flag & 64) { 
							
							map_count[0]++;
						}else if  (flag & 128){  
							
							map_count[1]++;
						}else {
							cout << "ERROR: Paring unknown... " << endl;
							cout <<  curr_read_group_it->data << endl;
							exit(2);
						}
						
						
					}else { 
						map_count[0]++;
					}
					
					curr_read_group_it++;
					
				}
			}
			
			
			bool non_unique = false;
			if (map_count[0]>1 || map_count[1]>1) {
				non_unique = true;
			}
			
			
			vector<posline_t>::iterator curr_read_group_it = curr_read_group.begin();
			while (curr_read_group_it != curr_read_group.end()) {
				
				update_jun_coverage(jundict_chr, curr_read_group_it->data, non_unique);
				curr_read_group_it++;
			}
			
			
			
			
			if (non_unique) {
				
				curr_read_group_it = curr_read_group.begin();
				while (curr_read_group_it != curr_read_group.end()) {
					string line_str = curr_read_group_it->data;
					
					size_t tabloc = line_str.find('\t');
					size_t tabloc2 = line_str.find('\t',tabloc+1); 
					
					int flag = atoi(curr_read_group_it->data.substr(tabloc+1,tabloc2 - tabloc - 1).c_str());
					
					tabloc = tabloc2;
					tabloc2 = line_str.find('\t',tabloc+1); 
					tabloc = tabloc2;
					tabloc2 = line_str.find('\t',tabloc+1); 
					tabloc = tabloc2;
					tabloc2 = line_str.find('\t',tabloc+1); 
					
					line_str = line_str.substr(0, tabloc) + "\t0\t" + line_str.substr(tabloc2+1); 
					
					if (flag & 1){ 
						if (flag & 64) { 
							line_str += "\tNH:i:" + IntToStr(map_count[0]);
						}else if (flag & 128){  
							line_str += "\tNH:i:" + IntToStr(map_count[1]);
						}else {
							cout << "ERROR: Paring unknown... " << endl;
							cout << line_str << endl;
							exit(2);
						}
						
						
					}else { 
						line_str += "\tNH:i:" + IntToStr(map_count[0]);
					}
					
					(*(outsam_list[curr_read_group_it->file_index])) << line_str << '\n';
					curr_read_group_it++;
				}
			}else {
				
				curr_read_group_it = curr_read_group.begin();
				while (curr_read_group_it != curr_read_group.end()) {
					(*(outsam_list[curr_read_group_it->file_index])) << curr_read_group_it->data << '\n';
					curr_read_group_it++;
				}
			}

			
		}
		
		
		
		
	}
	
	delete [] curr_line;
	
	
	cout << "Amalgamating SAM files" << endl;
	ofstream output;
	temp = outfile_name + ".sam";
	output.open(temp.c_str(), ios::out);
	
	if (!output.is_open()) {
		cout << "ERROR: Could not create output file -- " << temp << endl;
		exit(1);
	}
	
	output << SAM_header; 
	
	chr_list_it = chr_list.begin();
	while (chr_list_it != chr_list.end()) { 
		temp = chr_list_it->second + "_uniq";
		curr_SAM_file.open(temp.c_str(), ios::in);
		if (!curr_SAM_file.is_open()) {
			cout << "ERROR: could not open " << chr_list_it->second << "   :("<< endl;
			exit(2);
		}
		

		while (!curr_SAM_file.eof()) {
			getline(curr_SAM_file, line);
			trim2(line);
			if (line.length() == 0 || line[0] == '@') {
				continue;
			}
			
			output << line << '\n';
		}
		
		curr_SAM_file.close();
		curr_SAM_file.clear();
		
		chr_list_it++;
	}
	
	
	if(jundict_chr.size() > 0){
		cout << "Printing junctions... for " << jundict_chr.size() << " chromosomes" << endl;
		
		map<string,jundict_ij_t>::iterator jundict_chr_it = jundict_chr.begin();
		
		print_jun_dict_ij(jundict_chr_it->second,outfile_name+".bed",jundict_chr_it->first,"SpliceMap junctions");
		
		jundict_chr_it++;
		while(jundict_chr_it != jundict_chr.end()){
			
			
			print_jun_dict_ij(jundict_chr_it->second,outfile_name+".bed",jundict_chr_it->first,"");
			
			
			jundict_chr_it++;
		}
	}else {
		cout << "Warning: No junctions found" << endl;
	}
	
	
	

	gettimeofday(&tv, NULL);
	cout << "______________________" << endl;
	cout<<"SAM File Amalgamation Time: "<<diffclock(start_tv,tv)<<" s."<<endl;
	
	insam_list_it = insam_list.begin();
	while (insam_list_it != insam_list.end()) {
		(*insam_list_it)->close();
		(*insam_list_it)->clear();
		delete *insam_list_it;
		
		insam_list_it++;
	}
	
	outsam_list_it = outsam_list.begin();
	while (outsam_list_it != outsam_list.end()) {
		(*outsam_list_it)->close();
		(*outsam_list_it)->clear();
		delete *outsam_list_it;
		
		outsam_list_it++;
	}
	
	return 0;
}



void update_jun_coverage(map<string,jundict_ij_t> &jundict_chr, string line, bool multi_mapped)
{
	size_t strand_pos = line.rfind("XS:A:");
	
	if(strand_pos != string::npos){  
		
		bool strand = true;
		if (line[strand_pos+5] == '-') {
			strand = false;
		}
		
		int num_clip = 0;
		size_t clip_pos = line.rfind("XC:i:");
		if(clip_pos != string::npos){
			clip_pos = clip_pos+5;
			size_t clip_pos_end = line.find('\t',clip_pos);
			num_clip = atoi(line.substr(clip_pos, clip_pos_end-clip_pos).c_str());
		}
		
		vector<string> line_list = split(line);
		string chr_name = line_list[2];
		
		vector<uint_fast32_t> exon_start_list; 
		vector<uint_fast32_t> exon_end_list;
		uint_fast32_t front_clip = 0;
		
		
		
		
		string cigar = line_list[5];
		uint_fast32_t start_loc = atoi(line_list[3].c_str());
		size_t read_ptr = 0;
		
		while (read_ptr != string::npos) {
			size_t temp_read_ptr = cigar.find_first_of("SMN", read_ptr);
			
			if(temp_read_ptr != string::npos){
				if(cigar[temp_read_ptr] == 'M'){
					uint_fast32_t range = (uint_fast32_t)atoi(cigar.substr(read_ptr, temp_read_ptr-read_ptr).c_str());
					
					exon_start_list.push_back(start_loc);
					exon_end_list.push_back(start_loc+range-1);
					
					start_loc = start_loc + range;
				}else if(cigar[temp_read_ptr] == 'N'){ 
					uint_fast32_t range = (uint_fast32_t)atoi(cigar.substr(read_ptr, temp_read_ptr-read_ptr).c_str());
					
					start_loc = start_loc + range;
				} else if (read_ptr == 0) { 
					front_clip = (uint_fast32_t)atoi(cigar.substr(read_ptr, temp_read_ptr-read_ptr).c_str());
					
					
				}
				
				
				read_ptr = temp_read_ptr + 1;
				
				if (read_ptr == cigar.length()) {
					break;
				}
			}else {
				read_ptr = temp_read_ptr;
			}
			
			
		}
		
		
		
		
		
		if (exon_start_list.size() > 1) {
			
			pair<map<string,jundict_ij_t>::iterator, bool> jundict_pair 
			= jundict_chr.insert(pair<string,jundict_ij_t>(chr_name,jundict_ij_t()));
			
			vector<uint_fast32_t>::iterator start_list = exon_start_list.begin();
			vector<uint_fast32_t>::iterator end_list = exon_end_list.begin();
			
			start_list++; 
			
			while (start_list !=  exon_start_list.end()) {
				
				uint_fast32_t chr_pos = *end_list;  
				uint_fast32_t right_start = *start_list - 1;  
				
				bool direc_type;
				if (strand){
					direc_type = true; 
				} else {
					direc_type = false;
				}
				
				
				
				end_list++;
				
				
				
				
				
				jundict_ij_t::iterator jun_it = ((jundict_pair.first)->second).find(chr_pos);
				
				if (jun_it != ((jundict_pair.first)->second).end()) {
					jundict_b_ij_t::iterator jun_b_it = (jun_it->second).find(right_start);
					if (jun_b_it != (jun_it->second).end()) {
						
						
						
						
						
						if(multi_mapped){
							(jun_b_it->second).j++;  
						}else {
							(jun_b_it->second).i++; 
						}

						addnNR((jun_b_it->second).nNR,exon_start_list,exon_end_list,num_clip);
						
					}else {
						
						nNR_t nNR;
						
						addnNR(nNR,exon_start_list,exon_end_list,num_clip);
						
						
						jun_store_ij store = {direc_type,0,0,nNR};
						if(multi_mapped){
							store.j++;
						}else {
							store.i++;
						}
						
						
						
						(jun_it->second).insert(pair<int,jun_store_ij>(right_start,store));
						
					}
					
				}else {
					
					nNR_t nNR;
					
					addnNR(nNR,exon_start_list,exon_end_list,num_clip);

					jun_store_ij store = {direc_type,0,0,nNR};
					if(multi_mapped){
						store.j++;
					}else {
						store.i++;
					}

					
					
					jundict_b_ij_t jundict_b;
					jundict_b.insert(pair<int,jun_store_ij>(right_start,store));
					((jundict_pair.first)->second).insert(pair<int,jundict_b_ij_t>(chr_pos,jundict_b));
					
				}
				
				start_list++;
				
				
			}
		}
		
	}
	

}

void addnNR(nNR_t &store, vector<uint_fast32_t> exon_start_list, vector<uint_fast32_t> exon_end_list, int num_clip)
{
	
	
	nNR_t::iterator it = store.begin();
	while (it != store.end()) {
		bool same = true;
		if (!compareIntVector(it->chr_start,exon_start_list)) {
			same = false;
		}else if (!compareIntVector(it->chr_end,exon_end_list)) {
			same = false;
		}else if (it->num_clip != num_clip) {
			same = false;
		}
		
		if (same) {
			return;  
		}
		it++;
	}
	
	
	nNR_inside_t contents;
	contents.chr_start = exon_start_list;
	contents.chr_end = exon_end_list;
	contents.num_clip = num_clip;  
	
	store.push_back(contents);
}


bool compareIntVector(vector<uint_fast32_t> &a, vector<uint_fast32_t> &b)
{
	if (a.size() != b.size()) {
		return false;
	}
	
	vector<uint_fast32_t>::iterator it_a = a.begin();
	vector<uint_fast32_t>::iterator it_b = b.begin();
	
	while (it_a != a.end()) {
		if (*it_a != *it_b) {
			return false;
		}
		
		it_a++;
		it_b++;
	}
	
	return true;
}

void print_usage_and_exit()
{
	cout << "For internal SpliceMap usage only... it joins the individual chromosomeSAM files together" << endl;
	exit(1);
}






