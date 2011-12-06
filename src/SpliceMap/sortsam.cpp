


#include "sortsam.h"
using namespace std;








int main (int argc, char * const argv[]) 
{
	string line;
	string temp;
	int tabloc;
	int tabloc2;
	string chr_name;   
	uint_fast32_t curr_pos;
	int curr_idx;
	posline_t temp_data;
	
	int mode = 1; 
	
	string infile_name;
	string outfile_name;
	
	ifstream infile;
	ofstream outfile;
	
	
	vector<posline_t> full_data;
	
	if (argc > 4 || argc < 3) {
		print_usage_and_exit();
	}
	
	if (argc == 4) {
		temp = argv[1];
		
		
		
		if (temp.compare("-pos") == 0) {
			mode = 1;
		}else if (temp.compare("-idx") == 0){
			mode = 2;
		}else {
			mode = 0;
		}

		infile_name = argv[2];
		outfile_name = argv[3];
	}else {
		infile_name = argv[1];
		outfile_name = argv[2];
	}

	if (mode < 1 || mode > 2) {
		void print_usage_and_exit();
	}
	
	
	
	infile.open(infile_name.c_str(), ios::in);
	
	if (!infile.is_open()) {
		cout << "ERROR: I'm sorry, could not open the infile, " << infile_name << endl;
		exit(1);
	}
	
	outfile.open(outfile_name.c_str(),ios::out);
	
	

	while (!infile.eof()) {
		getline(infile, line);
		trim2(line);
		if (line.length() == 0) {
			continue;
		}
		
		if (line[0] == '@') {
			outfile << line + '\n';
		}else {
			break;
		}
	}

	outfile.flush();
	
	
	
	if (mode == 1) {  
		
		tabloc = (int) line.find('\t'); 
		tabloc = (int) line.find('\t',tabloc+1); 
		tabloc2 = (int) line.find('\t',tabloc+1); 
		chr_name = line.substr(tabloc+1, tabloc2 - tabloc - 1);
		
		cout << "Processing... " << chr_name << endl;
		
		tabloc = tabloc2;
		tabloc2 = (int) line.find('\t',tabloc+1); 
		curr_pos = (uint_fast32_t)atoi(line.substr(tabloc+1, tabloc2 - tabloc - 1).c_str());
		temp_data.pos = curr_pos;
		temp_data.data = line;
		full_data.push_back(temp_data);
		
		while (!infile.eof()) {
			getline(infile, line);
			if (line.length() < 3) {
				continue;
			}
			
			
			tabloc = (int) line.find('\t');
			tabloc = (int) line.find('\t',tabloc+1); 
			tabloc2 = (int) line.find('\t',tabloc+1); 
			string temp_chr_name = line.substr(tabloc+1, tabloc2 - tabloc - 1);
			tabloc = tabloc2;
			tabloc2 = (int) line.find('\t',tabloc+1); 
			curr_pos = (uint_fast32_t)atoi(line.substr(tabloc+1, tabloc2 - tabloc - 1).c_str());
			
			
			if (temp_chr_name.compare(chr_name) != 0) {
				
				cout << "Sorting... " << endl;
				sort(full_data.begin(), full_data.end(), compare_posline);
				
				
				
				vector<posline_t>::iterator line_it = full_data.begin();
				while (line_it != full_data.end()) {
					outfile << line_it->data << '\n';
					line_it++;
				}
				
				full_data.clear();
				chr_name = temp_chr_name;
				
				cout << "Processing... " << chr_name << endl;
			}
			
			
			temp_data.pos = curr_pos;
			temp_data.data = line;
			full_data.push_back(temp_data);
			
			
		}
		
		
		cout << "Sorting... " << endl;
		sort(full_data.begin(), full_data.end(), compare_posline);
		
		vector<posline_t>::iterator line_it = full_data.begin();
		while (line_it != full_data.end()) {
			outfile << line_it->data << '\n';
			line_it++;
		}
		
	}else if (mode == 2){  
		priority_queue<posline_t,vector<posline_t>, comparison_reverse_t> data_queue;
		int up_to_num = 0;
		
		
		tabloc = (int) line.find('\t'); 
		curr_idx = atoi(line.substr(0,tabloc).c_str());
		
		cout << "Processing... " << endl;
		
		
		
		tabloc2 = (int) line.find('\t',tabloc+1); 
		tabloc = tabloc2;
		tabloc2 = (int) line.find('\t',tabloc+1); 
		tabloc = tabloc2;
		tabloc2 = (int) line.find('\t',tabloc+1); 
		
		
		line = line.substr(0,tabloc2) + "\t" +line.substr(line.length()-1);
		
		
		temp_data.pos = curr_idx;
		temp_data.data = line;
		data_queue.push(temp_data);
		
		
		while (!infile.eof()) {
			int temp_curr_idx;
			
			getline(infile, line);
			if (line.length() < 3) {
				continue;
			}
			
			
			tabloc = (int) line.find('\t'); 
			temp_curr_idx = atoi(line.substr(0,tabloc).c_str());
			
			
			
			tabloc2 = (int) line.find('\t',tabloc+1); 
			tabloc = tabloc2;
			tabloc2 = (int) line.find('\t',tabloc+1); 
			tabloc = tabloc2;
			tabloc2 = (int) line.find('\t',tabloc+1); 
			
			
			line = line.substr(0,tabloc2) + "\t" +line.substr(line.length()-1);
			
			if (temp_curr_idx == curr_idx) {
				temp_data.pos = curr_idx;
				temp_data.data = line;
				data_queue.push(temp_data);
			}else {
				
				int top_index = data_queue.top().pos; 
				
				
				if (top_index == up_to_num) {
					
					
					while (top_index == up_to_num && data_queue.size()>0) {
						posline_t contents = data_queue.top();

						outfile << contents.data << "\n";
						
						data_queue.pop();
						
						top_index = data_queue.top().pos;
						if(top_index != up_to_num){
							up_to_num++;
						}
					}
					if (data_queue.size() == 0) {
						up_to_num++;
					}
					
					curr_idx = temp_curr_idx;
					
					temp_data.pos = curr_idx;
					temp_data.data = line;
					data_queue.push(temp_data);	
					
				}else {
					curr_idx = temp_curr_idx;
					
					temp_data.pos = curr_idx;
					temp_data.data = line;
					data_queue.push(temp_data);					
				}

			}

		}
		
		
		while (data_queue.size()>0) {
			posline_t contents = data_queue.top();
			
			outfile << contents.data << "\n";
			
			data_queue.pop();
		}
		

	}

	infile.close();
	infile.clear();
	outfile.close();
	outfile.clear();
	
	return 0;
}



void print_usage_and_exit()
{
	cout << "usage: ./sortsam [-pos|-idx] infile.sam outfile_sorted.sam" << endl;
	cout << "\t-pos -- Sort by chromosome position" << endl;
	cout << "\t-idx -- Sort by read_index position" << endl;
	cout << "Memory efficiently sorts a SAM file by coordinate or read index" << endl;
	cout << "WARNING: This is for internal SpliceMap use only, behaviour not guranteed" << endl;
	exit(2);
}
