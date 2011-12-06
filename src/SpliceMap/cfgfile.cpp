

#include "cfgfile.h"

cfgfile::cfgfile(string file_name)
{
	ifstream infile;
	string line;
	string temp;
	
	infile.open(file_name.c_str(), ios::in);
	
	if (!infile.is_open()) {
		cout << "ERROR: I'm sorry, could not open the configuration file, " << file_name << endl;
		exit(1);
	}
	
	while (!infile.eof()) {
		getline(infile, line);
		trim2(line);
		if (line.length() == 0) {
			continue;
		}
		
		if (line[0] == '#') {
			
		}else if(line[0] == '>'){
			string key = line.substr(1);
			trim2(key);
			
			vector<string> val;
			
			
			while (true) {
				getline(infile, line);
				trim2(line);
				if (line.length() == 0) {
					continue;
				}
				
				if (line.compare("<") == 0){
					break;
				}
				
				if (infile.eof()) {
					cout << "ERROR: I'm sorry, please check the format of the configuration file line, reached end of file " << endl;
					exit(1);
				}
				
				val.push_back(line);
			}
			
			pair<map<string,vector<string> >::iterator,bool> out = data.insert(pair<string,vector<string> >(key,val));
			
			if (!out.second) {
				cout << "ERROR: I'm sorry, there is a duplicate key,  " << key << endl;
				exit(1);
			}
			
		}else {
			size_t equal_pos = line.find('=');
			
			if (equal_pos == string::npos) {
				cout << "ERROR: I'm sorry, please check the format of the configuration file line, " << line << endl;
				exit(1);
			}else {
				string key = line.substr(0, equal_pos);
				trim2(key);
				vector<string> val;
				temp = line.substr(equal_pos+1);
				trim2(temp);
				val.push_back(temp);
				
				pair<map<string,vector<string> >::iterator,bool> out = data.insert(pair<string,vector<string> >(key,val));
				
				if (!out.second) {
					cout << "ERROR: I'm sorry, there is a duplicate key,  " << key << endl;
					exit(1);
				}
			}

		}

		
	}
	
	infile.close();
	infile.clear();

}



string cfgfile::getVal(string key)
{
	string result;
	
	map<string,vector<string> >::iterator data_it = data.find(key);
	
	if (data_it != data.end()) {
		result = data_it->second.front();
	}else {
		result = "";
	}

	
	return result;
}



vector<string> cfgfile::getList(string key)
{
	vector<string> result;
	
	map<string,vector<string> >::iterator data_it = data.find(key);
	
	if (data_it != data.end()) {
		result = data_it->second;
	}
	
	
	return result;
}






