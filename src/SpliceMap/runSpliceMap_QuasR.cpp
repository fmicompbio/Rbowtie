


#include "runSpliceMap_QuasR.h"


int main (int argc, char * const argv[]) {
    
    string bowtie_base_dir; 
    string ref_filename = "";
	
    string temp_path;
    vector<string> reads_filename[2];
    vector<string>::iterator readslist_it;
    string cfg_filename;
    string read_format;
	
    string temp;
    string line;
    string cmd;
	
    ofstream debug_out;

    pair<string,string> package_path_filename;
    pair<string,string> ref_path_filename;   
	
    map<string,reference_t> ref_map; 
    map<string,reference_t>::iterator ref_map_it;
	
    int num_pair = 0;
    int_fast32_t num_seed_mismatch = 1; 
    int_fast32_t num_read_mismatch = 2;  
    int_fast32_t max_clip_allowed = 40;  
	
    int fullread_length = -1;
    unsigned int head_clip_length = 0;
    bool print_sam = false; 
    bool cufflinks = false;
	
    string merge_opt = "";
    int max_intron = 400000;
    int min_intron = 20000;
    int min_extend_req = 0;
    int num_chromosome_together = 1;
    int quality_format = 0; 

    enum QusaR_Mode {GENERATE25MERS, INDEX25MERALIGNMENTS};
    int quasr_mode = 0;
	
    struct timeval tv;
    struct timeval start_tv;
    struct timeval temp_start_tv;
    struct timeval map_temp_start_tv;
	
    gettimeofday(&start_tv, NULL);
	
    if (argc == 3) {
	temp = argv[0];
	package_path_filename = get_path_and_filename(temp);
		
	ifstream temp_infile;
		
	cfg_filename = argv[1];
	cout << "Loading configuration file... " << cfg_filename << endl;
		
	cfgfile *run_cfg = new cfgfile(cfg_filename);
		
	/* QuasR mode */
	temp = argv[2];
	if (temp.compare("generate25mers") == 0) {
	    quasr_mode = GENERATE25MERS;
	}else if (temp.compare("index25merAlignments") == 0) {
	    quasr_mode = INDEX25MERALIGNMENTS;
	}else {
	    cout << "ERROR: The quasr_mode " << temp << " is not supported. " << endl;
	    debug_out << "ERROR: The quasr_mode " << temp << " is not supported. " << endl;
			
	    cout << "Supported modes are generate25mers and index25merAlignments" << endl;
	    exit(2);
	}

	temp_path = run_cfg->getVal("temp_path");
	if (temp_path.length() == 0) {
	    cout << "Warning: Temp directory location (temp_path) not specified, using default: temp/" << endl;
	    temp_path = "temp/";
	}else {
	    if (temp_path[temp_path.length()-1] != '/') {
		temp_path = temp_path + '/';
	    }
	}
		
	temp = temp_path+debug_path+"run_debug.log";
	debug_out.open(temp.c_str(), ios::out);
		
	if (quasr_mode == GENERATE25MERS) {
	    string genome_direc = run_cfg->getVal("genome_dir");
	    if (genome_direc.length() == 0) {
		cout << "ERROR: The location of the genome files should be specified in the configuration file tag \"genome_dir\". " << endl;
		debug_out << "ERROR: The location of the genome files should be specified in the configuration file tag \"genome_dir\". " << endl;
		exit(2);
	    }
		
	    if (genome_direc[genome_direc.length()-1] != '/' && genome_direc[genome_direc.length()-1] != '\\') {
		genome_direc = genome_direc+'/';  
	    }
		
	    vector<string> chr_file_list = run_cfg->getList("genome_files");
	    vector<string>::iterator chr_file_list_it = chr_file_list.begin();
		
	    cout << "Scaning genome: " << genome_direc << endl;
	    debug_out << "Scaning genome: " << genome_direc << endl;
		
	    if (!read_reference_map(genome_direc, ref_map, chr_file_list)) {
		cout << "ERROR: something wrong reading the genome files..." << endl;
		debug_out << "ERROR: something wrong reading the genome files..." << endl;
		exit(1);
	    }
		
	    if(ref_map.size() == 0){
		cout << "ERROR: No chromosomes found within the reference files" << endl;
		cout << "Please check that they are FASTA format" << endl;
		exit(1);
	    }
		
	    ofstream ref_list_file;
	    temp = temp_path+reference_filename;
	    ref_list_file.open(temp.c_str(), ios::out);
		
	    ref_map_it = ref_map.begin();
	    while (ref_map_it != ref_map.end()) {
			
		ref_list_file << ref_map_it->first << "\t" 
			      << ref_map_it->second.file_name << "\t" 
			      << ref_map_it->second.file_path << "\t" 
			      << ref_map_it->second.file_index_start << "\t"
			      << ref_map_it->second.file_index_end << "\n"; 
			
		ref_map_it++;
	    }
		
	    ref_list_file.close();
	    ref_list_file.clear();
	}
		
	cout << "List of chromosomes to be searched: " << endl;
	debug_out << "List of chromosomes to be searched: " << endl;
		
	ref_map_it = ref_map.begin();
	while (ref_map_it != ref_map.end()) {
	    cout << ref_map_it->first << " | " << ref_map_it->second.file_path << ref_map_it->second.file_name << " | pos:" << ref_map_it->second.file_index_start << " - " << ref_map_it->second.file_index_end << "\n";
	    debug_out << ref_map_it->first << " | " << ref_map_it->second.file_path << ref_map_it->second.file_name << " | pos:" << ref_map_it->second.file_index_start << " - " << ref_map_it->second.file_index_end << "\n";
			
	    ref_map_it++;
	}
		
	reads_filename[0] = run_cfg->getList("reads_list1");
	reads_filename[1] = run_cfg->getList("reads_list2");
		
	if(reads_filename[0].size() == 0){
	    cout << "ERROR: There should be at least one reads file in the first list (reads_list1). " << endl;
	    debug_out << "ERROR: There should be at least one reads file in the first list (reads_list1). " << endl;
	    exit(2);
	}
		
	if(reads_filename[1].size() == 0){
	    num_pair = 1;
	}else {
	    num_pair = 2;
			
	    if (reads_filename[0].size() != reads_filename[1].size()) {
		cout << "ERROR: Both read lists must be the same length. Please check reads_list1 and reads_list2" << endl;
		debug_out << "ERROR: Both read lists must be the same length. Please check reads_list1 and reads_list2" << endl;
		exit(2);
	    }
	}
		
	temp = run_cfg->getVal("full_read_length");
	if (temp.length() > 0) {
	    fullread_length = atoi(temp.c_str());
			
	    if (fullread_length == 0) {
		cout << "ERROR: Zero full_read_length" << endl;
		debug_out << "ERROR: Zero full_read_length" << endl;
		exit(2);
	    }
	}
		
	temp = run_cfg->getVal("head_clip_length");
	if (temp.length() > 0) {
	    head_clip_length = atoi(temp.c_str());
	    if (head_clip_length < 0) {
		cout << "Warning: Negative head_clip_length (reseting to 0): " << head_clip_length << endl;
		debug_out << "Warning: Negative head_clip_length (reseting to 0): " << head_clip_length << endl;
		
		head_clip_length = 0;
	    }
	    if (fullread_length>0 && fullread_length-head_clip_length < 50) {
		cout << "ERROR: fullread_length-head_clip_length < 50, the minimum read length accepted by SpliceMap is 50" << endl;
		debug_out << "ERROR: fullread_length-head_clip_length < 50, the minimum read length accepted by SpliceMap is 50" << endl;
				
		exit(2);
	    }
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
		
	ref_filename = run_cfg->getVal("annotations");
	ifstream temp_file;
	temp_file.open(ref_filename.c_str(), ios::in);
	if (ref_filename.length() > 0 && !temp_file.is_open()) {
	    cout << "ERROR: I'm sorry, cannot open the annotation file, " << ref_filename << endl;
	    debug_out << "ERROR: I'm sorry, cannot open the annotation file, " << ref_filename << endl;
	    
	    exit(2);
	}
	temp_file.close();
	temp_file.clear();
		
	ref_path_filename = get_path_and_filename(ref_filename);
		
	temp = run_cfg->getVal("sam_file");
	if (temp.compare("sam") == 0) {
	    print_sam = true;
	    merge_opt+=" -sam ";
	}else if(temp.compare("cuff") == 0){
	    cufflinks = true;
	    merge_opt+=" -cuff ";
	}else if(temp.length() == 0){
	    print_sam = false;
	    cufflinks = false;
	}else {
	    cout << "ERROR: The choices for SAM output are \"sam\", \"cuff\" ... not  " << temp << endl;
	    debug_out << "ERROR: The choices for SAM output are \"sam\", \"cuff\" ... not  " << temp << endl;
			
	    exit(2);
	}
		
	read_format = run_cfg->getVal("read_format");
	if (read_format.compare("FASTQ") == 0) {
			
	}else if (read_format.compare("FASTA") == 0) {
			
	}else if (read_format.compare("RAW") == 0) {
			
	}else {
	    cout << "ERROR: The read format " << read_format << " is not supported. " << endl;
	    cout << "Supported formats are FASTQ, FASTA and RAW" << endl;
			
	    debug_out << "ERROR: The read format " << read_format << " is not supported. " << endl;
	    debug_out << "Supported formats are FASTQ, FASTA and RAW" << endl;
			
	    exit(2);
	}
		
	temp = run_cfg->getVal("quality_format");
	if (temp.compare("phred-33") == 0) {
	    quality_format = 0;
	}else if (temp.compare("phred-64") == 0) {
	    quality_format = 1;
	}else if (temp.compare("solexa") == 0) {
	    quality_format = 2;
	}else if (temp.length() == 0){
	    
	}else {
	    cout << "ERROR: The quality format " << quality_format << " is not supported. " << endl;
	    debug_out << "ERROR: The quality format " << quality_format << " is not supported. " << endl;
			
	    cout << "Supported formats are phred-33(default), phred-64 and solexa" << endl;
	    cout << "Note that \"Illumnina 1.3+\" format is the same as phred-64" << endl;
	    exit(2);
	}
		
	int temp_intron_int;
	temp = run_cfg->getVal("max_intron");
	temp_intron_int = atoi(temp.c_str());
	if (temp_intron_int > 0) {
	    max_intron = temp_intron_int;
	}
		
	temp = run_cfg->getVal("min_intron");
	temp_intron_int = atoi(temp.c_str());
	if (temp_intron_int > 50) {
	    min_intron = temp_intron_int;
	}else if (temp_intron_int == 0) {
			
	}else {
	    cout << "ERROR: \"min_intron\" option should be greater than 50. " << endl;
	    cout << "Note that this is not the minimum search distance, SpliceMap will find small introns. :)" << endl;
	    
	    debug_out << "ERROR: \"min_intron\" option should be greater than 50. " << endl;
	    debug_out << "Note that this is not the minimum search distance, SpliceMap will find small introns. :)" << endl;
			
	    exit(2);
	}

	if (max_intron>0 && min_intron>0){
	    if (min_intron > max_intron) {
		cout << "Warning: minimum intron parameter greater than maximum intron parameter... " << endl;
		debug_out << "Warning: minimum intron parameter greater than maximum intron parameter... " << endl;
				
		max_intron = min_intron;
	    }
	}
		
	temp = run_cfg->getVal("min_extend_req");
	min_extend_req = atoi(temp.c_str()); 
		
	if (min_extend_req < 0 || min_extend_req > 25) {
	    cout << "ERROR: valid values for min_extend_req are 0-25 not " << temp << endl;
	    debug_out << "ERROR: valid values for min_extend_req are 0-25 not " << temp << endl;
			
	    exit(2);
	}
		
	temp = run_cfg->getVal("num_chromosome_together");
	num_chromosome_together = atoi(temp.c_str()); 
		
	if (num_chromosome_together <= 0)
	    num_chromosome_together = 1;
		
	delete run_cfg;
			
    } else {
	print_usage_and_exit();
    }
	
    if (cufflinks && print_sam) {
	cout << "ERROR: Please choose either -cuff or -sam. " << endl;
	debug_out << "ERROR: Please choose either -cuff or -sam. " << endl;
		
	cout << "-cuff will output a SAM file in a format compatible with cufflinks" << endl;
	cout << "-sam will output a proper SAM file" << endl;
	print_usage_and_exit();
    }
	
    cout << "__________" << endl;
    debug_out << "__________" << endl;
	
    cout << "Temp directory:   " << temp_path << endl;
    debug_out << "Temp directory:   " << temp_path << endl;
	
    cout << "Maximum number of mismatches allowed in 25-mer seed: " << num_seed_mismatch << endl;
    debug_out << "Maximum number of mismatches allowed in 25-mer seed: " << num_seed_mismatch << endl;
	
    cout << "Maximum number of mismatches allowed in full read: " << num_read_mismatch << endl;
    debug_out << "Maximum number of mismatches allowed in full read: " << num_read_mismatch << endl;
	
    cout << "Maximum number of bases SpliceMap is allowed to clip: " << max_clip_allowed << endl;
    debug_out << "Maximum number of bases SpliceMap is allowed to clip: " << max_clip_allowed << endl;
	
    cout << "Mapper used: ";
    debug_out << "Mapper used: ";
    cout << "bowtie" << endl;
    debug_out << "bowtie" << endl;
    
    cout << "(25th-percentile) intron size: " << min_intron << endl;
    cout << "(99th-percentile) intron size: " << max_intron << endl;
    debug_out << "(25th-percentile) intron size: " << min_intron << endl;
    debug_out << "(99th-percentile) intron size: " << max_intron << endl;
	
    cout << "Annotations path: " << ref_path_filename.first << " name: " << ref_path_filename.second << endl;
    cout << "Package path:     " << package_path_filename.first << " name: " << package_path_filename.second << endl;
    debug_out << "Annotations path: " << ref_path_filename.first << " name: " << ref_path_filename.second << endl;
    debug_out << "Package path:     " << package_path_filename.first << " name: " << package_path_filename.second << endl;
	
    cout << "Read format: " << read_format << endl;
    debug_out << "Read format: " << read_format << endl;
	
    if (read_format == "FASTQ") {
	cout << "Quality format: ";
	debug_out << "Quality format: ";
	if (quality_format == 0) {
	    cout << "phred-33" << endl;
	    debug_out << "phred-33" << endl;
	}else if (quality_format == 1){
	    cout << "phred-64" << endl;
	    debug_out << "phred-64" << endl;
	}else {
	    cout << "solexa" << endl;
	    debug_out << "solexa" << endl;
	}
    }
	
    cout << "Number of chromosomes to run together: " << num_chromosome_together << endl;
    debug_out << "Number of chromosomes to run together: " << num_chromosome_together << endl;
	
    if(print_sam){
	cout << "Will print regular SAM file" << endl;
	debug_out << "Will print regular SAM file" << endl;
    }
    if(cufflinks){
	cout << "Will print Cufflinks compatible SAM file" << endl;
	debug_out << "Will print Cufflinks compatible SAM file" << endl;
    }

    for(int j = 0;j<2;j++){
	cout << "Reads List " << (j+1) << ":" << endl;
	debug_out << "Reads List " << (j+1) << ":" << endl;
	for (unsigned int i = 0; i<reads_filename[j].size(); i++) {
	    cout << reads_filename[j][i] << endl;
	    debug_out << reads_filename[j][i] << endl;
	}
    }
	
    if (num_seed_mismatch<0 || num_seed_mismatch > 2){
	cout << "I'm sorry, SpliceMap " << VERSION << " supports a maximum of 2 mismatchs in each seed" << endl;
	debug_out << "I'm sorry, SpliceMap " << VERSION << " supports a maximum of 2 mismatchs in each seed" << endl;
	exit(1);
    }
	
    if (num_read_mismatch<0){
	cout << "ERROR: Read mismatch number negative: " << num_read_mismatch << endl;
	debug_out << "ERROR: Read mismatch number negative: " << num_read_mismatch << endl;
	exit(1);
    }
	
    if (max_clip_allowed<0){
	cout << "ERROR: Maximum clipping number negative: " << max_clip_allowed << endl;
	debug_out << "ERROR: Maximum clipping number negative: " << max_clip_allowed << endl;
	exit(1);
    }

    if (quasr_mode == GENERATE25MERS) {
		
	cout << "Preparing the reads!..." << endl;
	debug_out << "Preparing the reads!..." << endl;
		
	if(fullread_length > 0){
	    cout << "Using at most bases " << (head_clip_length+1) << " to " << fullread_length << " (inclusive) of each read for mapping." << endl ;
	    cout << "At most " << (fullread_length-head_clip_length) <<  " bases in total." << endl;
	    debug_out << "Using at most bases " << (head_clip_length+1) << " to " << fullread_length << " (inclusive) of each read for mapping." << endl ;
	    debug_out << "At most " << (fullread_length-head_clip_length) <<  " bases in total." << endl;
			
	}else { 
	    cout << "Bases removed from front: " << head_clip_length << endl;
	    debug_out << "Bases removed from front: " << head_clip_length << endl;
	    cout << "Using as many bases as possible." << endl;
	    debug_out << "Using as many bases as possible." << endl;
	}


	/*
	  the following block will create files:
	    temp_path/read_1_1       (sequence reads, where read_i_j refers to read j (e.g. j=2 for the second read in PE) in the i-th input file
	    temp_path/read_1_1.names (sequence names (FASTA/FASTQ) or "R" (RAW) with appended '[read_idx]', where read_idx is a sequential int value)
	    temp_path/read_1_1.qual  (sequence qualities in phred33 format, '0' for FASTA/RAW)

	  REMARK
	    create these files ourselves from R?
	*/
	if (read_format.compare("FASTQ") == 0) {
				
	    int read_count = 1;
	    vector<string>::iterator reads_list_it[2];
	    reads_list_it[0]= reads_filename[0].begin();
	    reads_list_it[1]= reads_filename[1].begin();
			
	    uint_fast32_t read_idx = 1;
			
	    while (reads_list_it[0] != reads_filename[0].end()) { 
		ofstream temp_outfile[2];
		ofstream temp_name_outfile[2];
		ofstream temp_quality_outfile[2];
				
		ifstream temp_infile[2];
				
		for (int k = 0; k<num_pair; k++) {
					
		    temp_infile[k].open(reads_list_it[k]->c_str(), ios::in);
					
		    string new_filename = temp_path+"read_"+IntToStr(read_count)+"_"+IntToStr(k+1);
		    temp_outfile[k].open(new_filename.c_str(), ios::out);
					
		    *(reads_list_it[k]) = new_filename;
					
		    temp = new_filename + ".names";
		    temp_name_outfile[k].open(temp.c_str(), ios::out);
		    temp = new_filename + ".quals";
		    temp_quality_outfile[k].open(temp.c_str(), ios::out);
		}
				
		uint_fast32_t line_num = 1;
				
		bool line_good[2] = {false,false};
		string line_name[2] = {"",""};
		string line_seq[2] = {"",""};
		string line_qual[2] = {"",""};
				
		while (!temp_infile[0].eof()) {  
					
		    bool end_reached = false;
					
		    for (int k = 0; k<num_pair; k++) {
						
			getline(temp_infile[k], line);
			trim2(line);
						
			if (line.length() == 0) {
			    end_reached = true;
							
			    continue;
			}
						
			if (line_num%4 == 1) {
			    line_good[k] = false; 
							
			    if (line[0] != '@') {
				cout << "Malformed FASTQ reads file: " << *(reads_list_it[k]) << endl;
				cout << "Line: " << line_num << endl;
				cout << line << endl;
				cout << "Should begin with an @" << endl;
				debug_out << "Malformed FASTQ reads file: " << *(reads_list_it[k]) << endl;
				debug_out << "Line: " << line_num << endl;
				debug_out << line << endl;
				debug_out << "Should begin with an @" << endl;
				exit(1);
			    }else {
				size_t name_end_loc = line.find_first_of(" \t/", 1);
				line_name[k] = line.substr(1, name_end_loc-1); 
								
			    }
							
			}else if(line_num%4 == 2){
			    if (line.length() > head_clip_length) {
				if (fullread_length > 0){
				    line = line.substr(head_clip_length, fullread_length-head_clip_length);   
				}else if(head_clip_length > 0 ){
				    line = line.substr(head_clip_length);   
				}else {
									
				}
								
				if (!(line[0] == 'A' || line[0] == 'C' || line[0] == 'T' || line[0] == 'G'
				      || line[0] == 'a' || line[0] == 'c' || line[0] == 't' || line[0] == 'g'
				      || line[0] == 'N' || line[0] == 'n')) {
				    cout << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
				    cout << "Ignoring line" << endl;
				    cout << line << endl;
				    debug_out << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
				    debug_out << "Ignoring line" << endl;
				    debug_out << line << endl;
									
				    line_good[k] = false;
				}else {
				    if (line.length() >= 50) {
					line_seq[k] = line;
										
					line_good[k] = true;
				    }else {
					cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
					debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
				    }
				}
								
			    }else {
				cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
				debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
			    }
							
			}else if (line_num%4 == 3) {
			    if (line[0] != '+') {
				cout << "Malformed FASTQ reads file: " << *(reads_list_it[k]) << endl;
				cout << "Line: " << line_num << endl;
				cout << line << endl;
				cout << "Should begin with a +" << endl;
				debug_out << "Malformed FASTQ reads file: " << *(reads_list_it[k]) << endl;
				debug_out << "Line: " << line_num << endl;
				debug_out << line << endl;
				debug_out << "Should begin with a +" << endl;
				exit(1);
			    }
			}else {
							
			    if (quality_format == 0) { 
								
			    }else if(quality_format == 1){
				phred642phred33(line);
			    }else { 
				solexa2phred33(line);
			    }
							
			    if (line.length() > head_clip_length) {
				if (fullread_length > 0){
				    line = line.substr(head_clip_length, fullread_length-head_clip_length);   
				}else if(head_clip_length > 0){
				    line = line.substr(head_clip_length);   
				}else {
				    
				}
			    }
							
			    line_qual[k] = line;  
							
			}
						
		    }

		    if (!end_reached && line_good[0] && (line_good[1]||num_pair==1) && line_num%4 == 0) { 
						
			for (int k = 0; k<num_pair; k++) {
			    temp_quality_outfile[k] << line_qual[k] << "\n"; 
							
			    temp_outfile[k] << line_seq[k] << "\n";  
			    temp_name_outfile[k] << line_name[k] <<"[" << read_idx << "]\n";
			}
		    }
					
		    line_num++;
		    read_idx++;
		}
				
		for (int k = 0; k<num_pair; k++) {
					
		    temp_infile[k].close();
		    temp_infile[k].clear();
		    temp_outfile[k].close();
		    temp_outfile[k].clear();
		    temp_name_outfile[k].close();
		    temp_name_outfile[k].clear();
		    temp_quality_outfile[k].close();
		    temp_quality_outfile[k].clear();
					
		    (reads_list_it[k])++;
		}
				
		read_count++;
				
	    }
			
	}else if (read_format.compare("FASTA") == 0) {
			
	    uint_fast32_t read_count = 1;
	    vector<string>::iterator reads_list_it[2];
	    reads_list_it[0]= reads_filename[0].begin();
	    reads_list_it[1]= reads_filename[1].begin();
			
	    uint_fast32_t read_idx = 1;
			
	    while (reads_list_it[0] != reads_filename[0].end()) { 
		ofstream temp_outfile[2];
		ofstream temp_name_outfile[2];
		ofstream temp_quality_outfile[2];
				
		ifstream temp_infile[2];
				
		for (int k = 0; k<num_pair; k++) {
					
		    temp_infile[k].open(reads_list_it[k]->c_str(), ios::in);
					
		    string new_filename = temp_path+"read_"+IntToStr(read_count)+"_"+IntToStr(k+1);
		    temp_outfile[k].open(new_filename.c_str(), ios::out);
					
		    *(reads_list_it[k]) = new_filename;
					
		    temp = new_filename + ".names";
		    temp_name_outfile[k].open(temp.c_str(), ios::out);
		    temp = new_filename + ".quals";
		    temp_quality_outfile[k].open(temp.c_str(), ios::out);
					
		}
				
		uint_fast32_t line_num = 1;
				
		bool line_good[2] = {false,false};
		string line_name[2] = {"",""};
		string line_seq[2] = {"",""};
				
		while (!temp_infile[0].eof()) {  
					
		    bool end_reached = false;
					
		    for (int k = 0; k<num_pair; k++) {
								
			getline(temp_infile[k], line);
			trim2(line);
						
			if (line.length() == 0) {
			    end_reached = true;
			    continue;
			}
						
			if (line_num%2 == 1) {
			    line_good[k] = false; 
							
			    if (line[0] != '>') {
				cout << "Malformed FASTA reads file: " << *(reads_list_it[k]) << endl;
				cout << "Line: " << line_num << endl;
				cout << line << endl;
				cout << "Should begin with an >" << endl;
				debug_out << "Malformed FASTA reads file: " << *(reads_list_it[k]) << endl;
				debug_out << "Line: " << line_num << endl;
				debug_out << line << endl;
				debug_out << "Should begin with an >" << endl;
				exit(1);
			    }else {
				size_t name_end_loc = line.find_first_of(" \t/", 1);
				line_name[k] = line.substr(1, name_end_loc-1); 
								
			    }
							
			}else if(line_num%2 == 0){
							
			    if (line.length() > head_clip_length) {
				if (fullread_length > 0){
				    line = line.substr(head_clip_length, fullread_length-head_clip_length);   
				}else if(head_clip_length > 0){
				    line = line.substr(head_clip_length);   
				}else {
									
				}
								
				if (!(line[0] == 'A' || line[0] == 'C' || line[0] == 'T' || line[0] == 'G'
				      || line[0] == 'a' || line[0] == 'c' || line[0] == 't' || line[0] == 'g'
				      || line[0] == 'N' || line[0] == 'n')) {
				    cout << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
				    cout << "Ignoring line" << endl;
				    cout << line << endl;
				    debug_out << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
				    debug_out << "Ignoring line" << endl;
				    debug_out << line << endl;
									
				    line_good[k] = false;
				}else {
									
				    if (line.length() >= 50) {
					line_seq[k] = line;
					
					line_good[k] = true;
				    }else {
					cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
					debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
				    }
				}
								
			    }else {
				cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
				debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
			    }
			    
			}

		    }
					
		    if (!end_reached && line_good[0] && (line_good[1]||num_pair==1) && line_num%2 == 0) { 
						
			for (int k = 0; k<num_pair; k++) {
			    temp_quality_outfile[k] << (char)0 << "\n"; 
							
			    temp_outfile[k] << line_seq[k] << "\n";  
			    temp_name_outfile[k] << line_name[k] << "[" << read_idx << "]\n";
			}
		    }
					
		    line_num++;
		    read_idx++;
		}
				
		for (int k = 0; k<num_pair; k++) {
					
		    temp_infile[k].close();
		    temp_infile[k].clear();
		    temp_outfile[k].close();
		    temp_outfile[k].clear();
		    temp_name_outfile[k].close();
		    temp_name_outfile[k].clear();
		    temp_quality_outfile[k].close();
		    temp_quality_outfile[k].clear();
		    
		    (reads_list_it[k])++;
		}
				
		read_count++;
				
	    }
			
	}else if (read_format.compare("RAW") == 0) {
				
	    uint_fast32_t read_count = 1;
	    vector<string>::iterator reads_list_it[2];
	    reads_list_it[0]= reads_filename[0].begin();
	    reads_list_it[1]= reads_filename[1].begin();
			
	    uint_fast32_t read_idx = 1;
			
	    while (reads_list_it[0] != reads_filename[0].end()) { 
		ofstream temp_outfile[2];
		ofstream temp_name_outfile[2];
		ofstream temp_quality_outfile[2];
		
		ifstream temp_infile[2];
				
		for (int k = 0; k<num_pair; k++) {
							
		    temp_infile[k].open(reads_list_it[k]->c_str(), ios::in);
					
		    string new_filename = temp_path+"read_"+IntToStr(read_count)+"_"+IntToStr(k+1);
		    temp_outfile[k].open(new_filename.c_str(), ios::out);
					
		    *(reads_list_it[k]) = new_filename;
							
		    temp = new_filename + ".names";
		    temp_name_outfile[k].open(temp.c_str(), ios::out);
		    temp = new_filename + ".quals";
		    temp_quality_outfile[k].open(temp.c_str(), ios::out);
					
		}
				
		uint_fast32_t line_num = 1;
				
		bool line_good[2] = {false,false};
				
		string line_seq[2] = {"",""};
				
		while (!temp_infile[0].eof()) {  
					
		    bool end_reached = false;
					
		    for (int k = 0; k<num_pair; k++) {

			getline(temp_infile[k], line);
			trim2(line);
						
			if (line.length() == 0) {
			    end_reached = true;
			    continue;
			}
						
			if (line.length() > head_clip_length) {
			    if (fullread_length > 0){
				line = line.substr(head_clip_length, fullread_length-head_clip_length);   
			    }else if(head_clip_length > 0){
				line = line.substr(head_clip_length);   
			    }else {
								
			    }
							
			    if (!(line[0] == 'A' || line[0] == 'C' || line[0] == 'T' || line[0] == 'G'
				  || line[0] == 'a' || line[0] == 'c' || line[0] == 't' || line[0] == 'g'
				  || line[0] == 'N' || line[0] == 'n')) {
				cout << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
				cout << "Ignoring line" << endl;
				cout << line << endl;
				debug_out << "Warning: Bad character in sequencer reads line " << *(reads_list_it[k]) << endl;
				debug_out << "Ignoring line" << endl;
				debug_out << line << endl;
								
				line_good[k] = false;
			    }else {
								
				if (line.length() >= 50) {
				    line_seq[k] = line;
									
				    line_good[k] = true;
				}else {
				    cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
				    debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(less than 50bp after clipping)" << endl; 
				}
			    }
			    
			}else {
			    cout << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
			    debug_out << "Warning: excluding line " << (reads_list_it[k])->c_str() << "|" <<  line_num << "(line shorter or equal to head_clip_length)" << endl;
			}
		    }
					
		    if (!end_reached && line_good[0] && (line_good[1]||num_pair==1)) { 
						
			for (int k = 0; k<num_pair; k++) {
			    temp_quality_outfile[k] << (char)0 << "\n"; 
							
			    temp_outfile[k] << line_seq[k] << "\n";  
			    temp_name_outfile[k] << "R["<< read_idx << "]\n";
			}
		    }
					
		    line_num++;
		    read_idx++;
		}
				
		for (int k = 0; k<num_pair; k++) {
					
		    temp_infile[k].close();
		    temp_infile[k].clear();
		    temp_outfile[k].close();
		    temp_outfile[k].clear();
		    temp_name_outfile[k].close();
		    temp_name_outfile[k].clear();
		    temp_quality_outfile[k].close();
		    temp_quality_outfile[k].clear();
		    
		    (reads_list_it[k])++;
		}
				
		read_count++;
	    }
	}

	
	/*
	  this block writes out all input files into the file:
	  temp_path/reads_list1
	  temp_path/reads_list2  (only if paired end)
	*/
	ofstream reads_list_temp_file;
	temp = temp_path+"reads_list1";
	reads_list_temp_file.open(temp.c_str(), ios::out);
	readslist_it = reads_filename[0].begin();
	while (readslist_it != reads_filename[0].end()) {
	    reads_list_temp_file << *readslist_it << endl;
	    readslist_it++;
	}
	reads_list_temp_file.close();
	reads_list_temp_file.clear();

	if(num_pair == 2){
	    temp = temp_path+"reads_list2";
	    reads_list_temp_file.open(temp.c_str(), ios::out);
	    readslist_it = reads_filename[1].begin();
	    while (readslist_it != reads_filename[1].end()) {
		reads_list_temp_file << *readslist_it << endl;
		readslist_it++;
	    }
	    reads_list_temp_file.close();
	    reads_list_temp_file.clear();
	}
		
	
	/*
	  extraction of 25-mers and creation of index files, output files are:
	    temp_path/map_filename  (map_filename is set to "25mers.map" in params.cpp)
	    temp_path/reads_1_1.index (same number of rows as "read_1_1" with sequences, with the length of each sequence)
 
	    25-mers from multiple input files, and for paired-end experiments, 25-mers from both reads are written to a single map_filename
	    work is done in output_index()

	    designsuffix() in SpliceMap_utils_QuasR.cpp defined the number and location of the 25-mers
	*/
		
	cout << "Extracting 25-mers... "  << endl;
	debug_out << "Extracting 25-mers... "  << endl;
		
	gettimeofday(&temp_start_tv, NULL);
			
	ofstream map_out;
		
	temp = temp_path+map_filename;
	map_out.open(temp.c_str(), ios::out);
			
	readslist_it = reads_filename[0].begin();
	while (readslist_it != reads_filename[0].end()) {

	    output_index(*readslist_it,map_out);
	    readslist_it++;
	}
		
	if(num_pair == 2){
	    readslist_it = reads_filename[1].begin();
	    while (readslist_it != reads_filename[1].end()) {
				
		output_index(*readslist_it,map_out);			
		readslist_it++;
	    }
	}
		
	map_out << flush;
	map_out.close();
	map_out.clear();
		
	gettimeofday(&tv, NULL);
	cout<<"Total 25-mer extraction section time: "<<diffclock(temp_start_tv,tv)<<" s."<<endl;
	debug_out<<"Total 25-mer extraction section time: "<<diffclock(temp_start_tv,tv)<<" s."<<endl;
	gettimeofday(&temp_start_tv, NULL);
    } // end quasr_mode == GENERATE25MERS



    if (quasr_mode == INDEX25MERALIGNMENTS) {

	/*
	  this block was added for modular running mode used by QuasR
	    -> it renames the reads_filename to the one in RAW format stored in temp_path
               (used to be done during conversion of FASTA/FASTQ to RAW above)
	*/
	int read_count = 1;
	vector<string>::iterator reads_list_it[2];
	reads_list_it[0]= reads_filename[0].begin();
	reads_list_it[1]= reads_filename[1].begin();
	while (reads_list_it[0] != reads_filename[0].end()) { 
	    for (int k = 0; k<num_pair; k++) {
		string new_filename = temp_path+"read_"+IntToStr(read_count)+"_"+IntToStr(k+1);
		*(reads_list_it[k]) = new_filename;
		(reads_list_it[k])++;
	    }
	    read_count++;
	}

	string output_file_path = temp_path+map_filename + ".out";
	
	/*
	  now, the extracted 25-mers (temp_path/25_mers.map) were mapped to the genome using the selected aligner and creates the file:
  	    temp_path/25_mers.map.out

	  the alignments are now generated by QuasR

	  bowtie alginments (using version 0.12.7) were done with parameters ("-y" only if 'tryhard'=='yes'):
	    "-y -S -k "+IntToStr(max_multi_hit)+" -m "+IntToStr(max_multi_hit)+" -v 2 -r -p "+IntToStr(num_threads)+" --best --strata"

	  initial bowtie alignments were written to temp_path/25mers.map.out_unsorted and then sorted using "sortsam":
	    sortsam -idx temp_path/25_mers.map.out_unsorted temp_path/25mers.map.out
	*/
	    

	/*
	  create a summary file for the alignments (the 'dot-t' file):
	    temp_path/25mers_map.out.t
	*/
	cout << "Generating .t file!..." << endl;
	debug_out << "Generating .t file!..." << endl;
	
	gettimeofday(&map_temp_start_tv, NULL);
		
	string infile = output_file_path;
	string outfile = output_file_path+".t";
	processbowtie(infile,outfile);
		
	gettimeofday(&tv, NULL);
	cout<<"Total .t creation section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
	debug_out<<"Total .t creation section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;

	    
	/*
	  this block reads the "dot-t" mapping summary file (temp_path/25mers.map.out.t) and the reads index file(s) (temp_path/read_1_1.index)
	  and combines them:
	    temp_path/read_1_1.index.t
	*/
	cout << "Creating mapping index!..." << endl;
	debug_out << "Creating mapping index!..." << endl;
	
	gettimeofday(&map_temp_start_tv, NULL);
	
	ifstream dot_t_file;
	temp = output_file_path+".t";
	dot_t_file.open(temp.c_str(), ios::in);
	
	if (!dot_t_file.is_open()) {
	    cout << "ERROR: Fatal error opening .t file: " << output_file_path+".t" << endl;
	    cout << "Have the reads been mapped?" << endl;
	    debug_out << "ERROR: Fatal error opening .t file: " << output_file_path+".t" << endl;
	    debug_out << "Have the reads been mapped?" << endl;
	    exit(2);
	}
	
	string dot_t_line;
	uint_fast32_t mapping_index = 1;
	getline(dot_t_file, dot_t_line);
	trim2(dot_t_line);
	if (dot_t_line.length() == 0) {
	    cout << "ERROR: There is something wrong with the mapping results... unexpected length" << endl;
	    exit(2);
	}
	
	for(int k = 0; k<num_pair;k++){

	    readslist_it = reads_filename[k].begin();
	    while (readslist_it != reads_filename[k].end()) {
		ifstream index_file;
		temp = *readslist_it + ".index";
		index_file.open(temp.c_str(), ios::in);
			
		if (!index_file.is_open()) {
		    cout << "ERROR: Fatal error opening index file ... "  << temp << endl;
		    debug_out << "ERROR: Fatal error opening index file ... "  << temp << endl;
		    exit(2);
		}
			
		ofstream out_index_file;
		temp = *readslist_it + ".index.t";
		out_index_file.open(temp.c_str(), ios::out);
			
		while (!index_file.eof()) {
		    getline(index_file, line);
		    trim2(line);
				
		    if (line.length() == 0)
			continue;
				
		    int curr_read_length = atoi(line.c_str());
		    out_index_file << curr_read_length << "\n";
				
		    vector<coord_t> suffix_list = designsuffix(curr_read_length);
				
		    for (unsigned int i = 0; i<suffix_list.size(); i++) {
							
			vector<dot_t_t> read_segment[2];
					
			for(int j = 0; j<2;j++){
			    bool reached_end = false; 
			    read_segment[j] = vector<dot_t_t>();
			    while (!reached_end && !dot_t_file.eof()) {
				
				vector<string> line_list = split(dot_t_line, '\t');
				unsigned int read_idx = (unsigned int)atoi(line_list[0].c_str());
							
				if (read_idx != mapping_index) {
				    reached_end = true;
				    mapping_index++;
								
				}else {

				    if (line_list.size() == 2) {
													
				    }else if (line_list.size() > 2){
									
					int curr_num_mismatch = -1;
					switch (line_list[1][1]) {  
					case '0':
					    curr_num_mismatch = 1;
					    break;
					case '1':
					    curr_num_mismatch = 2;
					    break;
					case '2':
					    curr_num_mismatch = 3;
					    break;
					default:
					    cout << "ERROR: num mismatch not matching... " << line << endl;
					    debug_out << "ERROR: num mismatch not matching... " << line << endl;
					    exit(2);
					    break;
					}
									
					int chr_pos = atoi(line_list[3].c_str());
					
					int direction = 0;
					switch (line_list[4][0]) {
					case 'F':
					    direction = 1;
					    break;
					case 'R':
					    direction = -1;
					    break;
					default:
					    cout << "ERROR: direction not matching... " << line << endl;
					    debug_out << "ERROR: direction not matching... " << line << endl;
					    exit(2);
					    break;
					}
									
					dot_t_t contents;
					contents.mismatch_dir = curr_num_mismatch*direction;
					contents.location = chr_pos;
					contents.chr_name = line_list[2];
														
					read_segment[j].push_back(contents);
									
				    }else {
					cout << "ERROR: Fatal error parsing .t file: " << line << endl;
					debug_out << "ERROR: Fatal error parsing .t file: " << line << endl;
					exit(2);
				    }
												
				    getline(dot_t_file, dot_t_line);
				    trim2(dot_t_line);
				}
			    }
			}
								
			for (int j = 0; j < 2; j++) {
						
			    out_index_file << read_segment[j].size() << "\n";
						
			    vector<dot_t_t>::iterator dot_t_it = read_segment[j].begin();
			    while (dot_t_it != read_segment[j].end()) {
				out_index_file << dot_t_it->chr_name << "\t" << dot_t_it->location << "\t" << (int)dot_t_it->mismatch_dir << "\n";
				    
				dot_t_it++;
			    }
			}
		    }
		}
					
		index_file.close();
		index_file.clear();
			
		out_index_file.close();
		out_index_file.clear();
			
		readslist_it++;
	    }
	}
	
	dot_t_file.close();
	dot_t_file.clear();
	
	gettimeofday(&tv, NULL);
	cout<<"Total mapping index creation section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
	debug_out<<"Total mapping index creation section execution time: "<<diffclock(map_temp_start_tv,tv)<<" s."<<endl;
    }

	
    /*
      here would start the calling of the SpliceMap binaries (as individual children threads created with fork())

      this is now done by QuasR

      there is one command for each chromosome, and the command lines are like:
        SpliceMap run.cfg chr19

      this creates one sam file for each chromosome:
        temp_path/chr19.fa_7.sam
    */
	

    /*
      merging of SpliceMap outputs (chromosomal sam files: temp_path/chr19.fa_7.sam) to create:
        temp_path/junction

      this is now done by QuasR

      this calls:
        amalgateSAM temp_path temp_path/junction

      amalgateSAM uses:
        temp_path/ref_list (names, paths and lengthes of chromosomal fasta files)
      and creates:
        temp_path/chr19.fa_7.sam_uniq
	temp_path/junction.sam (later sorted using "samsort", see below)
    */
	
	
    /*
      sorting of the amalgated junction sam file using "sortsam"
      no done by QuasR
    */

	
    /*
      create summary information from juction.sam file
      this is not used by QuasR anymore
	  
      work is done by:
        precipitateSAM temp_path/junction.sam temp_path/junction

	output files:
	  junction.bed (don't know where it comes from, did not find it in precipitateSAM.cpp)
	  junction_all.wig (later moved to outdir/coverage_all.wig) 
	  junction_up.wig (later moved to outdir/coverage_up.wig)
	  junction_down.wig (later moved to outdir/coverage_down.wig)
    */
	

    /*
      run additional post-processing (e.g. filtering of junctions, coloring, etc.)
      this is not used by QuasR anymore
    */
	

    debug_out.close();
    debug_out.clear();
	
	
    return 0;
}





inline int sam_char2int(char c)
{
    if (c == '0') {
	return 0;
    }else if (c == '1') {
	return 1;
    }else if(c == '2'){
	return 2;
    }
	
    cout << "Warning: Unexpected character in mismatch string, " << c << endl;
    return -1;
}


inline bool sam_is_mapped(unsigned int flag)
{
    return !(flag & 4);
}


inline string sam_get_direction(unsigned int flag)
{
    if (flag & 16){
	return "R";
    }else {
	return "F";
    }
}


inline void processbowtie(string &input_filename, string &output_filename)
{
    ifstream input;
    ofstream output;

    string temp;
    string line;
    list<sam_mismatch_t> sam_buf;
    sam_mismatch_t curr_sam;
    int curr_mismatch;
    int curr_read_idx;
    unsigned int curr_flag;
    int tabloc;
    int tabloc2;
    int min_mismatch;
	
    input.open(input_filename.c_str(), ios::in);
    if(!input.is_open()){
	cout << "ERROR: I'm sorry, cannot open the mapping output file, " << input_filename << endl;
	exit(1);
    }
	
    while (!input.eof()) {
	getline(input, line);
	trim2(line);
	if (line.length() == 0) {
	    continue;
	}
		
	if (line[0] == '@') {
			
	}else {
	    break;
	}
		
    }
	
    output.open(output_filename.c_str(), ios::out);
	
    tabloc = (int) line.find('\t');
    tabloc2 = (int) line.find('\t',tabloc+1);
    curr_flag = atoi(line.substr(tabloc+1, tabloc2-tabloc-1).c_str());
	
    curr_read_idx = atoi(line.substr(0,tabloc).c_str()) + 1;  
	
    if (sam_is_mapped(curr_flag)) {
		
	curr_mismatch = sam_char2int(line[line.length()-1]); 
		
	curr_sam.num_mismatch = curr_mismatch;
	curr_sam.data = line;
	sam_buf.push_back(curr_sam);
		
    }else {
	curr_mismatch = -1;
    }
	
	
    while (!input.eof()) {
	getline(input, line);
	trim2(line);
		
	if (line.length() == 0){
	    continue;  
	}
		
	tabloc = (int) line.find('\t');
	tabloc2 = (int) line.find('\t',tabloc+1);
	curr_flag = atoi(line.substr(tabloc+1, tabloc2-tabloc-1).c_str());
	int temp_curr_read_idx = atoi(line.substr(0,tabloc).c_str()) + 1;
		
	if (temp_curr_read_idx != curr_read_idx) {
				
	    if (curr_mismatch == -1) {
				
		output << curr_read_idx << "\t" << "NM" << endl;
				
	    }else{ 
		list<sam_mismatch_t>::iterator buf_it;
		
		buf_it = sam_buf.begin();
				
		bool curr_unique = (sam_buf.size()==1);
				
		while (buf_it != sam_buf.end()) {
		    temp = buf_it->data;
		    string curr_chr_name;
		    int curr_pos;
					
		    tabloc = (int) temp.find('\t'); 
		    tabloc2 = (int) temp.find('\t',tabloc+1); 
					
		    int temp_curr_flag = atoi(temp.substr(tabloc+1, tabloc2-tabloc-1).c_str());
					
		    tabloc = tabloc2;
		    tabloc2 = (int) temp.find('\t',tabloc+1); 
					
		    curr_chr_name = temp.substr(tabloc+1, tabloc2-tabloc-1).c_str();
					
		    tabloc = tabloc2;
		    tabloc2 = (int) temp.find('\t',tabloc+1); 
					
		    curr_pos = atoi(temp.substr(tabloc+1, tabloc2-tabloc-1).c_str());
					
		    int temp_curr_mismatch = sam_char2int(temp[temp.length()-1]); 
		    
		    if (curr_unique) {
			output << curr_read_idx << "\t" << "U" << temp_curr_mismatch << "\t" 
			       << curr_chr_name << "\t" << curr_pos << "\t" << sam_get_direction(temp_curr_flag) << endl;
		    } else{
			output << curr_read_idx << "\t" << "R" << temp_curr_mismatch << "\t" 
			       << curr_chr_name << "\t" << curr_pos << "\t" << sam_get_direction(temp_curr_flag) << endl;
						
		    }

		    buf_it++;
		}
	    }
			
	    sam_buf.clear();
	    min_mismatch = 2;
	    curr_read_idx = temp_curr_read_idx;
	}
		
	tabloc = (int) line.find('\t');
	tabloc2 = (int) line.find('\t',tabloc+1);
		
	curr_read_idx = atoi(line.substr(0,tabloc).c_str()) + 1;
		
	if (sam_is_mapped(curr_flag)) {
			
	    curr_mismatch = sam_char2int(line[line.length()-1]); 
			
	    curr_sam.num_mismatch = curr_mismatch;
	    curr_sam.data = line;
	    sam_buf.push_back(curr_sam);
			
	}else {
	    curr_mismatch = -1;
	}
    }
	
    if (curr_mismatch == -1) {
		
	output << curr_read_idx << "\t" << "NM" << endl;
		
    }else{ 
	list<sam_mismatch_t>::iterator buf_it;
		
	buf_it = sam_buf.begin();
		
	bool curr_unique = (sam_buf.size()==1);
		
	while (buf_it != sam_buf.end()) {
	    temp = buf_it->data;
	    string curr_chr_name;
	    int curr_pos;
			
	    tabloc = (int) temp.find('\t'); 
	    tabloc2 = (int) temp.find('\t',tabloc+1); 
			
	    int temp_curr_flag = atoi(temp.substr(tabloc+1, tabloc2-tabloc-1).c_str());
			
	    tabloc = tabloc2;
	    tabloc2 = (int) temp.find('\t',tabloc+1); 
			
	    curr_chr_name = temp.substr(tabloc+1, tabloc2-tabloc-1).c_str();
			
	    tabloc = tabloc2;
	    tabloc2 = (int) temp.find('\t',tabloc+1); 
			
	    curr_pos = atoi(temp.substr(tabloc+1, tabloc2-tabloc-1).c_str());
			
	    int temp_curr_mismatch = sam_char2int(temp[temp.length()-1]); 
			
	    if (curr_unique) {
		output << curr_read_idx << "\t" << "U" << temp_curr_mismatch << "\t" 
		       << curr_chr_name << "\t" << curr_pos << "\t" << sam_get_direction(temp_curr_flag) << endl;
	    } else{
		output << curr_read_idx << "\t" << "R" << temp_curr_mismatch << "\t" 
		       << curr_chr_name << "\t" << curr_pos << "\t" << sam_get_direction(temp_curr_flag) << endl;
	    }
			
	    buf_it++;
	}
    }
	
    input.close();
    input.clear();
    output.close();
    output.clear();
}


inline void print_usage_and_exit()
{
    cout << "usage: ./runSpliceMap run.cfg quasr_task" << endl;
    cout << "  run.cfg     --  Configuration options for this run, see comments in file for details" << endl;
    cout << "  quasr_task  --  task to perform as part of QuasR spliced alignment process" << endl;
    cout << "See website for further details" << endl;
    exit(0);
}


pair<string,string> get_path_and_filename(const string& path)
{
    string::size_type lastslash = path.rfind('/');
    if (lastslash == string::npos) {
	lastslash = path.rfind('\\');
    }
    if (lastslash == string::npos) {
	return pair<string,string>("",path);
    }
	
    return pair<string,string>(path.substr(0,lastslash+1),path.substr(lastslash+1,path.length()-lastslash));
}


void output_index(string reads_filename,ofstream &map_out)
{
    string temp;
    string line;
	
    ifstream reads_in;
	
    reads_in.open(reads_filename.c_str(), ios::in);
	
    if (!reads_in.is_open()) {
	cout << "ERROR: Fatal error opening reads file... " << reads_filename << endl;
	exit(1);
    }
	
    ofstream index_out;
    temp = reads_filename+".index";
    index_out.open(temp.c_str(), ios::out);
	
    while (!reads_in.eof()) {
	getline(reads_in, line);
		
	trim2(line);
		
	if (line.length() == 0)
	    continue;
		
	transform(line.begin(), line.end(), line.begin(), (int(*)(int))toupper); 
		
	index_out << line.length() << "\n";
		
	vector<coord_t> suffix = designsuffix((int)line.length());
		
	for (vector<coord_t>::iterator it = suffix.begin(); it < suffix.end(); it++) {
				
	    map_out << line.substr(it->first - 1,25) << '\n';    
	    map_out << line.substr(it->first+25 - 1,25) << '\n'; 
	}
    }
	
    index_out.close();
    index_out.clear();
    reads_in.close();
    reads_in.clear();
}

