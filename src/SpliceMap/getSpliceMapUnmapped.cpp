/*
  from intermediate output of SpliceMap, identify and output unmapped reads
  for paired-end experiments, unmapped reads are read pairs of which none of the reads were aligned

  Michael Stadler
  started on February 9, 2012
*/

#include <iostream>
#include <fstream>
#include <string>
#include <climits>
#include <cstdlib>

using namespace std;

class chrFile {
    int id, newid, start_pos, end_pos;
    const char *fname;
    ifstream *fh;
    string buffer;
public:
    chrFile(const char*);
    ~chrFile();
    const char* getFname() { return fname; }
    int getId() { return id; }
    int advance();
};

// constructor
chrFile::chrFile (const char* myfname) {
    //cout << "constructor called with " << fname << endl;
    // open file
    fname = myfname;
    fh = new ifstream;
    fh->open(fname, ifstream::in);
    id = -1;
    // advance
    advance();
}

// destructor
chrFile::~chrFile () {
    //cout << "destructor called" << endl;
    if(fh->is_open())
	fh->close();
    delete fh;
}

// read lines from fh until new id is found
int chrFile::advance() {
    if(id < INT_MAX) {
	newid = -1;
	while(newid<0 || newid==id) {
	    // read line(s)
	    getline (*fh, buffer, '\t'); // read to tab and store
	    fh->ignore(INT_MAX, '\n');   // ignore the rest of the line
	    if(fh->eof()) {
		newid = INT_MAX;
		break;
	    } else if(fh->fail() || fh->bad()) {
		cerr << "FATAL ERROR: could not read from file " << fname << endl;
		exit(1);
	    }

	    // skip header line
	    if(buffer[0] == '@')
		continue;

	    // extract id
	    start_pos = buffer.rfind('[')+1;
	    end_pos = buffer.rfind(']');
	    newid = atoi(buffer.substr(start_pos, end_pos-start_pos).c_str());
	}

	// set and return new id
	id = newid;
   }
    return id;
}

int extractNextId (const string &buffer) {
    static int start_pos;
    static int end_pos;

    start_pos = buffer.rfind('[')+1;
    end_pos = buffer.rfind(']');
    if(start_pos != string::npos && end_pos != string::npos)
	return atoi(buffer.substr(start_pos, end_pos-start_pos).c_str());
    else
	return -1;
}

int print_usage_and_exit () {
    cerr << "usage (single read experiments):" << endl
	 << "  getSpliceMapUnmapped S fname.seq fname.id fname.qual fname.out fname.chr1 [fname.chr2 ...]" << endl;
    cerr << "usage (paired-end experiments):" << endl
	 << "  getSpliceMapUnmapped P fname1.seq fname1.id fname1.qual fname2.seq fname2.id fname2.qual fname.out fname.chr1 [fname.chr2 ...]" << endl;
    exit(1);
}

int main(int argc, char** argv) {
    if(argc<7)
	print_usage_and_exit();

    // variable declaration, open files
    int i, id, nbIn = 0;
    bool alignmentFound;
    string seq_str1, id_str1, qual_str1, seq_str2, id_str2, qual_str2;
    char *fnameSeq1, *fnameId1, *fnameQual1, *fnameSeq2, *fnameId2, *fnameQual2, *fnameOut, **fnamesIn;
    chrFile **infiles;

    if(argv[1][0] == 'S') { // single read mode
	fnameSeq1  = argv[2];
	fnameId1   = argv[3];
	fnameQual1 = argv[4];
	fnameOut   = argv[5];
	fnamesIn   = argv + 6;
	nbIn       = argc-6;

	ifstream fhSeq1 (fnameSeq1, ifstream::in);
	if(! fhSeq1.good()) {
	    cerr << "FATAL ERROR opening " << fnameSeq1 << endl;
	    exit(1);
	}
	ifstream fhId1 (fnameId1, ifstream::in);
	if(! fhId1.good()) {
	    cerr << "FATAL ERROR opening " << fnameId1 << endl;
	    exit(1);
	}
	ifstream fhQual1 (fnameQual1, ifstream::in);
	if(! fhQual1.good()) {
	    cerr << "FATAL ERROR opening " << fnameQual1 << endl;
	    exit(1);
	}
	fstream fhOut (fnameOut, fstream::out);
	if(! fhOut.good()) {
	    cerr << "FATAL ERROR opening " << fnameOut << endl;
	    exit(1);
	}

	infiles = new chrFile*[nbIn];
	for(i=0; i<nbIn; i++)
	    infiles[i] = new chrFile(fnamesIn[i]);

	// iterate through input identifiers
	// read info for first sequence
	getline (fhSeq1, seq_str1);
	getline (fhId1,  id_str1);
	getline (fhQual1, qual_str1);
	while(fhSeq1.good() && fhId1.good() && fhQual1.good()) {
	    // extract identifier
	    id = extractNextId(id_str1);
	    if(id > 0) {
		// check if id is in one of the infiles
		alignmentFound = false;
		for(i=0; i<nbIn; i++)
		    if(infiles[i]->getId() == id) {
			alignmentFound = true;
			infiles[i]->advance();
		    }

		// output sequence if unmapped
		if(! alignmentFound) {
		    if( qual_str1[0] == 0 )
			qual_str1 = string(seq_str1.length(), 'I');
		    fhOut << id_str1 << "\t4\t*\t0\t0\t*\t*\t0\t0\t" << seq_str1 << "\t" << qual_str1 << endl;
		}

	    } else {
		cerr << "FATAL ERROR: id not found in " << id_str1 << endl;
		exit(1);
	    }

	    // read info for next sequence
	    getline (fhSeq1, seq_str1);
	    getline (fhId1,  id_str1);
	    getline (fhQual1, qual_str1);
	}

	// clean up
	fhSeq1.close();
	fhId1.close();
	fhQual1.close();
	fhOut.close();


    } else if(argv[1][0] == 'P') { // paired-end mode
	fnameSeq1  = argv[2];
	fnameId1   = argv[3];
	fnameQual1 = argv[4];
	fnameSeq2  = argv[5];
	fnameId2   = argv[6];
	fnameQual2 = argv[7];
	fnameOut   = argv[8];
	fnamesIn   = argv + 9;
	nbIn       = argc-9;

	ifstream fhSeq1 (fnameSeq1, ifstream::in);
	if(! fhSeq1.good()) {
	    cerr << "FATAL ERROR opening " << fnameSeq1 << endl;
	    exit(1);
	}
	ifstream fhId1 (fnameId1, ifstream::in);
	if(! fhId1.good()) {
	    cerr << "FATAL ERROR opening " << fnameId1 << endl;
	    exit(1);
	}
	ifstream fhQual1 (fnameQual1, ifstream::in);
	if(! fhQual1.good()) {
	    cerr << "FATAL ERROR opening " << fnameQual1 << endl;
	    exit(1);
	}
	ifstream fhSeq2 (fnameSeq2, ifstream::in);
	if(! fhSeq2.good()) {
	    cerr << "FATAL ERROR opening " << fnameSeq2 << endl;
	    exit(1);
	}
	ifstream fhId2 (fnameId2, ifstream::in);
	if(! fhId2.good()) {
	    cerr << "FATAL ERROR opening " << fnameId2 << endl;
	    exit(1);
	}
	ifstream fhQual2 (fnameQual2, ifstream::in);
	if(! fhQual2.good()) {
	    cerr << "FATAL ERROR opening " << fnameQual2 << endl;
	    exit(1);
	}
	fstream fhOut (fnameOut, fstream::out);
	if(! fhOut.good()) {
	    cerr << "FATAL ERROR opening " << fnameOut << endl;
	    exit(1);
	}

	infiles = new chrFile*[nbIn];
	for(i=0; i<nbIn; i++)
	    infiles[i] = new chrFile(fnamesIn[i]);

	// iterate through input identifiers
	// read info for first sequence
	getline (fhSeq1, seq_str1);
	getline (fhSeq2, seq_str2);
	getline (fhId1,  id_str1);
	getline (fhId2,  id_str2);
	getline (fhQual1, qual_str1);
	getline (fhQual2, qual_str2);
	while(fhSeq1.good() && fhId1.good() && fhQual1.good() && fhSeq2.good() && fhId2.good() && fhQual2.good()) {
	    // extract identifier
	    id = extractNextId(id_str1); // remark: could check if id_str2 contains identical id
	    if(id > 0) {
		// check if id is in one of the infiles
		alignmentFound = false;
		for(i=0; i<nbIn; i++)
		    if(infiles[i]->getId() == id) {
			alignmentFound = true;
			infiles[i]->advance();
		    }

		// output sequence if unmapped
		if(! alignmentFound) {
		    if( qual_str1[0] == 0 )
			qual_str1 = string(seq_str1.length(), 'I');
		    fhOut << id_str1 << "\t77\t*\t0\t0\t*\t*\t0\t0\t" << seq_str1 << "\t" << qual_str1 << endl;
		    if( qual_str2[0] == 0 )
			qual_str2 = string(seq_str2.length(), 'I');
		    fhOut << id_str2 << "\t141\t*\t0\t0\t*\t*\t0\t0\t" << seq_str2 << "\t" << qual_str2 << endl;
		}

	    } else {
		cerr << "FATAL ERROR: id not found in " << id_str1 << endl;
		exit(1);
	    }

	    // read info for next sequence
	    getline (fhSeq1, seq_str1);
	    getline (fhSeq2, seq_str2);
	    getline (fhId1,  id_str1);
	    getline (fhId2,  id_str2);
	    getline (fhQual1, qual_str1);
	    getline (fhQual2, qual_str2);
	}

	// clean up
	fhSeq1.close();
	fhId1.close();
	fhQual1.close();
	fhSeq2.close();
	fhId2.close();
	fhQual2.close();
	fhOut.close();


    } else {
	cerr << "ERROR: unknown mode '" << argv[1] << "', must be either 'S' or 'P'" << endl;
	print_usage_and_exit();
    }
	

    // clean up
    for(i=0; i<nbIn; i++)
	delete infiles[i];
    delete[] infiles;

    cout << "Finished." << endl;
    return 0;
}
