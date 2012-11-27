/*
  from intermediate outputs of getSpliceMapUnmapped ("unmapped.sam") and amalgamateSAM ("junction.sam", "unmapped2.sam"),
  produce "unmapped2.sam" the header from "junction.sam" and all context from the three sam files in order of the input reads

  Michael Stadler
  started on July 2, 2012
*/

#include <iostream>
#include <fstream>
#include <string>
#include <climits>
#include <cstdlib>
#include <queue>
//#include <vector>

using namespace std;


int print_usage_and_exit () {
    cerr << "usage:" << endl
	 << "  fuseReorder tmpdir" << endl;
    exit(1);
}


string extractChr(const string &line) {
    static size_t start_pos, end_pos;
    string chr("");

    start_pos = line.find('\t', 0);
    start_pos = line.find('\t', start_pos+1);
    if(start_pos != string::npos) {
        end_pos = line.find('\t', start_pos+1);
        if(end_pos != string::npos)
	    chr = line.substr(start_pos+1, end_pos-start_pos-1);
    }

    return chr;
}


class idLine { // stores alignments from file for a integer id
public:
    int id;        // aligment integer id
    int ifile;     // file index (origin of the alignment)
    string aln;    // SAM line(s) including line end characters
    idLine() {id=-1; ifile=-1; aln=""; }
    idLine(const int &newid, const int &newifile, const string &newaln) {id=newid; ifile=newifile; aln=newaln; }
    bool operator() (const idLine& lhs, const idLine&rhs) const {return (lhs.id>rhs.id);}
    void print() { cerr << id << " (" << ifile << ":\n" << aln << endl; }
};

class samFile { // ifstream to sam file (-region)
public:
    string fname;      // file name
    streampos start;   // start offset in file
    streampos end;     // end offset in file
    int ifile;         // file index
    idLine idline;     // aligment lines
    ifstream fh;       // input file handle
    string linebuffer; // buffer for getline()
    int linebufferId;  // id extracted from linebuffer
    samFile(const string &newfname, const streampos newstart, const streampos newend, const int newifile) {
	fname = newfname;
	start = newstart;
	end = newend;
	ifile = newifile;

	fh.open(fname.c_str(), ifstream::in | ifstream::binary);
	if(! fh.good()) {
	    cerr << "FATAL ERROR opening " << fname << endl;
	    exit(3);
	}

	if(end == (streampos)0) { // set end to the end of the file
	    fh.seekg(0, ios::end);
	    end = fh.tellg();
	}

	fh.seekg(start);
	this->readNextAlnIntoLinebuffer();
    }
    ~samFile() { fh.close(); }
    void readNextAlnIntoLinebuffer() {
	getline(fh, linebuffer);
	if(fh.good()) {
	    if(linebuffer[linebuffer.size()-1] == '\r')
		linebuffer.erase(linebuffer.size()-1,1);
	    this->extractAndRemoveId();
	} else {
	    linebuffer.clear();
	    linebufferId = INT_MAX;
	}
    }
    void extractAndRemoveId () {
	static size_t start_pos, end_pos;

	linebufferId = INT_MAX;

	start_pos = linebuffer.find('\t'); // first tab: end of identifier field
	if(start_pos != string::npos) {
	    end_pos = linebuffer.rfind(']', start_pos); // back to ']'

	    if(end_pos != string::npos) {
		start_pos = linebuffer.rfind('[', end_pos)+1; // back to '['

		if(start_pos != string::npos) {
		    linebufferId = atoi(linebuffer.substr(start_pos, end_pos-start_pos).c_str());
		    linebuffer.erase(start_pos-1, end_pos-start_pos+2);
		}
	    }
	}
    }
    bool loadAllAlignmentsForNextId() {
	int currId, newId;
	string alnlines;
	if(fh.good() && fh.tellg() < end) {

	    currId = linebufferId;
	    alnlines = linebuffer + '\n';
	    this->readNextAlnIntoLinebuffer();
	    newId = linebufferId;

	    while(currId == newId && fh.good() && fh.tellg() < end) {
		alnlines += (linebuffer + '\n');
		this->readNextAlnIntoLinebuffer();
		newId = linebufferId;
	    }

	    idline.id = currId;
	    idline.ifile = ifile;
	    idline.aln = alnlines;

	    return true;
	} else {
	    if(!linebuffer.empty()) {
		currId = linebufferId;
		idline.id = currId;
		idline.ifile = ifile;
		idline.aln = linebuffer + '\n';
		linebuffer.clear();
		return true;
	    } else {
		return false;
	    }
	}
    }
};


int main(int argc, char** argv) {

    if(argc<2)
	print_usage_and_exit();


    // variable declaration
    int i, nbChr;
    streampos currPos;
    vector<streampos> chrPos;
    vector<string> chrName;
    string currChr, newChr;
    string tmpdir(argv[1]);
    string fnameIn[3] = { tmpdir + "/junction.sam",
			  tmpdir + "/unmapped.sam",
			  tmpdir + "/unmapped2.sam" };
    string fnameOut = tmpdir + "/junction2.sam";
    string linebuffer;
    ifstream fhIn;
    fstream fhOut;
    priority_queue<idLine, vector<idLine>, idLine> queue; // stores alignments
    idLine topline;

    // open alignment and output file
    fhIn.open(fnameIn[0].c_str(),  ifstream::in | ifstream::binary);
    if(! fhIn.good()) {
	cerr << "FATAL ERROR opening " << fnameIn[0] << endl;
	exit(2);
    }
    fhOut.open(fnameOut.c_str(), fstream::out | fstream::binary);
    if(! fhOut.good()) {
	cerr << "FATAL ERROR opening " << fnameOut << endl;
	exit(2);
    }


    // copy header from junction.sam to junction2.sam
    while(fhIn.peek()=='@' && fhIn.good()) {
	getline(fhIn, linebuffer);
	if(linebuffer[linebuffer.size()-1] == '\r')
	    linebuffer.erase(linebuffer.size()-1,1);
	fhOut << linebuffer << endl;
    }


    // find chromosome block starts
    nbChr = 0;
    currChr = "";

    currPos = fhIn.tellg();
    getline(fhIn, linebuffer);
    while( fhIn.good() ) {
	newChr = extractChr(linebuffer);

	if(newChr != currChr) {
	    nbChr++;

	    chrPos.push_back(currPos);
	    chrName.push_back(newChr);
	    currChr = newChr;
	}

	currPos = fhIn.tellg();
	getline(fhIn, linebuffer);
    }
    fhIn.close();
    for(i=0; i<nbChr; i++) cout << chrName[i] << '\t' << chrPos[i] << endl;


    // open input streams
    samFile *samfiles[nbChr+2];
    for(i=0; i<nbChr; i++)
	samfiles[i] = new samFile(fnameIn[0], chrPos[i], i<nbChr-1 ? chrPos[i+1] : (streampos)0, i); // junction.sam (chromosome blocks)
    samfiles[nbChr]   = new samFile(fnameIn[1], 0, 0, nbChr);   // unmapped.sam
    samfiles[nbChr+1] = new samFile(fnameIn[2], 0, 0, nbChr+1); // unmapped2.sam


    // get next aligment from each file and initialize queue
    for(i=0; i<nbChr+2; i++)
	if(samfiles[i]->loadAllAlignmentsForNextId())
	    queue.push(samfiles[i]->idline);


    //    get next aligments, put on queue and output smallest
    while(!queue.empty()) {
	topline = queue.top();
	fhOut << topline.aln;
	i = topline.ifile;
	queue.pop();
	if(samfiles[i]->loadAllAlignmentsForNextId())
	    queue.push(samfiles[i]->idline);
    }


    // clean up
    fhOut.close();
    for(i=0; i<nbChr+2; i++)
	delete samfiles[i];

    cout << "Finished." << endl;
    return 0;
}
