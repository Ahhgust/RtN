#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <string>
#include <unordered_map>
#include <random>

#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/FastqReader.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"

// by default, pad 10 nucleotides 5' and 3'
// this approximates adding the primers back into the sequence
// so that the alignment can be made correctly (e.g.,
// insertions at the 5' and 3' ends can easily be misaligned).

#define PADDING 10

#define NUM_FASTQS 5

#define NO_IUPACS -1
#define MULTIPLE_IUPACS -2

using namespace SeqLib;
using namespace std;


struct Args {
  char *bedFile;
  char *refFasta;
  char *queryFasta;
  char *outdir;
};

char defaultDir[] = ".";

struct Amplicon {
  string chrom;
  int startPos;
  int stopPos;
  const string refSequence;
  string alignedSequence;
  int index;
  friend std::ostream& operator<<(std::ostream& o, const Amplicon &s) {
    o << s.chrom << "\t" << s.startPos << "\t" << s.stopPos << "\t" << s.alignedSequence << "\t" << s.index;
    return o;
  }
  
};

typedef std::vector<Amplicon> AmpliconVector;

bool
parseArgs(char **argv, int argc, Args &args) {

  if (argc < 4) {
    return false;
  }
  args.bedFile = argv[1];
  args.refFasta = argv[2];
  args.queryFasta = argv[3];
  args.outdir = defaultDir;
  if (argc < 5) {
    return true;
  }
  
  args.outdir = argv[4];
  return true;
}

/*
  Parses a bed file-- it only uses the first 3 columns of which
  and it fills out an Amplicon vector with:
  the bed coordinates as well as the reference sequence of the locus...
 */

bool
parseBed(const char* filename, AmpliconVector &vec, const RefGenome &ref) {
  ifstream file;
  file.open(filename);
  if (! file.is_open() )
    return false;
  string line;
  string chrom;
  string startPos;
  string stopPos;

  int startI, stopI;
  
  while (getline(file, line)) {
    if (line.length() == 0 || line[0] == '#')
      continue;

    // parse the first three records from the bed file
    istringstream is( line );
    is >> chrom >> startPos >> stopPos;
    
    //    vec.push_back( GenomicRegion(chrom, startPos, stopPos, hdr) );

    startI = atoi(startPos.c_str() );
    stopI = atoi(stopPos.c_str() );

    if (startI >= stopI) {
      cerr << "Malformed bed record " << endl << line << endl;
      return false;
    }
    
    string sequence = ref.QueryRegion(chrom, startI, stopI-1);
    vec.push_back( {chrom, startI, stopI, sequence, ""} );
  }
  return true;
}

// borrowed froM: https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
// tells you if a file exists.
// really, we want it to exist and be a directory...  but close enough!
bool
file_exists (char * name) {
  struct stat buffer;   
  return (stat (name, &buffer) == 0); 
}

/*
  Takes a string input (in)
  which has a 2-base ambiguity code at index idx
  and creates two strings (out1 and out2)
  which contain the 2-base versions
 */
bool
iupacSplit(const string &in, string &out1, string &out2, int idx) {

  char ambig = in[idx];
  out1 = out2 = in;
  if (ambig == 'R') {
    out1[idx] = 'A';
    out2[idx] = 'G';
  } else if (ambig == 'Y') {
    out1[idx] = 'T';
    out2[idx] = 'C';
  } else if (ambig == 'K') {
    out1[idx] = 'G';
    out2[idx] = 'C';
  } else if (ambig == 'M') {
    out1[idx] = 'A';
    out2[idx] = 'C';
  } else if (ambig == 'S') {
    out1[idx] = 'G';
    out2[idx] = 'C';
  } else if (ambig == 'W') {
    out1[idx] = 'A';
    out2[idx] = 'T';
  } else {
    out1[idx] = 'N';
    out2[idx] = 'N';
    return false;
  }
  return true;
}

// to string method for ints... joy
string
toString(int i) {

  std::string out_string;
  std::stringstream ss;
  ss << i;
  out_string = ss.str();
  return out_string;
}

void
writeFastq(Amplicon amp, char *outdir) {
  string dirname(outdir);
  string filename("/amp" + toString(amp.index) + ".fastq");
  
  ofstream outfile;
  outfile.open(dirname + filename);
  for (int i = 0; i < NUM_FASTQS; ++i) {
    outfile << "@amp" <<
      toString(amp.index) <<
      "." <<
      amp.chrom << ":" <<
      toString(amp.startPos) << "-" <<
      toString(amp.stopPos) << endl <<
      amp.alignedSequence <<
      endl << "+" << endl <<
      string( amp.alignedSequence.size(), '~') << endl;
  }
  outfile.close();

}


int
hasIUPAC(const string &stringy) {

  int iupacIndex;

  const char *beg, *s;
  iupacIndex = NO_IUPACS;

  beg = s = stringy.c_str();
  
  while (*s) {
    if (*s == 'R' ||
        *s == 'Y' ||
        *s == 'K' ||
        *s == 'M' ||
        *s == 'S' ||
        *s == 'W') {

      if (iupacIndex == NO_IUPACS)
        iupacIndex = s - beg;
      else
        return MULTIPLE_IUPACS; // more than 1 iupac code. can't split the haplotypes...
    }
    ++s;
  }
  
  return iupacIndex;

}


int
main(int argc, char** argv) {

  Args args;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> uniform(0.0,1.0);

  
  if (! parseArgs(argv, argc, args)) {

    cerr << "Usage: " << endl <<
      argv[0] << " bedFile referenceFasta queryFasta" << endl;
    return 1;
  }

  if (! file_exists(args.outdir)) {
    cerr << "The specified directory " << endl <<
      args.outdir << endl <<
      "doesn't exist. Please make it for me!" << endl;
    return 1;
  }
  
  RefGenome ref;
  if (!ref.LoadIndex(args.refFasta) ) {
    cerr << "Failed to load " << args.refFasta << endl;
    return 1;
  }

  AmpliconVector amps;
  
  if (!parseBed(args.bedFile, amps, ref)) {
    cerr << "Failed to parse the bed file" << endl;
    return 1;
  }

  int uniqCount=0;
  unordered_map<string, int> uniquify;

  // supports fastq or fasta...
  FastqReader fa( args.queryFasta);


  BamRecordVector results;
  UnalignedSequence seq;
  while (fa.GetNextSequence(seq) ) {

    BWAWrapper bwa;
    
    UnalignedSequenceVector usv;
    usv.push_back({ seq.Name, seq.Seq });
    
    bwa.ConstructIndex(usv);

    bool gotAll=true;
    for( std::vector<Amplicon>::iterator it = amps.begin();
         it != amps.end();
         ++it) {

      results.clear();
      bwa.AlignSequence(it->refSequence,
                        "foo",
                        results,
                        false,
                        0.99,
                        1);
      
      if (results.size()) {
        BamRecord hit = results[0];
        int startCoord = hit.Position();
        startCoord -= PADDING;
        if (startCoord < 0)
          startCoord = 0;
        // pad twice; once to make up for the deficit of -10, and again to add 10 additional bases
        // couldbe a little cleaner...
        int len = hit.GetCigar().NumReferenceConsumed() + PADDING + PADDING;

        string querySeq;
        if ((uint) len + (uint)startCoord <  seq.Seq.size() ) {
          querySeq = seq.Seq.substr(startCoord, len);
        } else {
          querySeq = seq.Seq.substr(startCoord);
        }

        int iupacIndex = hasIUPAC(querySeq);
        if (iupacIndex != NO_IUPACS) {
          if (iupacIndex == MULTIPLE_IUPACS) {
            gotAll=false; // we don't know the phase of IUPAC codes; 1 is fine, 2 is ambiguous as to what the haplotypes are (with 2 iupacs: 3 are possible; 4 are computationally available)
            break;
          }
          string s1, s2;
          iupacSplit(querySeq, s1, s2, iupacIndex);
          
          if (uniform( generator ) < 0.5) {
            querySeq = s1;
          } else {
            querySeq = s2;
          }
        }
        
        
        if (hit.ReverseFlag() ) {
          reverse(querySeq.begin(), querySeq.end());
          for (int i =0; i< len; ++i) {
            if (querySeq[i] == 'A')
              querySeq[i] = 'T';
            else if (querySeq[i] == 'T')
              querySeq[i] = 'A';
            else if (querySeq[i] == 'G')
              querySeq[i] = 'C';
            else if (querySeq[i] == 'C')
              querySeq[i] = 'G';
            else
              querySeq[i] = 'N';
          }
          
        }

        it->alignedSequence = querySeq;

        // new sequence
        // add it to the dictionary
        // and write the fastq
        if ( uniquify.count(querySeq) == 0) {
          uniquify[querySeq] = uniqCount;
          it->index = uniqCount;

          writeFastq(*it, args.outdir);
          
          ++uniqCount;
        } else {
          it->index = uniquify.at(querySeq);
        }
        
      } else {
        gotAll=false;
        break;
      }
    }

    if (gotAll) {
      for( std::vector<Amplicon>::iterator it = amps.begin();
           it != amps.end();
           ++it) {

        cout  << *it << "\t" << seq.Name  << endl;
      }
      
    }
  }
  
  
  return 0;
}


