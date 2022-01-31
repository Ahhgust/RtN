#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <random>
#include <algorithm>
#include <unordered_map>

#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/FastqReader.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"



using namespace SeqLib;
using namespace std;

// used when categorizing amplicons
// reads should neither be within amps (fully) or span multiple amps
#define WITHIN_AMP -2
#define MULTIPLE_AMPS -1

// when the alignment returns nothing to the Numt database
// this gives the default distance
#define NUMT_NONN_PLACEHOLDER 99

#define MITOLEN 16561

// taken from: https://stackoverflow.com/questions/5590381/easiest-way-to-convert-int-to-string-in-c
#define INT2STRING( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

char DEFAULT_CHROM[] = "chrM"; 

// tags for the BAM files...
std::string DEFAULT_ZTAG_HUMAN = "ZH";
std::string DEFAULT_ZTAG_NUMT  = "ZN";

#define DEFAULT_TRAINING_SIZE 100000

// BWA mapping parameters
#define DEFAULT_MISMATCH 4
#define DEFAULT_GAPOPEN 4
#define DEFAULT_GAPEXTEND 4
#define DEFAULT_BANDWIDTH 1000
#define DEFAULT_THREEPRIMECLIP 100000
#define DEFAULT_FIVEPRIMECLIP 100000

// 1% sequencing error is assumed if the qual string is empty
#define DEFAULT_SEQ_ERROR 0.01
// '5' in the ascii table corresponds to a phred score of 20 which is a 1% chance of error
#define DEFAULT_SEQ_CHAR 53

#define DEFAULT_MIN_MAP_QUALITY 4

// measured in bases of the genome
#define DEFAULT_MIN_READ_LEN 20

#define DEFAULT_MIN_LIKELIHOOD 1e-5


struct Options {
  char *humanFastaDbFilename;
  char *nonhumanFastaDbFilename;
  char *bamFilename;
  char *chrom;
  char *bedFilename;
  int gapOpen;
  int gapExtend;
  int mismatch;
  int bandwidth;
  int threePrimePenalty;
  int fivePrimePenalty;
  int minMappingQuality;
  bool ignoreIndels;
  bool verbose;
  bool train;
  bool writeOffTarget;
  bool scaleLikelihoodByReadlen;
  bool filterReadPairs;
  bool indexJump;
  bool removeDups;
  double minLikelihood;
  int minReadSize;
};

struct SummaryStat {
  string readName;
  int numMismatches; // # of mismatches to best alignment in HmtDB
  double readLikelihood; // from MAQ: product of mismatched bases' read qualities
  int genomeMidpoint; // where in the genome this read corresponds to
  double averageMismatchQuality; // similar to likelihood, but just the mean
  int numQ20Mismatches; // how many Q20+ mismatches were there?
  double meanPhred; // mean phred score
  bool isReverse; // + or - strand
  int numBases; // number of ALIGNED and MATCHING/MISMATCHING bases (no indels)
  int numIndels;
  int numQ20Indels;
  double meanIndelQuality;
  double readLikelihoodWithIndels;
  
  friend std::ostream& operator<<(std::ostream& o, const SummaryStat &s) {
    o << s.readName <<"\t" << s.numMismatches << "\t" << s.readLikelihood << "\t" << s.genomeMidpoint << "\t" <<
      s.averageMismatchQuality << "\t" << s.numQ20Mismatches << "\t" <<
      s.meanPhred << "\t" << s.isReverse << "\t" << s.numBases <<
      "\t" << s.numIndels << "\t" << s.numQ20Indels << "\t" << s.meanIndelQuality << "\t" << s.readLikelihoodWithIndels;
    return o;
  }

};


/*
  Perl-inspired print stderr + exit + usage statement
 */
void
die(const char * message) {
  if (message != NULL) {
    cerr << message << endl;

  }

  cerr << "How to run this program:" << endl <<
    "Required " << endl
       << "\t-h humanFasta" << endl
       << "\t-n numtFasta" << endl
       << "\t-b bam" << endl
    << endl << "Optional (with default)" << endl << endl
    << "\t-c chromosome (" << DEFAULT_CHROM << ") (ie, what we call the mitochondrial contig)" << endl
       << "\t-m minMappingQuality (" << DEFAULT_MIN_MAP_QUALITY << ")" << endl
    << "\t-t (output training data; summary statistics on the reads)" << endl
    << "\t-s bedFile (soft-clip reads to the amplicons specified in the bedFile)" << endl
    << "\t-w (writes off-target reads; defaults to FALSE)" << endl
    << "\t-l length ((" << DEFAULT_MIN_READ_LEN << ") the minimum read length, as mapped to the genome)" << endl
    << "\t-L likelihood (" << DEFAULT_MIN_LIKELIHOOD << ") the minimum read likelihood)" << endl
    << "\t-S (this scales the read's likelihood (log(likelihood)/read length); with -L must reflect this (e.g., the threshold must be negative; -0.05 == 1e-5 with 100bp reads))" << endl
    << "\t-p (this filters on read pairs; if either read fails the likelihood requirement, they both do)" << endl
    << "\t-i (default: FALSE, ignores indels in the likelhood function. Recommended for Ion sequencing)" << endl
    << "\t-j (default: FALSE, index jumps to the mito in the bam. Recommended for WGS data )" << endl
    << "\t-d (default: FALSE, Removes PCR/optical duplicates. Recommended for WGS data )" << endl
    
    << endl << "Read mapping parameters..." << endl
    << "\t-o gapOpen (" << DEFAULT_GAPOPEN << ")" << endl
    << "\t-e gapExtend (" << DEFAULT_GAPEXTEND << ")" << endl
    << "\t-w bandWidth (" << DEFAULT_BANDWIDTH << ")" << endl
    << "\t-x mismatchPenalty (" << DEFAULT_MISMATCH << ")" << endl

    << "\t-3 3-primeclipping (" << DEFAULT_THREEPRIMECLIP << ")" << endl
       << "\t-5 5-primeclipping (" << DEFAULT_FIVEPRIMECLIP << ")" << endl << endl;
  
  exit(1);
}

/*
  Boilerplate command-line interaction
 */
Options
parseOptions(char **argv) {

  bool gotMinlFlag=false;
  char f, *arg;
  Options opt = {NULL, NULL, NULL};
  opt.chrom=DEFAULT_CHROM;
  opt.bedFilename=NULL;
  opt.indexJump=false; // do we seek to the mito
  opt.removeDups=false;
  opt.filterReadPairs=false;
  opt.mismatch = DEFAULT_MISMATCH;
  opt.gapOpen = DEFAULT_GAPOPEN;
  opt.gapExtend = DEFAULT_GAPEXTEND;
  opt.bandwidth = DEFAULT_BANDWIDTH;
  opt.threePrimePenalty = DEFAULT_THREEPRIMECLIP;
  opt.fivePrimePenalty = DEFAULT_FIVEPRIMECLIP;
  opt.verbose = opt.ignoreIndels=false;
  opt.minMappingQuality = DEFAULT_MIN_MAP_QUALITY;
  opt.writeOffTarget=false;
  opt.minReadSize = DEFAULT_MIN_READ_LEN;
  opt.minLikelihood = DEFAULT_MIN_LIKELIHOOD;
  opt.scaleLikelihoodByReadlen = false;
  
  for (; *argv != NULL; ++argv) {
    arg = *argv;
    if (*arg != '-') {
      if (*(argv-1) == NULL)
        cerr << "Unexpected argument: (you specified a flag that required an argument... yet you gave no argument!) "<< endl;
      else
        cerr << "Unexpected argument: " << arg << endl;
      
      die(NULL);
    }
    // all flags have arguments
    f = arg[1]; // flag
    ++argv; // and its argument

    
    if (f == 'h') 
      opt.humanFastaDbFilename = *argv;
    else if (f == 'n')
      opt.nonhumanFastaDbFilename = *argv;
    else if (f == 'c')
      opt.chrom = *argv;
    else if (f == 'b')
      opt.bamFilename = *argv;
    else if (f == 's')
      opt.bedFilename = *argv;
    else if (f == 'o')
      opt.gapOpen = atoi(*argv);
    else if (f == 'e')
      opt.gapExtend = atoi(*argv);
    else if (f == 'w')
      opt.bandwidth = atoi(*argv);
    else if (f == '3')
      opt.threePrimePenalty = atoi(*argv);
    else if (f == '5')
      opt.fivePrimePenalty = atoi(*argv);
    else if (f == 'm')
      opt.minMappingQuality = atoi(*argv);
    else if (f == 'l')
      opt.minReadSize = atoi(*argv);
    else if (f == 'x')
      opt.mismatch = atoi(*argv);
    else if (f == 'L') {
      opt.minLikelihood = atof(*argv);
      gotMinlFlag=true;
    } else if (f == 'w') {
      opt.writeOffTarget = true;
      --argv; // only a flag; no argument.
    } else if (f == 'i') {
      opt.ignoreIndels = true;
      --argv; // only a flag; no argument.
    } else if (f == 'p') {
      opt.filterReadPairs=true;
      --argv;
    } else if (f == 't') {
      opt.train = true;
      --argv; // only a flag; no argument.
    } else if (f == 'v') {
      opt.verbose = true;
      --argv; // only a flag; no argument.
    } else if (f == 'j') {
      opt.indexJump = true;
      --argv; // only a flag; no argument.
    } else if (f == 'd') {
      opt.removeDups=true;
      --argv;
    } else if (f == 'S') {
      opt.scaleLikelihoodByReadlen = true;
      --argv;
      if (gotMinlFlag==false) // reset the default. 1e-5 was designed for 100bp reads. maps into -0.05
        opt.minLikelihood = log10(DEFAULT_MIN_LIKELIHOOD)/100.0;
    } else {
      cerr << "Unknown flag: " << f << endl;
      die(NULL);
    }
  }

  if (opt.humanFastaDbFilename == NULL) 
    die("I need a human fasta database of mitochondrial genomes! (-h humanFasta)!");
  
  if (opt.nonhumanFastaDbFilename == NULL) 
    die( "I need a human fasta database of mitochondrial genomes! (-n humanFasta!)");
  
  if (opt.nonhumanFastaDbFilename == NULL)
    die("I need a bam file! (-b bam)");

  if (opt.scaleLikelihoodByReadlen && opt.minLikelihood >0) {
    die("Incompatible parameters: if you want to scale the likelihood by the read length, the threshold is based on log10(likelihood)/read length; thus the threshold must be NEGATIVE");
  }

  if (!opt.scaleLikelihoodByReadlen && opt.minLikelihood <0) {
    die("Incompatible parameters: the likelihood threshold must be non-negative. Did you forget to include -S?");
  }

  return opt;
}

/*
  pure C style implementation used to find the substring postion of the last . in a filename...
  and it's smart; only uses the one characters since the last file separator (assumed to be either a dos-separator
  or a unix style. which isn't 100% correct, but it's close enough)
  
 */
int
findExtensionSubstringPosition(const char *c) {

  int lastDotPos=-1;
  int i = 0;
  while (*c) {
    if (*c == '.')
      lastDotPos = i;
    else if (*c == '/' || *c == '\\') { // previous '.' was in a directory name. assume either dir symbol and that people aren't jerks (putting a / in a filename)
      lastDotPos = -1;
    }
    ++i;
    ++c;
  }
  // no file extension. the substring wanted is the whole string
  if (lastDotPos==-1)
    return i;
  
  return lastDotPos;
}

/*
  Parses a bed file-- it only uses the first 3 columns of which
  this will fill out a vector vec
  GenomicRegionVectors are a typedef of a vector of GenomicRegions (from SeqLib)
 */

bool
parseBed(const char* filename, GenomicRegionVector &vec, const BamHeader &hdr) {
  ifstream file;
  file.open(filename);
  if (! file.is_open() )
    return false;
  string line;
  string chrom;
  string startPos;
  string stopPos;
  
  while (getline(file, line)) {
    if (line.length() == 0 || line[0] == '#')
      continue;

    // parse the first three records from the bed file
    istringstream is( line );
    is >> chrom >> startPos >> stopPos;
    
    vec.push_back( GenomicRegion(chrom, startPos, stopPos, hdr) );

  }
  return true;
}

/*
gets the start position in the read
i.e., the first nucleotide after all soft-clipping
*/
int
getStartSubstr(Cigar *ciggy) {

  vector<CigarField>::iterator itr = ciggy->begin();
  int startPositionRead=0;

  char type;
  uint32_t len;
  // infer the start/stop positions in the read
  // i.e., trim off the soft clipping

  // start with soft-clip 5'
  while (itr != ciggy->end() ) {
    type = itr->Type();
    len = itr->Length();
    
    if (type == 'S')
      startPositionRead += len;
    else if (type != 'H')
      break;

    ++itr;
  }
  return startPositionRead;
}


/*
Returns the number of nucleotides to be trimmed off the end
*/
int
getLastNClipped(Cigar *ciggy) {
    // and read the cigar backwards...
  int nOps = ciggy->size();
  int clipLastN=0;
  for (int i = nOps-1; i >= 0; --i) {
    CigarField f = (*ciggy)[i];
    if (f.Type() == 'S') {
      clipLastN += f.Length();
    } else if (f.Type() != 'H')
      break;
  }
  return clipLastN;
}

int
getStopSubstr(string &seq, Cigar *ciggy) {
  int lastN = getLastNClipped(ciggy);
  return seq.size() - lastN;
}


/*
  This computes the number of mismatches between two sequences;
  the number as taken as the sum of the number of substitution differences +
  the number of indel differences.
  Both have a cost of 1
  further, all indels have a flat cost of 1 (regardless of length)
  TODO:
  homopolymers!
 */

bool
getMismatches(BamRecord &r, const char* seq, const char *ref, int *mismatches, int *tot, bool ignoreIndels) {

  

  Cigar ciggy = r.GetCigar();
  *mismatches = *tot = 0;
  
  int i, len;
  Cigar::const_iterator c = ciggy.begin();
  if (c == ciggy.end()) // empty cigar?! I don't think this is legitimate, but you never know!
    return false;
  char t = c->Type();
  
  
  // scroll past the front-most soft/hard clipping tag
  // don't need to adjust indexing b/c the input sequence are pre-adjusted for soft clipping
  if (t == 'S' || t == 'H') 
    ++c;
  
  if (c == ciggy.end()) // a cigar that's only clipping...? also weird, but okay I suppose.
    return false;
  
  for ( ; c != ciggy.end(); ++c) {
    t = c->Type();
    len = c->Length();
    
    if (c->ConsumesReference() ) {
      if (c->ConsumesQuery() ) {
        for ( i = 0; i < len; ++i, ++ref, ++seq) {
          if (*ref != *seq) {
            ++(*mismatches);
          }
          ++(*tot);
        }
      } else { // consumes ref, not query ( deletion in query of clipping)
        if (! ignoreIndels && t == 'D') {
          ++(*mismatches);
          ++(*tot);
        }
        ref += len;
      }
    } else { // consumes query, not reference (insertion in query or clipping)
      if (! ignoreIndels && t == 'I') {
        ++(*mismatches);
        ++(*tot);
      }
      seq += len;
    }
  }


  return true;
}


/*
  This computes the number of mismatches (substitution only)
  based on the alignment
  AND 
  it fills out the mismatches array of their quality scores
  THe array should be made to tbe the size of the quals string.
 */

int
getWeightedMismatches(BamRecord &r, const char* seq, const char *ref, const char *quals, SummaryStat &stat, char qMax, Options &opt) {


  Cigar ciggy = r.GetCigar();

  int i, len;

  const char *qualsInit = quals;
  
  stat.readLikelihood= stat.readLikelihoodWithIndels=1.;
  stat.numBases=stat.numIndels=stat.numQ20Indels=0;
  stat.numMismatches= stat.numQ20Mismatches=0;
  stat.meanIndelQuality=stat.averageMismatchQuality = 0.;
  
  Cigar::const_iterator c = ciggy.begin();
  if (c == ciggy.end()) // empty cigar?! I don't think this is legitimate, but you never know!
    return -1;

  char t = c->Type();
  
  int mismatchQualsum=0;
  int sumPhred=0;
  
  // hard clipping... skip!
  if (t == 'H') 
    ++c;

  
  if (c == ciggy.end()) // a cigar that's only clipping...? also weird, but okay I suppose.
    return -1;
  

  for ( ; c != ciggy.end(); ++c) {
    t = c->Type();
    
    len = c->Length();
    
    if (c->ConsumesReference() ) {
      if (c->ConsumesQuery() ) {

        
        stat.numBases += len;
          
        for ( i = 0; i < len; ++i, ++ref, ++seq, ++quals) {
          sumPhred += (int)(*quals-33);
          
          if (*ref != *seq) {
            ++stat.numMismatches;
            // quals is encoded either as 0 or -1 by seqlib if the quality string is empty.
            if (*quals > 33) {
              if (qMax <= *quals-33) { // high base-quality mismatches
                mismatchQualsum += ((int)qMax);
                ++stat.numQ20Mismatches;
                stat.averageMismatchQuality += ((int)qMax);
              } else {
                mismatchQualsum += ((int)*quals)-33;
                stat.averageMismatchQuality += ((int)*quals)-33;
              }
              
            } else // quality string is empty... defaults to -1 in seqlib
              mismatchQualsum += qMax;
          }
        }
      } else { // consumes ref, not query ( deletion in query of clipping)
        if (! opt.ignoreIndels && t == 'D') {
          int qual = qMax; // defaults to Q20
          /* average of the two adjacent quality scores, iff they exist */
          if (quals >qualsInit) {
            if (quals[-1] > 33)
              qual = quals[-1]-33;
              
            if ( (c+1) < ciggy.end() ) {
              if (*quals > 33) {
                qual += *quals-33;
                qual /= 2;
              }
            } 
          } else if ( (c+1) < ciggy.end() ) {

            if (*quals > 33) 
              qual = *quals-33;
          }

          ++stat.numIndels;
          if (qual >= qMax) {
            ++stat.numQ20Indels;
            qual = qMax;
          }
          stat.meanIndelQuality += qual;
        }

        ref += len;
      }
    } else if ( c->ConsumesQuery() ) { // consumes query, not reference (insertion in query or clipping)

      
      // treat soft-clipping as mismatches
      if (t == 'S') {
        for ( i = 0; i < len; ++i, ++seq, ++quals) {
          ++stat.numMismatches;
          // quals is encoded either as 0 or -1 by seqlib if the quality string is empty.
          if (*quals > 33) {
            if (qMax <= *quals-33) { // high base-quality mismatches
              mismatchQualsum += ((int)qMax);
              ++stat.numQ20Mismatches;
              stat.averageMismatchQuality += ((int)qMax);
            } else {
              mismatchQualsum += ((int)*quals)-33;
              stat.averageMismatchQuality += ((int)*quals)-33;
            }
            
          } else // quality string is empty... defaults to -1 in seqlib
            mismatchQualsum += qMax;
        }

      } else if (! opt.ignoreIndels && t == 'I') {

        int thisSumPhred=0;
        for ( i = 0; i < len; ++i, ++seq, ++quals) {
          if (*quals < 33) { // empty quality string, treated as Q20
            sumPhred += qMax;
            thisSumPhred += qMax;
          } else {
            sumPhred += (int)(*quals-33);
            thisSumPhred += (int)(*quals-33);
          }
        }
        thisSumPhred /= i;
        
        ++stat.numIndels;
        if (thisSumPhred >= qMax) {
          ++stat.numQ20Indels;
          thisSumPhred = qMax;
        }
        
        stat.meanIndelQuality += thisSumPhred;
          
      } else {
        seq += len;
        quals += len;
      }
    }
  }

  if (stat.numMismatches) {
    stat.readLikelihood = pow(10, -mismatchQualsum/10.0); // taken from Heng Li's MAQ paper
    stat.averageMismatchQuality /= stat.numMismatches;
  }
  stat.meanPhred = (double)sumPhred / stat.numBases;
  if (stat.numIndels) {
    // stat.meanIndelQuality is the SUM at this point.
    stat.readLikelihoodWithIndels = pow(10, -(mismatchQualsum+stat.meanIndelQuality)/10.0); // taken from Heng Li's MAQ paper
    stat.meanIndelQuality = (double) stat.meanIndelQuality / stat.numIndels; // and now it's the mean
  } else {
    stat.readLikelihoodWithIndels = stat.readLikelihood;
  }
  
  return stat.numMismatches;
}

bool
realignIt(std::string seq, std::string qSeq, BWAWrapper &bwa, RefGenome &ref,  BamRecordVector &results, SummaryStat &stat, Options &opt) {
  
  results.clear();  
  bwa.AlignSequence(seq, "foo", results, false, 0.99, 100);
  vector<BamRecord>::iterator it;

  string refSubstr, chromName;
  int seqlen = seq.size();
  
  bool gotOne=false;

  SummaryStat thisStat = stat;
  
  // take a look at the best hits from BWA
  for (it = results.begin(); it != results.end(); ++it) {
    //int startPositionRead = it->AlignmentPosition();    
    int genomeStart = it->Position();
    int genomeStop =it->PositionEnd();
    chromName = bwa.ChrIDToName( it->ChrID()  );    
    refSubstr = ref.QueryRegion( chromName, genomeStart, genomeStop-1);

    // Ensure you get the sequences from the bam record; not seq and qSeq
    // (this accounts for the strand; bam alignments are always on the + strand)
    string thisqSeq = it->QualitySequence();
    const char *qseqRaw = thisqSeq.c_str();

    string thisSeq = it->Sequence();
    const char *seqRaw = thisSeq.c_str();
    

    // and compute summary statistics for this alignment
    int retval = getWeightedMismatches(*it,
                                       //&(seqRaw[ startPositionRead ]),
                                       seqRaw,
                                       refSubstr.c_str(),
                                       qseqRaw,
                                       //&(qseqRaw[ startPositionRead ]),
                                       //&(seqRaw[ thisSeq.size() ]),
                                       thisStat, DEFAULT_SEQ_CHAR-33, opt);

    if (retval < 0)
      continue;

    
    //    if ( thisTot > nTot || // aligned more bases
      //         (thisTot==nTot && thisProb > *readProb )){ // aligned same # of bases, but better posterior

    // and pick the alignment that's "best" 
    if (! gotOne ||
        thisStat.numBases > stat.numBases || // aligned more bases
        (thisStat.numBases == stat.numBases && thisStat.readLikelihood > stat.readLikelihood) // better alignment
        ) {
      
      stat = thisStat;
      // perfect hit. No need to iterate through the subsequent alignments
      if (opt.ignoreIndels) {
        if (stat.numBases == seqlen && stat.numMismatches==0)
          return true;
      } else {
        if (stat.numBases == seqlen && stat.readLikelihoodWithIndels==1.)
          return true;
      }
    }
    gotOne=true;
        
  }

  return gotOne;
}


bool
getSummaryStats(BamRecord &r, BWAWrapper &bwa, RefGenome &ref,  BamRecordVector &results, SummaryStat &stat, Options &opt, Cigar *ciggy) {


  // Cannot use the cigar field associated with the read (yay seqlib memory error)
  //  int startPositionRead = r.AlignmentPosition();
  //  int stopPositionRead = r.AlignmentEndPosition();


  int startPositionRead = getStartSubstr(ciggy);
  int genomeStart = r.Position();
  
  string qualities = r.Qualities();
  string seq =r.Sequence();
  int stopPositionRead = getStopSubstr(seq, ciggy);

  stat.genomeMidpoint = genomeStart + (int)(startPositionRead/2. + stopPositionRead/2.);
  stat.isReverse=r.ReverseFlag();
  stat.readName = r.Qname();

  
  bool hasNN = realignIt(
                         seq.substr( startPositionRead, stopPositionRead - startPositionRead),
                         qualities.substr( startPositionRead, stopPositionRead - startPositionRead),
                         bwa, ref, results, stat, opt);

  if (! hasNN)
    return false;


  return true;
}

void
softclipEverything(BamRecord &r) {
  string s = INT2STRING(r.Length() ) + "S";
  r.SetCigar(s);
}

/*
  This takes in two coordinate values-- start and stop, which are the positions of the read
  in the genome.
  and it takes a vector of GenomicRegions (regions that coorespond to amplicons)
  and an index i in amps in which we start looking
  AND
  it takes in a boolean; true means that the read is on the + strand, false - strand.
  
  If this returns a negative number then the read has no clear corresponding amplicon
  otherwise it returns the index in amps that corresponds to the amp that this read belongs to.

  If the read is on the + strand it must begin before the amplicon begins (<=)
  and span exactly 1 amplicon
  (stop < the stop coordinate of the next amp)

  If it's on the negative strand it must begin strictly after an amp ends (>=)
  and start after the start coordinate of the previous amplicon.

  To run this, first try the clipped coordinates (0/1 based half indexing. e.g., start coords are fine, but add 1 to the stop)
  if it fails (<0)
  then try the unclipped
 */
int
getAmpIndex(int32_t start, int32_t stop, const GenomicRegionVector &amps, unsigned j, bool forward) {

  int32_t a = -1;
  int32_t b = -1;
  int32_t initJ = j;
  int32_t ampSize = (int32_t) amps.size();
  // grab the left-most amplicon (by index) that this read can belong to if it spanned the whole amp
  while (j < amps.size()) {
    if (start > amps[j].pos1 ) {
      ++j;
    } else {
      a = j;
      break;
    }
  }

  j = initJ;

  // ditto for the right-most amp
  while (j < amps.size()) {
    if (stop >= amps[j].pos2 ) {
      b = j;
    } else
      break;
    
    ++j;
  }

  
  // lazy evaluation. the easy case.
  if (a==b && a != -1)
    return a;
  
  // rules for the read on the FORWARD strand
  if (forward) {
    
    if (a == -1) { // the read may be strictly within an amp. 
      return WITHIN_AMP; // WITHIN amp.
    } else if (a != b) { // the read doesn't span, but does begin strictly before one amp and dies out before the amp boundary
      
      // and ensure that it actually covers at least 1 base in the amp it's suppsoed to belong to.
      while (a < ampSize && amps[a].pos1 <= stop) {
        if (start <= amps[a].pos1 &&  stop > amps[a].pos1) {
          return a;
        }
        ++a;
      }

      
      return WITHIN_AMP;
    } 

  } else { // reverse strand. the read is moving from right to left in the coordinate system
    
    if (b == -1) { // strictly within
      return WITHIN_AMP;
    } else if(
              b != a
              ) {
      
      while (b < ampSize && stop >= amps[b].pos2) {
        if (start < amps[b].pos2)  // and again; at least 1 base is covered.
          return b;
        
        ++b;

      }


      return WITHIN_AMP;
    }


  }
    // the read spans TWO amplicons. no good!
    return MULTIPLE_AMPS;

}


/*
  This trims a read (bam record)
  to the boundaries of an amplicon
  If in doing so the read is of length 0 it is dropped  (function returns true)
  (can be written to disk optionally, but the entire read is set to soft-clipped)
  the trimming is done using soft-clipping, and thus &r is modified, as is &i
  which is used to update the index in the amps array
 */

bool
trimToAmpBoundaries(BamRecord &r, Cigar &cigar, Cigar &newCigar, const GenomicRegionVector &amps, unsigned &i) {


  
  // 0-based coordinates (Genome)
  int32_t startPosition = r.Position();
  int32_t startPositionBeforeClip = r.PositionWithSClips();

  
  // does the current amp STOP after this read STARTS
  // if so the current amp index needs adjusting.
  while (i < amps.size() && amps[i].pos2 <= startPositionBeforeClip) {
    ++i; // note the &i
  }
  
  // the read is AFTER all amps.
  // or the bam's unsorted...?
  if (i >= amps.size() ) {
    if (startPositionBeforeClip < amps[ amps.size() -1 ].pos2) {
      cerr << r << " vs " << amps[amps.size() -1 ] << endl;      
      cerr << "Your bam does not appear to be sorted. I need it to be..." << endl;
      die(NULL);
    }
    
    softclipEverything(r);
    return true;
  }

  // these use 1-based coordinates. together they are 0-based half open coordinates (google it). just like the bed file.
  int32_t stopPosition = startPosition + cigar.NumReferenceConsumed();

  // I tried adding soft-clipped positions back in to infer the amplicon
  // this is a mistake; this procedure can only increase the amount of soft-clipping (5' and 3')
  // so I cannot/should not try and undo this (soft-clipping is performed by the aligner or by up-stream analyses)
  //int32_t stopPositionBeforeClip = stopPosition;

  
  // first try the coordinates after soft-clipping
  int myAmp = getAmpIndex(startPosition, stopPosition, amps, i, ! r.ReverseFlag() );

  if (myAmp == MULTIPLE_AMPS) { // the coords after clipping are smaller. if they're still too big, kick it out!
    softclipEverything(r);
    return true;
  } else if (myAmp == WITHIN_AMP) {

    /*
    // if the read is too short, try re-adding the soft clipping back in and see if that changes the inference.
    if (r.ReverseFlag()) {
      myAmp = getAmpIndex(startPosition, stopPositionBeforeClip, amps, i, ! r.ReverseFlag() );
    } else {
      myAmp = getAmpIndex(startPositionBeforeClip, stopPosition, amps, i, ! r.ReverseFlag() );
    }
    */
    //if (myAmp < 0) {
    softclipEverything(r);
    return true;
      //}
  }

  

  int nextPos = r.Position(); // work in 1-based indexing (same indexing as amps[myAmp].pos2)
  
  int qBasesConsumed=0;
  int queryLength = r.Length(); // trimming can be inferred from the sequence length - the number of query bases consumed.
  
  bool softclipLeft = false;
  if (startPosition < amps[myAmp].pos1)
    softclipLeft = true;

  bool softclipRight = false;
  if (stopPosition > amps[myAmp].pos2)
    softclipRight = true;


  
  uint32_t softclipLeftSize=0;
  
  // adjust the soft-clipping 5'
  // this sets the start position (new read)
  if (startPosition < amps[myAmp].pos1) {
    r.SetPosition( amps[myAmp].pos1 );
  }
  
  vector<CigarField>::iterator itr = cigar.begin();


  char type;
  int len;
  while (itr != cigar.end()) {
    
    type = itr->Type();
    len = itr->Length();

    if (itr->ConsumesQuery()) {
      qBasesConsumed += len;
    }

    if (itr->ConsumesReference()) {
      nextPos += len;
    }

    if (type == 'H') { // added no matter what...
      newCigar.add(*itr); // and doesn't impact the positions as the sequence isn't present
      ++itr;
      continue;
      
    } else if (softclipLeft) {

      // note that in the ielse if
      // the current cigar OP is added in the if OUTSIDE of this else if code brace
      if (nextPos == amps[myAmp].pos1) { // the current tag doesn't need to be split (it's consumed entirely). But softclipping must occur
        newCigar.add( CigarField('S', qBasesConsumed));
        softclipLeft = false;
        softclipLeftSize = qBasesConsumed;
        ++itr;
        continue;
        
        // len is correct, so don't update it.
      } else if (nextPos > amps[myAmp].pos1) { // harder case. the tag *does* need to be split
        // consumesReference must be true here (otherwise nextPos would've been evaluated on the last round

        // the number of bases into the amplicon 
        int diff = nextPos -  amps[myAmp].pos1;
        
        softclipLeftSize = qBasesConsumed;
        if (itr->ConsumesQuery()) 
          softclipLeftSize -= diff;
        
        newCigar.add( CigarField('S', softclipLeftSize));
        len = diff; // effectively trims the op to the boundary of the amp (even deletions)
        // note: no continue!
        softclipLeft = false;
      } else { // whole OP is consumed in soft clipping...
        ++itr;
        continue;
      }
      
    }


    // whole op is within the amp...
    if (nextPos <= amps[myAmp].pos2) {
      newCigar.add( CigarField(type, len));
    } else if (softclipRight) {

      
      int diff = nextPos - amps[myAmp].pos2; // the number of bases that need to be decreased in the current op
      if (diff < 0) {
        cerr << "should never happen" << endl;
      }


      // qBasesConsumed includes len as well. This may need adjusting
      int rightclipSize= queryLength - qBasesConsumed;
      if (diff < len) // avoid a 0-len OP
        newCigar.add( CigarField(type, len-diff));
      
      if (itr->ConsumesQuery() ) {
        rightclipSize += diff;
      }

      newCigar.add(CigarField('S', rightclipSize ));
      
      // search for more hard clipping, and add that OP as well...
      // this OP does NOT impact the query or the ref
      while (++itr != cigar.end()) {
        if (itr->Type()=='H')
          newCigar.add(*itr);

      }
      
      break;
    }
    
    ++itr;
  }

  
  int numQConsumedNew = newCigar.NumQueryConsumed();
  int numQConsumedOld = cigar.NumQueryConsumed();
  if (numQConsumedOld != numQConsumedNew) {
    cout << "Should not happen!\n" <<
      numQConsumedNew << " " << numQConsumedOld << endl <<
      qBasesConsumed << "\t" << queryLength << endl <<
      newCigar << " " <<
      cigar << endl <<
      r << endl <<
      amps[myAmp].pos1 << " "  <<     amps[myAmp].pos2 << endl <<
      r.Position() << " " <<  nextPos << endl;
    exit(1);
  }
  
  return false;
}

inline void
add2Cache(BamRecord &r, Cigar *ciggy, unordered_map<string, pair<int, int>> &cache, int numtMatches) {
  string seq = r.Sequence();
  
  int startPositionRead = getStartSubstr(ciggy);
  int stopPositionRead = getStopSubstr(seq, ciggy);
  cache[seq.substr( startPositionRead, stopPositionRead - startPositionRead)] = std::make_pair(r.Position(), numtMatches);
}

inline bool
isInReadCache(BamRecord &r, Cigar *ciggy, unordered_map<string, pair<int, int>> &cache, int &genomePos, int &numtMatches) {

  //  int startPositionRead = r.AlignmentPosition();
  //int stopPositionRead = r.AlignmentEndPosition();
  string seq = r.Sequence();
  int startPositionRead = getStartSubstr(ciggy);
  int stopPositionRead = getStopSubstr(seq, ciggy);
  if (cache.find( seq.substr( startPositionRead, stopPositionRead - startPositionRead)) == cache.end() )
    return false;

  std::pair<int, int> p = cache[ seq.substr( startPositionRead, stopPositionRead - startPositionRead) ];
  genomePos = p.first;
  numtMatches = p.second;
  return true;
}


int
main(int argc, char** argv) {

  Options opt = parseOptions(&(argv[1]));
  
  RefGenome ref;
  RefGenome numtRef;
  
  SummaryStat stat;
  SummaryStat numtStat;
  
  vector<SummaryStat> stats;
  vector<SummaryStat> stats2numts;

  // associates a read sequence (after softclipping)
  // which has a perfect hit to the human mito (likelihood is 0)
  // to the genome position of this perfect hit as well as the mismatch distance to the nearest
  // numt
  unordered_map<string, pair<int, int>> perfectHmtHits;

  unordered_map<string, int> readPairs;
  
  int trainingSize=0;
  int numReadsSeen=0;

  std::default_random_engine generator;
  std::uniform_real_distribution<double> uniform(0.0,1.0);
  
  if (opt.train) {
    trainingSize = DEFAULT_TRAINING_SIZE;
    stats.resize(trainingSize+1);
    stats2numts.resize(trainingSize+1);
  }
  
  ref.LoadIndex(opt.humanFastaDbFilename);

  /* set the mapping parameters. affine gaps are a bad idea here... */
  BWAWrapper bwa2hum;
  if (!bwa2hum.LoadIndex(opt.humanFastaDbFilename) ) {
    cerr << endl << "Failed to load the index files for "  << endl << "\t" << opt.humanFastaDbFilename << endl << endl ;
    die(NULL);
  }
 
  bwa2hum.SetGapOpen(opt.gapOpen);
  bwa2hum.SetGapExtension(opt.gapExtend);
  bwa2hum.SetBandwidth(opt.bandwidth);
  
  bwa2hum.Set3primeClippingPenalty(opt.threePrimePenalty);
  bwa2hum.Set5primeClippingPenalty(opt.fivePrimePenalty);

  
  BWAWrapper bwa2nonh;

  if (!bwa2nonh.LoadIndex(opt.nonhumanFastaDbFilename)) {
    cerr << endl << "Failed to load the index files for "  << endl << "\t" << opt.humanFastaDbFilename << endl << endl ;
    die(NULL);
  }
  
  bwa2nonh.SetGapOpen(opt.gapOpen);
  bwa2nonh.SetGapExtension(opt.gapExtend);
  bwa2nonh.SetBandwidth(opt.bandwidth);
  bwa2nonh.Set3primeClippingPenalty(opt.threePrimePenalty);
  bwa2nonh.Set5primeClippingPenalty(opt.fivePrimePenalty);

  numtRef.LoadIndex(opt.nonhumanFastaDbFilename);
  
  SeqLib::BamReader br;
  if (!br.Open( opt.bamFilename)) {
    cerr << endl << endl << "Failed to open " << opt.bamFilename << " for reading " << endl << endl;
    die(NULL);
  }
  
  
  
  BamWriter bw(SeqLib::BAM);


  if (! opt.train) {
    string outFilename(opt.bamFilename, findExtensionSubstringPosition( opt.bamFilename) );
    outFilename += ".rtn.bam";
    if (!bw.Open(outFilename)) {
      cerr << endl << endl << "Failed to open " << outFilename << " for writing!!" << endl << endl;
      die(NULL);
    }
    
    bw.SetHeader( br.Header() );
    bw.WriteHeader();
  }
  
  // Get the chromosome name (as an int)
  // in the INPUT bam
  // of the chromosome we want
  int chromIndex=-1;
  string chromosomeName(opt.chrom);
  chromIndex = br.Header().Name2ID(chromosomeName);
  if (chromIndex < 0) {
    cerr << endl << "Failed to find chromosome/contig" << endl << "\t'" <<
      opt.chrom << "'" << endl<< "In your bam... I need this to go on..." << endl << endl;
    die(NULL);
  }

  // testing: for whole genome data, this lets us skip to the mito
  if (opt.indexJump) {
    GenomicRegion mito(chromIndex, 1, MITOLEN);
    if (!br.SetRegion(mito)) {
      cerr << endl << "Failed to seek to the mitochondrial genome...? Is your bam file not indexed?" << endl;
      die(NULL);
    }
  }
  
  BamRecord r;
  r.init();
  
  int chrID;

  BamRecordVector results;

  BamRecordVector out;
  
  GenomicRegionVector amps;
  unsigned ampIndex=0;
  if (opt.bedFilename != NULL) {
    if (!parseBed(opt.bedFilename, amps, br.Header())) {
      die("Failed to parse the bed file");
    }
  }

  
  
  while (br.GetNextRecord(r) ) {

    if (opt.removeDups && r.DuplicateFlag())
      continue;

    if ( r.DuplicateFlag() || ! r.MappedFlag() || r.SecondaryFlag() || r.SupplementaryFlag() ||
         r.MapQuality() < opt.minMappingQuality ||
         (r.PairedFlag() && ! r.ProperPair())
         ) {
      if (! opt.train && opt.writeOffTarget) {
        out.push_back(r);
      }
      continue;
    }
  
    chrID = r.ChrID();
    if (chrID != chromIndex) {
      if (! opt.train && opt.writeOffTarget) {
        out.push_back(r);
      }
      continue;
    }

    
    Cigar orig = r.GetCigar();
    Cigar newCig;
    Cigar *cig2use= &orig;
    
    /* 
       trim the read to the amplicon boundaries... 
    */
    if (opt.bedFilename != NULL) {
      
      bool readDropped = trimToAmpBoundaries( r, orig, newCig, amps, ampIndex);
     
      if (readDropped) {
        if (! opt.train && opt.minReadSize >= 0) {
          out.push_back(r);
        }
        continue;
      }
      
      cig2use = &newCig;
    } 
    
    if (cig2use->NumReferenceConsumed() < opt.minReadSize)
      continue;

    // todo: pass cig2use to all constructors!
    if (opt.train) {
      bool hasNN=true;
      bool hasNumtNN=true;
      // resevoir sample.
      if (numReadsSeen < trainingSize) { // haven't filled the resv.
        hasNN = getSummaryStats(r, bwa2hum, ref, results, stats[numReadsSeen], opt, cig2use);
        hasNumtNN = getSummaryStats(r, bwa2nonh, numtRef, results, stats2numts[numReadsSeen], opt, cig2use);
      } else { // randomly replace an index with prob resSize/nReadsSeen
        int index = uniform( generator ) * numReadsSeen;
        if (index < trainingSize) {
          // if hasNN is false then stats[index] is unchanged.
          hasNN = getSummaryStats(r, bwa2hum, ref, results, stats[index], opt, cig2use);
          hasNumtNN = getSummaryStats(r, bwa2nonh, numtRef, results, stats2numts[index], opt, cig2use);
        }
      }

      if (hasNN && hasNumtNN)
        ++numReadsSeen;
      
    } else {

      // check to see if we've seen this string sequence before
      // (after soft clipping)
      // if we have AND it was a perfect match to some human sequence
      // then the read passes, regardless of the filtering strategy used.
      // (i.e., the likelihood of the read is 1, which is not less than your threshold)
      int newgenomePos=-1;
      int numtDistance=-1; // only used in the if below...
      if (isInReadCache(r, cig2use, perfectHmtHits, newgenomePos, numtDistance)) {

        out.push_back(r);
        out.back().AddIntTag(DEFAULT_ZTAG_HUMAN,  0);  // definitional. that's the only case that we add it to the cache.
        out.back().AddIntTag(DEFAULT_ZTAG_NUMT,  numtDistance); 
        if (cig2use != &orig) 
          out.back().SetCigar(*cig2use);
        
        continue;
      }

      
      bool hasNN=getSummaryStats(r, bwa2hum, ref, results, stat, opt, cig2use);

      bool hasNumtNN = getSummaryStats(r, bwa2nonh, numtRef, results, numtStat, opt, cig2use);


      bool mqToZero=false;
      
      if (hasNN) {
        // try just a simple filter on the read likelihood
        if (opt.ignoreIndels) {

          double like = stat.readLikelihood;
          if (opt.scaleLikelihoodByReadlen) {
            like = log10(like)/stat.numBases;
          }

          if (opt.minLikelihood > like) {
            r.SetMapQuality(0);
            mqToZero=true;
          }

        } else {

          double like = stat.readLikelihoodWithIndels;
          if (opt.scaleLikelihoodByReadlen)
            like = log10(like)/stat.numBases;

          if (opt.minLikelihood > like) {
            r.SetMapQuality(0);
            mqToZero=true;
          }
          
        }

        out.push_back(r);
        
        out.back().AddIntTag(DEFAULT_ZTAG_HUMAN,  stat.numMismatches);
        if (hasNumtNN)
          out.back().AddIntTag(DEFAULT_ZTAG_NUMT,  numtStat.numMismatches);
        else
          out.back().AddIntTag(DEFAULT_ZTAG_NUMT,  NUMT_NONN_PLACEHOLDER);



        // a bug in seqlib occurs if you modify the bam *after* you modify the cigar. modifying the cigar must occur last.
        if (cig2use != &orig) {
          out.back().SetCigar(*cig2use);
        }

        // again try the singleton design pattern.
        // that is, record the perfect hits to HmtDB. If you see another perfect hit it should be in your
        // final bam, regardless of the q-score string.
        // It only works on a subset of cases (really, numMismatches==0 -> likelihood == 1,
        // which will be included no MATTER what.
        // otherwise, you have to consider the quality scores.
        if (stat.numMismatches==0) {
          if (hasNumtNN)
            add2Cache(r, cig2use, perfectHmtHits, numtStat.numMismatches);
          else
            add2Cache(r, cig2use, perfectHmtHits, NUMT_NONN_PLACEHOLDER);
        }
        
      } else { // alignment to HmtDB returned... nothing?! 
        r.SetMapQuality(0);
        out.push_back(r);
        mqToZero=true;
        out.back().AddIntTag(DEFAULT_ZTAG_HUMAN,  NUMT_NONN_PLACEHOLDER);
        if (hasNumtNN)
          out.back().AddIntTag(DEFAULT_ZTAG_NUMT,  numtStat.numMismatches);
        else
          out.back().AddIntTag(DEFAULT_ZTAG_NUMT,  NUMT_NONN_PLACEHOLDER);

      }

      if (opt.filterReadPairs && mqToZero) {
        // todo: fix read-mate-pair info
        // only true if trimming required (to amplicons)
        readPairs[ r.Qname() ] = 1;
      }
      
    }
  }
  

  
  if (opt.train) {
    cout << "Readname\tNumMismatches\tReadLikelihood\tGenomePosition\tAverageMismatchQuality\tNumQ20Mismatches\tMeanPhred\tIsReverse\tNumAlignedBases\t" <<
      "NumIndels\tNumQ20Indels\tMeanIndelQuality\tReadLikelihoodWithIndels\t" <<
      "Readname2\tNumMismatchesNUMT\tReadLikelihoodNUMT\tGenomePositionNUMT\tAverageMismatchQualityNUMT\tNumQ20MismatchesNUMT\tMeanPhredNUMT\tIsReverseNUMT\tNumAlignedBasesNUMT\t" <<
      "NumIndelsNUMT\tNumQ20IndelsNUMT\tMeanIndelQualityNUMT\tReadLikelihoodWithIndelsNUMT" << endl;
    int nMin = min(trainingSize, numReadsSeen);
    for (int i = 0; i < nMin; ++i)
      cout << stats[i] << "\t" << stats2numts[i] << endl;
    
  } else {
    // in-place sort
    BamRecordSort::ByReadPosition sorter;
    std::sort(out.begin(), out.end(), sorter);

    // write it to the bam
    vector<BamRecord>::iterator it;
    
    for (it = out.begin(); it != out.end(); ++it) {

      if (opt.filterReadPairs &&
          readPairs.find( it->Qname() ) != readPairs.end()) {

        it->SetMapQuality(0);
      }
      
      bw.WriteRecord(*it);

    }
    bw.Close();
    bw.BuildIndex(); // and be a good person. make an index too.
  }

  return 0;
}


