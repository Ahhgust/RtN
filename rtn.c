#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <random>
#include <algorithm>

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


// taken from: https://stackoverflow.com/questions/5590381/easiest-way-to-convert-int-to-string-in-c
#define INT2STRING( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

char DEFAULT_CHROM[] = "chrM"; 

#define DEFAULT_TRAINING_SIZE 100000

#define DEFAULT_GAPOPEN 5
#define DEFAULT_GAPEXTEND 5
#define DEFAULT_BANDWIDTH 1000
#define DEFAULT_THREEPRIMECLIP 1000
#define DEFAULT_FIVEPRIMECLIP 1000
// how many PCR errors are expected (per bp in a given read)
#define PCR_ERRORS_PER_BP 0.02
// 1% sequencing error is assumed if the qual string is empty
#define DEFAULT_SEQ_ERROR 0.01
// '5' in the ascii table corresponds to a phred score of 20 which is a 1% chance of error
#define DEFAULT_SEQ_CHAR 53

#define DEFAULT_MIN_MAP_QUALITY 4

// measured in bases of the genome
#define DEFAULT_MIN_READ_LEN 20

#define DEFAULT_MIN_LIKELIHOOD 1e-5

// The old ways are best.
// good old macro to convert a phred quality score in ASCII format
// into a probability
#define PHRED2PROB(c) ( (c>32) ? (pow(10, ( ((int)(c-33))/-10.0))) : DEFAULT_SEQ_ERROR )

// the mito genome is circularized...
// thus 16569*2 bases in length
#define MITO_LENGTH 33118

struct Options {
  char *humanFastaDbFilename;
  char *nonhumanFastaDbFilename;
  char *bamFilename;
  char *chrom;
  char *bedFilename;
  int gapOpen;
  int gapExtend;
  int bandwidth;
  int threePrimePenalty;
  int fivePrimePenalty;
  int minMappingQuality;
  double pcrErrorPerBp;
  bool ignoreIndels;
  bool verbose;
  bool train;
  bool writeOffTarget;
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
       << "\t-n nonhumanFasta" << endl
       << "\t-b bam" << endl
    << endl << "Optional (with default)" << endl
    << endl << "The naming (of the mito genome)..." << endl
    << "\t-c chromosome (" << DEFAULT_CHROM << ")" << endl
       << "\t-m minMappingQuality (" << DEFAULT_MIN_MAP_QUALITY << ")" << endl
    << "\t-t (train an SVM)" << endl
    << "\t-s bedFile (soft-clip reads to the amplicons specified in the bedFile)" << endl
    << "\t-w (writes off-target reads; defaults to FALSE)" << endl
    << "\t-l length (" << DEFAULT_MIN_READ_LEN << ") the minimum read length, as mapped to the genome)" << endl
    
    << endl << "The error estimate..." << endl
    << "\t-p pcrErrorPerReadPerBP (" << PCR_ERRORS_PER_BP << ")" << endl
    << "\t-i ignoreIndels ( false )" << endl
    
    << endl << "Read mapping parameters..." << endl
    << "\t-o gapOpen (" << DEFAULT_GAPOPEN << ")" << endl
    << "\t-e gapExtend (" << DEFAULT_GAPEXTEND << ")" << endl
    << "\t-w bandWidth (" << DEFAULT_BANDWIDTH << ")" << endl

    << "\t-3 3-primeclipping (" << DEFAULT_THREEPRIMECLIP << ")" << endl
       << "\t-5 5-primeclipping (" << DEFAULT_FIVEPRIMECLIP << ")" << endl << endl;
  
  exit(1);
}

/*
  Boilerplate command-line interaction
 */
Options
parseOptions(char **argv) {

  char f, *arg;
  Options opt = {NULL, NULL, NULL};
  opt.chrom=DEFAULT_CHROM;
  opt.bedFilename=NULL;
  
  opt.gapOpen = DEFAULT_GAPOPEN;
  opt.gapExtend = DEFAULT_GAPEXTEND;
  opt.bandwidth = DEFAULT_BANDWIDTH;
  opt.threePrimePenalty = DEFAULT_THREEPRIMECLIP;
  opt.fivePrimePenalty = DEFAULT_FIVEPRIMECLIP;
  opt.pcrErrorPerBp = PCR_ERRORS_PER_BP;
  opt.verbose = opt.ignoreIndels=false;
  opt.minMappingQuality = DEFAULT_MIN_MAP_QUALITY;
  opt.writeOffTarget=false;
  opt.minReadSize = DEFAULT_MIN_READ_LEN;
  opt.minLikelihood = DEFAULT_MIN_LIKELIHOOD;
  
  for (; *argv != NULL; ++argv) {
    arg = *argv;
    if (*arg != '-') {
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
    else if (f == 'p')
      opt.pcrErrorPerBp = atoi(*argv);
    else if (f == 'm')
      opt.minMappingQuality = atoi(*argv);
    else if (f == 'l')
      opt.minReadSize = atoi(*argv);
    else if (f == 'L')
      opt.minLikelihood = atof(*argv);
    else if (f == 'w') {
      opt.writeOffTarget = true;
      --argv; // only a flag; no argument.
    } else if (f == 'i') {
      opt.ignoreIndels = true;
      --argv; // only a flag; no argument.
    } else if (f == 't') {
      opt.train = true;
      --argv; // only a flag; no argument.
    } else if (f == 'v') {
      opt.verbose = true;
      --argv; // only a flag; no argument.
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
  
  return opt;
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
  The expected number of sequencing errors in the read is estimated as the sum of the error probabilities.
  (this is textbook)

  Special accomodations are made if we're using Ion sequencing, in which case homopolymers have a high rate of error
  BUT 
  only when try to say how long the homopolymer stretch is. 
  This function optionally compresses these stretches (effecitvely) 
  by taking the min probability of error WITHIN the stretch
  e.g.
  GAAATCC
  Would consider this sequence to be
  GATC
  and use the MIN of the 3 A probabilities
  and use the MIN of the trailing 2 C probabilties would be used

 */
int
getExpectedNumberOfSequencingErrors(const string &seq, const string &quals, int startPos, int stopPos, bool compressHomopolymers) {
  int i;
  
  double probError, prevProb=0.;
  double sumOfError=0.;
  int qualInt;
  
  for( i = startPos; i < stopPos; ++i) {

    // BAM format...
    // convert quality to probability of error...
    qualInt = ((int)(quals[i]-33));
    if (qualInt >= 0) {
      probError = pow(10,    -(qualInt /10.0 ))  ;
    } else {
      probError = DEFAULT_SEQ_ERROR; // TODO: Make this an argument!
    }
    
    //    cout << seq[i] << "\t" << quals[i] << "\t"  << ((int)(quals[i]-33)) << "\t" << probError << endl;
    
    if (compressHomopolymers) {
      // we're in a homopolymer
      // use the min error for the homopolymer as the estimator of the total error.
      if (i && seq[i-1] == seq[i]) {
        if (probError < prevProb)
          prevProb = probError;
      } else { // done w/ the homopolymer run...
        sumOfError += prevProb;
      }
    } else {
      sumOfError += probError;
    }
    prevProb = probError;
  }
  
  // trailing case; last base used is a homopolymer.
  if (compressHomopolymers && i && seq[i-1] == seq[i]) {
    sumOfError += prevProb;
  }
  return ceil(sumOfError);
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
    return 0;

  char t = c->Type();
  
  int mismatchQualsum=0;
  int sumPhred=0;
  // scroll past the front-most soft/hard clipping tag
  // don't need to adjust indexing b/c the input sequence are pre-adjusted for soft clipping
  if (t == 'S' || t == 'H') 
    ++c;

  
  if (c == ciggy.end()) // a cigar that's only clipping...? also weird, but okay I suppose.
    return 0;
  
  
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
              if (qMax <= *quals) { // high base-quality mismatches
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
    } else { // consumes query, not reference (insertion in query or clipping)
      /* TODO; qualities of indels (take 2)! */
      if (! opt.ignoreIndels && t == 'I') {

        int thisSumPhred=0;
        for ( i = 0; i < len; ++i, ++seq, ++quals) {
          sumPhred += (int)(*quals-33);
          thisSumPhred += (int)(*quals-33);
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
  bwa.AlignSequence(seq, "foo", results, false, 0.99, results.capacity());
  vector<BamRecord>::iterator it;

  const char *qseqRaw = qSeq.c_str();
  const char *seqRaw = seq.c_str();
  
  string refSubstr, chromName;

  bool gotOne=false;

  SummaryStat thisStat;

  // take a look at the best hits from BWA
  for (it = results.begin(); it != results.end(); ++it) {
    int startPositionRead = it->AlignmentPosition();    
    int genomeStart = it->Position();
    int genomeStop =it->PositionEnd();
    chromName = bwa.ChrIDToName( it->ChrID()  );    
    refSubstr = ref.QueryRegion( chromName, genomeStart, genomeStop-1);

    // and compute summary statistics for this alignment
    getWeightedMismatches(*it,
                          &(seqRaw[ startPositionRead ]),
                          refSubstr.c_str(),
                          &(qseqRaw[ startPositionRead ]), thisStat, DEFAULT_SEQ_CHAR-33, opt);
    
    //    if ( thisTot > nTot || // aligned more bases
      //         (thisTot==nTot && thisProb > *readProb )){ // aligned same # of bases, but better posterior

    // and pick the alignment that's "best" 
    if (! gotOne ||
        thisStat.numBases > stat.numBases || // aligned more bases
        (thisStat.numBases == stat.numBases && thisStat.readLikelihood > stat.readLikelihood) // better alignment
        ) {
      
      stat = thisStat;
    }
    gotOne=true;
        
  }

  return gotOne;
}


bool
getSummaryStats(BamRecord &r, BWAWrapper &bwa, RefGenome &ref,  BamRecordVector &results, SummaryStat &stat, Options &opt) {



  int startPositionRead = r.AlignmentPosition();
  int stopPositionRead = r.AlignmentEndPosition();
  int genomeStart = r.Position();
  
  string qualities = r.Qualities();
  string seq =r.Sequence();

  string ciggy = r.CigarString();
  if (ciggy == "*") {
    cerr << "Should never happen." << endl << r;
    exit(1);
  }
  
  bool hasNN = realignIt(
                         seq.substr( startPositionRead, stopPositionRead - startPositionRead),
                         qualities.substr( startPositionRead, stopPositionRead - startPositionRead),
                         bwa, ref, results, stat, opt);

  if (! hasNN)
    return false;

  stat.genomeMidpoint = genomeStart + (int)(startPositionRead/2. + stopPositionRead/2.);
  stat.isReverse=r.ReverseFlag();
  stat.readName = r.Qname();


  return true;
}





bool
realignItOld(std::string seq, BWAWrapper &bwa, RefGenome &ref,  BamRecordVector &results, int nMax, double targetDistance, Options &opt) {

  results.clear();  
  bwa.AlignSequence(seq, "foo", results, false, 0.8, nMax);

  vector<BamRecord>::iterator it;

  int matches, tot;
  
  string refSubstr, chromName;
  
  for (it = results.begin(); it != results.end(); ++it) {
    string seq  = it->Sequence();
    const char* rawSeq  =seq.c_str();
    int startPositionRead = it->AlignmentPosition();

    int genomeStart = it->Position();
    int genomeStop =it->PositionEnd();
    chromName = bwa.ChrIDToName( it->ChrID()  );
    
    refSubstr = ref.QueryRegion( chromName, genomeStart, genomeStop-1);

    if (getMismatches(*it, &(rawSeq[ startPositionRead ]), refSubstr.c_str(), &matches, &tot, opt.ignoreIndels) ) {
      //  cout << matches << "\t" << tot << "\t" << seqlen << endl;
      if ((double)matches/tot < targetDistance) {
        //        cout << "BOO\n";
        return true;
      }
    }
    
  }
  return false;
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
        if (start <= amps[a].pos1 &&  stop > amps[a].pos1)
          return a;
        
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
trimToAmpBoundaries(BamRecord &r, const GenomicRegionVector &amps, unsigned &i) {


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
  int32_t stopPosition = startPosition + r.GetCigar().NumReferenceConsumed(); 

  Cigar origCiggy = r.GetCigar();
  int32_t stopPositionBeforeClip = stopPosition;
  // re-add the soft-clipping to the stop-index (this is approximate)
  if (origCiggy.back().Type() == 'S')
    stopPosition += origCiggy.back().Length();
  

  
  // first try the coordinates after soft-clipping
  int myAmp = getAmpIndex(startPosition, stopPosition, amps, i, ! r.ReverseFlag() );

  if (myAmp == MULTIPLE_AMPS) { // the coords after clipping are smaller. if they're still too big, kick it out!
    softclipEverything(r);
    return true;
  } else if (myAmp == WITHIN_AMP) {

    // if the read is too short, try re-adding the soft clipping back in and see if that changes the inference.
    if (r.ReverseFlag()) {
      myAmp = getAmpIndex(startPosition, stopPositionBeforeClip, amps, i, ! r.ReverseFlag() );
    } else {
      myAmp = getAmpIndex(startPositionBeforeClip, stopPosition, amps, i, ! r.ReverseFlag() );
    }
    
    if (myAmp < 0) {
      softclipEverything(r);
      return true;
    }
  }

  // adjust the soft-clipping 5'
  if (startPosition < amps[myAmp].pos1) {
    Cigar orig = r.GetCigar();
    Cigar newCig;
    r.SetPosition( amps[myAmp].pos1 );


    // the number of reference bases that need to be consumed 5'...
    int deficitRef = amps[myAmp].pos1 - startPosition;
    int softclipSize = 0;
    
    for(vector<CigarField>::iterator itr = orig.begin(); itr != orig.end(); ++itr) {
      if (deficitRef <= 0 || itr->Type() == 'H') { // the soft-clipping has been accounted for. just copy the rest!
        newCig.add( *itr );
      } else { // these bases are to be (partially) soft-clipped.

        
        if (itr->ConsumesQuery()) {
          softclipSize += itr->Length();
        }

        if (itr->ConsumesReference()) {
          deficitRef -= itr->Length();

          if (deficitRef == 0) { // clean break; I don't need to fracture the current cigar field
            newCig.add( CigarField('S', softclipSize));
            softclipSize=0;
          } else if (deficitRef < 0) { // cut the current cigar-op in two. 

            if (itr->ConsumesQuery() )
              softclipSize += deficitRef; // consumes both Q and R,
            
            newCig.add( CigarField('S', softclipSize));             
            newCig.add( CigarField( itr->Type(), -deficitRef));
            softclipSize=0;
          }
        }
      }
      
    }
    // iterated through the read. never wrote the  soft-clipping
    if (softclipSize) {
      newCig.add( CigarField('S', softclipSize));             
    }

    // sanity test
    if (newCig.NumQueryConsumed() !=  r.Length()) {
      cerr << "Deficit " << deficitRef << " " << startPosition << "\t" << amps[myAmp].pos1 << endl;
      cerr << "new "  << newCig << endl;
      cerr << "old " << r.GetCigar( ) << endl;
      cerr << "Amp " << amps[myAmp] << endl;
      cerr << r.Position() << "\t" << r.Position() + newCig.NumReferenceConsumed()  << "\tVERSUS\t" << amps[myAmp].pos1 << "\t" << amps[myAmp].pos2 << endl;
      cerr << r << endl;
      exit(1);
    }

    r.SetCigar(newCig);
    
  }

  if (stopPosition > amps[myAmp].pos2+1) {
    Cigar orig = r.GetCigar();
    Cigar newCig;
    int posGenome = r.Position(); // work in 1-based indexing (same indexing as amps[myAmp].pos2)
    int nextPos = posGenome; 
    int qBasesConsumed=0;
    int queryLength = r.Length(); // trimming can be inferred from the sequence length - the number of query bases consumed.
    
    // write all of the cigar operations up until the end of the read
    for(vector<CigarField>::iterator itr = orig.begin(); itr != orig.end(); ++itr) {
      if (itr->ConsumesReference() ) {
        nextPos = posGenome + itr->Length();
      }
      
      if (itr->ConsumesQuery() ) {
        qBasesConsumed += itr->Length();
      }

      //      cout  << posGenome << "\t" << nextPos << "\t" << itr->Type() << "\t" << itr->Length() << endl;
      
      if (nextPos <= amps[myAmp].pos2) {
        newCig.add( *itr );
      } else { // we walked forward enough in the read where now we can start soft-clipping
        // by construction, ConsumesReference is true here. (otherwise the if statement above never flips from TRUE to FALSE)
        
        int diff = nextPos - amps[myAmp].pos2;
        
        // minimum; if the current cigar element spans the boundary, this needs to be adjusted
        int softclipSize = queryLength - qBasesConsumed;
        
        //        cerr << itr->Type() << "\t" << itr->Length() << "\t" << diff << endl;
        if (itr->Length() - diff > 0) {
          newCig.add( CigarField( itr->Type(), itr->Length() - diff) );

          if (itr->ConsumesQuery() )
            softclipSize += diff;
          
        } else if (itr->ConsumesQuery() )
          softclipSize += itr->Length();

        newCig.add( CigarField( 'S', softclipSize ));
        
        // and maintain hard clipping, but only if it was written in a way that's sane.
        if (orig.back().Type() == 'H') {
          newCig.add( orig.back() );
        }
        break;
      }
      
      posGenome = nextPos;
    }


    if (newCig.NumQueryConsumed() !=  r.Length()) {
      cerr << "Also new "  << newCig << endl;
      cerr << "Also old " << r.GetCigar( ) << endl; 
      cerr << r.Position() << "\t" << r.Position() + newCig.NumReferenceConsumed()  << "\tVERSUS\t" << amps[myAmp].pos1 << "\t" << amps[myAmp].pos2 << endl;
      cerr << r << endl;
      exit(1);
    }
    
    r.SetCigar(newCig);    
            
  }


  return false;
}


int
main(int argc, char** argv) {

  Options opt = parseOptions(&(argv[1]));
  
  RefGenome ref;
  RefGenome numtRef;
  
  SummaryStat stat;
  vector<SummaryStat> stats;
  vector<SummaryStat> stats2numts;
  
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
  br.Open( opt.bamFilename);

  BamWriter bw(SeqLib::BAM);


  if (! opt.train) {
    string outFilename(opt.bamFilename);
    outFilename += ".out.bam";
    bw.Open(outFilename);
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

  
  //  string refSequence = ref.QueryRegion(chromosomeName, 0, MITO_LENGTH);
  //  const char *refAsCharStar = refSequence.c_str();
  
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

    if ( r.DuplicateFlag() || ! r.MappedFlag() || r.SecondaryFlag() || r.SupplementaryFlag() ||
         r.MapQuality() < opt.minMappingQuality ||
         (r.PairedFlag() && ! r.ProperPair())
         ) {
      if (! opt.train && opt.writeOffTarget) {
        //        bw.WriteRecord(r);
        out.push_back(r);
      }
      continue;
    }
  
    chrID = r.ChrID();
    if (chrID != chromIndex) {
      if (! opt.train && opt.writeOffTarget) {
        //bw.WriteRecord(r);
        out.push_back(r);
      }
      continue;
    }

    /* 
       trim the read to the amplicon boundaries... 
    */
    if (opt.bedFilename != NULL) {

      
      bool readDropped = trimToAmpBoundaries( r, amps, ampIndex);
     
      if (readDropped) {
        if (! opt.train && opt.minReadSize >= 0) {
          // bw.WriteRecord(r);
          out.push_back(r);
        }
        continue;
      }       
    }
    
    if (r.GetCigar().NumReferenceConsumed() < opt.minReadSize)
      continue;

    
    if (opt.train) {
      bool hasNN=true;
      bool hasNumtNN=true;
      // resevoir sample.
      if (numReadsSeen < trainingSize) { // haven't filled the resv.
        hasNN = getSummaryStats(r, bwa2hum, ref, results, stats[numReadsSeen], opt);
        hasNumtNN = getSummaryStats(r, bwa2nonh, numtRef, results, stats2numts[numReadsSeen], opt);
      } else { // randomly replace an index with prob resSize/nReadsSeen
        int index = uniform( generator ) * numReadsSeen;
        if (index < trainingSize) {
          // if hasNN is false then stats[index] is unchanged.
          hasNN = getSummaryStats(r, bwa2hum, ref, results, stats[index], opt);
          hasNumtNN = getSummaryStats(r, bwa2nonh, numtRef, results, stats2numts[numReadsSeen], opt);
        }
      }

      if (hasNN && hasNumtNN)
        ++numReadsSeen;
      
    } else { 

      bool hasNN=getSummaryStats(r, bwa2hum, ref, results, stat, opt);

      if (hasNN) {
        // try just a simple filter on the read likelihood
        if (opt.ignoreIndels) {
          
          if (opt.minLikelihood <= stat.readLikelihood) {
            //bw.WriteRecord(r);
            out.push_back(r);     
          }
        } else if (opt.minLikelihood <= stat.readLikelihoodWithIndels) {
          //bw.WriteRecord(r);
          out.push_back(r);
        }
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
    
    for (it = out.begin(); it != out.end(); ++it)
      bw.WriteRecord(*it);
    
    bw.Close();
    bw.BuildIndex(); // and be a good person. make an index too.
  }
  
  return 0;
}


