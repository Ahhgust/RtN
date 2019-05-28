#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>

#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/FastqReader.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BamReader.h"

using namespace SeqLib;
using namespace std;

char DEFAULT_CHROM[] = "chrM"; 

#define DEFAULT_GAPOPEN 5
#define DEFAULT_GAPEXTEND 5
#define DEFAULT_BANDWIDTH 1000
#define DEFAULT_THREEPRIMECLIP 100
#define DEFAULT_FIVEPRIMECLIP 100
// how many PCR errors are expected (per bp in a given read)
#define PCR_ERRORS_PER_BP 0.02
// 1% sequencing error is assumed if the qual string is empty
#define DEFAULT_SEQ_ERROR 0.01
// '5' in the ascii table corresponds to a phred score of 20 which is a 1% chance of error
#define DEAULT_SEQ_CHAR 53

// The old ways are best.
// good old macro to convert a phred quality score in ASCII format
// into a probability
#define PHRED2PROB(c) ( (c>=32) ? (pow(10, ( ((int)(c-33))/-10.0))) : DEFAULT_SEQ_ERROR )

// the mito genome is circularized...
// thus 16569*2 bases in length
#define MITO_LENGTH 33118

struct Options {
  char *humanFastaDbFilename;
  char *nonhumanFastaDbFilename;
  char *bamFilename;
  char *chrom;
  int gapOpen;
  int gapExtend;
  int bandwidth;
  int threePrimePenalty;
  int fivePrimePenalty;
  double pcrErrorPerBp;
  bool ignoreIndels;
  bool verbose;
};


/*
  Perl-inspired print stderr + exit + usage statement
 */
void
die(const char * message) {
  if (message != NULL) {
    cerr << message << endl;

  }

  cerr << "How to run this program!" << endl <<
    "Required " << endl
       << "\t-h humanFasta" << endl
       << "\t-n nonhumanFasta" << endl
       << "\t-b bam" << endl
    << endl << "Optional (with default)" << endl
    << endl << "The naming (of the mito genome)..." << endl
    << "\t-c chromosome (" << DEFAULT_CHROM << ")" << endl
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
  opt.gapOpen = DEFAULT_GAPOPEN;
  opt.gapExtend = DEFAULT_GAPEXTEND;
  opt.bandwidth = DEFAULT_BANDWIDTH;
  opt.threePrimePenalty = DEFAULT_THREEPRIMECLIP;
  opt.fivePrimePenalty = DEFAULT_FIVEPRIMECLIP;
  opt.pcrErrorPerBp = PCR_ERRORS_PER_BP;
  opt.verbose = opt.ignoreIndels=false;
  
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
    else if (f == 'i') {
      opt.ignoreIndels = true;
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
  
  double probError, prevProb;
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
    sumOfError = prevProb;
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
getWeightedMismatches(BamRecord &r, const char* seq, const char *ref, const char *quals, char *mismatches) {


  Cigar ciggy = r.GetCigar();
  
  int i, len, nMismatch=0;
  Cigar::const_iterator c = ciggy.begin();
  if (c == ciggy.end()) // empty cigar?! I don't think this is legitimate, but you never know!
    return 0;
  char t = c->Type();
  
  
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
        for ( i = 0; i < len; ++i, ++ref, ++seq, ++quals) {
          if (*ref != *seq) {
            mismatches[ nMismatch++] += *quals;
          }
        }
      } else { // consumes ref, not query ( deletion in query of clipping)
        /* TODO; qualities of indels?!
        if (! ignoreIndels && t == 'D') {

          
        }
        */
        ref += len;
      }
    } else { // consumes query, not reference (insertion in query or clipping)
      /* TODO; qualities of indels (take 2)!
      if (! ignoreIndels && t == 'I') {
        ++(*mismatches);
        ++(*tot);
      }
      */
      seq += len;
    }
  }


  return nMismatch;
}



bool
realignIt(std::string seq, BWAWrapper &bwa, RefGenome &ref,  BamRecordVector &results, int nMax, double targetDistance, Options &opt) {

  results.clear();  
  bwa.AlignSequence(seq, "foo", results, false, 0.8, nMax);

  vector<BamRecord>::iterator it;

  int matches, tot;
  int seqlen = seq.size();
  
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


int
main(int argc, char** argv) {

  Options opt = parseOptions(&(argv[1]));
  
  RefGenome ref;

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

  int numHumanSeqs =  bwa2hum.NumSequences();

  
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

  
  SeqLib::BamReader br;
  br.Open( opt.bamFilename);

  BamWriter bw(SeqLib::BAM);
  
  string outFilename(opt.bamFilename);
  outFilename += ".out.bam";
  bw.Open(outFilename);
  bw.SetHeader( br.Header() );
  bw.WriteHeader();
  
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

  string refSequence = ref.QueryRegion(chromosomeName, 0, MITO_LENGTH);
  const char *refAsCharStar = refSequence.c_str();
  
  BamRecord r;

  int startPositionRead, stopPositionRead, chrID;
  int genomeStart, mismatches, tot;
  string alignmentHit;

  BamRecordVector results;

  while (br.GetNextRecord(r) ) {

    if ( r.DuplicateFlag() || ! r.MappedFlag() ) {
      bw.WriteRecord(r);
      continue;
    }
    chrID = r.ChrID();
    if (chrID != chromIndex) {
      bw.WriteRecord(r);
      continue;
    }
    startPositionRead = r.AlignmentPosition();
    stopPositionRead = r.AlignmentEndPosition();

    genomeStart = r.Position();
    
    string qualities = r.Qualities();
    string seq =r.Sequence();
    const char *seqRaw = seq.c_str();

    // # of errors as estimated from the quality scores (from sequencing)
    int numErrors = getExpectedNumberOfSequencingErrors(seq, qualities, startPositionRead, stopPositionRead, true);
    //    cout << numErrors << endl;
    if ( getMismatches(r, &seqRaw[startPositionRead], &refAsCharStar[genomeStart], &mismatches, &tot, opt.ignoreIndels) ) {

      // and how many PCR AND sequencing errors might we expect?
      numErrors += ceil(opt.pcrErrorPerBp * tot);
      //      cout << numErrors << endl;
      
      if (numErrors >= mismatches) {
        bw.WriteRecord(r);        
      } else{
        //  cout << "here\n";
        bool hasNN = realignIt(
                           seq.substr( startPositionRead, stopPositionRead - startPositionRead),
                           bwa2hum, ref, results, numHumanSeqs, (double)numErrors/tot, opt);

        if (hasNN) {
          bw.WriteRecord(r);
        } else {
          //          cout << "flagged\n";
          
          r.SetMapQuality(1);
          bw.WriteRecord(r);
          
        }
      }
    }
    
    
  }
  bw.Close();
  bw.BuildIndex();
  
  /*

  
  string headerString = //"@HD\tVN:1.6\tSO:coordinate\n" +
    bwa.HeaderFromIndex().AsString();

  
  vector< string > seqNames;
  UnalignedSequence s;
  getReads(argv[2], seqNames); // read the fastq file, get all of the Names

  string biggerHeaderString = "";
  for (auto &name : seqNames) { // construct a read group header
    biggerHeaderString += getReadgroup(name);

  }

  // before we do the alignment, add the read groups in
  BamWriter writer;

  if (argc > 4) {
    writer.Open(argv[4]);
  } else {
    writer.Open("out1.bam");
  }

  writer.SetHeader( headerString  + biggerHeaderString);
  writer.WriteHeader();

  
  FastqReader fq(argv[2]);

  int numReps = DEFAULT_DUP;
  if (argc >3) {
    numReps = atoi(argv[3]);
    if (numReps < 1) {
      cerr << "Must choose a number of reps > 0!\n";
      return 1;
    }
  }
  */
  /*  
  while (fq.GetNextSequence(s) ) {
    BamRecordVector results;
    bwa.AlignSequence(s.Seq, s.Name, results, false, 0.9, 1);
    
    if (results.size() > 0) {


      BamRecord rec = results[0];
      rec.SmartAddTag("RG", s.Name);
      rec.SetQualities(s.Qual, 33);
      
      for (int i =1; i <= numReps; ++i) {
        
        rec.SetQname( s.Name + "_" + to_string(i) );
        string s = rec.QualitySequence();

        writer.WriteRecord(rec);

      }
    }
  }
  */
  
  /*
  BamRecordVector v;
  while (fq.GetNextSequence(s) ) {
    BamRecordVector results;
    bwa.AlignSequence(s.Seq, s.Name, results, false, 0.9, 1);
    
    if (results.size() > 0) {
      BamRecord rec = results[0];
      rec.AddZTag("RG", s.Name); // was SmartAddTag. 
      rec.SetQualities(s.Qual, 33);
      v.push_back(rec);
    }
  }
  cout << v[0].Qname() << endl << v[1].Qname() << endl;

  BamRecordSort::ByReadPosition sorter;
  std::sort(v.begin(), v.end(), sorter);

  for (auto rec = v.begin(); rec < v.end(); rec++) {
    //    std::string name = rec->Qname().c_str(); // deep copy of the name of the read/individual
    
    for (int i =1; i <= numReps; ++i) { // now write it a bunch of times

      //rec->SetQname( name + to_string(i) );
      writer.WriteRecord(*rec);
      
    }
  }
  
  writer.Close();
  writer.BuildIndex();
  */
  
  return 0;
}


