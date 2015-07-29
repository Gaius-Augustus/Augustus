// Convert mRNA-to-genome alignments in BAM format into a hint file for AUGUSTUS in gff format.
// Sebastian Adler, 2011-09-27
// Some stuff patched by Tonatiuh Pena-Centeno,  2-December-2012

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h> // getopt_long
#include <list>     // data structure for all hints
#include <set>      // data structure for processed target names
#include <string.h> // strcmp
#include <algorithm>
#include <math.h>

// BAMTools
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAlignment.h>
#include <api/BamAux.h>


using namespace std;
using namespace BamTools;

// type that holds hint information
struct hint_t {
   int start;   // 1-based begin coordinate on the reference
   int end;     // 1-based last coordinate on the reference
   char strand; // one of +-.
   unsigned short int mult;    // multiplicity for intron hints
   // unsigned short is sufficient, as the number of alignments at one location is bounded by the (maximum) alignment coverage (unsigned short int)
   // constructor
   hint_t(){
      strand = '.';
      mult = 1;
   };
};

// type to label the hint lists with their IDs (see RefNameByID)
struct hintListLabel_t {
   list<hint_t>* ref;
   char* label; // e.g. exonpart, intron
   // constructor
   hintListLabel_t(list<hint_t>& r, const char* l){
      ref = &r;
      label = (char*) l;
   }
};


/*
 *  Global Variables
 */


// option parameters and default settings
char* InFileName;
char* OutFileName;
int Pri = 4;
int MaxGapLen = 14;
int MinIntLen = 32;
int MaxIntLen = 350000;
int MaxQGapLen = 5;
int EpCutoff = 10;
int MinEndBlockLen = 8;
const char* Source = "E";   //causes deprecation warning without 'const'
bool IntOnly = true;
bool Mult = true;
bool RemRed = false;
unsigned short int MaxCov = 0; // was 3000 before (from PSL and EST times)
bool SSOn = false;
bool TrunkSS = false;
double Score = 0;
char* CloneFileName;
char* TermFileName;
int MaxGeneLen = 400000;
bool Help = false;


const char* PrgSrc = "b2h"; // entry for gff field #2
BamAlignment* pal = new BamAlignment; // pointer to a single alignment, marks space for the new alignment at beginning of the while-loop
int TargetID = -2; // identifier of the currently processed reference sequence
int OldTargetID; // identifier of the formerly processed reference sequence
char* TargetName; // name of the processed reference sequence, used for printing and checking sortedness
vector<char*> RefNameByID; // array of reference sequence names, accessible with their BAM index
vector<int> RefLengthByID; // array of reference sequence lengths, accessible with their BAM index

// hint lists
list<hint_t> eplist;   // list of exonpart hints, ordered iff RemRed
list<hint_t> intronlist; // list of intron hints
list<hint_t> exonlist;   // list of exon hints
list<hint_t> dsslist;   // list of DSS hints
list<hint_t> asslist;   // list of ASS hints

vector<hintListLabel_t> hintList; // labels of the hint lists, defines printing order


/*
 *  Subroutines
 */


// gives value at 2^bit position of binary form of variable
int getbit (int bit, int variable)
{
  return (((1<<bit) & variable)>>bit);
}

// check and insert new hint into the list of exonpart hints
void addExonpartHint(hint_t* ph)
{
  if(!RemRed || eplist.empty()) // no need to sort
  {
    // simply append the hint to an (empty) list
    eplist.push_back(*ph);
  }
  else // check hint redundancy and sort into non-empty list
  {
    // search for nested hints
    bool redundant = false; // whether the hints information is already present
    list<hint_t>::iterator listIter = eplist.end();
    do
    {
      listIter--;
      if(listIter->start <= ph->start && listIter->end >= ph->end && listIter->strand == ph->strand)
      {
	// found including hint, drop the new hint
	redundant = true;
      }
      else if(listIter->start >= ph->start && listIter->end <= ph->end && listIter->strand == ph->strand)
      {
	// found included hint, delete it from the list
	eplist.erase(listIter);
      }
    }
    while(listIter != eplist.begin() && listIter->start > ph->start - 10000 && !redundant);
    // TODO: find better radius than fixed 10000 positions range (e.g. max exponpart length so far)

    // do sorting
    if(!redundant)
    {
      // search insertion site
      listIter = eplist.end();
      while(listIter != eplist.begin())
      {
	listIter--;
	if(listIter->start <= ph->start)
	{
	  // set insertion site AFTER the found position and quit the while loop
	  listIter++;
	  break;
	}
      }
      // insert
      eplist.insert(listIter, *ph);
    }
  }
} // end of addExonpartHint

// add new hint into the list of intron hints
void addIntronHint(hint_t* ph)
{
  // simply append hint to list
  intronlist.push_back(*ph);
}

// add new hint into the list of exon hints
void addExonHint(hint_t* ph)
{
  // simply append hint to list
  exonlist.push_back(*ph);
}

// add new hint into the list of DSS hints
void addDSS_Hint(hint_t* ph)
{
  // simply append hint to list
  dsslist.push_back(*ph);
}

// add new hint into the list of ASS hints
void addASS_Hint(hint_t* ph)
{
  // simply append hint to list
  asslist.push_back(*ph);
}

// output all exonpart hints, for test issues
void printExonpartList()
{
  printf("exonpart hints:\n");
  for(list<hint_t>::iterator listIter = eplist.begin(); listIter != eplist.end(); listIter++)
    {
      printf("%10i %10i\n", listIter->start, listIter->end);
    }
}

// output all intron hints, for test issues
void printIntronList()
{
  printf("intron hints:\n");
  for(list<hint_t>::iterator listIter = intronlist.begin(); listIter != intronlist.end(); listIter++)
    {
      printf("%10i %10i %i\n", listIter->start, listIter->end, listIter->mult);
    }
}

// output all exon hints, for test issues
void printExonList()
{
  printf("exon hints:\n");
  for(list<hint_t>::iterator listIter = exonlist.begin(); listIter != exonlist.end(); listIter++)
    {
      printf("%10i %10i\n", listIter->start, listIter->end);
    }
}


// ordering of hints used in list sorting
//  decide if hint1 is strictly before hint2
//  i.e. hint1 starts earlier or they start together and hint1 ends earlier
bool cmpIntervalHints(hint_t hint1, hint_t hint2)
{
  if(hint1.start < hint2.start || (hint1.start == hint2.start && hint1.end < hint2.end))
    return true;
  else
    return false;
}

// ordering of positional hints used in list sorting
//  decide if hint1 is strictly before hint2
//  i.e. hint1 starts earlier than hint2
bool cmpSiteHints(hint_t hint1, hint_t hint2)
{
  return (hint1.start < hint2.start);
}

// sort the whole list of intron hints
void sortIntronList()
{
  intronlist.sort(cmpIntervalHints);
}

// sort the whole list of exon hints
void sortExonList()
{
  exonlist.sort(cmpIntervalHints);
}

// sort the whole list of DSS hints
void sortDSS_List()
{
  dsslist.sort(cmpSiteHints);
}

// sort the whole list of DSS hints
void sortASS_List()
{
  asslist.sort(cmpSiteHints);
}


// print all hints (default) or the inalterable ones ("filter") to the outfile given with OUT
void printHints(FILE* OUT, const char* tag = "")
{
  vector<hintListLabel_t>::iterator hintListIter; // pointing to current hint list
  list<hint_t>::iterator listIter;                // pointing to current element in a list

  // preprocess certain lists
  sortIntronList();
  if(Mult && !intronlist.empty())
  {
    // compress similar intron hints
    listIter = intronlist.begin();
    list<hint_t>::iterator protoHint = listIter; // pointing to a prototype, i.e. the first in a row of similar hints
    listIter++; // start comparing to the second element

    // go through the intron hints
    while(listIter != intronlist.end())
    {
      // check similarity (equality of the intron boundaries)
      if(protoHint->start == listIter->start && protoHint->end == listIter->end)
      {
	// count the duplicate
	protoHint->mult++;
	// TODO: after filtered printing and sorting there may be similar hints with mult > 1 anywhere, then do   protoHint->mult += listIter->mult;

	// delete the duplicate and move to the returned following element
	listIter = intronlist.erase(listIter);
      }
      else
      {
	// get new prototype
	protoHint = listIter;
	// go to the next element
	listIter++;
      }
    }
  }

  sortExonList();
  sortDSS_List();
  sortASS_List();

  // TODO: drop exonpart hints contained in exon hints

  //list<hint_t>* hintList [] = { &eplist, &intronlist, &exonlist }; // references to the hint lists
  //const char* typeList [] = { "ep", "intron", "exon" }; // corresponding gff field #3 (type labels)

  //bool filter = false; // whether to print only finalized hints (e.g. when reaching the memory capacity)

  if( strcmp(tag, "filter") == 0 )
  {
    // filter before printing
    cout<<"filtering started\n";
    printExonpartList();
    printIntronList();
    printExonList();

    int pos = pal->Position; // hints beginning before here won't multiply or change and can thus be printed
    // TODO: also incorporate filtered hint printing

    for(hintListIter = hintList.begin(); hintListIter != hintList.end(); hintListIter++)
    {
      // traverse the (ordered) list until reaching the current position
      for(listIter = hintListIter->ref->begin(); listIter->start < pos && listIter != hintListIter->ref->end(); listIter++)
      {
	// print gff fields #1-#8
	fprintf(OUT, "%s\t%s\t%s\t%i\t%i\t%g\t%c\t.\t", TargetName, PrgSrc, hintListIter->label, listIter->start, listIter->end, Score, listIter->strand);
	// print gff field #9
	fprintf(OUT, "Key=Value\n");

	// remove the printed element
	hintListIter->ref->erase(listIter);
      }
    }
    cout<<"filtering done\n";
    printExonpartList();
    printIntronList();
    printExonList();
  }
  else
  {
    // print all remaining elements of each list

    for(hintListIter = hintList.begin(); hintListIter != hintList.end(); hintListIter++)
    {
      // output one gff line for every element
      for(listIter = hintListIter->ref->begin(); listIter != hintListIter->ref->end(); listIter++)
      {
	// gff fields #1-#8
	fprintf(OUT, "%s\t%s\t%s\t%i\t%i\t%g\t%c\t.\t", TargetName, PrgSrc, hintListIter->label, listIter->start, listIter->end, Score, listIter->strand);

	// gff field #9
	// TODO: print "grp" and "cdna" data for clonefiles
	// print multiplicity (of intron hints)
	if(listIter->mult > 1)
	{
	  fprintf(OUT, "mult=%i;", listIter->mult);
	}
	// AUGUSTUS priority and source values
	fprintf(OUT, "pri=%i;src=%s\n", Pri, Source);
      }

      // reset the list
      hintListIter->ref->clear();
    }
  } // end if (filtering)
}

// fetch and print fields of the alignment                                                                                                                                                                                                 
void showAlnLine(BamAlignment* pal)
{
  /*
    string Name = al.Name;
    int32_t Length = al.Length;
    string QueryBases = al.QueryBases;
    string AlignedBases = al.AlignedBases;
    string Qualities = al.Qualities;
    string TagData = al.TagData;
    int32_t RefID = al.RefID;
    int32_t Position = al.Position;
    uint16_t Bin = al.Bin;
    uint16_t MapQuality = al.MapQuality;
    uint32_t AlignmentFlag = al.AlignmentFlag;
    vector<CigarOp> CigarData = al.CigarData;
    int32_t MateRefID = al.MateRefID;
    int32_t MatePosition = al.MatePosition;
    int32_t InsertSize = al.InsertSize;
  // accessing tag data
  string tagvalue;
  bool tagfound = false;
  bool tagtypefound = -1;
  string strtag = "MD";
  const char* tag = strtag.c_str();
  char tagtype = 'x';
  tagtypefound = pal->GetTagType(strtag, tagtype); // bug fixed on GitHub Oct 16
  if(tagtypefound)    {      tagfound = pal->GetTag(tag, tagvalue);    }
  cout << tag << ":" << tagtypefound << " " << tagtype << "," << tagfound << " " << tagvalue << endl;
  */

  // Displaying the public attributes of the alignment structure

  cout << pal->Name << ";" << pal->Length << ";" << pal->QueryBases << ";" << pal->AlignedBases << ";" << pal->Qualities << ";" << pal->TagData << ";" << pal->RefID << ";" << pal->Position << ";" << pal->Bin << ";" << pal->MapQuality << ";" << pal->AlignmentFlag << ";";
  for( vector<CigarOp>::iterator CIGARiter = pal->CigarData.begin(); CIGARiter != pal->CigarData.end(); CIGARiter++ )
  {
    cout << (*CIGARiter).Length << CIGARiter->Type;
  }
  cout << pal->MateRefID << ";" << pal->MatePosition << ";" << pal->InsertSize << endl;
}

// return on/off value for a flag
const char* OnOff(bool flag)
{
  return flag ? "1=On" : "0=Off";
}


/*
 *  Begin of Execution
 */


int main(int argc, char* argv[])
{
  int opt;
  int option_index = -1;

  // define full name, argument necessity, output address, short form of the long options
  static option long_options[] = {
    {"in", 1, 0, 'i'},
    {"out", 1, 0, 'o'},
    {"priority", 1, 0, 'p'},
    {"maxgaplen", 1, 0, 'g'},
    {"minintronlen", 1, 0, 'm'},
    {"maxintronlen", 1, 0, 'M'},
    {"MinEndBlockLen", 1, 0, 'b'},
    {"maxQgaplen", 1, 0, 'q'},
    {"exonhints", 0, 0, 'x'},
    {"ep_cutoff", 1, 0, 'e'},
    {"source", 1, 0, 's'},
    {"intronsonly", 0, 0, 'I'},
    {"nomult",0, 0, 'n'},
    {"remove_redundant",0, 0, 'r'},
    {"maxcoverage", 1, 0, 'C'},
    {"ssOn",0, 0, 'S'},
    {"trunkSS", 0, 0, 'T'},
// option coloffset superfluous
    {"score", 1, 0, 'v'},
    {"clonefile", 1, 0, 'c'},
    {"terminusfile", 1, 0, 't'},
    {"maxgenelen", 1, 0, 'G'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };
  //printf("output of long options\n"); for(int i=0; i<5; i++){ printf("long_options[%i]: %s %i %p %c ,", i, long_options[i].name, long_options[i].has_arg, long_options[i].flag, long_options[i].val);}

  // receive the input options

  while((opt = getopt_long(argc, argv, "i:o:p:g:m:M:q:xe:s:InrC:STv:c:t:G:b:h", long_options, &option_index)) != -1)
  {
    switch(opt)
    {
      // TODO: error checking after string-to-int conversion ("-m=16" gives MinIntLen of 0, as "=" isn't expected)
    case 'i': InFileName    = optarg;       break;
    case 'o': OutFileName   = optarg;       break;
    case 'p': Pri           = atoi(optarg); break;
    case 'g': MaxGapLen     = atoi(optarg); break;
    case 'm': MinIntLen     = atoi(optarg); break;
    case 'M': MaxIntLen     = atoi(optarg); break;
    case 'b': MinEndBlockLen = atoi(optarg); break;
    case 'q': MaxQGapLen    = atoi(optarg); break;
    case 'e': EpCutoff      = atoi(optarg); break;
    case 's': Source        = optarg;       break;
    case 'x': IntOnly       = false;        break;
    case 'I': IntOnly       = true;         break; // if contradicting options -x and -I are specified, the later is preferred
    case 'n': Mult          = false;        break;
    case 'r': RemRed        = true;         break;
    case 'C': MaxCov        = atoi(optarg); break;
    case 'S': SSOn          = true;         break;
    case 'T': TrunkSS       = true;         break;
    case 'v': Score         = atof(optarg); break;
    case 'c': CloneFileName = optarg;       break;
    case 't': TermFileName  = optarg;       break;
    case 'G': MaxGeneLen    = atoi(optarg); break;
    case 'h': Help          = true;         break;
    }
  }

  // consistency check of supplied options
  if(MaxGapLen >= MinIntLen)
  {
    cerr << "Need to have maxgaplen < minintronlen\n";
    return -1;
  }

  // help text
  if(Help)
  {
    cout << "bam2hints -- Convert mRNA-to-genome alignments in BAM format into a hint file for AUGUSTUS in gff format.\n"
	 << "\n"
	 << "Usage:   bam2hints --in=example.bam --out=hints.gff\n"
	 << "  PREREQUISITE: input BAM file must be sorted by target (=genome) sequence names\n"
         << "                and within the sequences by begin coordinates\n"
      // TODO: add example of sorting on command line
         << "\n"
         << "  Options:\n"
         << "  --priority=n       -p   priority of hint group (set to " << Pri << ")\n"
         << "  --maxgaplen=n      -g   gaps at most this length are simply closed (set to " << MaxGapLen << ")\n"
         << "  --minintronlen=n   -m   alignments with gaps shorter than this and longer than maxgaplen are discarded (set to " << MinIntLen << ")\n"
         << "  --maxintronlen=n   -M   alignments with longer gaps are discarded (set to " << MaxIntLen << ")\n"
         << "  --MinEndBlockLen=n -b   minimum length of a 'dangling' exon (set to " << MinEndBlockLen << ")\n"
         << "  --maxQgaplen=n     -q   maximum length of gap in query (cDNA) sequence (set to " << MaxQGapLen << ")\n"
       	 << "  --exonhints        -x   compute exonpart, exon and splice site hints in addition to intron hints (set to " << OnOff(!IntOnly) << ")\n"
	 << "                          You should generate exonpart hints from RNA-Seq using wiggle (.wig) input to wig2hints.\n"
         << "  --ep_cutoff=n      -e   this many bp are cut off of each exonpart hint at end of alignment (set to " << EpCutoff << ")\n"
         << "  --source=s         -s   source identifier (set to '" << Source << "')\n"
         << "  --intronsonly      -I   only retreive intron hints (e.g. because the exon(part) hints are retreived by converting to a wig track, set to " << OnOff(IntOnly) << ")\n"
	 << "                          deprecated as this is the default now\n"
         << "  --nomult           -n   do not summarize multiple identical intron hints to a single one (set to " << OnOff(!Mult) << ")\n"
         << "  --remove_redundant -r   only keep the strongest hint for a region (set to " << OnOff(RemRed) << ")\n"
         << "  --maxcoverage=n    -C   maximal number of hints at a given position (0: filtering deactivated). A high value causes long running time of\n"
         << "                          AUGUSTUS in regions with thousands of cDNA alignments. (set to " << MaxCov << ")\n"
         << "  --ssOn             -S   include splice site (dss, ass) hints in output (set to " << OnOff(SSOn) << ")\n"
         << "  --trunkSS          -T   include splice sites hints from the ends of a truncated alignment (contig too short, set to " << OnOff(TrunkSS) << ")\n"
         << "  --score=f          -s   fill this number in in the score column (set to " << Score << ")\n"
      /*
         << "  --clonefile=s      -c   provide a file with clone names so close alignments from the same clone can be grouped.\n"
         << "                          AUGUSTUS will try to put those hints into a single transcripts even if different ends of\n"
         << "                          the clones were sequenced. File format (tab delimited):\n"
         << "                          cloneA\tread1\tread2\n"
         << "                          cloneA\tread3\n"
         << "                          cloneB\tread4\tread5\n"
         << "  --terminusfile=s   -t   provide a file with EST terminus information to infer tss/tts hints.\n"
         << "                          AUGUSTUS will use tss/tts hints to predict transcription start/termination sites\n"
         << "                          File format (tab delimited):\n"
         << "                          # ESTname    EstDir    Type FrontTerminus  EndTerminus\n"
         << "                          CACW5781.b1     5       A2      5TSS       Unknown\n"
         << "                          CACW6759.g1     3       F23     5TNS       3TNS\n"
         << "                          CACW14459.g2    3       D2      Unknown    3TNS\n"
         << "                          CACW21662.g1    3       C2      5TNS       Unknown\n"
         << "                          CACW25491.g1    3       F21     5TNS       3TNS-NP\n"
      */
         << "  --maxgenelen=n     -G   alignments of the same clone are considered to be of the same gene if not separeted by more than this (set to " << MaxGeneLen << ")\n"
         << "                          Alignments that span more than this are ignored, but better filter long introns through an alignment program.\n"
         << "  --help             -h   show this help text\n"
         << "\n";
    return 0;
  } 


  // open the input BAM file
  BamReader BAM;

  if( InFileName == NULL || !BAM.Open(InFileName) )
	{
	  cerr << "Could not open input BAM file: " << (InFileName ? InFileName : "No input file") << endl;
	  return -1;
	}

  // Estimating the right value of the arrays:  PSLb;  PSLq;  PSLt; BlockBegins; BlockEnds and FolOK;
  // This is done by sweeping through all the alignments and calculating the maxBlock size.
  cout << "Wait a moment, calculating maximum block size that needs to be allocated... ";
  int alignmentBlock, maxBlock=0;
  while (BAM.GetNextAlignment(*pal)) 
	{ 	
	  alignmentBlock=0;
	  // Retrieving maximum number of "blocks" in the BAM file
	  for( vector<CigarOp>::iterator tempIter = pal->CigarData.begin(); tempIter != pal->CigarData.end(); tempIter++ )
		{
		  if(tempIter->Type == 'M' || tempIter->Type == 'X' || tempIter->Type == '=')
			{
			  alignmentBlock++;
			}
		} 

	  if (alignmentBlock >= maxBlock) maxBlock = alignmentBlock;
	} // end while
  cout << ".. done" << endl;

  // closing and opening handle of BAM file
  BAM.Close();
  BAM.Open(InFileName);
 

  // check sortedness according to BAM
  SamHeader header = BAM.GetHeader();
  if(header.HasSortOrder() && header.SortOrder == "unsorted" && IntOnly && Mult)
  {
    cout << "\nBAM file MUST be sorted by target sequence names when 'intronsonly' and 'mult' options are active\n";
    return -1;
  }


  // open the output gff file
  FILE* GFF = fopen(OutFileName, "w");

  if( GFF == NULL )
  {
    cerr << "Could not open output file: " << (OutFileName ? OutFileName : "No output file") << endl;
    return -1;
  }


  /*
  // prints the SORTED(!) SAM header lines
  printf("%s", BAM.GetHeaderText().data());  //printf("%s", (char*) BAM.GetHeaderText().c_str());
  */

  // get a mapping from reference sequence IDs to their names and lengths, respectively

  RefVector RefSeq = BAM.GetReferenceData(); // implicit ID-to-name listing
  RefNameByID.resize(RefSeq.size()); // set the map vector's size to the number of reference sequence entries
  RefLengthByID.resize(RefSeq.size());

  // Obtaining the maximum reference sequence length
  int maxRefLen = *max_element(RefLengthByID.begin(),RefLengthByID.end());
  int maxCovBins = ceil(maxRefLen/10+1.5);
  // printf("\nmaxCovBins=%d\n", maxCovBins);

  // initialize the labelling of hint lists
  hintList.push_back(hintListLabel_t(eplist, "exonpart"));
  hintList.push_back(hintListLabel_t(intronlist, "intron"));
  hintList.push_back(hintListLabel_t(exonlist, "exon"));
  hintList.push_back(hintListLabel_t(dsslist, "DSS"));
  hintList.push_back(hintListLabel_t(asslist, "ASS"));

  
  // // print reference mapping
  // printf("RefNameByID:\n");
  // for(int rmit = 0; rmit < RefSeq.size(); rmit++)
  // {
  //   printf("%2i %s\n", rmit, RefNameByID[rmit]);
  // }


  // iterate through the alignments and get hints for the proper ones

  /*
  // save alignment in a map (for multi-line alignments)
  multimap<string,BamAlignment*> mal;
  mal[pal->Name] = pal; //   OR   mal.insert(make_pair(pal->Name,pal));   OR   mal.insert(pair<string,BamAlignment*>(pal->Name,pal));
  */

  long int Line = 0;    // line count
  int QOffset, TOffset; // 1-based start coordinate of the actual segment in query/target

  // PSL-like alignment data
  // TODO: ensure sufficient array length / throw overflow warning
  int block;    // index of next matching block, holds the element count of the "PSL?" arrays
  int PSLb[maxBlock];
  int PSLq[maxBlock];
  int PSLt[maxBlock]; // may need to be 'long int' if refseq longer than 400 Mbp

  // filtered block data
  int blockNew;        // index of next filtered block, holds the element count of the following arrays
  int BlockBegins[maxBlock]; // 1-based start coordinates of filtered alignment blocks
  int BlockEnds[maxBlock];   // 1-based end coordinates of filtered alignment blocks
  bool FolIntOK[maxBlock];   // whether the gap following a block is considered an intron

  set<char*> seenRefSet; // list of already encountered reference sequences to check sortedness
  bool badAlignment;     // alignment quality flag
  int BlockIter;         // index of current block
  unsigned short int * alnCoverage = new unsigned short int [maxCovBins]; // alignment coverage data of the current reference sequence
    // as alternative use STL container vector<unsigned short int>
  int CovIter;           // index of current bin of alignment coverage
  int GapLen;            // length of the gap preceding the current block on the target sequence
  hint_t* hint = new hint_t; // pointer to a new hint, passed to the "add...Hint" routines

  //  cout << "size of BamAlignment: " << sizeof(BamAlignment) << endl;
  

  while(BAM.GetNextAlignment(*pal))
  {
    // check and derive hints from a single alignment line


    // increase line count
    Line++;

    /*
    // print processing information every 1000 lines
    if(Line % 1000 == 1)
    {
      if(Line > 1){ cout << "\r"; }
      cout << "Processed line " << Line - 1 ;
      cout.flush(); // to really output the message
      // print a newline after ending while
    }
    */

    /*
    // for multi-line alignments
    // TASK: check completeness of actual alignment
        // TASK: BAM segment flag
        // TASK: hash previous parts
            // TASK: search for key "query+target"
            // TASK: append/prepend/insert actual alignment into previous segments using BAM's NEXT info
    // TASK: if complete or eof: (for now assuming complete alignments in every single line)
    */

    // deduction of alignment notation similar to PSL
    // i.e. block lengths "PSLb", 1-BASED(!) start positions in query "PSLq" and start positions in target "PSLt"

    badAlignment = false; // whether this alignment should be dropped
    block = 0; // reset PSL block count
    QOffset = 1; // refers to the first base in the alignment, not necessary the first base of the read itself!
    TOffset = pal->Position + 1; // transform 0-based alignment start to 1-based coordinate

    /*
    // replace every CIGAR string with a test dummy
    vector<CigarOp> CIGARtest;
    CIGARtest.push_back(CigarOp('H',1));
    CIGARtest.push_back(CigarOp('M',2));
    CIGARtest.push_back(CigarOp('I',4));
    CIGARtest.push_back(CigarOp('M',8));
    CIGARtest.push_back(CigarOp('D',16));
    CIGARtest.push_back(CigarOp('M',32));
    CIGARtest.push_back(CigarOp('X',64));
    CIGARtest.push_back(CigarOp('P',4096));
    CIGARtest.push_back(CigarOp('=',128));
    CIGARtest.push_back(CigarOp('N',256));
    CIGARtest.push_back(CigarOp('P',512));
    CIGARtest.push_back(CigarOp('X',1024));
    CIGARtest.push_back(CigarOp('H',2048));
    pal->CigarData = CIGARtest;
    */

    // parse the CIGAR string

    for( vector<CigarOp>::iterator CIGARiter = pal->CigarData.begin(); CIGARiter != pal->CigarData.end(); CIGARiter++ )
    {
      // decide by CIGAR type (one of "MIDNSHP=X")

      // TODO: incorporate 'B' type?

      if(CIGARiter->Type == 'M' || CIGARiter->Type == 'X' || CIGARiter->Type == '=')
      {
  	// if extension of the previous block is necessary...
  	if(block > 0 &&  PSLt[block-1] + PSLb[block-1] == TOffset && PSLq[block-1] + PSLb[block-1] == QOffset)
  	{
  	  // add this segment to the existing block, update the block length
          PSLb[block-1] += CIGARiter->Length;
  	}
  	else
  	{
  	  // create a new block
  	  PSLb[block] = CIGARiter->Length;
  	  PSLq[block] = QOffset;
  	  PSLt[block] = TOffset;
  	  block++;
  	}

  	// set positions in query and target to right after that block
  	QOffset += CIGARiter->Length;
  	TOffset += CIGARiter->Length;
      }
      else if(CIGARiter->Type == 'H' || CIGARiter->Type == 'S' || CIGARiter->Type == 'P')
      {
  	// ignore clipping (unaligned rad sequence) and padding (gap-to-gap-match used in multiple sequence alignments)
      }
      else if(CIGARiter->Type == 'D' || CIGARiter->Type == 'N')
      {
  	// update position in target
  	TOffset += CIGARiter->Length;
      }
      else if(CIGARiter->Type == 'I')
      {
  	// update position in query
  	QOffset += CIGARiter->Length;
      }
      else
      {
  	printf("Found unknown CIGAR type %c\n", CIGARiter->Type);
  	// drop alignment and suppress further filtering of blocks
  	badAlignment = true;
  	block = -1;
      } // end if (recognizing CIGAR type)
    } // end for (parsing of CIGAR string)

    /*
    // print infos for parsing of CIGAR
    for( vector<CigarOp>::iterator CIGARiter = pal->CigarData.begin(); CIGARiter != pal->CigarData.end(); CIGARiter++ )
    {
      printf("%i%c",CIGARiter->Length,CIGARiter->Type);
    }
    int it;
    printf("\n    b:");
    for(it=0; it<block; it++){ printf(" %i",PSLb[it]); }
    printf("\n    q:");
    for(it=0; it<block; it++){ printf(" %i",PSLq[it]); }
    printf("\n    t:");
    for(it=0; it<block; it++){ printf(" %i",PSLt[it]); }
    printf("\n");
    */

    // drop the alignment if its length on the target exceeds the theoretical maximum length of a gene

    if(badAlignment || PSLt[block-1] + PSLb[block-1] - PSLt[0] > MaxGeneLen)
    {
      /*
      printf("Alignment dropped: exceeding max length: ");
      for( vector<CigarOp>::iterator CIGARiter = pal->CigarData.begin(); CIGARiter != pal->CigarData.end(); CIGARiter++ )
      {
  	printf("%i%c",CIGARiter->Length,CIGARiter->Type);
      }
      printf("\n");
      */
      continue;
    }


    // get the reference sequence IDs
    OldTargetID = TargetID;
    TargetID = pal->RefID;

    // check alteration of the reference sequence, if so print hints
    if(TargetID != OldTargetID)
    {
      // check if there is a proper old target
      if(OldTargetID >= 0) // -1 <-> "*"(unknown reference), -2 <-> uninitialized
      {
  	// check previous occurence of that target
  	if(seenRefSet.find(TargetName) != seenRefSet.end())
  	{
  	  if(IntOnly && Mult)
  	  {
  	    // require sorting and abort
  	    cout << "\nBAM file MUST be sorted by target sequence names when 'intronsonly' and 'mult' options are active\n";
  	    return -1;
  	  }
  	}
  	seenRefSet.insert(TargetName);

  	// print remaining hints using the remaining TargetName
  	printHints(GFF, "");
      }

      // reset data structures for the new target
      if(TargetID >= 0)
      {
		TargetName = strdup(RefSeq.at(TargetID).RefName.c_str()); // update target name

  	// free the alignment coverage array
  	delete [] alnCoverage;

  	int CovBinCount = RefSeq.at(TargetID).RefLength/10 + 1; // needed number of entries
  	// cout << "CovBinCount=" << CovBinCount << endl;
  	// allocate alignment coverage array
  	alnCoverage = new(nothrow) unsigned short int [CovBinCount]; // disable exceptions for failures
  	// handle failed allocation
  	if(alnCoverage == NULL)
  	{
  	  cout << "Could not allocate memory for " << TargetName << "\n"
  	       << "Aborting!\n";
  	  return -1;
  	}

  	// initialize the coverage with zeros
  	for(CovIter = 0; CovIter < CovBinCount; CovIter++)
  	{
  	  alnCoverage[CovIter] = 0;
  	}
      }
    }


    // apply a coverage threshold
    // check each 10bp bin for too high abundance of alignments
    for(CovIter = PSLt[0]/10; CovIter <= (PSLt[block-1] + PSLb[block-1] - 1)/10 - 1; CovIter++)
    {
      if(MaxCov > 0 && alnCoverage[CovIter] >= MaxCov)
      {
  	// stop scanning and ...
  	badAlignment = true;
  	break;
      }
    }
  	// cout << "CovIter=" << CovIter << endl;

    if(badAlignment)
    {
      //cerr<<"reached coverage at "<<TargetName<<", position "<<CovIter*10<<"\n";
      // ... drop the alignment
      continue;
    }

    // update the coverage data with the accepted alignment
    for(CovIter = PSLt[0]/10; CovIter <= (PSLt[block-1] + PSLb[block-1] - 1)/10 - 1; CovIter++)
    {
       if (CovIter < maxCovBins)
	  alnCoverage[CovIter]++;
       // there is a bug here because above range check is not always satisfied
       // however, this maxCov filtering is not active by default since May25th, 2015
    }


    // filter the blocks as in "blat2hints.pl"

    blockNew = 0;

    for(BlockIter = 0; BlockIter < block; BlockIter++)
    {
      // blat2hints notation: mstart = PSLt[BlockIter] , mend = PSLt[BlockIter] + PSLb[BlockIter] - 1
      // get the length of the gap to the former block
      if(blockNew == 0)
      {
  	GapLen = MinIntLen;
      }
      else
      {
  	GapLen = PSLt[BlockIter] - BlockEnds[blockNew-1] - 1;
      }
      // printf("        GapLen=%i blockNew=%i MinIntLen=%i MaxIntLen=%i MaxGapLen=%i MaxQGapLen=%i\n", GapLen, blockNew, MinIntLen, MaxIntLen, MaxGapLen, MaxQGapLen);

      // decide by gap length
      if(MinIntLen <= GapLen && GapLen <= MaxIntLen)
      {
  	// gap represents an intron, add new block
  	BlockBegins[blockNew] = PSLt[BlockIter];
        BlockEnds[blockNew] = PSLt[BlockIter] + PSLb[BlockIter] - 1;
  	if(BlockIter < block - 1 && PSLq[BlockIter+1] - PSLq[BlockIter] - PSLb[BlockIter] <= MaxQGapLen)
  	{
  	  FolIntOK[blockNew] = true;
  	}
  	else
  	{
  	  FolIntOK[blockNew] = false;
  	}
  	blockNew++;
      }
      else if(GapLen <= MaxGapLen)
      {
  	// gap represents a deletion in the read, expand former block
  	BlockEnds[blockNew-1] = PSLt[BlockIter] + PSLb[BlockIter] - 1;
        if(BlockIter < block - 1 && PSLq[BlockIter+1] - PSLq[BlockIter] - PSLb[BlockIter] <= MaxQGapLen)
  	{ 
  	  FolIntOK[blockNew-1] = true;
  	}
        else
  	{ 
  	  FolIntOK[blockNew-1] = false;
  	}
      }
      else
      {
  	// gap is not feasible, drop alignment
  	badAlignment = true;
      }
    } // end for

    /*
    // print info on filtered blocks
    int it2;
    printf("    BlockBegins:");
    for(it2=0; it2<blockNew; it2++){ printf(" %i",BlockBegins[it2]); }
    printf("\n    BlockEnds:");
    for(it2=0; it2<blockNew; it2++){ printf(" %i",BlockEnds[it2]); }
    printf("\n    FolIntOK:");
    for(it2=0; it2<blockNew; it2++){ printf(" %i",FolIntOK[it2]); }
    printf("\n");
    */

    if(!badAlignment)
    {
      // derive hints
      //printf("    Deriving hints\n");
      for(BlockIter=0; BlockIter<blockNew; BlockIter++)
      {
  	if(BlockIter == 0)
  	{
  	  // first block of an alignment
  	  if(blockNew == 1 && !IntOnly)
  	  {
  	    // one-block-alignment, only exonpart hint derivable
  	    if(BlockEnds[BlockIter] - BlockBegins[BlockIter] >= 2*EpCutoff)
  	    {
  	      //printf("        Exonpart %i %i\n", BlockBegins[BlockIter]+EpCutoff, BlockEnds[BlockIter]-EpCutoff);
  	      hint->start = BlockBegins[BlockIter]+EpCutoff;
  	      hint->end = BlockEnds[BlockIter]-EpCutoff;
  	      addExonpartHint(hint);
  	    }
  	  }
  	  else if(BlockEnds[BlockIter] - BlockBegins[BlockIter] + 1 >= MinEndBlockLen)
  	  {
  	    // first block of multi-block-alignment and having a minimum length
  	    if(BlockEnds[BlockIter] - BlockBegins[BlockIter] >= EpCutoff && !IntOnly)
  	    {
  	      // exonpart hint
              //printf("        Exonpart %i %i\n", BlockBegins[BlockIter]+EpCutoff, BlockEnds[BlockIter]);
              hint->start = BlockBegins[BlockIter]+EpCutoff;
  	      hint->end = BlockEnds[BlockIter];
              addExonpartHint(hint);
  	    }
  	    if(SSOn && !IntOnly)
  	    {
  	      // SS hints for the following intron (all 4 possibilities)
  	      //printf("        DSS/ASS %i %i\n", BlockEnds[BlockIter]+1, BlockEnds[BlockIter]+1);
              //printf("        DSS/ASS %i %i\n", BlockBegins[BlockIter+1]-1, BlockBegins[BlockIter+1]-1);
              hint->start = BlockEnds[BlockIter] + 1;
              hint->end = BlockEnds[BlockIter] + 1;
  	      addDSS_Hint(hint);
              addASS_Hint(hint);
              hint->start = BlockBegins[BlockIter+1] - 1;
              hint->end = BlockBegins[BlockIter+1] - 1;
              addDSS_Hint(hint);
              addASS_Hint(hint);
  	    }
  	    if(FolIntOK[BlockIter] && (BlockIter < blockNew - 2 || BlockEnds[BlockIter+1] - BlockBegins[BlockIter+1] + 1 >= MinEndBlockLen))
  	    {
              //printf("        Intron %i %i\n", BlockEnds[BlockIter]+1, BlockBegins[BlockIter+1]-1);
  	      hint->start = BlockEnds[BlockIter] + 1;
              hint->end = BlockBegins[BlockIter+1] - 1;
              addIntronHint(hint);
  	    }
  	  }
  	}
  	else if(BlockIter == blockNew - 1 && !IntOnly)
  	{
  	  // last block of multi-block-alignment
  	  if(BlockEnds[BlockIter] - BlockBegins[BlockIter] + 1 >= MinEndBlockLen && BlockEnds[BlockIter] - BlockBegins[BlockIter] >= EpCutoff)
  	  {
  	    //printf("        Exonpart %i %i\n", BlockBegins[BlockIter], BlockEnds[BlockIter]-EpCutoff);
  	    hint->start = BlockBegins[BlockIter];
  	    hint->end = BlockEnds[BlockIter]-EpCutoff;
  	    addExonpartHint(hint);
  	  }
  	}
  	else
  	{
  	  // inner block of an alignment
  	  if(!IntOnly)
  	  {
  	    //printf("        Exon %i %i\n", BlockBegins[BlockIter], BlockEnds[BlockIter]);
  	    hint->start = BlockBegins[BlockIter];
            hint->end = BlockEnds[BlockIter];
            addExonHint(hint);
  	  }
  	  if(FolIntOK[BlockIter] && (BlockIter < blockNew - 2 || BlockEnds[BlockIter+1] - BlockBegins[BlockIter+1] + 1 >= MinEndBlockLen))
  	  {
  	    //printf("        Intron %i %i\n", BlockEnds[BlockIter]+1, BlockBegins[BlockIter+1]-1);
  	    hint->start = BlockEnds[BlockIter] + 1;
  	    hint->end = BlockBegins[BlockIter+1] - 1;
  	    addIntronHint(hint);
  	    if(SSOn && !IntOnly)
  	    {
  	      //printf("        DSS/ASS %i %i\n", BlockEnds[BlockIter]+1, BlockEnds[BlockIter]+1);
              //printf("        DSS/ASS %i %i\n", BlockBegins[BlockIter+1]-1, BlockBegins[BlockIter+1]-1);
              hint->start = BlockEnds[BlockIter] + 1;
              hint->end = BlockEnds[BlockIter] + 1;
              addDSS_Hint(hint);
              addASS_Hint(hint);
              hint->start = BlockBegins[BlockIter+1] - 1;
              hint->end = BlockBegins[BlockIter+1] - 1;
              addDSS_Hint(hint);
              addASS_Hint(hint);
  	    }
  	  }
  	} // end if (figuring out the blocks position)
      } // end for (looping through blocks)
    } // end if (deriving hints)


    // TASK: if hint buffers are full:
        // TASK: scan for finished hints
        // TASK: print to output file
    // TASK: insert new hints
        // TASK: free memory, remember one site for "pal"
    // TASK: else store alignment pointer into map
    // TASK: set new or reuse "pal"

  } // end while (parsing through bam aligment lines)


  // newline after display of the line count
  cout << "\n";

  // print the last hints
  printHints(GFF, "");

  // close input and output files
  BAM.Close();
  fclose(GFF);

  return 0;

} // end of main
