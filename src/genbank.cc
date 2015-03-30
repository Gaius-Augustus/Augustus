/**********************************************************************
 * file:    genbank.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  sequence and gene struture input from file
 * authors: Mario Stanke, mario@gobics.de, Stafilarakis
 *
 * date     |   author      |  changes
 * ---------|---------------|------------------------------------------
 * 03.01.07 | Mario Stanke  | fixed bug with clength (was 3 too long for every exon except the last)
 * 05.04.07 | Mario Stanke  | Do not read out incomplete mRNA if it does not have a 
 *          |               | matching but an overlapping CDS. Usefull for training with 3'incompete genes.
 **********************************************************************/

#include "genbank.hh"
#include "fasta.hh"

// project includes
#include "properties.hh"
#include "gene.hh"
#include "projectio.hh"  // for goto_line_after

// standard C/C++ includes
#include <iostream>

#ifdef ZIPINPUT
using boost::iostreams::filtering_istream;
using boost::iostreams::gzip_decompressor;
#endif

GBProcessor::GBProcessor(string filename) :
    gbs(filename)
{
    try{
	gbVerbosity = Properties::getIntProperty("/genbank/verbosity");
    } catch (...){
	gbVerbosity = 3;
    }
}

GBPositions* GBProcessor::nextPosition( ) throw( GBError ){
    return gbs.nextData( );
}

AnnoSequence* GBProcessor::getAnnoSequence( GBPositions* pos ){
    if( !pos )
        return (AnnoSequence*) 0;
    bool stopCodonExcludedFromCDS = false;
    try {
	stopCodonExcludedFromCDS = Properties::getBoolProperty("stopCodonExcludedFromCDS");
    } catch (ProjectError e) {
	stopCodonExcludedFromCDS = false;
    }
    bool withUTR;
    try {
      string utr = Properties::getProperty("UTR");
      withUTR = (utr == "both" || utr == "1" || utr == "on" || utr == "5");
    } catch (ProjectError e) {
      withUTR = false;
    }

    AnnoSequence *annoseq = new AnnoSequence();
    if (!annoseq)
	return NULL;
    annoseq->anno = new Annotation;
    annoseq->sequence = getSequence( *pos );
    annoseq->length = pos->seqlength;
    annoseq->bc.addSequence(annoseq->sequence, strlen(annoseq->sequence) );
    char* idtmp = strstr(pos->buffer, "LOCUS");
    if (idtmp == NULL) {
	annoseq->seqname = newstrcpy("unknown");
    } else {
	annoseq->seqname = new char[100];
	idtmp += strlen("LOCUS");
	while( isspace(*idtmp) )
	    idtmp++;
	int i;
	for( i = 0; i < 99 && !isspace(*idtmp); i++, idtmp++ ){
	    annoseq->seqname[i] = *idtmp;
	}
	annoseq->seqname[i] = '\0';
    }
    int numGenes = pos->CDS.size();
    if (numGenes == 0 && gbVerbosity > 0) {
      cerr << "No 'CDS' found in sequence " << annoseq->seqname << endl;
    }
    /*
     * make a list of genes from the data in pos
     *
     */
    list<Gene*> allgenes, *allgenesSorted;
    int curGeneNr = 1;
    list<GBFeature>::iterator fit, mRNAit;
    int lastgeneend=-1;
    Gene *gene;

    pos->CDS.sort();
    pos->mRNA.sort();

    // make genes from the CDS first
    
    for (fit = pos->CDS.begin(); fit != pos->CDS.end(); fit++) {
      bool cont=true;
      if (fit->begin <= lastgeneend) {
	cerr << "CDS overlaps with other CDS (" << (fit->begin+1) << " <= " << (lastgeneend+1) << ") in sequence " << annoseq->seqname
	     << ". Ignoring this gene." << endl;
	cont = false;
      } 
      if (cont && gbVerbosity > 0 && (!fit->complete_l || !fit->complete_r)){
	  cerr << "Coding sequence incomplete in sequence " << annoseq->seqname << "." << endl;
	  //cont = false;
      }
      if (cont) {
	gene = new Gene;
	// set the id, geneid, seqname and source
	gene->id = "1";
	
	if (fit->geneid != "") {
	  gene->geneid = fit->geneid;
	} else {
	    gene->geneid = string(annoseq->seqname) + "-" + itoa(curGeneNr);
	    
	  // char *temp = new char[strlen(annoseq->seqname)+6];
	  // sprintf(temp, "%s-%d",annoseq->seqname, curGeneNr );
	  // gene->geneid = temp;
	}
	gene->seqname = annoseq->seqname;

	gene->source = "database";

	gene->length   = pos->seqlength;
	gene->clength  = 0;
	gene->codingstart = gene->codingend = -1;
	gene->strand   = fit->strand;
	gene->bc       = annoseq->bc;
	gene->exons    = fit->ranges;
	/*
	 * construct the coding exons and introns
	 */
	State *exon_old=NULL, *intron_old=NULL, *intron=NULL, *st = gene->exons;
	while( st ){
	  gene->clength += st->length();
	  if (st == gene->exons)
	    if (st->next == NULL)
	      st->type =  (fit->strand == plusstrand)? singleG : rsingleG;
	    else
	      if (fit->strand == plusstrand)
		st->type = initialExon(gene->clength);
	      else 
		st->type = rterminalExon(2 - gene->clength);
	  else if (st->next == NULL)
	    if (fit->strand == plusstrand)
	      st->type = terminal;
	    else 
	      st->type = rinitial;
	  else 
	    if (fit->strand == plusstrand)
	      st->type = internalExon(gene->clength);
	    else 
	      st->type = rinternalExon(2 - gene->clength);
	  
	  if (stopCodonExcludedFromCDS) {
	      if (st->type == singleG || st->type == terminal){
		  	    gene->clength += 3;
			    st->end += 3;
	      }
	      if (st->type == rsingleG || isRTerminalExon(st->type)){
		  gene->clength += 3;
		  st->begin -= 3;
	      }
	      if (st->begin < 0 || st->end >= annoseq->length)
		  throw GBError("Stop codon out of sequence bounds. Ignoring sequence.");
	  }

	  if (gene->codingstart < 0 || gene->codingstart > st->begin)
	    gene->codingstart = st->begin;
	  if (gene->codingend < st->end)
	    gene->codingend = st->end;
	  if( exon_old ){
	    intron = new State();
	    intron->begin = exon_old->end+1;
	    intron->end   = st->begin-1;
	    if (intron->begin > intron->end){
		cerr << "Error: In sequence " << annoseq->seqname << ": One CDS exon does not begin properly after the previous CDS exon."
		     << exon_old->end << " >= " << st->begin << endl;  
		throw GBError("Intron has non-positive length.");
	    }
		
	    intron->type = intron_type;
	    intron->next  = 0;
	    if( !gene->introns )
	      gene->introns = intron;
	    if( intron_old )
	      intron_old->next = intron;
	    intron_old = intron;
	  }
	  exon_old=st;
	  st = st->next;
	}
	if (withUTR){
	  /*
	   * look if there is any fitting mRNA: positions compatible (and same geneid)
	   * take the first fitting but delete all fitting mRNAs
	   * construct the utr
	   */
	  bool mrnafound = false;
	  mRNAit = pos->mRNA.begin();
	  GBFeature mrnaFeature;
	  while(mRNAit != pos->mRNA.end()){
	    if (fit->matches(*mRNAit)){
	      if (!mrnafound)
		mrnaFeature = *mRNAit;
	      mrnafound = true;
	      mRNAit = pos->mRNA.erase(mRNAit);
	    } else {
	      mRNAit++;
	    }
	  }
	  if (mrnafound){
	    // construct the utr
	    gene->addUTR(mrnaFeature.ranges, mrnaFeature.complete_l, mrnaFeature.complete_r);
	  }
	}

	// gene->computeBC(annoseq->sequence);
	allgenes.push_back(gene);
	lastgeneend = fit->end;
	curGeneNr++;
      }
    }

    if (withUTR){
	/*
	 * make incomplete genes from the remaining mRNAs
	 */
	for (mRNAit = pos->mRNA.begin(); mRNAit != pos->mRNA.end(); mRNAit++) {
	    gene = new Gene;
	    // set the id, geneid, seqname and source
	    gene->id = "1";
	    if (mRNAit->geneid != "") {
		gene->geneid = mRNAit->geneid;
	    } else {
		gene->geneid = string(annoseq->seqname) + "-" + itoa(curGeneNr);

		// char *temp = new char[strlen(annoseq->seqname)+6];
		// sprintf(temp, "%s-%d",annoseq->seqname, curGeneNr );
		// gene->geneid = temp;
	    }
	    gene->seqname = annoseq->seqname;
	  
	    gene->source = "database";
	  
	    gene->length   = pos->seqlength;
	    gene->clength  = 0;
	    gene->codingstart = gene->codingend = -1;
	    gene->strand   = mRNAit->strand;
	    gene->bc       = annoseq->bc;
	    gene->exons    = NULL;
	    // add the UTR
	    gene->addUTR(mRNAit->ranges, mRNAit->complete_l,  mRNAit->complete_r);
	    // check whether it overlaps a CDS
	    // If yes, don't use it.
	    bool overlapsCDS = false;
	    for (list<Gene*>::iterator git = allgenes.begin(); git != allgenes.end(); git++)
		if ((gene->transstart <= (*git)->codingstart && gene->transend >= (*git)->codingstart) ||
		    (gene->transstart >= (*git)->codingstart && gene->transstart <= (*git)->codingend)){
		    overlapsCDS = true;
		}

	    if ((gene->utr5exons || gene->utr3exons) && !overlapsCDS){
		// gene->computeBC(annoseq->sequence);
		allgenes.push_back(gene);
		curGeneNr++;
	    }
	}
    }
    if (allgenes.empty()) {
      cerr << "No correct 'CDS'"; 
      if (withUTR)
	cerr << "or 'mRNA'";
      cerr << " found in sequence " << annoseq->seqname << endl;
    }
    allgenesSorted = sortGenePtrList(allgenes);
    for (list<Gene*>::iterator git = allgenesSorted->begin(); git != allgenesSorted->end(); git++)
      annoseq->anno->appendGene(*git);
    //annoseq->printGFF();
    return annoseq;
}

/*
 * GBProcessor::getAnnoSequenceList()
 *  
 * gets the list of all annotated sequences from 
 *
 */
AnnoSequence* GBProcessor::getAnnoSequenceList(){
    AnnoSequence *first = NULL, *annoseq, *last = NULL;
    bool notend = true;
    int annocount = 0;
    while (notend){
	try {
	    annoseq = getAnnoSequence(nextPosition( ));
	    if ((notend = (annoseq != NULL))){
		annoseq->next = NULL;
		if (last)
		    last->next = annoseq;
		else
		    first = annoseq;
		last = annoseq;
		annocount++;
	    }
	} catch (GBError e) {
	  cerr << "GBProcessor::getGeneList(): " << e.getMessage() << endl;
	  cerr << "Encountered error after reading " << annocount << " annotations." << endl;
 	}
    }
    if (annocount==0) 
	throw ProjectError("No genbank sequences found.");
    if (gbVerbosity)
	cout << "# Read in " << annocount << " genbank sequences." << endl;
    
    return first;
}

AnnoSequence* GBProcessor::getSequenceList(){
    AnnoSequence *seqlist = NULL, *seq, *last = NULL;
    while((seq = gbs.getNextFASTASequence())){
	seq->next = NULL;
	if (last)
	    last->next = seq;
	else
	    seqlist = seq;
	last = seq;
    };
    gbs.clear();
    if (seqlist == NULL)
	throw ProjectError("No sequences found.");
    return seqlist;
}

char* GBProcessor::getSequence( GBPositions& pos ) throw( GBError) {
    char* seq = new char[pos.seqlength+1];
    if( seq ){
        char* seqp = seq;
	istringstream isstrm( pos.seqbegin );
        char c;
        int  i;
	int count = 0;
        isstrm.getline(seq, pos.seqlength+1);
        while( isstrm ){
            isstrm >> i;  // Zahl zu Begin der Zeile lesen!
            while( isstrm.get( c ) && c != '\n' )
                if( isalpha( c ) ){
		    if (count++ < pos.seqlength) 
			*seqp++ = c; // tolower now postponed to after softmasking detection
		    else {
			char buffer [140];
			sprintf (buffer, "Sequence was longer than the expected %d bp. Please compare the line 'source 1..X' with the actual sequence length.", pos.seqlength);
			throw GBError(buffer);
		    }
		}
	}
	if (count < pos.seqlength) {
	    cerr << "Warning: Sequence was shorter than expected frome the source line of the Genbank entry ( expected " << pos.seqlength 
		 << " bp). Please compare the line 'source 1..X' with the actual sequence length (" << count << ")." << endl;
	  pos.seqlength = count;
	} 
	seq[pos.seqlength]='\0'; 
    }
    return seq;
}


GBFeature::GBFeature(const char *pos) throw( GBError ){
  begin = end = -1;
  ranges= NULL;
  complete_r=complete_l=true;
  istringstream isstrm( pos );
  char c;
  char tmp[GBMAXLINELEN] = { '\0' }; // need more than 4200 chars for human muscle gene TTN
  char* join = (char*)0;
  int i = 0;

  isstrm >> fkey >> ws;    // read 'CDS','mRNA' or signal FKEY
  c = isstrm.peek( );
  if (c == 'c') { // "complement"
    strand = minusstrand;
    while( isstrm.get( c ) && c != '(' && c != '\n')
      ;
  } else
    strand = plusstrand;
  c = isstrm.peek( );
  if( isalpha( c ) ){
    while( isstrm.get( c ) && c != '(' && c != '\n')
      ;
    if ( c == '\n') {
      cerr << "Feature key = " << fkey << ", c = " << c << endl;
      throw GBError( "GBProcessor::getJoin( ):  failed!!!" );
    }
    while( isstrm.get( c ) && c != ')' )
      if( !isspace( c ) ) {
	if (isdigit(c) || c =='.' || c == ',' || c=='<' || c=='>' )
	  tmp[i++] = c;
	else {
	  cerr << fkey << " contains character " << c << endl;
	  throw GBError( "GBProcessor::getJoin( ):  failed!!!" );
	}
      }
    while( isstrm.get( c ) && c != '\n' )
      ;         // read up to the end of the line
  } else if( isdigit( c ) || c == '<'){
    while( isstrm.get( c ) && c != ')' && c != '\n' )
      tmp[i++] = c;
    while (isstrm && (c != '\n'))  // read until end of line
      isstrm.get( c );
  } else {  
    cerr << "Feature key = " << fkey << ", c = " << c << endl;
    throw GBError( "GBProcessor::getJoin( ):  failed!!!" );
  }
  
  join = newstrcpy(tmp,i);

  // try to determine the name of the gene
  char *genename;
  char buf[GBMAXLINELEN];
  bool searchFinished = false;
  char *gtag;
  while( isstrm && !searchFinished){
    isstrm.getline( buf, GBMAXLINELEN-1 );
    if (strncmp(buf, "                     ", 21) != 0) {  // name of the gene not found
      geneid="";
      searchFinished = true;
    }
    else {
      gtag = strstr(buf, "/gene=");
      if (gtag){// found the name of the gene in the annotation
	int cplen = strcspn(gtag + 7, "\"");
	genename = newstrcpy(gtag + 7, cplen);
	geneid = genename;
	searchFinished = true; 
      }
    }    
  }
 
  /*
   * Now construct the ranges state sequence
   */
  try {
    State *exon = 0, *exon_old = 0;
    istringstream joinstrm( join );
    while( joinstrm ){
      char c;
      int ebegin=-1, eend=-1;
      exon = new State();

      c = joinstrm.peek( );
      if (c=='<'){
	complete_l = false;
	joinstrm.get( c );
      }
    
      joinstrm >> ebegin >> c >> c;
      c = joinstrm.peek( );
      if (c=='>'){
	complete_r = false;
	joinstrm.get( c );
      }
      
      joinstrm >> eend >> c;  // format:  'Numer..Number,'
      if (ebegin < 1 || eend < 1)
	throw ProjectError(string("Wrong format for coordinates: ") + join);
      if (ebegin > eend)
	throw ProjectError(string("Feature begins after it ends: ") + join);
      exon->begin = ebegin - 1; // correct for the fact that indices
      exon->end   = eend - 1;   // start with 0 in the c++ sequence
      exon->next = 0;

      if( !ranges )
	ranges = exon;
      if (begin == -1 || exon->begin < begin)
	begin = exon->begin;
      if (end == -1 || exon->end > end)
	end = exon->end;
    
      if( exon_old ){
	exon_old->next = exon;
      }    
      exon_old = exon;
    }
  } catch (ProjectError e) {
    cerr << "Constructing GenBank feature: " << e.getMessage() << endl;
    throw GBError("GBFeature constructor:Format error when reading genbank format.");
  }
}

bool GBFeature::checkRange(int len){
  return (begin >= 0 && end >= 0 && begin<len && end<len);
}

/*
 * GBFeature::matches
 * checks whether other extends this, as mRNA should extend CDS
 * in the range covered by this, the ranges should be identical
 */

bool GBFeature::matches(GBFeature &other){
  // if genenames exist they should be identical
  if (geneid != "" && other.geneid != "" && geneid != other.geneid)
    return false;
  if (strand != other.strand)
    return false;
  State *st=ranges, *otherst=other.ranges;
  if (!st || !otherst)
    return false;
  while (otherst && otherst->end < st->end)
    otherst = otherst->next;
  if (!(otherst && (otherst->begin <= st->begin) && (otherst->end == st->end || (st->next == NULL && otherst->end >= st->end))))
    return false;
  while (otherst->next && st->next){
    st=st->next;
    otherst=otherst->next;
    if (!(st->next == NULL || st == ranges) // internal exon
	&& !(st->begin == otherst->begin && st->end == otherst->end))// doesn't match exactly
      return false;
  }
  if (st->next)
    return false;
  if (! (otherst->end >= st->end && (otherst->begin == st->begin || (st == ranges && otherst->begin <= st->begin))))
    return false;  
  return true;
}

GBSplitter::GBSplitter( string fname ) : ftype(unknown) {
    if (fname != "-") {
	ifstrm.open(fname.c_str());    
	if( !ifstrm )
	    throw GBError("Could not open input file \"" + fname + "\"!");
    } else { // read from standard input
	throw GBError("Input from STDIN not supported anymore since the introcution of .gzipped input (with version 3.0).");
	//ifstrm.ios::rdbuf(cin.rdbuf());
	//ftype = fasta;
    }
#ifdef ZIPINPUT
    // deflate if gzipped
    boost::iostreams::filtering_istream zin;
    try {
	zin.push(gzip_decompressor());
	zin.push(ifstrm);
	zin.peek();
	if (!zin)
	    throw("Could not read first character assuming gzip format.");
	cout << "# Looks like " << ((fname != "-")? fname: "STDIN") 
	     << " is in gzip format. Deflating..." << endl;
    } catch (...) { // boost::iostreams::gzip_error& 
	// not a gzip file or ill-formatted
	zin.reset();
	ifstrm.seekg(0);
	zin.push(ifstrm);
    }
    boost::iostreams::copy(zin, sin);

#else
    // only normal file input available
    sin << ifstrm.rdbuf();
#endif
    determineFileType();
}

GBSplitter::~GBSplitter( ){
    if (ifstrm.is_open())
        ifstrm.close();
}



void GBSplitter::determineFileType(){
    std::stringstream csin(sin.str());
    char wrongChar=' ';
    ftype = unknown;
    csin >> ws;

    // check whether it could be genbank format
    // search for at least one LOCUS and ORIGIN at the beginning of a line
    csin >> goto_line_after( "LOCUS" );
    csin >> goto_line_after( "ORIGIN" );
    if (csin)
	ftype = genbank;
    else {
	// check, whether it is FASTA or plain sequence
	// only letters allowed in the sequence part, except for the end of the line
	bool whiteSpaceSeen;
	csin.clear();
	csin.seekg(0, ios::beg);
	string line;
	bool haveWrongChar=false;
	while (csin && !haveWrongChar) {
	    getline(csin, line);
	    if (line[0]!='>') {
		whiteSpaceSeen=false;
		for (int i=0; i < line.length() && !haveWrongChar; i++) 
		    if (isspace(line[i]))
			whiteSpaceSeen = true;
		    else if (!isalpha(line[i]) && whiteSpaceSeen){
			haveWrongChar = true;
			wrongChar = line[i];
		    }
	    }
	}
	if (!haveWrongChar)
	    ftype = fasta;
	
	csin.clear();
	csin.seekg(0, ios::beg);
    }
    if (ftype == unknown) {
     	string errmsg = "GBProcessor::determineFileType(): Couldn't determine input file type. Found bad character '";
	errmsg += wrongChar;
     	errmsg += "'";
     	throw GBError(errmsg);
    }

    //    csin.clear();
    //    csin.seekg(0, ios::beg);
}


GBPositions* GBSplitter::nextData( ) throw( GBError ){
    GBPositions* pos = new GBPositions;
    if( !findPositions( *pos ) ){     // No other data available!
        delete pos;
        return (GBPositions*)0;
    }
    return pos;
}

AnnoSequence *GBSplitter::getNextFASTASequence( ) throw( GBError ){
    char *sequence = NULL, *name = NULL;
    int length;
    readOneFastaSeq(sin, sequence, name, length);
    if (sequence == NULL || length == 0)
	return NULL;

    AnnoSequence *seq = new AnnoSequence();
    seq->seqname = name;
    seq->sequence = sequence;
    seq->length = length;
    return seq;
}


Boolean GBSplitter::gotoEnd( ){
    char buf[GBMAXLINELEN];
    do{
        int i = 0;
        sin.getline( buf, GBMAXLINELEN );
	if ((sin.rdstate() & ios_base::failbit) && !(sin.rdstate() & ios_base::eofbit))
	    sin.clear(sin.rdstate() & ~ios_base::failbit);
        while( i < GBMAXLINELEN-1 && isspace(buf[i]))
            i++;
        if( buf[i] == '/' && buf[i+1] == '/' )
            return true;
    } while (sin);
    return false;
}


Boolean GBSplitter::findPositions( GBPositions& pos ) throw( GBError ){
    int fposb, fpose;
    fposb = sin.tellg();
    if( !gotoEnd( ) )
        return false;
    fpose = sin.tellg();
    sin.seekg( fposb );

    pos.length = fpose-fposb /*+1*/;         // Without the '\0'!!
    pos.buffer = new char[pos.length+1];
    sin.read( pos.buffer, pos.length );
    pos.buffer[pos.length] = '\0';
    pos.length++;                       // Now with the '\0'!!!
    pos.seqlength = 0;
    istringstream isstrm( pos.buffer );
    char buf[GBMAXLINELEN];
    while( isstrm ){
        int curpos = isstrm.tellg();
        isstrm >> ws;
        isstrm.getline( buf, GBMAXLINELEN-1 );
	if ((!sin.eof() && (sin.rdstate() & ios_base::failbit)) || strlen(buf) >= GBMAXLINELEN-2){
	    throw GBError(string("Could not read the following line in Genbank file.\n") + buf 
			  + "\nMaximum line length is \n" + itoa(GBMAXLINELEN-2) + ".\n");
	}
        char *src;
        char *rna;
        char *cds;
        char *seq;

        src = strstr( buf, "source" );
        if( src && src == buf ){
            istringstream strm(buf);
            string str;
            int a, b;
            char c;
            strm >> str >> a >> c >> c >> b;
	    if (b-a+1 > pos.seqlength)
		pos.seqlength = b-a+1;
            continue;
        }
        rna = strstr( buf, "mRNA         " );

        if( rna && rna == buf ){
	    GBFeature mrnaFeature(pos.buffer+curpos);
	    if (mrnaFeature.checkRange(pos.seqlength))
		pos.mRNA.push_back(mrnaFeature);
	    else {
		cout << "mRNA " << mrnaFeature.geneid << " out of range. Ignoring it." << endl;
		cout << "end= " << mrnaFeature.end << " len= " << pos.seqlength << endl;
	    }
        }
	cds = strstr( buf, "CDS             " );
	if( cds && cds == buf ){
	    GBFeature cdsFeature( pos.buffer+curpos );
	    if(cdsFeature.checkRange(pos.seqlength))
		pos.CDS.push_back( cdsFeature );
	    else
	        cout << "CDS " << cdsFeature.geneid << " out of range. Ignoring it." << endl;
	}
	seq = strstr( buf, "ORIGIN" );
	if( seq && seq == buf ){
	    if (curpos >= pos.length)
		cerr << "In gene " << pos.CDS.front().geneid << " a coordinate points after the sequence." << endl;
	  pos.seqbegin = pos.buffer+curpos;
	  break;
	}
    }  // while
    if (pos.seqlength == 0) {
	throw GBError("Sequence has 0 length. Maybe 'source' Feature missing?");
    }
    return true;
}
