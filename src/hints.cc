/**********************************************************************
 * file:    hints.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  hints, user constraints
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 20.10.06| Mario Stanke  | creation of the file by splitting extrinsicinfo.cc
 **********************************************************************/

#include "hints.hh"
#include "properties.hh"
#include "projectio.hh"

#include <iostream>
#include <cstdlib>
#include <set>

const char* featureTypeNames[NUM_FEATURE_TYPES]= {"start", "stop", "ass", "dss", "tss", "tts",
						  "exonpart", "exon", "intronpart", "intron",
						  "irpart", "CDS", "CDSpart", "UTR", "UTRpart", "nonexonpart", "genicpart"};

bool isSignalType(FeatureType type){
    return (type == startF || type == stopF || type == assF || type == dssF || type == tssF || type == ttsF);
}

bool isGFF(ifstream &ifstrm){
    ifstrm.seekg(0);
    // skip comments  
   ifstrm >> comment >> ws;
    if(!(ifstrm))
        return false;
    streampos spos = ifstrm.tellg();
    string line;
    getline(ifstrm,line);
    size_t pos = line.find('\t');
    if (pos != string::npos) {
        ifstrm.clear();
        ifstrm.seekg(spos);
        return true;
    }
    // line not tap separated   
    return false;
}

ostream& operator<<(ostream&out, Feature& feature){
    out << feature.seqname << "\t"
	<< feature.source << "\t"
	<< feature.feature << "\t"
	<< Feature::offset + feature.start + 1 << "\t" 
	<< Feature::offset + feature.end + 1 << "\t"
	<< feature.score << "\t";
    switch (feature.strand) {
	case plusstrand: out << "+"; break;
	case minusstrand: out << "-"; break;
	default: out << ".";
    }
    out << "\t";
    if (feature.frame != -1)
	out << feature.frame;
    else 
	out << ".";
    out << "\t";
    if (feature.groupname != "")
	out << "grp=" << feature.groupname << ";";
    if (feature.mult > 1)
	out << "mult=" << feature.mult << ";";
    if (feature.priority >= 0)
	out << "pri=" << feature.priority << ";";
    out	<< "src=" << feature.esource
	//<< feature.attributes
	<< " \"" << feature.bonus << ";" << feature.malus << ";" << feature.numSupporting << ":" << feature.numContradicting <<"\"";
    return out;
}

istream& operator>>( istream& in, Feature& feature ){
    char buff[1024]; 
    char copybuff[1024];
    char *temp;
    const char *spos;
    try {
	in.getline(buff, 1024);
	strncpy(copybuff, buff, 1014);

	if (strstr(buff, "\t")==NULL) {
	    throw ProjectError("Line not tab separated.");
	}
	temp = strtok(buff, "\t");
	if (temp)
	    feature.seqname = temp;
	else 
	    throw ProjectError("Could not read sequence name.");
	
	temp = strtok(NULL, "\t");
	if (temp)
	    feature.source = temp;
	else
	    throw ProjectError("Could not read second column.");
	temp = strtok(NULL, "\t");
	if (temp)
	    feature.feature = temp;
	else 
	    throw ProjectError("Could not read feature type.");
	temp = strtok(NULL, "\t");
	if (temp)
	    feature.start = atoi(temp);
	else 
	    throw ProjectError("Could not read start position.");
	temp = strtok(NULL, "\t");
	if (temp)
	    feature.end = atoi(temp);
	else 
	    throw ProjectError("Could not read end position.");
	temp = strtok(NULL, "\t");
	if (temp) 
	    feature.score = atof(temp);
	else 
	    throw ProjectError("Could not read score.");
	temp = strtok(NULL, "\t"); // stand
	if (!temp)
	    throw ProjectError("Could not read strand.");
	if (strcmp(temp, "+") == 0)
	    feature.strand = plusstrand;
	else if (strcmp(temp, "-") == 0)
	    feature.strand = minusstrand;
	else 
	    feature.strand = STRAND_UNKNOWN;
	temp = strtok(NULL, "\t"); // frame
	if (!temp)
	    throw ProjectError("Could not read frame.");
	if (strcmp(temp, "0") == 0)
	    feature.frame = 0;
	else if (strcmp(temp, "1") == 0)
	    feature.frame = 1;
	else if (strcmp(temp, "2") == 0)
	    feature.frame = 2;
	else 
	    feature.frame = -1;
	
	temp = strtok(NULL, "\t");
	if (temp){
	    feature.attributes = temp;
	    //replace dos carriage return by space
	    for (int i=0; i < feature.attributes.length(); i++) 
		if (feature.attributes[i]=='\r')
		    feature.attributes[i]=' ';
	}
	else 
	    throw ProjectError("Could not read last column.");
	
	/*
	 * find groupname of hint, specified in gff as: group=xxx; or grp=xxx;
	 */
	spos = strstr(feature.attributes.c_str(), "group=");
	if (spos)
	    spos += 6;
	if (!spos) {
	    spos = strstr(feature.attributes.c_str(), "grp=");
	    if (spos)
		spos += 4;
	}
	if (spos){
	    int skeylen=1;
	    while (*(spos+skeylen) != ';' && *(spos+skeylen) != ' ' && *(spos+skeylen) != '\0'){
		skeylen++;
	    }
	    feature.groupname = string(spos, skeylen);
	} else {
	    feature.groupname = "";
	}
	
	/*
	 * find priority of hint, specified in gff as: priority=N; or pri=N;
	 * higher number means higher priority
	 */
	spos = strstr(feature.attributes.c_str(), "priority=");
	if (spos)
	    spos += 9;
	if (!spos) {
	    spos = strstr(feature.attributes.c_str(), "pri=");
	    if (spos)
		spos += 4;
	}
	if (spos){
	    feature.priority = atoi(spos);
	} else {
	    feature.priority = -1;
	}
	/*
	 * find multiplicity of hint, specified in gff as: mult=N;
	 */
	spos = strstr(feature.attributes.c_str(), "mult=");
	if (spos)
	    spos += 5;
	if (spos)
	    feature.mult = atoi(spos);
	else
	    feature.mult = 1;
	
	/*
	 * find source of extrinsic info, specified in gff as: source=X or src=X
	 */
	spos = strstr(feature.attributes.c_str(), "source=");
	if (spos)
	    spos += 7;
	if (!spos) {
	    spos = strstr(feature.attributes.c_str(), "src=");
	    if (spos)
		spos += 4;
	}
	if (spos && isalpha(*spos)){
	    int skeylen=1;
	    while (isalpha(*(spos+skeylen))){
		skeylen++;
	    }
	    feature.esource = string(spos, skeylen);
	} else {
	    cerr << "Error in hint line: " << copybuff << endl;
	    cerr << "No source specified (e.g. by source=M in the last column)" << endl;
	}
	feature.gradeclass = 0;// is set by class FeatureCollection
	
	feature.start--; // shift by -1, cause I start indexing with 0
	feature.end--;   // shift by -1, cause I start indexing with 0
	
	feature.type = Feature::getFeatureType(feature.feature); // may fail: type = -1
    } catch (ProjectError e) {
	cerr << "Error in hint line: " << copybuff << endl;
	cerr << e.getMessage() << endl;
	cerr << "Maybe you used spaces instead of tabulators?" << endl;
	throw e;
    }
    return in;
}

int Feature::offset = 0;

FeatureType Feature::getFeatureType(string typestring){
    FeatureType type;
    if (typestring == "dss")
	type = dssF;
    else if (typestring == "ass")
	type = assF;
    else if (typestring == "stop")
	type = stopF;
    else if (typestring == "start")
	type = startF;
    else if (typestring == "exonpart" || typestring == "ep")
	type = exonpartF;
    else if (typestring == "exon")
	type = exonF;
    else if (typestring == "intronpart" || typestring == "ip")
	type = intronpartF;  
    else if (typestring == "intron")
	type = intronF;
    else if (typestring == "tss")
	type = tssF;
    else if (typestring == "tts")
	type = ttsF;
    else if (typestring == "irpart")
	type = irpartF;
    else if (typestring == "CDS")
	type = CDSF;
    else if (typestring == "CDSpart" || typestring == "cp") 
	type = CDSpartF;
    else if (typestring == "UTR")
	type = UTRF;
    else if (typestring == "UTRpart" || typestring == "up")
	type = UTRpartF;
    else if (typestring == "nonexonpart" || typestring == "nep")
	type = nonexonpartF;
    else if (typestring == "nonirpart" || typestring == "genicpart")
	type = nonirpartF;
    else {
      cerr << "Unknown hint feature '" << typestring << "'. Ignoring it. For a list of the 17 allowed features ";
      cerr << "see the first column of the table in config/extrinsic/extrinsic.cfg." << endl;
      type = (FeatureType) -1;
    }
    return type;
}

FeatureType Feature::getFeatureType(int typeint){
    if(typeint < NUM_FEATURE_TYPES && typeint >=0)
	return (FeatureType)typeint;
    return (FeatureType) -1;
}

/*
 * Feature::compatibleWith
 * True if a gene structure could exist that is compatible with both features.
 */
bool Feature::compatibleWith(Feature &other){
    const int transcript_fuzzy_margin = 50; // a tss/tts at most this far from the end of an exon/UTR end is still considered compatible
    const int term3_M = 1000; // two tts at most this far apart but not overlapping and not suggesting tail to tail genes are imcompatible
    const int term5_M = 0; // two tss at most this far apart but not overlapping and not suggesting head to head genes are imcompatible
    
    if (start > other.end || end < other.start){ // no overlap
	if ((type == tssF && other.type == tssF) &&
	    (abs((end+start)-(other.end+other.start))/2 <= term5_M) && // close
	    (strand == other.strand)){
		return false; // tss on same strand and close but not overlapping (macro-heterogeneity)
	}
	if ((type == ttsF && other.type == ttsF) &&
	    (abs((end+start)-(other.end+other.start))/2 <= term3_M) && // close
	    (strand == other.strand)){
		return false; // tts on same strand and close but not overlapping (macro-heterogeneity)
	}
	return true;
    }
    
    // overlapping case follows
    // opposite strands are incompatible unless we have at least one signal type
    if ((strand == minusstrand && other.strand == plusstrand) || (strand == plusstrand && other.strand == minusstrand)) {
	// this is an approximation, should take into account that min length of start and stop is 3
	if (isSignalType(type) && isSignalType(other.type))
	    return true;
	if (isSignalType(type) && (start < other.start || end > other.end))
	    return true;
	if (isSignalType(other.type) && (start > other.start || end < other.end))
	    return true;
	return false;
    }
    // same strand (or bothstrands) follows
    if (type == other.type) { // same type
	if (start == other.start && end == other.end)
	    return true;
	if (type == exonF || type == intronF || type == CDSF || type == UTRF)
	    return false;
	return true;
    }
    Feature *f1=this, *f2=&other;
    // different types, wlg let f1 have the smaller index, otherwise swap
    if ((int) f1->type > (int) f2->type){
	f1 = &other;
	f2 = this;
    }
    if (f1->type == startF || f1->type == stopF){
	if ((f2->type == intronpartF || f2->type == intronF || f2->type == irpartF || f2->type == nonexonpartF || f2->type == UTRF || f2->type == UTRpartF) &&
	    f1->start>= f2->start && f1->end <= f2->end) // start/stop interval included in non-CDS interval
	    return false;
	if (f2->type == CDSpartF || f2->type == CDSF){
	    if (f1->start > f2->start && f1->end < f2->end) // start/stop interval strictly included in CDS interval
		return false;
	    // rest: true overlap or start/stop is included in CDS/CDSpart and right at the end
	    // must not be at the wrong end
	    if (strand == plusstrand && ((f1->type == startF && f1->start > f2->start) || (f1->type == stopF && f1->end < f2->end)))
		return false;
	    if (strand == minusstrand && ((f1->type == startF && f1->end < f2->end) || (f1->type == stopF && f1->start > f2->start)))
		return false;
	    if (f2->type == CDSF && (f1->end < f2->start + 2 || f1->start > f2->end-2)) // overlap less than 3
		return false;
	    return true;
	}
	return true;
    }
    if (f1->type == assF || f1->type == dssF){
	if ((f2->type == irpartF || f2->type == UTRF || f2->type == UTRpartF || f2->type == exonpartF || f2->type == exonF || f2->type == CDSF || f2->type == CDSpartF) &&
	    f1->start>= f2->start && f1->end <= f2->end) // splice site interval included in non-CDS interval
	    return false;
	if (f2->type == intronF || f2->type == intronpartF){
	    if (f1->start > f2->start && f1->end < f2->end) // splice site interval strictly included in intron interval
		return false;
	    // rest: true overlap or splice site is included in intron/intronpart and right at an end
	    // must not be at the wrong end
	    if (strand == plusstrand && ((f1->type == dssF && f1->start > f2->start) || (f1->type == assF && f1->end < f2->end)))
		return false;
	    if (strand == minusstrand && ((f1->type == dssF && f1->end < f2->end) || (f1->type == assF && f1->start > f2->start)))
		return false;
	    return true;
	}
	return true;
    }
    if (f1->type == tssF || f1->type == ttsF){
	if ((f2->type == irpartF || f2->type == intronF || f2->type == intronpartF || f2->type == nonexonpartF || f2->type == CDSF || f2->type == CDSpartF) &&
	    f1->start >= f2->start && f1->end <= f2->end) // transcript start/end interval included in non-UTR interval
	    return false;
	if (f2->type == UTRF || f2->type == UTRpartF || f2->type == exonF || f2->type == exonpartF){
	    if (f1->start > f2->start + transcript_fuzzy_margin && f1->end < f2->end - transcript_fuzzy_margin) // transcript start/end interval strictly included in UTR/exon interval
		return false;
	    // rest: true overlap or transcript start/end site is included in UTR/exon and right at an end
	    // must not be at the wrong end
	    if (strand == plusstrand && ((f1->type == tssF && f1->start > f2->start + transcript_fuzzy_margin) || (f1->type == ttsF && f1->end < f2->end - transcript_fuzzy_margin)))
		return false;
	    if (strand == minusstrand && ((f1->type == tssF && f1->end < f2->end - transcript_fuzzy_margin) || (f1->type == ttsF && f1->start > f2->start + transcript_fuzzy_margin)))
		return false;
	    return true;
	}
	return true;
    }
    if (f1->type == exonpartF) {
	if (f2->type == intronpartF || f2->type == intronF || f2->type == irpartF || f2->type == nonexonpartF) // types that are simply not compatible with exonpart
	    return false;
	if (f2->type == exonF && (f1->start < f2->start || f1->end > f2->end)) // exonpart not contained in exon
	    return false;
	if (f2->type == UTRF && (f1->start < f2->start && f1->end > f2->end)) // exonpart properly contains UTR
	    return false;
	return true;
    }
    if (f1->type == exonF) {
	if (f2->type == intronpartF || f2->type == intronF || f2->type == irpartF || f2->type == nonexonpartF) // types that are simple not compatible with exonpart
	    return false;
	if (f2->type == CDSF && !(f1->start <= f2->start && f1->end >= f2->end)) // exon must contain CDS
	    return false;
	if (f2->type == CDSpartF && (f1->start > f2->start || f1->end < f2->end)) // CDSpart covers something that exon does not cover
	    return false;
	if (f2->type == UTRF &&
	    !((f1->start == f2->start && f1->end >= f2->end) || (f1->end == f2->end && f1->start <= f2->end))) 
	    return false;
	if (f2->type == UTRpartF && (f1->start > f2->start || f1->end < f2->end))
	    return false;
	return true;
    }
    if (f1->type == intronpartF) {
	if (f2->type == intronF && (f1->start < f2->start || f1->end > f2->end))
	    return false;
	if (f2->type == irpartF || f2->type == CDSF || f2->type == CDSpartF || f2->type == UTRF || f2->type == UTRpartF)
	    return false;
	return true;
    }
    if (f1->type == intronF) {
	if (f2->type == irpartF || f2->type == CDSF || f2->type == CDSpartF || f2->type == UTRF || f2->type == UTRpartF)
	    return false;
	return true;
    }
    if (f1->type == irpartF) {
	if (f2->type == nonexonpartF)
	    return true;
	return false;
    }
    if (f1->type == CDSF) {
	if (f2->type == CDSpartF && (f1->start <= f2->start && f1->end >= f2->end)) // CDS contains CDSpart
	    return true;
	return false;
    }
    if (f1->type == CDSpartF) 
	return false; // incompatible with UTRF, UTRpartF, nonexonpartF
    if (f1->type == UTRF){
	if (f2->type == UTRpartF && (f1->start <= f2->start && f1->end >= f2->end))
	    return true;
	return false;
    }
    return false;
}

/*
 * Feature::weakerThan
 * True if every gene structure that is compatible with 'other' also is compatible with 'this'.
 * Strictly is true, if there could be gene structures that are only compatible with 'this' and not with other.
 * Usually, when a **partF inverval is larger in other. 'strictly' is only correctly set when the return value is true.
 */
bool Feature::weakerThan(Feature &other, bool &strictly){
    strictly = false;
    if (other.end < start || other.start > end)
	return false;
    if (type == other.type && start == other.start && end == other.end)
	return true;
    if (start != other.start || end != other.end)
	strictly = true;
    if (type == other.type && (type == startF || type == stopF || type == assF || type == dssF || type == tssF || type == ttsF) && (start <= other.start && end >= other.end))
	return true;
    if (type == exonpartF && (other.type == exonF || other.type == exonpartF) && (start >= other.start && end <= other.end))
	return true;
    if (type == intronpartF && (other.type == intronF || other.type == intronpartF) && (start >= other.start && end <= other.end))
	return true;
    if (type == irpartF && other.type == irpartF && (start >= other.start && end <= other.end))
	return true;
    if (type == CDSpartF && (other.type == CDSF || other.type == CDSpartF) && (start >= other.start && end <= other.end))
	return true;
    if (type == UTRpartF && (other.type == UTRF || other.type == UTRpartF) && (start >= other.start && end <= other.end))
	return true;
    if (type == nonexonpartF && other.type == nonexonpartF && (start >= other.start && end <= other.end))
	return true;
    if (type == nonirpartF && other.type != irpartF &&  (start >= other.start && end <= other.end))
	return true;
    return false;
}

/*
 * conformance
 * between 0 and 1. Close to 1 if lots of support and no contradiction.
 */
double Feature::conformance(){
    int pseudoSupporting = 5, pseudoContradicting = 5;
    return (double) (pseudoSupporting + numSupporting) / (pseudoContradicting + pseudoSupporting + numSupporting + numContradicting);
}

/*
 * shift the coordinates of a feature by -seqStart
 * if seqStrand is minusstrand, the coordinates and strand of the
 * feature are set relative to the minusstrand 
 */
void Feature::shiftCoordinates(int seqStart, int seqEnd, bool rc){
    if (!rc) {
	end -= seqStart;
	start -= seqStart;
    } else {
	int temp = end;
	end = seqEnd - start;
	start = seqEnd - temp;
	if (strand == plusstrand)
	    strand = minusstrand;
	else if (strand == minusstrand)
	    strand = plusstrand;
    }
}

void Feature::setFrame(string f){
    if(f == "0")
	frame=0;
    else if (f == "1")
	frame=1;
    else if (f == "2")
	frame=2;
    else
	frame=-1;
}

void Feature::setStrand(string s){
    if (s == "+")
	strand=plusstrand;
    else if (s == "-")
	strand=minusstrand;
    else
	strand=STRAND_UNKNOWN;
}

/*
 * operator<
 * is needed for the set operations in printAccuracyForSequenceSet
 * sorting is by increasing end positions
 */

bool operator<(const Feature& f1, const Feature& f2){
    if (f1.end < f2.end)
	return true;
    else if (f1.end > f2.end)
	return false;
    if (f1.start < f2.start)
	return true;
    else if (f1.start > f2.start)
	return false;
    if (f1.strand < f2.strand)
	return true;
    else if (f1.strand > f2.strand)
	return false;
    if (f1.frame < f2.frame)
	return true;
    else if (f1.frame > f2.frame)
	return false;
    if (f1.esource < f2.esource && f2.esource != "annotrain" && f1.esource != "annotrain")
	return true;
    return false;
}

bool operator==(const Feature& f1, const Feature& f2){
    return ((f1.type == f2.type) && (f1.end == f2.end) && (f1.start == f2.start) && (f1.strand == f2.strand) && (f1.frame == f2.frame));
}

double Feature::distance_faded_bonus(int pos){
    double result;
    if (pos < start || pos > end)
	result=1.0;
    else {
	double delta = 2.0*(pos - (end+start)/2.0)/(end-start+1); 
	if (delta < 0)
	    delta = -delta; // delta between 0 and 1, delta=0 when pos in the middle of hint interval, delta=1 when pos at an end of hint interval
	if (delta == 0.0)
	    result=bonus; // shortcut to save the expensive log and exp operations
	else 
	    result = exp(log(bonus) * (1-delta));
    }
    return result;
}

/*
 * HintGroup::addFeature
 * put the feature at the right place: sorted increasingly by end positions
 */
void HintGroup::addFeature(Feature *hint){
    if (hints == NULL) {
	hints = new list<Feature*>;
	name = hint->groupname;
    }
// TODO. This could be more efficient when the list is traversed from the end
    list<Feature*>::iterator fit = hints->begin();
    while(fit != hints->end() && (*fit)->end < hint->end)
	fit++;
    hints->insert(fit, hint);
    // if feature is genic then reset begin and end if neccessary
    // nongenic features: irpartF, nonexonpartF
    if (begin < 0 || begin > hint->start)
	begin = hint->start;
    if (end < 0 || end < hint->end)
	end = hint->end;
    if (hint->type != irpartF && hint->type != nonexonpartF) {
	if (geneBegin <0 || geneBegin > hint->start)
	    geneBegin = hint->start;
	if (geneEnd <0 || geneEnd < hint->end)
	    geneEnd = hint->end;
    }
    if (hint->priority > priority)
	priority = hint->priority;
    if (hint->mult > copynumber)
      copynumber = hint->mult; 
}

/*
 * HintGroup::compatibleWith
 * Two HintGroups g1, g2 are compatible if a gene structure could exist that is compatible with both HintGroups.
 * (Could exist, because I don't check for ORFs and such.)
 * All pairs of features f1 in g2 and f2 in g2 that overlap must be compatible.
 * rascal1 and rascal2 are the hints that caused the incompatibility if any
 * weakerThan is true if this is strictly weaker than 'other'.
 */
bool HintGroup::compatibleWith(HintGroup &other, Feature *&rascal1, Feature *&rascal2, bool &weakerThan){
    bool compatible = true;
    weakerThan = false;
    bool strictly = false;
    rascal1 = rascal2 = NULL;
    if (begin>other.end || end<other.begin)
	return compatible;
    // go through both feature lists and check all pairs of overlapping features
    if (!hints || !other.hints) 
	return true;
    bool featureWeakerThan, wt, sly;
    weakerThan = true;
    for (list<Feature*>::iterator f1it = hints->begin(); f1it != hints->end() && compatible; f1it++) {
	featureWeakerThan = false;
	for (list<Feature*>::iterator f2it = other.hints->begin(); f2it != other.hints->end() && compatible; f2it++) {
	    compatible &= (*f1it)->compatibleWith(**f2it);
	    if (!compatible){
		rascal1 = *f1it;
		rascal2 = *f2it;
	    }
	    wt = (*f1it)->weakerThan(**f2it, sly);
	    featureWeakerThan |= wt;
	    if (wt)
		strictly |= sly;
	}
	weakerThan &= featureWeakerThan;
    }
    if (weakerThan && !strictly) {
	// check if proper 'weaker than' relation comes from a subset of features (as opposed to a single feature being weaker)
	// make a simple approximation. If a feature of other lies completely outside the range of this and this is weaker than other
	// then it is also strictly weaker
	for (list<Feature*>::iterator f2it = other.hints->begin(); f2it != other.hints->end() && compatible; f2it++)
	    if ((*f2it)->end < begin || (*f2it)->start > end) 
		strictly = true;
    }
    weakerThan &= strictly;
    if (!compatible)
	weakerThan = false;
    return compatible;
}


/*
 * HintGroup::updateFeatureConformance
 * Updates conformance of the Features of this group to the other group.
 * For each Feature f of 'this', numSupporting is increased by the multiplicity of 'other'
 * if there is one Feature in 'other' that supports f. numContradicting is increased for f if 
 * there is at least one hint in 'other' that contradicts f.
 */
void HintGroup::updateFeatureConformance(HintGroup &other){
    if (end < other.begin || begin > other.end || hints == NULL || other.hints == NULL)
	return;
    // If this is uncommented then minor "dirt" exonpart hints in a well-supported intron are taken
    // seriously.
    //if (nestedGenePossible(other))
    //return;
    if (this == &other) {
	for(list<Feature*>::iterator f = hints->begin(); f != hints->end(); f++)
	    (*f)->numSupporting += copynumber - 1;
	return;
    }
    bool supporting, contradicting, strictly, only_ep_conflicting_intron;
    bool lowerpriority = (other.priority < priority && other.priority >=0);
    float fractContra = 1.0;
    for(list<Feature*>::iterator f = hints->begin(); f != hints->end(); f++){
	supporting = false;
	contradicting = false;
	only_ep_conflicting_intron = true;
	for(list<Feature*>::iterator otherF = other.hints->begin(); otherF != other.hints->end(); otherF++){
	    if (!lowerpriority && !(*f)->compatibleWith(**otherF)){
	         contradicting = true;
		 /*
		  * special case: the amount an exonpart contradicts an intron depends on the lengths of the two
		  * Reason: Assume constant ep length and constant splice frequencies. Then a longer intron
		  * will have more ep hints that contradict it.
		  */
		 if ((*f)->type == intronF && ((*otherF)->type == exonpartF || 
					       (*otherF)->type == CDSpartF || 
					       (*otherF)->type == UTRpartF)){
		   int ilen = (*f)->length();
		   if (ilen > 2000) // exonpart hints covering a range of a full long exon can fully contradict an intron
		     ilen = 2000;
		   if (ilen < 1)
		     ilen = 1;
		   int eplen = (*otherF)->length();// todo: just consider overlap length
		   if (eplen > ilen)
		     eplen = ilen;
		   fractContra = (float) eplen/ilen;
		 } else {
		   only_ep_conflicting_intron = false;
		 }
	    }
		
	    if ((*f)->weakerThan(**otherF, strictly))
		supporting = true;
	}
	if (supporting && !contradicting)
	    (*f)->numSupporting += other.copynumber;
	else if (contradicting){
  	    if (!only_ep_conflicting_intron)
	        fractContra = 1.0;
	    (*f)->numContradicting += fractContra * other.copynumber;
	}
    }
}

/*
 * HintGroup::nestedGenePossible
 * Return true if one hint group is contained in an intron of the other hintgroup.
 * Require minimal distance from intron of 50. This could be tightened by 
 * checking that the total mRNA of the inner gene must have enough space.
 */
bool HintGroup::nestedGenePossible(HintGroup &other){
    for(list<Feature*>::iterator otherF = other.hints->begin(); otherF != other.hints->end(); otherF++)
	if ((*otherF)->type == intronF && (*otherF)->start < begin -50 && (*otherF)->end > end +50)
	    return true;
    for(list<Feature*>::iterator thisF = hints->begin(); thisF != hints->end(); thisF++)
	if ((*thisF)->type == intronF && (*thisF)->start < other.begin -50 && (*thisF)->end > other.end +50)
	    return true;
    return false;
}

/*
 * HintGroup::isTrashy
 * Mark the group to be discarded iff any feature of the group contradicts a large number of features of other groups
 * and is not supported by a sufficient number of other groups.
 */
bool HintGroup::isTrashy(){
    if (hints == NULL)
	return false;
    float c;
    int s;
    list<Feature*>::iterator hit = hints->begin();
    while(hit != hints->end() && !trashy) {
	c = (*hit)->numContradicting;
	s = (*hit)->numSupporting;
	if (c/(s+1) >= Constant::max_contra_supp_ratio)
	    trashy = true;
	/* was
	   (s == 0 && c>=10) ||
	   (s == 1 && c>=15) ||
	   (s == 2 && c>=25) ||
	   (s == 3 && c>=45) ||
	   (s == 4 && c>=60) ||
	   (s >= 5 && c/s >= 20))
	*/
	hit++;
    }
    if (trashy) {
	//cout << "This group is trashy:"<<endl;
	//print(cout, true);
	setDiscardFlag(true);
    }
    return trashy;
}

/*
 * HintGroup::canCauseAltSplice
 * optional TODO: also implement a minimal multiplicity, e.g. intron:3,exon:1
 */
bool HintGroup::canCauseAltSplice() {
  static set<FeatureType> *causingTypes = NULL;
  if (causingTypes == NULL){ // initialize on first call
    causingTypes = new set<FeatureType>;
    char *canCauseAltSpliceList = NULL;
    try {
      canCauseAltSpliceList = newstrcpy(Properties::getProperty( "canCauseAltSplice" ));
    } catch(...) {
      canCauseAltSpliceList = newstrcpy("intron,exon,tss,tts,start,stop,ass,dss,ip,CDS,CDSpart,UTR,UTRpart"); // default
    }
    char *typestr = strtok(canCauseAltSpliceList, ",");
    while (typestr != NULL) {
      FeatureType type = Feature::getFeatureType(typestr);
      if (type >= 0)
	causingTypes->insert(type);
      typestr = strtok (NULL, ",");
    }
    if (canCauseAltSpliceList)
      delete [] canCauseAltSpliceList;
  }
  bool ret = false;
  if (hints){
    list<Feature*>::iterator fit = hints->begin();
    while(ret == false && fit != hints->end()){
      ret |= (causingTypes->find((*fit)->type) != causingTypes->end()); // can cause altsplice if ANY of the types in the group match
      fit++;
    }
  }
  // was before a hack (for rGASP):
  // return hints && (hints->size()>1 || (hints->front()->type != exonpartF));
  return ret;
}

/*
 * HintGroup::setActiveFlag
 */
void HintGroup::setActiveFlag(bool active){
    if (!hints)
	return;
    for (list<Feature*>::iterator hit = hints->begin(); hit != hints->end(); hit++)
	(*hit)->active = active;
}

/*
 * HintGroup::setDiscardFlag
 */
void HintGroup::setDiscardFlag(bool discard){
    if (!hints)
	return;
    for (list<Feature*>::iterator hit = hints->begin(); hit != hints->end(); hit++){
	(*hit)->discard = discard;
    }
}

/*
 * HintGroup::addIncompGroup
 */
void HintGroup::addIncompGroup(HintGroup *otherGroup){
    if (!otherGroup)
	return;
    if (otherGroup == this) {
	cerr << "addIncompGroup:Group incompatible to itself" << endl;
	return;
    }
    if (!incompGroups)
	incompGroups = new list<HintGroup*>;
    incompGroups->push_back(otherGroup);
}

/*
 * HintGroup::addStrongerGroup
 */
void HintGroup::addStrongerGroup(HintGroup *otherGroup){
    if (!otherGroup)
	return;
    if (otherGroup == this) {
	cerr << "addStrongerGroup:Group stronger than itself" << endl;
	return;
    }
    if (!strongerGroups)
	strongerGroups = new list<HintGroup*>;
    strongerGroups->push_back(otherGroup);
}

/*
 * HintGroup::operator<
 * for sorting incrementally according to begin positions
 */
bool operator<(const HintGroup& g1, const HintGroup& g2){    
    if (g1.begin < g2.begin)
	return true;
    if (g1.begin > g2.begin)
	return false;
    // both groups begin at the same position
    // go through the groups feature by feature and return true if the first differing feature is < in g1
    if (!g1.hints || !g2.hints)
	return false;
    list<Feature*>::iterator fit1 = g1.hints->begin(), fit2 = g2.hints->begin();
    while(fit1 != g1.hints->end() && fit2 != g2.hints->end()) {
	if ((*fit1)->start < (*fit2)->start || ((*fit1)->start == (*fit2)->start && (*fit1)->end < (*fit2)->end))
	    return true;
	if (!((*fit1)->start == (*fit2)->start && (*fit1)->end == (*fit2)->end))
	    return false;
	fit1++;
	fit2++;
    }
    return false;
}

/*
 * HintGroup::operator==
 */
bool operator==(const HintGroup& g1, const HintGroup& g2){
    if (g1.begin != g2.begin || g1.end != g2.end)
	return false;
    if (g1.hints == NULL || g2.hints == NULL || g1.hints->size() != g2.hints->size())
	return false;
    list<Feature*>::iterator fit1, fit2;
    for (fit1 = g1.hints->begin(), fit2 = g2.hints->begin(); fit1 != g1.hints->end() && fit2 != g2.hints->end(); fit1++, fit2++) {
	if (!((**fit1) == (**fit2)) || ((*fit1)->bonus != (*fit2)->bonus))
	    return false;
    }
    if (fit1 == g1.hints->end() && fit2 == g2.hints->end())
	return true;
    else 
	return false;
}



/*
 * HintGroup::print
 */
void HintGroup::print(ostream& out, bool withHints){
  if (hints == NULL) {
    out << "HintGroup " << name << ", " << begin << "-" << end << ", mult= " << copynumber << ", priority= " << priority << " no features" << endl;
  } else {
    out << "HintGroup " << name << ", " << begin << "-" << end << ", mult= " << copynumber << ", priority= " << priority << " " << hints->size() << " features" << endl;
    if (withHints)
      for (list<Feature*>::iterator fit = hints->begin(); fit!= hints->end(); fit++)
	out << "\t" << **fit << endl;
  }
}

