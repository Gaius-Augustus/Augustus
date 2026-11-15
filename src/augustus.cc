/**********************************************************************
 * file:    augustus.cc
 * license: Artistic License, see file LICENSE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 *
 * authors: Mario Stanke, mario.stanke@uni-greifswald.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 23.01.07| Mario Stanke  | enabled redirection of cout to outfile and cerr to errfile
 **********************************************************************/

// project includes
#include "types.hh"
#include "gene.hh"
#include "genbank.hh"
#include "namgene.hh"
#include "evaluation.hh"
#include "statemodel.hh"
#ifdef COMPGENEPRED
#include "compgenepred.hh"
#endif
// standard C/C++ includes
#include <fstream>
#include <sys/stat.h>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <set>
#include <vector>
#ifndef _WIN32
#include <cerrno>
#include <csignal>
#include <cstring>
#include <unistd.h>
#include <sys/wait.h>
#endif

int verbosity=1;
bool mea_prediction = false;
Boolean checkExAcc = false;
Boolean noprediction = false;
ofstream outputfile, errorfile;
string outputfilename, errorfilename;


/*
 * evaluateOnTestSet
 *
 * takes a list of annotated sequences, predicts the genes using the extrinsic information and 
 * evaluates the accuracy of the prediction
 * gene: a list of genes (annotated sequences)
 * namgene: object for the prediction algorithms
 * extrinsicFeatures: extrinsic info, is NULL if there is none
 * props: contains the parameters
 */ 
void evaluateOnTestSet(AnnoSequence *annoseq, NAMGene &namgene, FeatureCollection &extrinsicFeatures, 
		       Strand strand);

/*
 * predictOnInputSequences
 * 
 * takes a list of raw sequences and predicts the genes using the extrinsic information
 * seq: a list of sequences
 * namgene: object for the prediction algorithms
 * extrinsicFeatures: extrinsic info, is NULL if there is none
 * props: contains the parameters
 */
void predictOnInputSequences(AnnoSequence *seq, NAMGene &namgene, FeatureCollection &extrinsicFeatures, 
			     Strand strand);


/*
 * Set the global variables according to command line options or configuration file
 */

void setParameters();

/*
 * checkExtrinsicAccuracy
 *
 * takes a list of annotated sequences and determines the bonus malus of the extrinsic information
 * evaluates the bonus of the hints without cleaning for redundancies
 * the malus is not correctly computed because 
 * This is a first step to decide of which source the hints are of higher quality.
 * namgene: object for the prediction algorithms
 * extrinsicFeatures: extrinsic info
 * props: contains the parameters
 */ 
void checkExtrinsicAccuracy(AnnoSequence *annoseq, NAMGene &namgene, FeatureCollection &extrinsicFeatures);


/*
 * cutRelevantPiece
 * If predictionStart and predictionEnd are set this function cuts
 * out the piece from predictionStart to predictionEnd, and stores the offset predictionStart in the returned
 * AnnoSequence. The memory of the parameter AnnoSequence is deleted in this case, so the whole chromosome doesn't sit
 * in memory when we actually need only a small part.
 */
void cutRelevantPiece(AnnoSequence *annoseq);


/*
 * main
 */
int main( int argc, char* argv[] ){
    string     configfile;
    string     commandline;
    Strand     strand = bothstrands; // default
    int        errorcode = 0;

    LLDouble::setOutputPrecision(3);

    for (int i=0; i<argc; i++){
	commandline += argv[i];
	if (i<argc-1) 
	    commandline += " ";
    }
    
    // determination of the configuration file
    if( argc <= 1 ){
	cout << GREETING << endl << endl;
	cout << HELP_USAGE << endl;
	return 0;
    }
    try{
	Properties::init( argc, argv );
	Constant::init();
	Gene::init();
	GeneticCode::init();
	setParameters(); // NOTE: need Constant and GeneticCode to be initialized first
	StateModel::init();   // set global parameters of state models	  

	if(Properties::hasProperty("trainFeatureFile"))
            Constant::MultSpeciesMode = true;

	if (Constant::MultSpeciesMode){
#ifdef COMPGENEPRED
	    clock_t start;
	    start = clock();
	    CompGenePred cgp;
	    cgp.start();
	    cout << "# total time: " << (double) (clock()-start) / CLOCKS_PER_SEC << "s" << endl;
#else
	    throw ProjectError("Comparative gene prediction not possible with this compiled version. Please recompile with flag COMPGENEPRED = true set in common.mk.");
#endif	  
	} else { // single species mode, default
            // check query filename
            if (!Properties::hasProperty(INPUTFILE_KEY)){
                throw ProjectError("No query file specified. Type \"augustus\" for help.");
            }
            string filename =  Properties::getProperty(INPUTFILE_KEY);
            GBProcessor gbank(filename);
            if (Gene::gff3)
                cout << "##gff-version 3" << endl;
	  
            cout << PREAMBLE << endl;
	  
            /*
             * check for extrinsic information and initialize when existent
             */
            FeatureCollection extrinsicFeatures;
            const char *extrinsicfilename;
            try {
                extrinsicfilename =  Properties::getProperty("hintsfile");
            } catch (...){
                extrinsicfilename = NULL;
                if (verbosity) 
                    cout << "# No extrinsic information on sequences given." << endl;
            }
            if (extrinsicfilename || Constant::softmasking) {
                extrinsicFeatures.readExtrinsicCFGFile();
                if (extrinsicfilename) {
                    if (verbosity)
                        cout << "# reading in the file " << extrinsicfilename << " ..." << endl;
                    extrinsicFeatures.readGFFFile(extrinsicfilename);
                    if (verbosity) 
                        cout << "# Have extrinsic information about " << extrinsicFeatures.getNumSeqsWithInfo()
                             << " sequences (in the specified range). " << endl;
                }
            }
	  
            if (verbosity > 1) 
                cout << "# Initializing the parameters using config directory " << Properties::getProperty(CFGPATH_KEY) << " ..." << endl;
            BaseCount::init();
            PP::initConstants();
            NAMGene namgene; // creates and initializes the states
            StateModel::readAllParameters(); // read in the parameter files: species_{igenic,exon,intron,utr}_probs.pbl
            try{
                string strandstr = Properties::getProperty("strand");
                if (strandstr == "forward" || strandstr == "Forward" || strandstr == "plus" || strandstr == "Plus" 
                    || strandstr == "+" || strandstr == "Watson" || strandstr == "watson" || strandstr == "w" )
                    strand = plusstrand;
                else if (strandstr == "backward" || strandstr == "Backward" || strandstr == "minus" 
                         || strandstr == "Minus" || strandstr == "-" || strandstr == "Crick" || strandstr == "crick" 
                         || strandstr == "c" || strandstr == "reverse" || strandstr == "Reverse")
                    strand = minusstrand;
                else if (strandstr == "both")
                    strand = bothstrands;
                else if (!(strandstr == ""))
                    cerr << "# Unknown option for strand: " << strandstr << endl;
            } catch (...){} // take default strand
	
            if (mea_prediction)
                cout <<"# Using MEA approach (Maximizing expected accuracy)."<<endl;

            if (gbank.fileType() == fasta) {
                /*
                 * Just predict the genes for every sequence in the file.
                 */
                if (verbosity>2) {
                    if (filename == "-")
                        cout << "# Reading sequences from standard input. Assuming fasta format." << endl;
                    else
                        cout << "# Looks like " << filename << " is in fasta format." << endl;
                }
                AnnoSequence *testsequence = gbank.getSequenceList();
                warnAllLowerCase(testsequence);
                cutRelevantPiece(testsequence);
                predictOnInputSequences(testsequence, namgene, extrinsicFeatures, strand);
                AnnoSequence::deleteSequence(testsequence);
            } else if (gbank.fileType() == genbank) {
                /*
                 * Sequences were already annotated. Predict and also check the accuracy.
                 */
                if (verbosity>2)
                    cout << "# Looks like " << filename << " is in genbank format. " 
                         << "Augustus uses the annotation for evaluation of accuracy." << endl;
                AnnoSequence *annoseq = gbank.getAnnoSequenceList();
                if (!Constant::softmasking_explicitly_requested)
                    Constant::softmasking = false; // default is false only for .gb files
                warnAllLowerCase(annoseq);
                cutRelevantPiece(annoseq);
                if (!checkExAcc)
                    evaluateOnTestSet(annoseq, namgene, extrinsicFeatures, strand);
                else { // do not predict just check the accuracy of the extrinsic information
                    // without deleting for redundancies
                    checkExtrinsicAccuracy(annoseq, namgene, extrinsicFeatures);
                }
                AnnoSequence::deleteSequence(annoseq);
            } else {
                throw ProjectError("File format of " + filename + " not recognized.");
            }
	} // single species mode
	
	//	if (verbosity>2)
	cout << "# command line:" << endl << "# " << commandline << endl;

    } catch( ProjectError& err ){
        cerr << "\n" <<  argv[0] << ": ERROR\n\t" << err.getMessage( ) << "\n\n";
        errorcode=1;
    } catch ( HelpException &help ) {
    	cerr << help.message << endl;
    }
    if (outputfile.is_open())
	outputfile.close();
    if (errorfile.is_open())
	errorfile.close();
    return errorcode;
}

/*
 * evaluateOnTestSet
 */

void evaluateOnTestSet(AnnoSequence *annoseq, NAMGene &namgene, FeatureCollection &extrinsicFeatures, 
		       Strand strand){
    int dnaproben = 0;
    Double quotient; 
    Evaluation eval;
    Transcript *genes;

    // TEMP Zeitmessung
    double total =0;

    /*
     * How good was the extrinsic information.
     */
    if (extrinsicFeatures.getNumSeqsWithInfo()){
	int numC = extrinsicFeatures.getNumCommonSeqs(annoseq);
	if (numC == 0) 
	    cout << "# WARNING: extrinsic information given but not on any of the sequences in the input set!" << endl
		 << "# Maybe different sequence names?" << endl;
	if (verbosity > 3){
	    cout << "# We have hints for " << numC << " of the sequences in the input set." << endl;
	    cout << "# Reliability of extrinsic information: " << endl;
	    extrinsicFeatures.printAccuracyForSequenceSet(annoseq, true);
	}
    }
    // to output the annotated amino acid sequence:
    if (noprediction){
        namgene.setPathAndProb(annoseq, extrinsicFeatures);
	while( annoseq ){
	    //cout << ">" << annoseq->seqname << endl;
	    printGeneList(annoseq->anno->genes, annoseq, Constant::codSeqOutput, Constant::proteinOutput);
	    cout << "# joint probability of gene structure and sequence = " << annoseq->anno->emiProb << endl;
	    annoseq = annoseq->next;
	}
    }

    while( annoseq ){
	dnaproben++;
	if (verbosity) {
	    cout << "#\n# ----- sequence number " << dnaproben << " (length = "
		 << annoseq->length << ", name = "
		 << annoseq->seqname << ") -----" << endl << "#" << endl;
	    cout << "\n# annotation: " << endl; 
	    annoseq->printGFF();
	}
	clock_t anfang, ende;
	anfang = clock();
       
	/*
	 * check for extrinsic information about this sequence
	 */
	SequenceFeatureCollection& sfc = extrinsicFeatures.getSequenceFeatureCollection(annoseq->seqname);
	sfc.prepare(annoseq, Constant::printHints && !(Gene::gff3));
	
	bool singlestrand = false; // use not the shadow states
	try {
            singlestrand = Properties::getBoolProperty("singlestrand");
	} catch (...) {}

	cout << "# Predicted genes for sequence number " << dnaproben <<  " on ";
	if (strand==plusstrand)
	    cout << "forward strand" << endl;
	else if (strand==minusstrand)
	    cout << "reverse strand" << endl;
	else 
	    cout << "both strands" << endl;
	if (singlestrand)
	    cout << "# Overlapping genes on opposite strand were allowed." << endl;
        unsigned int num_pieces;
        genes = namgene.doViterbiPiecewise(sfc, annoseq, strand, num_pieces);
	ende = clock();
	total += (double) (ende-anfang) / CLOCKS_PER_SEC;
	//cout << "time " << (double) (ende-anfang) / CLOCKS_PER_SEC << ", seqlen=" << annoseq->length << endl;

	eval.addToEvaluation(genes, annoseq->anno->genes, strand, quotient);
	/*
	 * clean up
	 */
	Transcript::destroyGeneSequence(genes);
	if (annoseq->anno->path){
	    delete annoseq->anno->path;
	    annoseq->anno->path = NULL;           // so annoseq destructor doesn't crash
	}
	if (annoseq->anno->forwardPath){
	    delete annoseq->anno->forwardPath;
	    annoseq->anno->forwardPath = NULL;           // so annoseq destructor doesn't crash
	}	
	if (annoseq->anno->backwardPath){
	    delete annoseq->anno->backwardPath;
	    annoseq->anno->backwardPath = NULL;           // so annoseq destructor doesn't crash
	}	
	if (annoseq->anno->condensedPath){
	    delete annoseq->anno->condensedPath;
	    annoseq->anno->condensedPath = NULL; // so annoseq destructor doesn't crash
	}
	if (annoseq->anno->condensedForwardPath){
	    delete annoseq->anno->condensedForwardPath;
	    annoseq->anno->condensedForwardPath = NULL; // so annoseq destructor doesn't crash
	}
	if (annoseq->anno->condensedBackwardPath){
	    delete annoseq->anno->condensedBackwardPath;
	    annoseq->anno->condensedBackwardPath = NULL; // so annoseq destructor doesn't crash
	}
	annoseq = annoseq->next;
    }
    if (!noprediction) {
        eval.finishEvaluation();
        eval.printQuotients();
        eval.print();
    }
    cout << "# total time: " << total << endl;
}


/*
 * predictOneSequence
 * predict and print the genes of one sequence; throws ProjectError on failure
 * coutFormat: the format of cout before the first sequence
 */
static void predictOneSequence(AnnoSequence *curseq, int dnaproben, NAMGene &namgene,
			       FeatureCollection &extrinsicFeatures, Strand strand, const ios &coutFormat){
    // Reseed the random generator used for sampling, so that the result for this sequence does
    // not depend on the sequences predicted before it (nor on --cpus). srand(1) is the initial
    // state of rand(), so the output for the first or only sequence is unchanged.
    srand(1);
    // Likewise for the format of cout: the gene output sets e.g. setprecision(3) and keeps it,
    // which formerly changed the printed hints (--printHints) of all but the first sequence.
    cout.copyfmt(coutFormat);
    if (verbosity)
	cout << "#\n# ----- prediction on sequence number " << dnaproben << " (length = "
	     << strlen(curseq->sequence) << ", name = "
	     << curseq->seqname << ") -----" << endl << "#" << endl;
    /*
     * check for extrinsic information about this sequence
     */
    SequenceFeatureCollection& sfc = extrinsicFeatures.getSequenceFeatureCollection(curseq->seqname);
    sfc.prepare(curseq, Constant::printHints && !(Gene::gff3));
    bool singlestrand = false; // use not the shadow states
    try {
	singlestrand = Properties::getBoolProperty("singlestrand");
    } catch (...) {}

    cout << "# Predicted genes for sequence number " << dnaproben <<  " on ";
    if (strand==plusstrand)
	cout << "forward strand" << endl;
    else if (strand==minusstrand)
	cout << "reverse strand" << endl;
    else 
	cout << "both strands" << endl;
    if (singlestrand)
	cout << "# Overlapping genes on opposite strand are allowed." << endl;
    unsigned int num_pieces;
    Transcript *genes = namgene.doViterbiPiecewise(sfc, curseq, strand, num_pieces);

    try {
	if (Properties::getBoolProperty("emiprobs")){ // get emission probs (special request from Ingo Ebersberger)
	    try {
		if (num_pieces > 1){
		    int maxDNAPieceSize = Properties::getIntProperty( "maxDNAPieceSize" );
		    cerr << "Error: --emiprobs option works only for input sequences of length at most maxDNAPieceSize (currently set to " <<
			maxDNAPieceSize << "). Use --maxDNAPieceSize to readjust." << endl;
		} else {
		    Annotation *a = new Annotation(), *olda = curseq->anno;
		    a->genes = genes;
		    curseq->anno = a;
		    namgene.setPathAndProb(curseq, extrinsicFeatures);
		    cout << "# joint probability of gene structure and sequence in " <<
			Properties::getProperty(SPECIES_KEY) << " model: " << curseq->anno->emiProb << endl;
		    curseq->anno = olda;
		}
	    } catch (ProjectError &e){
		cerr << e.getMessage() << endl << "Error: Could not compute the emission probabilities (--emiprobs)" << endl;
	    }
	}
    } catch (...) {}

    Transcript::destroyGeneSequence(genes); // don't need them anymore after they are printed
}

#ifndef _WIN32
/*
 * shiftGeneId
 * If the identifier that starts at position p of s begins with a gene identifier <prefix><n>,
 * e.g. g12 or g12.t1, replace it by <prefix><n+offset> and return n, otherwise return 0.
 * The identifier of a protein match, pp.<transcript id>.<block name>, is shifted in the
 * transcript id only.
 */
static long shiftGeneId(string &s, size_t p, const string &prefix, long offset){
    if (s.compare(p, prefix.size(), prefix) != 0 && s.compare(p, 3, "pp.") == 0)
	p += 3;
    if (s.compare(p, prefix.size(), prefix) != 0)
	return 0;
    size_t d = p + prefix.size(), e = d;
    while (e < s.size() && isdigit((unsigned char) s[e]))
	e++;
    if (e == d || !(e == s.size() || s[e] == '.' || s[e] == '"' || s[e] == ';'))
	return 0;
    long id = atol(s.substr(d, e-d).c_str());
    s.replace(d, e-d, itoa(id + offset));
    return id;
}

/*
 * appendRenumbered
 * Copy the output of a child process to cout, shifting its gene numbers (which start at 1)
 * by offset. Only the gene identifiers in the gene part of the output (after the printed hints)
 * are changed: in the "# start gene" / "# end gene" lines and, in the GFF lines of the
 * predictions (source AUGUSTUS), the identifier in the last column of the gene and transcript
 * lines (GTF) and the values of ID, Parent, transcript_id and gene_id. Other text, e.g. hint
 * groups or the block names of protein profiles, is copied unchanged.
 * Returns the number of genes of the file.
 */
static long appendRenumbered(const string &filename, const char *seqname, long offset){
    string prefix = "g";
    if (Constant::uniqueGeneId) {
	char buf[38];
	sprintf(buf, "%.30s.g", seqname); // as in NAMGene::doViterbiPiecewise
	prefix = buf;
    }
    static const string keys[] = {"ID=", "Parent=", "transcript_id \"", "gene_id \""};
    ifstream in(filename.c_str());
    string line;
    long numGenes = 0;
    bool genepart = false;
    while (getline(in, line)) {
	if (!genepart) {
	    genepart = (line.compare(0, 37, "# Predicted genes for sequence number") == 0);
	} else if (line.compare(0, 13, "# start gene ") == 0) {
	    numGenes = max(numGenes, shiftGeneId(line, 13, prefix, offset));
	} else if (line.compare(0, 11, "# end gene ") == 0) {
	    shiftGeneId(line, 11, prefix, offset);
	} else {
	    size_t p = line.find('\t');
	    if (p != string::npos && line.compare(p+1, 9, "AUGUSTUS\t") == 0) {
		int tabs = 1;
		while (tabs < 8 && (p = line.find('\t', p+1)) != string::npos)
		    tabs++;
		if (tabs == 8) {
		    size_t col9 = p+1;
		    shiftGeneId(line, col9, prefix, offset); // gene and transcript lines (GTF)
		    for (const string &key : keys) {
			for (size_t q = col9; (q = line.find(key, q)) != string::npos; q += key.size())
			    if (q == col9 || line[q-1] == ';' || line[q-1] == ' ') // start of an attribute
				shiftGeneId(line, q + key.size(), prefix, offset);
		    }
		}
	    }
	}
	cout << line << '\n';
    }
    return numGenes;
}

/*
 * predictInParallel
 * --cpus=n: predict up to n sequences at a time in forked child processes. The parameters are
 * shared copy-on-write. Each child writes its output to a temporary file; the parent writes the
 * outputs in input order and numbers the genes consecutively, as a run with --cpus=1 would.
 * Errors are handled as with --cpus=1: an error on the first sequence ends the run with an error,
 * an error on a later sequence is reported and the run continues. A worker that crashes or is
 * killed (e.g. by the OOM killer) ends the run with an error after the outputs before it.
 */
static void predictInParallel(vector<AnnoSequence*> &seqs, unsigned cpus, NAMGene &namgene,
			      FeatureCollection &extrinsicFeatures, Strand strand, const ios &coutFormat){
    const char *tmpenv = getenv("TMPDIR");
    string tmpl = string((tmpenv && *tmpenv) ? tmpenv : "/tmp") + "/augustus.XXXXXX";
    vector<char> tbuf(tmpl.begin(), tmpl.end());
    tbuf.push_back('\0');
    if (!mkdtemp(tbuf.data()))
	throw ProjectError("--cpus: could not create a temporary directory " + tmpl);
    string tmpdir(tbuf.data());

    size_t n = seqs.size();
    vector<string> outnames(n);
    // -1: not finished, 0: success, 1: error (ProjectError), 2: other error, crash or killed
    vector<int> status(n, -1);
    // start the longest sequences first for a better load balance; the output order is unchanged
    vector<size_t> order(n);
    for (size_t i = 0; i < n; i++)
	order[i] = i;
    stable_sort(order.begin(), order.end(), [&seqs](size_t a, size_t b) {
	return seqs[a]->length > seqs[b]->length;
    });

    cout.flush();
    cerr.flush();
    fflush(NULL); // nothing buffered may be inherited by the children
    map<pid_t, size_t> running;
    size_t next = 0, nextToWrite = 0;
    long geneOffset = 0;
    int successfull = 0;
    try {
	while (next < n || !running.empty()) {
	    while (next < n && running.size() < cpus) {
		size_t i = order[next];
		outnames[i] = tmpdir + "/seq" + itoa(i) + ".out";
		pid_t pid = fork();
		if (pid < 0)
		    throw ProjectError(string("--cpus: fork failed: ") + strerror(errno));
		if (pid == 0) { // child: predict sequence i into its own file
		    // No exception may leave this block, the child would continue in the code of the parent.
		    int code = 0;
		    ofstream out(outnames[i].c_str());
		    streambuf *oldbuf = cout.rdbuf(out.rdbuf());
		    try {
			predictOneSequence(seqs[i], i+1, namgene, extrinsicFeatures, strand, coutFormat);
		    } catch (ProjectError &err) {
			cerr << "\n augustus: ERROR\n\t" << err.getMessage( ) << "\n\n";
			code = 1;
		    } catch (std::exception &e) {
			cerr << "\n augustus: ERROR\n\t" << e.what() << "\n\n";
			code = 2;
		    } catch (...) {
			cerr << "\n augustus: ERROR\n\tunknown exception\n\n";
			code = 2;
		    }
		    cout.flush();
		    cout.rdbuf(oldbuf);
		    out.close();
		    cerr.flush();
		    _exit(out.fail() ? 2 : code); // no destructors, no atexit handlers of the parent
		}
		next++;
		running[pid] = i;
	    }
	    int st;
	    pid_t pid = wait(&st);
	    if (pid < 0) {
		if (errno == EINTR)
		    continue;
		throw ProjectError(string("--cpus: wait failed: ") + strerror(errno));
	    }
	    map<pid_t, size_t>::iterator it = running.find(pid);
	    if (it == running.end())
		continue;
	    size_t i = it->second;
	    running.erase(it);
	    if (WIFEXITED(st))
		status[i] = min(WEXITSTATUS(st), 2);
	    else {
		status[i] = 2;
		if (WIFSIGNALED(st))
		    cerr << "augustus: the prediction on sequence number " << i+1 << " (" << seqs[i]->seqname
			 << ") was terminated by signal " << WTERMSIG(st) << " (" << strsignal(WTERMSIG(st)) << ")." << endl;
	    }
	    // write all finished outputs that are next in input order
	    while (nextToWrite < n && status[nextToWrite] >= 0) {
		size_t k = nextToWrite++;
		geneOffset += appendRenumbered(outnames[k], seqs[k]->seqname, geneOffset);
		cout.flush();
		unlink(outnames[k].c_str());
		outnames[k] = "";
		if (status[k] == 0)
		    successfull++;
		else if (status[k] == 2 || successfull < 1) // as with --cpus=1
		    throw ProjectError("--cpus: the prediction on sequence number " + itoa(k+1) + " ("
				       + seqs[k]->seqname + ") failed.");
		// else the child has reported the error, continue as with --cpus=1
	    }
	}
    } catch (...) {
	// stop the running workers, wait for them and remove the temporary files
	for (map<pid_t, size_t>::iterator it = running.begin(); it != running.end(); ++it)
	    kill(it->first, SIGTERM);
	for (map<pid_t, size_t>::iterator it = running.begin(); it != running.end(); ++it)
	    while (waitpid(it->first, NULL, 0) < 0 && errno == EINTR);
	for (size_t i = 0; i < n; i++)
	    if (!outnames[i].empty())
		unlink(outnames[i].c_str());
	rmdir(tmpdir.c_str());
	throw;
    }
    rmdir(tmpdir.c_str());
}
#endif

/*
 * predictOnInputSequences
 */

void predictOnInputSequences(AnnoSequence *seq, NAMGene &namgene, FeatureCollection &extrinsicFeatures, 
			     Strand strand){
    int dnaproben = 0;
    int successfull = 0;
    AnnoSequence *curseq;

    int numC = extrinsicFeatures.getNumCommonSeqs(seq);
    if (numC == 0 && extrinsicFeatures.getNumSeqsWithInfo() > 0) 
	cout << "# WARNING: extrinsic information given but not on any of the sequences in the input set!" << endl
	     << "# Maybe different sequence names?" << endl;
    if (verbosity>0) {
	int anz = extrinsicFeatures.getNumSeqsWithInfo();
	cout << "# We have hints for " << anz << " sequence" << (anz!=1? "s ": " ") << "and for " 
	     << numC << " of the sequences in the input set." << endl;
    }

    ios coutFormat(NULL);
    coutFormat.copyfmt(cout);
    int cpus = 1;
    // only if set: the parameter list in a config directory of an older version lacks cpus, and
    // getIntProperty would print an error message
    if (Properties::hasProperty("cpus"))
	cpus = Properties::getIntProperty("cpus");
#ifndef _WIN32
    if (cpus > 1 && seq && seq->next) {
	// Sequences with the same name share the hints, which are changed when a sequence is
	// predicted (e.g. softmasking hints are added). Then the result for a sequence can depend on
	// the sequences before it, and only the prediction one after the other gives that result.
	set<string> names;
	for (AnnoSequence *s = seq; s && cpus > 1; s = s->next)
	    if (!names.insert(s->seqname).second) {
		cerr << "Warning: The input contains several sequences with the name " << s->seqname
		     << ". They are predicted one after the other, --cpus=" << cpus << " is ignored." << endl;
		cpus = 1;
	    }
    }
    if (cpus > 1 && seq && seq->next) {
	vector<AnnoSequence*> seqs;
	for (AnnoSequence *s = seq; s; s = s->next)
	    seqs.push_back(s);
	for (size_t i = 0; i < seqs.size(); i++)
	    seqs[i]->next = NULL; // each child works on one sequence
	try {
	    predictInParallel(seqs, (unsigned) cpus, namgene, extrinsicFeatures, strand, coutFormat);
	} catch (...) {
	    for (size_t i = 0; i+1 < seqs.size(); i++)
		seqs[i]->next = seqs[i+1];
	    throw;
	}
	for (size_t i = 0; i+1 < seqs.size(); i++)
	    seqs[i]->next = seqs[i+1]; // restore the list for the caller
	return;
    }
#endif

    while( seq ){
	dnaproben++;
	curseq = seq;
	seq = seq->next;  
	curseq->next = NULL;
	try {
	    predictOneSequence(curseq, dnaproben, namgene, extrinsicFeatures, strand, coutFormat);
	    successfull++;
	} catch (ProjectError& err ){
	    if (successfull < 1)
		throw err;
	    else 
		cerr << "\n augustus: ERROR\n\t" << err.getMessage( ) << "\n\n";
	}
	curseq->next = seq;
    }
}


/*
 * setParameters
 */

void setParameters(){
    try { 
	outputfilename = Properties::getProperty("outfile");
    } catch (...){}
    try { 
	errorfilename = Properties::getProperty("errfile");
    } catch (...){}
    try {
	verbosity = Properties::getIntProperty("/augustus/verbosity");
    } catch (...) {} 
    try {
	checkExAcc = Properties::getBoolProperty("checkExAcc");
    } catch (ProjectError &e) {
	cerr << e.getMessage();
    }
    try {
        noprediction = Properties::getBoolProperty("noprediction");
    } catch (...) {}   
    try {
        mea_prediction = Properties::getBoolProperty("mea");
    } catch (...) {}

    // directory for output files of comparative gene prediction
    // create this before outfile because outfile could be in this directory
    try {
        string outdir = Properties::getProperty("/CompPred/outdir");
	if (outdir != ""){
	    outdir = expandHome(outdir);         // replace "~" by "$HOME"
	    if(outdir[outdir.size()-1] != '/')
		outdir += '/';                   // append slash if necessary 
	    
	    Properties::addProperty("/CompPred/outdir", outdir);
	    // the directory actually exists? 
	    struct stat buffer;
	    if( stat(outdir.c_str(), &buffer) == -1 || !S_ISDIR(buffer.st_mode)){
		int err = mkdir(outdir.c_str(),0755); // create it, if not
		if (err != 0)
		    throw ProjectError("Could not create output directory " + outdir);
	    }
	}
    } catch (...) { }

    if (outputfilename != "") {
	outputfile.open(outputfilename.c_str());
	if (outputfile){
	    streambuf * fobuf = outputfile.rdbuf();
	    // the original stdout buffer could be stored at this point and restored later
	    cout.rdbuf(fobuf);
	} else {
	    cerr << "Could not open output file " << outputfilename << ". Will use stdout instead." << endl;
	}
    }
    if (errorfilename != "") {
	errorfile.open(errorfilename.c_str());
	if (errorfile){
	    cerr.rdbuf(errorfile.rdbuf());
	} else {
	    cerr << "Could not open error file " << errorfilename << ". Will use stderr instead." << endl;
	}
    }
    if (Gene::gff3 && Gene::stopCodonExcludedFromCDS)
	cerr << "Warning: The gff3 standard requires that the stop codon is included in the CDS. "
	     << "Unless this is your intention, set stopCodonExcludedFromCDS to false in your "
	     << "species' configuration file or on the command line." << endl;
     
    if (Properties::hasProperty("translation_table")){
        GeneticCode::chooseTranslationTable(Properties::getIntProperty("translation_table"));
    }
}


void checkExtrinsicAccuracy(AnnoSequence *annoseq, NAMGene &namgene, FeatureCollection &extrinsicFeatures){
    /*
     * How good was the extrinsic information.
     */
    if (extrinsicFeatures.getNumSeqsWithInfo()){
	if (verbosity > 1){
	    cout << "Reliability of extrinsic information: " << endl;
	    extrinsicFeatures.printAccuracyForSequenceSet(annoseq, false);
	}
    }
}


/*
 * cutRelevantPiece
 * If predictionStart and predictionEnd are set this function cuts
 * out the piece from predictionStart to predictionEnd, and stores the offset predictionStart in the returned
 * AnnoSequence. The memory of the complete sequence of AnnoSequence is deleted in this case, so the whole chromosome doesn't sit
 * in memory when we actually need only a small part.
 */
void cutRelevantPiece(AnnoSequence *annoseq){
    long predictionStart, predictionEnd;
    long seqlen = annoseq->length;

    try {
	predictionStart = Properties::getIntProperty( "predictionStart" ) - 1;
    } catch (...) {
	predictionStart = 0;
    }
    try {
	predictionEnd = Properties::getIntProperty( "predictionEnd" ) - 1;
    } catch (...) {
	predictionEnd = seqlen - 1;
    }

    if ((predictionStart != 0 || predictionEnd != seqlen - 1) 
	&& !(predictionEnd < 0 && predictionStart < 0)) {
	if (predictionStart < 0)
	    predictionStart = 0;
	if (predictionEnd > seqlen - 1)
	    predictionEnd = seqlen -1;
	if (predictionStart >= seqlen) 
	    throw ProjectError("predictionStart (" + itoa(predictionStart + 1) 
			       + ") is larger than sequence length (" + itoa(seqlen) 
			       + "). No predictions made.");
	if (predictionEnd < predictionStart)
	    throw ProjectError("predictionEnd (" + itoa(predictionEnd + 1) + 
			       ") is smaller than predictionStart (" + itoa(predictionStart + 1)
			       + "). No predictions made!");
	if (annoseq->next) {
	    cerr << "Warning: predictionStart or predictionEnd set but input consists of more than one sequence." << endl
		 << "Prediction will be made only on first sequence." << endl;
	    annoseq->next = 0;
	}

	annoseq->length = predictionEnd - predictionStart + 1;
	char *seq = newstrcpy(annoseq->sequence + predictionStart, annoseq->length);
	delete [] annoseq->sequence;
	annoseq->sequence = seq;
	annoseq->offset = predictionStart;
    } else if (predictionStart < 0 && predictionEnd < 0 && predictionEnd == predictionStart){ 
        // this is intended for applications in which a sequence fragment is cut out of a chromosome for input to augustus
        // and the hints have still original coordinates
        // => use the whole input sequence and just left-shift hints and right-shift output coordinates
        annoseq->offset = -predictionStart - 1;
    }
}
