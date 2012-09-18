/**********************************************************************
 * file: compgenepred.cc licence: Artistic Licence, see file
 * LICENCE.TXT or
 * http://www.opensource.org/licenses/artistic-license.php descr.:
 * comparative gene prediction on multiple species authors: Mario
 * Stanke, Alexander Gebauer, Stefanie König
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 07.03.12| Mario Stanke  | creation of the file
 **********************************************************************/

#include "compgenepred.hh"
#include "orthograph.hh"
#include "mea.hh"
#include "genomicMSA.hh"
#include "geneMSA.hh"
#include "orthoexon.hh"
#include "namgene.hh"


CompGenePred::CompGenePred(){
    if (Constant::Constant::dbaccess.empty()) { // give priority to database in case both exist
        rsa = new MemSeqAccess();
    } else {
        rsa = new DbSeqAccess();
    }
}

void CompGenePred::start(){
    // read in alignment, determine orthologous sequence fragments


#ifdef DEBUG
    cout << "reading in the phylogenetic tree" << endl;
#endif
    PhyloTree *tree = new PhyloTree(Constant::treefile);  //has to be initialized before OrthoGraph

    OrthoGraph::tree = tree;
    GeneMSA::tree = tree;
    OrthoGraph::numSpecies = OrthoGraph::tree->species.size();


#ifdef DEBUG
    OrthoGraph::tree->printWithGraphviz("tree.dot");
    cout << "-------------------------------\nparameters phylogenetic model\n-------------------------------" << endl;
    cout << "rate exon loss:\t" << OrthoGraph::tree->evo.getMu() << endl;
    cout << "rate exon gain:\t" << OrthoGraph::tree->evo.getLambda() << endl;
    cout << "phylo factor:\t" << OrthoGraph::tree->evo.getPhyloFactor() <<  "\n-------------------------------" << endl;
#endif
  
    //initialize output files of initial gene prediction and optimized gene prediction
    vector<ofstream*> baseGenes = initOutputFiles(".base"); // equivalent to MEA prediction
    vector<int> base_geneid(OrthoGraph::numSpecies, 1); // gene numbering
    vector<ofstream*> initGenes = initOutputFiles(".init"); // score added to all orthologous exons and penalty added to all non orthologous exons, then global path search repeated
    vector<int> init_geneid(OrthoGraph::numSpecies, 1);
    vector<ofstream*> optGenes = initOutputFiles();  //optimized gene prediction by applying majority rule move
    vector<int> opt_geneid(OrthoGraph::numSpecies, 1);
    vector<ofstream*> sampledExons = initOutputFiles(".sampled_ECs");

    BaseCount::init();
    PP::initConstants();
    NAMGene namgene; // creates and initializes the states
    FeatureCollection extrinsicFeatures; // hints, empty for now, will later read in hints for sequence ranges from database
    SequenceFeatureCollection sfc(&extrinsicFeatures);
    StateModel::readAllParameters(); // read in the parameter files: species_{igenic,exon,intron,utr}_probs.pbl

    vector<string> speciesNames = OrthoGraph::tree->species;
    GenomicMSA msa;
    msa.readAlignment(speciesNames);
    msa.prepareExons();
    vector<int> offsets;

    // determine object that holds a sequence range for each species
    // loop over species
    GeneMSA::openOutputFiles();
    while (GeneMSA *geneRange = msa.getNextGene()) {
        OrthoGraph orthograph;
        for (int s = 0; s < speciesNames.size(); s++) {
            string seqID = geneRange->getSeqID(s);
            if (!seqID.empty()) {
                AnnoSequence *seqRange = rsa->getSeq(speciesNames[s], seqID, geneRange->getStart(s), geneRange->getEnd(s), geneRange->getStrand(s));
#ifdef DEBUG
                cout << "retrieving sequence:\t" << speciesNames[s] << ":" << geneRange->getSeqID(s) << "\t" << geneRange->getStart(s) << "-" << geneRange->getEnd(s) << "\t";
                if( geneRange->getStrand(s) == plusstrand )
                    cout << "+\t";
                else
                    cout << "-\t";
                cout << "(" <<geneRange->getEnd(s) - geneRange->getStart(s) + 1 << "bp)" << endl;
#endif
                if (seqRange==NULL) {
                    cerr << "random sequence access failed on " << speciesNames[s] << ", " << geneRange->getSeqID(s) << ", " << geneRange->getStart(s) << ", " <<  geneRange->getEnd(s) << ", " << endl;
                    break;
                } else {
                    namgene.getPrepareModels(seqRange->sequence, seqRange->length); // is needed for IntronModel::dssProb in GenomicMSA::createExonCands
                    if (geneRange->getStrand(s)==plusstrand) {
                        offsets.push_back(geneRange->getStart(s));
                    } else {
                        offsets.push_back(geneRange->getSeqIDLength(s) - (geneRange->getEnd(s)) - 1);
                    }
                    geneRange->createExonCands(seqRange->sequence); // ToDo:make this reasonable after experience with the data
                    list<ExonCandidate*> additionalExons = *(geneRange->getExonCands(s));

                    namgene.doViterbiPiecewise(sfc, seqRange, bothstrands); // sampling
                    list<Gene> *alltranscripts = namgene.getAllTranscripts();
                    if(alltranscripts){
                        cout << "building Graph for " << speciesNames[s] << endl;
                        /* build datastructure for graph representation
                         * @stlist : list of all sampled states
                         */
                        list<Status> stlist;
                        if(!alltranscripts->empty()){
                            buildStatusList(alltranscripts, false, stlist);
                        }
                        //build graph
                        orthograph.graphs[s] = new SpeciesGraph(&stlist, seqRange, additionalExons, speciesNames[s], geneRange->getStrand(s), sampledExons[s]);
                        orthograph.graphs[s]->buildGraph();

                        orthograph.ptrs_to_alltranscripts[s] = alltranscripts; //save pointers to transcripts and delete them after gene list is build
                    }
                }
            } else {
                offsets.push_back(0);
                geneRange->exoncands.push_back(NULL);
                geneRange->existingCandidates.push_back(NULL);
                cout<< speciesNames[s] << " doesn't exist in this part of the alignment."<< endl;
            }
        }
        geneRange->printGeneRanges();
        geneRange->printExonCands(offsets);
        geneRange->createOrthoExons(offsets);
        geneRange->printOrthoExons(offsets);
        orthograph.all_orthoex = geneRange->getOrthoExons();
	
	orthograph.outputGenes(baseGenes,base_geneid);
	//add score for selective pressure of orthoexons
	orthograph.addScoreSelectivePressure();
	//determine initial path
	orthograph.globalPathSearch();
	orthograph.outputGenes(initGenes,init_geneid);

        if(!orthograph.all_orthoex.empty()){
	    orthograph.pruningAlgor();
	    orthograph.printCache();
	    // Iterative optimization of labelings in graphs
	    orthograph.optimize();
	}

	// transfer max weight paths to genes + filter + ouput
	orthograph.outputGenes(optGenes,opt_geneid);

        offsets.clear();
        delete geneRange;
    }

    GeneMSA::closeOutputFiles();

    closeOutputFiles(baseGenes);
    closeOutputFiles(initGenes);
    closeOutputFiles(optGenes);

    // free memory space of tree
    delete OrthoGraph::tree;

}
