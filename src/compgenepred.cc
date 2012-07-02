/**********************************************************************
 * file:    compgenepred.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  comparative gene prediction on multiple species
 * authors: Mario Stanke, Alexander Gebauer, Stefanie König
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
    if (!Constant::speciesfilenames.empty()) {
        rsa = new MemSeqAccess();
    } else {
        rsa = new DbSeqAccess();
    }
}

void CompGenePred::start(){
    // read in alignment, determine orthologous sequence fragments
    OrthoGraph::tree = new PhyloTree(Constant::treefile);  //has to be initialized before OrthoGraph

    OrthoGraph::numSpecies = OrthoGraph::tree->species.size();

    //OrthoGraph orthograph;
    //OrthoGraph::initOutputFiles();

    NAMGene namgene; // creates and initializes the states
    FeatureCollection extrinsicFeatures; // hints, empty for now, will later read in hints for sequence ranges from database
    SequenceFeatureCollection sfc(&extrinsicFeatures);
    StateModel::readAllParameters(); // read in the parameter files: species_{igenic,exon,intron,utr}_probs.pbl

    GenomicMSA msa;
    msa.readAlignment();
    msa.prepareExons();
    vector<string> speciesname = OrthoGraph::tree->species;
    GeneMSA *geneRange;

    /*//temp:
		orthograph.all_orthoex = readOrthoExons(Constant::orthoexons);
		vector<string> chrName;
		chrName.push_back("chr21");
		chrName.push_back("chr17");
		chrName.push_back("chr31");
		int start[] = {44836600, 31983200, 39789600} ;
		int end[] = {44846200, 31992100, 39798300};
		Strand strand = minusstrand;*/

    // determine object that holds a sequence range for each species
    // loop over species
    geneRange->openOutputFiles();
    while (geneRange = msa.getNextGene()) {
        for (int s = 0; s < speciesname.size(); s++) {
            if (!geneRange->getChr(s).empty()) {
                //AnnoSequence *seqRange = rsa->getSeq(speciesname[s], chrName[s], start[s], end[s]);
                //orthograph.orthoSeqRanges[s] = seqRange;
                AnnoSequence *seqRange = rsa->getSeq(speciesname[s], geneRange->getChr(s), geneRange->getStart(s), geneRange->getEnd(s)/*, geneRange->getStrand(s)*/);
                // maybe do not need getStrand because nameGene generates the reverse complement of an already reversed complement
                if (seqRange==NULL) {
                    cerr << "random sequence access failed on " << speciesname[s] << ", " << geneRange->getChr(s) << ", " << geneRange->getStart(s) << ", " <<  geneRange->getEnd(s) << ", " << endl;
                    break;
                } else {
                    namgene.getPrepareModels(seqRange->sequence, seqRange->length); // is needed for IntronModel::dssProb in GenomicMSA::createExonCands
                    geneRange->createExonCands(seqRange, 0.05, 0.23, 0.23); // ToDo:make this reasonable after experience with the data

                    /*
                     * build list of additional exoncandidates, which are inserted in the graph
                     */
                    /*list<ExonCandidate*> additionalExons;
					for(list<OrthoExon>::iterator it = orthograph.all_orthoex.begin(); it !=  orthograph.all_orthoex.end(); it++){
						if(it->orthoex[s] != NULL){
							it->orthoex[s]->begin -= orthograph.orthoSeqRanges[s]->offset;
							it->orthoex[s]->end -= orthograph.orthoSeqRanges[s]->offset;
							additionalExons.push_back(it->orthoex[s]);
						}
					}*/

                    //namgene.doViterbiPiecewise(sfc, orthograph.orthoSeqRanges[s], strand); // builds graph for each species

                    /*list<Gene> *alltranscripts = namgene.getAllTranscripts();
					if(alltranscripts){
						cout << "building Graph for " << speciesname[s] << endl;
						if(!alltranscripts->empty()){*
							/
							/* build datastructure for graph representation
                     * @stlist : list of all sampled states
                     */
                    /*list<Status> stlist;
							buildStatusList(alltranscripts, false, stlist);
							//build graph
							SpeciesGraph *singleGraph = new SpeciesGraph(&stlist, orthograph.orthoSeqRanges[s]->length, additionalExons, speciesname[s]);
							singleGraph->buildGraph();
							//find correct position in vector and add graph for species to OrthoGraph
							size_t pos = OrthoGraph::tree->getVectorPositionSpecies(speciesname[s]);
							if (pos < OrthoGraph::numSpecies){
								orthograph.graphs[pos] = singleGraph;
							} else {
								cerr << "species names in Orthograph and OrthoExon don't match" << endl;
							}
							orthograph.storePtrsToAlltranscripts(alltranscripts); //save pointers to transcripts and delete them after gene list is build
						}*/
                }
            } else {
                geneRange->exoncands.push_back(NULL);
                geneRange->existingCandidates.push_back(NULL);
                cout<< speciesname[s] << " doesn't exist in this part of the alignment."<< endl;
            }
        }
        geneRange->printExonCands();
        geneRange->createOrthoExons();
        geneRange->printOrthoExons();
        //list<OrthoExon> orthoExons = geneRange->getOrthoExons();

        //orthograph.pruningAlgor();

        // iterative optimization
        //orthograph.optimize();

        // transfer longest paths to genes + filter + ouput
        //orthograph.outputGenes(minusstrand);
    }
    geneRange->closeOutputFiles();
    //OrthoGraph::closeOutputFiles();

    // free memory space of tree
    delete OrthoGraph::tree;
}
