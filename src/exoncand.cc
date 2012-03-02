/**********************************************************************
 * file:    exoncand.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  Generation of exon candidates
 * author:  Alexander Gebauer
 *
 * date    |   author           |  changes
 * --------|--------------------|------------------------------------------
 * 03.11.11| Mario Stanke       | creation of the file
 * 06.12.11| Alexander Gebauer  | implementation of getExonCands
 **********************************************************************/

#include "exoncand.hh"
#include "intronmodel.hh"
#include "geneticcode.hh"
#include <fstream>

void getExonCands(const char* dna, float assqthresh, float dssqthresh){
  int n = strlen(dna);
  int max_exon_length = 12000;
  int frame;
  int single_fw=0;
  int single_rv=0;
  int initial_fw=0;
  int initial_rv=0;
  int internal_fw=0;
  int internal_rv=0;
  int terminal_fw=0;
  int terminal_rv=0;
  list<int> exonStart;
  list<int> exonRCStop;
  list<int> exonASS;
  list<int> exonRDSS;

  OpenReadingFrame orf(dna, max_exon_length, n);
  
  for (int i=0; i<=n-1; i++) {
	// Liste der Positionen aller Startcodons "atg"
	if (onStart(dna+i)) {
	    	exonStart.push_back(i+1);
	}
	// Liste der Positonen aller ASSs "ag"
    if (onASS(dna+i)) {
    	exonASS.push_back(i);
    }
    // Liste der Positionen aller reversen DSS "ac"
    if (onRDSS(dna+i)) {
    	exonRDSS.push_back(i);
    }
    // Liste der Positionen aller reversen komplementären Stopcodons "cta, tta, tca"
    if (onRCStopcodon(dna+i)) {
    	exonRCStop.push_back(i);
    }
}
list<int>::reverse_iterator ritStart = exonStart.rbegin();
list<int>::reverse_iterator ritStart_cur=ritStart;
list<int>::reverse_iterator ritASS = exonASS.rbegin();
list<int>::reverse_iterator ritASS_cur=ritASS;
list<int>::reverse_iterator ritRDSS = exonRDSS.rbegin();
list<int>::reverse_iterator ritRDSS_cur=ritRDSS;
list<int>::reverse_iterator ritRCStop = exonRCStop.rbegin();
list<int>::reverse_iterator ritRCStop_cur=ritRCStop;

for (int i=n-1; i>=2; i--) {
	// Berechnung der Single Genes auf dem Vorwärtsstrang mit mindestens Startcodon plus Codon
	if (ochre(dna+i) || opal(dna+i) ||amber(dna+i)) {
		ritStart=ritStart_cur;
		while ((i<*ritStart)&&(ritStart!=exonStart.rend())){
			ritStart++;
		}
		ritStart_cur=ritStart;
		int lmb = orf.leftmostExonBegin(0,i,true);
		while ((lmb<=*ritStart)&&(i-*ritStart<=max_exon_length)&&(ritStart!=exonStart.rend())) {
			if ((i-*ritStart>=4)&&((i-*ritStart+1)%3==0)) {
				cout << "exon candidate single gene forward  " << *ritStart << "   " << i <<endl;
		    	single_fw++;
		   	}
		   	ritStart++;
		};
	}

	// Berechnung der Initials auf dem Vorwärtsstrang mit mindestens Startcodon plus Base
	if (onDSS(dna+i)) {
		for (frame=0; frame<=2; frame++) {
			ritStart=ritStart_cur;
			while ((i<*ritStart)&&(ritStart!=exonStart.rend())){
				ritStart++;
			}
			ritStart_cur=ritStart;
			int lmb = orf.leftmostExonBegin(frame,i,true);
			while((lmb<=*ritStart)&&(i-*ritStart<=max_exon_length)&&(ritStart!=exonStart.rend())) {
				if ((i-*ritStart>=2)&&((i-*ritStart+1)%3==frame)) {
					cout << "exon candidate initial"<< frame <<" forward  " << *ritStart << "   " << i <<endl;
					initial_fw++;
				}
				ritStart++;
			};
		}

    	// Berechnung der Internals auf dem Vorwärtsstrang mit mindestens einem Codon
    	for (frame=0; frame<=2; frame++) {
    		ritASS=ritASS_cur;
    		while ((i<*ritASS)&&(ritASS!=exonASS.rend())){
    			ritASS++;
    		}
    		ritASS_cur=ritASS;
    		int lmb = orf.leftmostExonBegin(frame,i,true);
    		while((lmb<=*ritASS)&&(i-*ritASS<=max_exon_length)&&(ritASS!=exonASS.rend())) {
    			if (i-*ritASS>=5) {
    				cout << "exon candidate internal"<< frame <<" forward  " << *ritASS+3 << "   " << i <<endl;
    				internal_fw++;
    			}
    			ritASS++;
    		};
      	}
	}

	// Berechnung der Terminals auf dem Vorwärtstrang mit mindestens einer Base exklusive Stopcodon
	if (ochre(dna+i) || opal(dna+i) ||amber(dna+i)) {
		for (frame=0; frame<=2; frame++) {
	    	ritASS=ritASS_cur;
	    	while ((i<*ritASS)&&(ritASS!=exonASS.rend())){
	    		ritASS++;
	    	}
	    	ritASS_cur=ritASS;
	    	while ((i-*ritASS<=max_exon_length)&&(ritASS!=exonASS.rend())) {
	    		if ((i-*ritASS>=3)&&((i-*ritASS+1)%3==frame)&&(*ritASS>=orf.leftmostExonBegin(0,i,true))) {
	    	    	cout << "exon candidate terminal"<< frame <<" forward  " << *ritASS+3 << "   " << i <<endl;
	    	    	terminal_fw++;
	    		}
	    		ritASS++;
	    	};
		}
	}

	// Berechnung der Single Genes auf dem Rückwärtsstrang mit mindestens Startcodon+Codon exklusive Stopcodon
	if (onRStart(dna+i)) {
		ritRCStop=ritRCStop_cur;
		while ((i<*ritRCStop)&&(ritRCStop!=exonRCStop.rend())){
			ritRCStop++;
		}
		ritRCStop_cur=ritRCStop;
		while ((i-*ritRCStop<=max_exon_length)&&(ritRCStop!=exonRCStop.rend())) {
			if ((i-*ritRCStop>=6)&&((i-*ritRCStop)%3==0)) {
				cout << "exon candidate single gene reverse  " << *ritRCStop+4 << "   " << i+3 <<endl;
				single_rv++;
				break;
			} else {
				ritRCStop++;
			}
	   	};
	}

	// Berechnung der Initials auf dem Rückwärtstrang mit mindestens Startcodon plus Base
	if (onRStart(dna+i)) {
		for (frame=0; frame<=2; frame++) {
			ritRDSS=ritRDSS_cur;
			while ((i<*ritRDSS)&&(ritRDSS!=exonRDSS.rend())){
			    ritRDSS++;
			}
			ritRDSS_cur=ritRDSS;
			int lmb = orf.leftmostExonBegin(2,i,false);
			while((lmb<=*ritRDSS+2)&&(i-*ritRDSS<=max_exon_length)&&(ritRDSS!=exonRDSS.rend())) {
				if ((i-*ritRDSS>=2)&&((i+1-*ritRDSS)%3==frame)) {
					cout << "exon candidate Initial"<< frame <<" reverse  " << *ritRDSS+3 << "   " << i+3 <<endl;
					initial_rv++;
				}
				ritRDSS++;
			};
		}
	}

	// Berechnung der Internals auf dem Rückwärtstrang mit mindestens einem Codon
	if (onRASS(dna+i)) {
    	for (frame=0; frame<=2; frame++) {
    		ritRDSS=ritRDSS_cur;
          	while ((i<*ritRDSS)&&(ritRDSS!=exonRDSS.rend())){
          		  ritRDSS++;
          	}
          	ritRDSS_cur=ritRDSS;
          	int lmb = orf.leftmostExonBegin(frame,i,false);
          	while((lmb<=*ritRDSS)&&(i-*ritRDSS<=max_exon_length)&&(ritRDSS!=exonRDSS.rend())) {
    			if (i-*ritRDSS>=5) {
    			cout << "exon candidate internal"<< frame <<" reverse  " << *ritRDSS+3 << "   " << i <<endl;
    				internal_rv++;
    			}
    			ritRDSS++;
          	};
       	}
    }

	// Berechnung der Terminals auf dem Rückwärtsstrang mit mindestens einer Base exklusive Stopcodon
		if (onRASS(dna+i)) {
			for (frame=0; frame<=2; frame++) {
		    	ritRCStop=ritRCStop_cur;
		    	while ((i<*ritRCStop)&&(ritRCStop!=exonRCStop.rend())){
		    		ritRCStop++;
		    	}
		    	ritRCStop_cur=ritRCStop;
		    	while ((i-*ritRCStop<=max_exon_length)&&(ritRCStop!=exonRCStop.rend())) {
		    		if ((i-*ritRCStop>=4)&&((i-*ritRCStop)%3==frame)) {
		    	    	cout << "exon candidate terminal"<< frame <<" reverse  " << *ritRCStop+4 << "   " << i <<endl;
		    	    	terminal_rv++;
		    	    	break;
		    		} else {
		    			 ritRCStop++;
		    		}
		    	};
			}
		}
}

// int anzahl=single_fw+single_rv+initial_fw+initial_rv+internal_fw+internal_rv+terminal_fw+terminal_rv;

cout<< "Single Gene forward: "<<single_fw++<<endl;
cout<< "Single Gene reverse: "<<single_rv++<<endl;
cout <<"Initial forward: " <<initial_fw<<endl;
cout <<"Initial reverse: " <<initial_rv<<endl;
cout <<"Internal forward: " <<internal_fw<<endl;
cout <<"Internal reverse: " <<internal_rv<<endl;
cout <<"Terminal forward: " <<terminal_fw<<endl;
cout <<"Terminal reverse: " <<terminal_rv<<endl;
cout << "strand length: " <<n<<endl;
  
 // assqthresh, dssqthresh 
 Double assminprob = IntronModel::assBinProbs.getMinProb(assqthresh);
 Double dssminprob = IntronModel::dssBinProbs.getMinProb(dssqthresh);
 cout << "thresholds: dssminprob=" << dssminprob << " assminprob=" << assminprob << endl;
 // demo of splice site scores
 Double p;
 for (int i=0; i<200; i++){ // limit for testing and demo of code
   if (i + Constant::ass_whole_size() + Constant::ass_upwindow_size < n){ // see class IntronModel for boundaries
     p = IntronModel::aSSProb(i, true); // second argument: plusstrand?
     if (p>0)
       cout <<  ((p >= assminprob)? "strong" : "weak") << " intron may end at pos " 
	    << i+Constant::ass_upwindow_size+Constant::ass_start+ASS_MIDDLE-1
	    << " (0-based) on plusstrand, " << endl;
   }
   if (i + Constant::dss_whole_size() < n){ // see class IntronModel for boundaries
     p = IntronModel::dSSProb(i, true);
     if (p>0)
       cout << ((p >= dssminprob)? "strong" : "weak") << " intron may start at pos " 
	    << i+Constant::dss_start  << " (0-based) on plusstrand, "  << endl;
   }
 }
}
