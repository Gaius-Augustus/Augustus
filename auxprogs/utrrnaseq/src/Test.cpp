/**
 * \file Test.cpp
 */

#include "Genomic_Data.hpp"
#include "Supporting_Methods.hpp"
#include "Error.hpp"
#include "Global.hpp"
#include "Compute_UTRs.hpp"
#include "UTRs.hpp"
#include "Test.hpp"

#include <string>

#include <iostream>
#include <fstream>

#include <boost/assign/std/vector.hpp>
#include <boost/foreach.hpp>

using namespace std;
using namespace boost::assign;


#ifdef TEST_MODE
	#define BOOST_TEST_DYN_LINK
	#define BOOST_TEST_MODULE general

	#include <boost/test/included/unit_test.hpp>
	#include <boost/test/floating_point_comparison.hpp>

	struct Constructed_Test_Sets {
		Constructed_Test_Sets() {
			unsigned readLength = 5;
			unsigned smoothingWindowSize = 15;
			unsigned limit = 100;
			unsigned dropWindowSize = 5;
			unsigned minLength = 1;
			unsigned minCov = 1;
			double pWin = 0.9;
			double pInt = 0.9;
			double pMult = 0.1;
			bool ZeroCov = false;

			UTRs::clear();

			UTRs::parameters().window_size = smoothingWindowSize;
			UTRs::parameters().read_length = readLength;
			UTRs::parameters().limit = limit;
			UTRs::parameters().drop_window_size = dropWindowSize;
			UTRs::parameters().min_length = minLength;
			UTRs::parameters().min_average_cov = minCov;
			UTRs::parameters().p_win = pWin;
			UTRs::parameters().p_int = pInt;
			UTRs::parameters().p_mult = pMult;
			UTRs::parameters().zero_cov = ZeroCov;
		}
	};

	BOOST_FIXTURE_TEST_SUITE( Constructed_Test_Suites, Constructed_Test_Sets )

		#ifdef TEST_CASES

			BOOST_AUTO_TEST_CASE( Standard_Test_Case )
			{
				//2 genes, 5' and 3' utr for each gene, and intron in one 5' utr and one 3' utr
				//no repeats
				//no overlapping introns

				string scaffold_fname = "for-test-cases/standard-test-set/test.fa";
				string coding_region_fname = "for-test-cases/standard-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/standard-test-set/test-introns.hints";
				string repeat_fname = "for-test-cases/standard-test-set/test-repeat.hints";
				string wiggle_fname = "for-test-cases/standard-test-set/test.wig";
				string output_fname = "standard-test-case-utr.gff";
				string comparison_output_fname = "for-test-cases/comparison-files/comparison-standard-test-case.gff";
				vector< pair <string, string> > splice_sites;

				ERR_STRM.str(""); //clear stringstream

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);

				UTRs::parameters().output_fname = output_fname;
				compute_UTRs(wiggle_fname);
				UTRs::output();

				string str = ERR_STRM.str();
				BOOST_CHECK(str.empty());
				BOOST_CHECK(UTRs::size() == 4);
				BOOST_CHECK(files_identical( output_fname, comparison_output_fname ));
			}

			BOOST_AUTO_TEST_CASE( Repeat_Test_Case )
			{
				//2 genes, 5' and 3' utr for each gene, and intron in one 5' utr and one 3' utr
				//all utrs are discarded due to a repeat in the utr
				//no introns
				//no overlapping introns

				string scaffold_fname = "for-test-cases/repeat-test-set/test.fa";
				string coding_region_fname = "for-test-cases/repeat-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/repeat-test-set/test-introns.hints";
				string repeat_fname = "for-test-cases/repeat-test-set/test-repeat.hints";
				string wiggle_fname = "for-test-cases/repeat-test-set/test.wig";
				string output_fname = "repeat-test-case-utr.gff";
				vector< pair <string, string> > splice_sites;

				ERR_STRM.str("");

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);

				UTRs::parameters().output_fname = output_fname;
				compute_UTRs(wiggle_fname);
				UTRs::output();

				string str = ERR_STRM.str();
				BOOST_CHECK(str.empty());
				BOOST_CHECK(UTRs::size() == 0);
			}

			BOOST_AUTO_TEST_CASE( Repeat_Corner_Test_Case )
			{
				//2 genes, 5' and 3' utr for each gene, and intron in one 5' utr and one 3' utr
				//only the first and last position in the scaffold are repeats
				//no introns
				//no overlapping introns

				string scaffold_fname = "for-test-cases/repeat-border-test-set/test.fa";
				string coding_region_fname = "for-test-cases/repeat-border-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/repeat-border-test-set/test-introns.hints";
				string repeat_fname = "for-test-cases/repeat-border-test-set/test-repeat.hints";
				string wiggle_fname = "for-test-cases/repeat-border-test-set/test.wig";
				string output_fname = "repeat-border-test-case-utr.gff";
				vector< pair <string, string> > splice_sites;

				ERR_STRM.str("");
				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);

				UTRs::parameters().output_fname = output_fname;
				compute_UTRs(wiggle_fname);
				UTRs::output();

				string str = ERR_STRM.str();
				BOOST_CHECK(str.empty());
				BOOST_CHECK(UTRs::size() > 0);
			}

			BOOST_AUTO_TEST_CASE( Softmasking_Test_Case )
			{
				// Initialize data once with files from "Repeat_Test_Case" and once with empty repeat file but
				// softmasked genome sequences. Compare then whether the stored repeats are the same.

				string scaffold_fname = "for-test-cases/repeat-test-set/test.fa";
				string coding_region_fname = "for-test-cases/repeat-test-set/test-start-stop.gtf";
				string repeat_fname = "for-test-cases/repeat-test-set/test-repeat.hints";
				string intron_hints_fname = "for-test-cases/repeat-test-set/test-introns.hints";
				string wiggle_fname = "for-test-cases/repeat-test-set/test.wig";
				string output_fname = "repeat-test-case-utr.gff";
				vector< pair <string, string> > splice_sites;

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);
				vector<Genomic_Data::Repeat> reps_by_file = Genomic_Data::get_all_genes()[0].repeats;

				scaffold_fname = "for-test-cases/softmasking-test-set/test.fa";
				coding_region_fname = "for-test-cases/softmasking-test-set/test-start-stop.gtf";
				intron_hints_fname = "for-test-cases/softmasking-test-set/test-introns.hints";
				wiggle_fname = "for-test-cases/softmasking-test-set/test.wig";
				output_fname = "softmasking-test-case-utr.gff";

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, "", splice_sites, false);
				vector<Genomic_Data::Repeat> reps_by_mask = Genomic_Data::get_all_genes()[0].repeats;

				cout << "reps_by_file" << endl;
				BOOST_FOREACH(Genomic_Data::Repeat r, reps_by_file)
					cout << r << endl;
				cout << "reps_by_mask" << endl;
				BOOST_FOREACH(Genomic_Data::Repeat r, reps_by_mask)
					cout << r << endl;

				BOOST_CHECK(reps_by_file == reps_by_mask);
			}

			BOOST_AUTO_TEST_CASE( Softmasking_Corner_Test_Case )
			{
				// Initialize data once with files from "Repeat_Border_Test_Case" and once with empty repeat file but
				// softmasked genome sequences. Compare then whether the stored repeats are the same.

				string scaffold_fname = "for-test-cases/repeat-border-test-set/test.fa";
				string coding_region_fname = "for-test-cases/repeat-border-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/repeat-border-test-set/test-introns.hints";
				string repeat_fname = "for-test-cases/repeat-border-test-set/test-repeat.hints";
				string wiggle_fname = "for-test-cases/repeat-border-test-set/test.wig";
				string output_fname = "repeat-border-test-case-utr.gff";
				vector< pair <string, string> > splice_sites;

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);
				vector<Genomic_Data::Repeat> reps_by_file = Genomic_Data::get_all_genes()[0].repeats;

				scaffold_fname = "for-test-cases/softmasking-border-test-set/test.fa";
				coding_region_fname = "for-test-cases/softmasking-border-test-set/test-start-stop.gtf";
				intron_hints_fname = "for-test-cases/softmasking-border-test-set/test-introns.hints";
				repeat_fname = "for-test-cases/softmasking-border-test-set/test-repeat.hints";
				wiggle_fname = "for-test-cases/softmasking-border-test-set/test.wig";
				output_fname = "softmasking-border-test-case-utr.gff";

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, "", splice_sites, false);
				vector<Genomic_Data::Repeat> reps_by_mask = Genomic_Data::get_all_genes()[0].repeats;

				BOOST_CHECK(reps_by_file == reps_by_mask);
			}

			BOOST_AUTO_TEST_CASE( Overlap_Introns_Test_Case )
			{
				//3 (identical) scaffolds
				//2 genes, 5' and 3' utr for each gene, and intron in one 5' utr and one 3' utr
				//overlapping introns (only difference is the overlap situation between two introns)
				//no repeats

				string scaffold_fname = "for-test-cases/overlap-introns-test-set/test.fa";
				string coding_region_fname = "for-test-cases/overlap-introns-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/overlap-introns-test-set/test-introns-overlap.hints";
				string repeat_fname = "for-test-cases/overlap-introns-test-set/test-repeat.hints";
				string wiggle_fname = "for-test-cases/overlap-introns-test-set/test.wig";
				string output_fname = "overlap-introns-test-case-utr.gff";
				string comparison_output_fname = "for-test-cases/comparison-files/comparison-intron-overlap-test-case.gff";
				vector< pair <string, string> > splice_sites;

				ERR_STRM.str(""); //clear stringstream

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);

				UTRs::parameters().output_fname = output_fname;
				compute_UTRs(wiggle_fname);
				UTRs::output();

				string str = ERR_STRM.str();
				BOOST_CHECK(str.empty());
				BOOST_CHECK(UTRs::size() == 12);
				BOOST_CHECK(files_identical( output_fname, comparison_output_fname ));

			}

		#endif


		#ifdef TEST_ERRORS

			BOOST_AUTO_TEST_CASE( Scaffold_Errors_Test_Case )
			{
				//output of every error message in void Genomic_Data::read_scaffold_file(string scaffold_fname)
				//except error message with abort() as a result
				//no introns necessary
				//no repeats necessary
				//no genes

				string scaffold_fname = "for-test-cases/scaffold-errors-test-set/test.fa";
				string coding_region_fname = "for-test-cases/scaffold-errors-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/scaffold-errors-test-set/test-introns.hints";
				string repeat_fname = "for-test-cases/scaffold-errors-test-set/test-repeat.hints";
				string wiggle_fname = "for-test-cases/scaffold-errors-test-set/test.wig";
				string output_fname = "errors_utrs.gff"; //utrs not important
				vector< pair <string, string> > splice_sites;

				ERR_STRM.str(""); //clear stringstream

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);

				UTRs::parameters().output_fname = output_fname;
				compute_UTRs(wiggle_fname);
				UTRs::output();

				string test = "Error in scaffold file (line 1): Empty scaffold found in 'scaffold_empty'!\n"
				"Error in scaffold file (line 2): > in 'scaffold_>' is not first character in its line!\n"
				"Error in scaffold file (line 5): Empty scaffold found in 'scaffold_empty'!\n"
				"Error in scaffold file (line 10): Empty scaffold found in 'scaffold_empty'!\n";

				string str = ERR_STRM.str();
				BOOST_CHECK(str == test);
				BOOST_CHECK(UTRs::size() == 0);

			}

			BOOST_AUTO_TEST_CASE( Start_Stop_Errors_Test_Case )
			{
				//output of every error message in void Genomic_Data::read_coding_region_boundary_file(string coding_region_boundary_fname)
				//except error message with abort() as a result
				//no introns necessary
				//no repeats necessary
				//no errors in scaffold file

				string scaffold_fname = "for-test-cases/start-stop-errors-test-set/test.fa";
				string coding_region_fname = "for-test-cases/start-stop-errors-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/start-stop-errors-test-set/test-introns.hints";
				string repeat_fname = "for-test-cases/start-stop-errors-test-set/test-repeat.hints";
				string wiggle_fname = "for-test-cases/start-stop-errors-test-set/test.wig";
				string output_fname = "errors_utrs.gff"; //UTRs not important
				vector< pair <string, string> > splice_sites;

				ERR_STRM.str(""); //clear stringstream

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);

				UTRs::parameters().output_fname = output_fname;
				compute_UTRs(wiggle_fname);
				UTRs::output();

				string test = "Error in coding region file (line 1): 'scaffold	AUGUSTUS	start_codon	60	62	.	+	0' does not have the right number of columns!\n"
				"Error in coding region file (line 2): 'scaffold1	AUGUSTUS	stop_codon	168	170	.	+	0	transcript_id \"au13.g1.t1\"; gene_id \"au13.g1\";' contains a sequence name, which is not in the scaffold file!\n"
				"Error in coding region file (line 3): In 'scaffold	AUGUSTUS	intron	530	532	.	-	0	transcript_id \"au13.g2.t1\"; gene_id \"au13.g2\";' only start_codon or stop_codon are allowed as features!\n"
				"Error in coding region file (line 4): 'scaffold	AUGUSTUS	start_codon	760	762.1	.	-	0	transcript_id \"au13.g2.t1\"; gene_id \"au13.g2\";' contains a codon position which is not an integer!\n"
				"Error in coding region file (line 5): 'scaffold	AUGUSTUS	stop_codon	368	370	.	.	0	transcript_id \"au13.g1.t1\"; gene_id \"au13.g1\";' contains a strand other than '+' and '-'\n";

				string str = ERR_STRM.str();
				BOOST_CHECK(str == test);
				BOOST_CHECK(UTRs::size() == 0);
			}

			BOOST_AUTO_TEST_CASE( Intron_Errors_Test_Case )
			{
				//output of every error message in void Genomic_Data::read_intron_hints_file(string intron_hints_fname)
				//except error message with abort() as a result
				//no repeats necessary
				//no errors in scaffold file
				//no errors in file with start and stop codons

				string scaffold_fname = "for-test-cases/intron-errors-test-set/test.fa";
				string coding_region_fname = "for-test-cases/intron-errors-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/intron-errors-test-set/test-introns.hints";
				string repeat_fname = "for-test-cases/intron-errors-test-set/test-repeat.hints";
				string wiggle_fname = "for-test-cases/intron-errors-test-set/test.wig";
				string output_fname = "errors_utrs.gff"; //utrs not important
				vector< pair <string, string> > splice_sites;

				ERR_STRM.str(""); //clear stringstream

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);

				UTRs::parameters().output_fname = output_fname;
				compute_UTRs(wiggle_fname);
				UTRs::output();

				string test = "Error in intron file (line 1): 'scaffold	b2h_lib1	intron	775	780	0	-	.' does not have the right number of columns!\n"
				"Error in intron file (line 2): 'scaffold1	b2h_lib1	intron	175	178	0	+	.	mult=10_GTAG' contains a sequence name, which is not in the scaffold file!\n"
				"Error in intron file (line 3): Feature in line 'scaffold	b2h_lib1	exon	75	112	0	.	.	mult=10_GTAG' does not have the feature type 'intron'!\n"
				"Error in intron file (line 4): 'scaffold	b2h_lib1	intron	323	356.1	0	+	.	mult=10_GTAG' contains a coordinate which is not an integer!\n"
				"Error in intron file (line 5): 'scaffold	b2h_lib1	intron	823.5	886	0	+	.	mult=10_GTAG' contains a coordinate which is not an integer!\n"
				"Error in intron file (line 6): 'scaffold	b2h_lib1	intron	A	B	0	+	.	mult=10_GTAG' contains a coordinate which is not an integer!\n"
				"Error in intron file (line 7): 'scaffold	b2h_lib1	intron	175	178	0	+	.	mult=' contains a multiplicity which is not an integer!\n"
				"Error in intron file (line 8): 'scaffold	b2h_lib1	intron	175	178	0	+	.	mult=Eins_GTAG' contains a multiplicity which is not an integer!\n"
				"Error in intron file (line 9): 'scaffold	b2h_lib1	intron	175	178	0	+	.	leereSpalte' contains a multiplicity which is not an integer!\n";

				string str = ERR_STRM.str();
				BOOST_CHECK (str == test);
				BOOST_CHECK(UTRs::size() == 4);
			}

			BOOST_AUTO_TEST_CASE( Repeat_Errors_Test_Case )
			{
				//output of every error message in void Genomic_Data::read_intron_hints_file(string intron_hints_fname)
				//except error message with abort() as a result
				//no introns necessary
				//no errors in scaffold file
				//no errors in file with start and stop codons
				//no errors in intron file

				string scaffold_fname = "for-test-cases/repeat-errors-test-set/test.fa";
				string coding_region_fname = "for-test-cases/repeat-errors-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/repeat-errors-test-set/test-introns.hints";
				string repeat_fname = "for-test-cases/repeat-errors-test-set/test-repeat.hints";
				string wiggle_fname = "for-test-cases/repeat-errors-test-set/test.wig";
				string output_fname = "errors_utrs.gff"; //utrs not important
				vector< pair <string, string> > splice_sites;

				ERR_STRM.str(""); //clear stringstream

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);

				UTRs::parameters().output_fname = output_fname;
				compute_UTRs(wiggle_fname);
				UTRs::output();

				string test = "Error in repeat file (line 1): 'scaffold	repmask	nonexonpart	53	55	0	.	.' does not have the right number of columns!\n"
				"Error in repeat file (line 2): 'scaffold1	repmask	nonexonpart	175	178	0	.	.	src=RM' contains a sequence name, which is not in the scaffold file!\n"
				"Error in repeat (line 3): Feature in line 'scaffold	repmask	exonpart	520	525	0	.	.	src=RM' is not 'nonexonpart'!\n"
				"Error in repeat file (line 4): 'scaffold	repmask	nonexonpart	700	790.7	0	.	.	src=RM' contains a coordinate which is not an integer!\n"
				"Error in repeat file (line 5): 'scaffold	repmask	nonexonpart	960.0	999	0	.	.	src=RM' contains a coordinate which is not an integer!\n"
				"Error in repeat file (line 6): 'scaffold	repmask	nonexonpart	start	end	0	.	.	src=RM' contains a coordinate which is not an integer!\n";

				string str = ERR_STRM.str();
				BOOST_CHECK(str == test);
				BOOST_CHECK(UTRs::size() == 4);
			}

			BOOST_AUTO_TEST_CASE( Wiggle_Errors_Test_Case )
			{
				//output of every error message in class "Process_Wiggle"
				//except error message with abort() as a result
				//Problem most of the error message end with abort()
				//no genes
				//no introns necessary
				//no repeats necessary
				//no errors in scaffold file
				//no errors in file with start and stop codons
				//no errors in intron file
				//no errors in repeat file

				string scaffold_fname = "for-test-cases/wiggle-errors-test-set/test.fa";
				string coding_region_fname = "for-test-cases/wiggle-errors-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/wiggle-errors-test-set/test-introns.hints";
				string repeat_fname = "for-test-cases/wiggle-errors-test-set/test-repeat.hints";
				string output_fname = "errors_utrs.gff"; //utrs not important
				vector< pair <string, string> > splice_sites;

				vector<string> wiggle_fnames;
				wiggle_fnames += "for-test-cases/wiggle-errors-test-set/test1.wig",
								 "for-test-cases/wiggle-errors-test-set/test2.wig",
								 "for-test-cases/wiggle-errors-test-set/test3.wig",
								 "for-test-cases/wiggle-errors-test-set/test4.wig",
								 "for-test-cases/wiggle-errors-test-set/test5.wig",
								 "for-test-cases/wiggle-errors-test-set/test6.wig",
								 "for-test-cases/wiggle-errors-test-set/test7.wig",
								 "for-test-cases/wiggle-errors-test-set/test8.wig",
								 "for-test-cases/wiggle-errors-test-set/test9.wig";

				vector<string> error_messages;
				error_messages += "Error in wiggle file (line 2): Coverage Value in '1 keinDouble' is not a double! Row will be ignored!\n"
						          "Error in wiggle file (line 3): Position in '2.5 12' is not an integer! Row will be ignored!\n",
				                  "Error in wiggle file (line 1): 'variableStep chrom=scaffold span=5' includes an optional span=WindowSize!\n",
				                  "Error in wiggle file! Scaffold 'scaffold1' in wiggle file does not have a sequence in the fasta file! Further information regarding line oder row not available!\n",
				                  "Error in wiggle file (line 1): 'fixedStep chrom=scaffold start=1 step=5 span=5' includes an optional span=WindowSize!\n",
				                  "Error in wiggle file (line 1): Scaffold Position in 'fixedStep chrom=scaffold start=1.1 step=5' is not an integer!\n",
				                  "Error in wiggle file (line 1): Step in 'fixedStep chrom=scaffold start=1 step=5.1' is not an integer!\n",
				                  "Error in wiggle file (line 1): Neither variableStep nor fixedStep, a track line or a comment!\n",
				                  "Error in wiggle file (line 2): Coverage Value '1/2' is not a double!\n",
				                  "Error in wiggle file (line 1): A data line occured before a declaration line!\n";

				ERR_STRM.str(""); //clear stringstream

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);

				UTRs::parameters().output_fname = output_fname;

				for (unsigned i = 0; i < wiggle_fnames.size(); i++) {
					string test;
					string str;
					ERR_STRM.str(""); //clear stringstream

					try {
						compute_UTRs(wiggle_fnames[i]);
					}
					catch (Error e) {
						test = error_messages[i];

						str = ERR_STRM.str();
						BOOST_CHECK (str == test);
						continue;
					}

					UTRs::output();
					test = error_messages[i];
					str = ERR_STRM.str();
					BOOST_CHECK (str == test);
					BOOST_CHECK(UTRs::size() == 0);
				}
			}


			BOOST_AUTO_TEST_CASE( Intron_Splice_Site_Error_Test_Case )
			{
				//testing error, when strands of intron undefined and no splice sites added for splice site filtering
				//no repeats necessary
				//no errors in scaffold file
				//no errors in file with start and stop codons

				string scaffold_fname = "for-test-cases/intron-splice-sites-errors-test-set/test.fa";
				string coding_region_fname = "for-test-cases/intron-splice-sites-errors-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/intron-splice-sites-errors-test-set/test-introns.hints";
				string repeat_fname = "for-test-cases/intron-splice-sites-errors-test-set/test-repeat.hints";
				string wiggle_fname = "for-test-cases/intron-splice-sites-errors-test-set/test.wig";
				vector< pair <string, string> > splice_sites;

				ERR_STRM.str(""); //clear stringstream

				try {
					Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);
				}
				catch (Error e) {
					string str = ERR_STRM.str();
					string test = "Error in intron file (line 1): 'scaffold	b2h_lib1	intron	775	780	0	.	.	mult=10_GTAG' contains an undefined strand '.'. In this case at least one splice site must be given to check for splice sites and define the strand.\n";
					BOOST_CHECK (str == test);
				}
			}

		#endif

		#ifdef TEST_BUGS

			BOOST_AUTO_TEST_CASE( Minimal_Coverage_Test_Case )
			{
				//testing if minimal coverage determination method works (otherwise change method)
				//no repeats
				//no introns
				//no errors in all files
				//currently not in use, should later extract the average coverage of the UTRs

				string scaffold_fname = "for-test-cases/minimal-coverage-test-set/test.fa";
				string coding_region_fname = "for-test-cases/minimal-coverage-test-set/test-start-stop.gtf";
				string intron_hints_fname = "for-test-cases/minimal-coverage-test-set/test-introns.hints";
				string repeat_fname = "for-test-cases/minimal-coverage-test-set/test-repeat.hints";
				string wiggle_fname = "for-test-cases/minimal-coverage-test-set/test.wig";
				string output_fname = "minimal_coverage_utrs.gff";
				vector< pair <string, string> > splice_sites;

				ERR_STRM.str(""); //clear stringstream

				Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites, true);

				UTRs::parameters().output_fname = output_fname;
				compute_UTRs(wiggle_fname);
				UTRs::output();

			}

		#endif

	BOOST_AUTO_TEST_SUITE_END()


	BOOST_AUTO_TEST_SUITE( Real_Data_Test_Suites )

		#ifdef TEST_HUMAN19

			#ifdef TEST_WITHOUT_Z

				BOOST_AUTO_TEST_CASE( Human19_Test_Case_Without_Z )
				{
					//real data test case
					//human chromosome 19
					//End of UTR through maximal coverage drop

					string scaffold_fname = "input/human-chr19/chr19.fa";
					string coding_region_fname = "input/human-chr19/start_stop.gff";
					string intron_hints_fname = "input/human-chr19/trimmed05.blat.sf.gtag.gff";
					string repeat_fname = "input/human-chr19/repeats.gff";
					string wiggle_fname = "input/human-chr19/trimmed-wig/trimmed.sf.wig";
					string comparison_output_fname = "for-test-cases/comparison-files/comparison-human19-test-case-without-z.gff";
					vector< pair <string, string> > splice_sites;
					splice_sites += make_pair("gt","ag");

					string output_fname = "human19_utrs.gff";
					unsigned readLength = 150;
					unsigned smoothingWindowSize = 150;
					unsigned limit = 5000;
					unsigned dropWindowSize = 50;
					unsigned minLength = 2;
					unsigned minCov = 10;
					double pWin = 0.6;
					double pInt = 0.5;
					double pMult = 0.1;
					bool ZeroCov = false;

					Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites);
					UTRs::clear();

					ERR_STRM.str("");

					UTRs::parameters().output_fname = output_fname;
					UTRs::parameters().window_size = smoothingWindowSize;
					UTRs::parameters().read_length = readLength;
					UTRs::parameters().limit = limit;
					UTRs::parameters().drop_window_size = dropWindowSize;
					UTRs::parameters().min_length = minLength;
					UTRs::parameters().min_average_cov = minCov;
					UTRs::parameters().p_win = pWin;
					UTRs::parameters().p_int = pInt;
					UTRs::parameters().p_mult = pMult;
					UTRs::parameters().zero_cov = ZeroCov;

					compute_UTRs(wiggle_fname);
					UTRs::output();

					BOOST_CHECK(UTRs::size() == 481);
				}

			#endif

			#ifdef TEST_WITH_Z

				BOOST_AUTO_TEST_CASE( Human19_Test_Case_with_Z )
				{
					//real data test case
					//human chromosome 19
					//UTR end when coverage drops down to zero the first time

					string scaffold_fname = "input/human-chr19/chr19.fa";
					string coding_region_fname = "input/human-chr19/start_stop.gff";
					string intron_hints_fname = "input/human-chr19/trimmed05.blat.sf.gtag.gff";
					string repeat_fname = "input/human-chr19/repeats.gff";
					string wiggle_fname = "input/human-chr19/trimmed-wig/trimmed.sf.wig";
					string comparison_output_fname = "for-test-cases/comparison-files/comparison-human19-test-case-with-z.gff";
					vector< pair <string, string> > splice_sites;
					splice_sites += make_pair("gt","ag");

					string output_fname = "human19_utrs.gff";
					unsigned readLength = 150;
					unsigned smoothingWindowSize = 1;
					unsigned limit = 5000;
					unsigned dropWindowSize = 10;
					unsigned minLength = 2;
					unsigned minCov = 1;
					double pWin = 0.1;
					double pInt = 0.5;
					double pMult = 0.1;
					bool ZeroCov = true;

					Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites);
					UTRs::clear();

					ERR_STRM.str("");

					UTRs::parameters().output_fname = output_fname;
					UTRs::parameters().window_size = smoothingWindowSize;
					UTRs::parameters().read_length = readLength;
					UTRs::parameters().limit = limit;
					UTRs::parameters().drop_window_size = dropWindowSize;
					UTRs::parameters().min_length = minLength;
					UTRs::parameters().min_average_cov = minCov;
					UTRs::parameters().p_win = pWin;
					UTRs::parameters().p_int = pInt;
					UTRs::parameters().p_mult = pMult;
					UTRs::parameters().zero_cov = ZeroCov;

					compute_UTRs(wiggle_fname);
					UTRs::output();

					BOOST_CHECK(UTRs::size() == 1597);
				}

			#endif

			#ifdef TEST_SPLICE_SITE_FILTER

				BOOST_AUTO_TEST_CASE( Human19_Test_Case_Splice_Site_Filter )
				{
					//real data test case
					//human chromosome 19
					//testing if splice site filtering works correctly
					//TODO: not working as it should, because the introns after splice site filtering are not easily extractable

					string scaffold_fname = "input/human-chr19/chr19.fa";
					string coding_region_fname = "input/human-chr19/start_stop.gff";
					string intron_hints_fname = "input/human-chr19/trimmed05.blat.sf.gff";
					string repeat_fname = "input/human-chr19/repeats.gff";
					string wiggle_fname = "input/human-chr19/trimmed-wig/trimmed.sf.wig";
					string comparison_output_fname = "for-test-cases/comparison-files/comparison-human19-test-case-without-z.gff";
					vector< pair <string, string> > splice_sites;
					splice_sites += make_pair("gt","ag");

					string output_fname = "human19_utrs.gff";
					unsigned readLength = 150;
					unsigned smoothingWindowSize = 150;
					unsigned limit = 5000;
					unsigned dropWindowSize = 50;
					unsigned minLength = 2;
					unsigned minCov = 10;
					double pWin = 0.6;
					double pInt = 0.5;
					double pMult = 0.1;
					bool ZeroCov = false;

					Genomic_Data::initialize(scaffold_fname, coding_region_fname, intron_hints_fname, repeat_fname, splice_sites);
					UTRs::clear();

					ERR_STRM.str("");

					UTRs::parameters().output_fname = output_fname;
					UTRs::parameters().window_size = smoothingWindowSize;
					UTRs::parameters().read_length = readLength;
					UTRs::parameters().limit = limit;
					UTRs::parameters().drop_window_size = dropWindowSize;
					UTRs::parameters().min_length = minLength;
					UTRs::parameters().min_average_cov = minCov;
					UTRs::parameters().p_win = pWin;
					UTRs::parameters().p_int = pInt;
					UTRs::parameters().p_mult = pMult;
					UTRs::parameters().zero_cov = ZeroCov;

					compute_UTRs(wiggle_fname);
					UTRs::output();

					BOOST_CHECK(files_identical( output_fname, comparison_output_fname ));

				}

			#endif

		#endif


	BOOST_AUTO_TEST_SUITE_END()

#endif




