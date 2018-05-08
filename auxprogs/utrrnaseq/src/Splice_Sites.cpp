/**
 * \file Splice_Sites.cpp
 */

#include "Splice_Sites.hpp"
#include "Supporting_Methods.hpp"
#include "Global.hpp"

#include <string>
#include <vector>

#include <utility>
#include <ctype.h>
#include <iostream>

#include <boost/assign/std/vector.hpp>
#include <boost/algorithm/string/case_conv.hpp>

using namespace std;
using namespace boost::assign;
using namespace boost::algorithm;



pair<bool, string> check_sp_sites (const vector< pair <string, string> >& sp_sites, const string& sequence,
		unsigned start, unsigned end, bool strand_defined, string strand) {
	string don_site_scaff;
	string acc_site_scaff;
	unsigned sp_site_no = 0;
	unsigned length_don_site;
	unsigned length_acc_site;

	while ( (sp_site_no < sp_sites.size()) ) {
		length_don_site = sp_sites[sp_site_no].first.size();
		length_acc_site = sp_sites[sp_site_no].second.size();

		if ( ((start + length_don_site - 1) < end) && ((end - length_acc_site + 1) > start) ) {

			if (!strand_defined || strand == "+" ) {//possible plus strand or plus strand
				don_site_scaff = sequence.substr(start-1, length_don_site);
				acc_site_scaff = sequence.substr(end-length_acc_site, length_acc_site);

				to_lower(don_site_scaff); //to lower case, so unified format, donor and acceptor sites in lower case
				to_lower(acc_site_scaff);

				if ( (sp_sites[sp_site_no].first == don_site_scaff)
						&& (sp_sites[sp_site_no].second == acc_site_scaff) ) {
					return make_pair(true, "+");
				}
			}

			if (!strand_defined || strand == "-") {//possible minus strand or minus strand
				don_site_scaff = sp_site_conversion(sequence, end-1, end-length_don_site);
				acc_site_scaff = sp_site_conversion(sequence, start+length_acc_site-2, start-1);

				to_lower(don_site_scaff); //to lower case, so unified format, donor and acceptor sites in lower case
				to_lower(acc_site_scaff);

				if ( (sp_sites[sp_site_no].first == don_site_scaff)
						&& (sp_sites[sp_site_no].second == acc_site_scaff) ) {
					return make_pair(true, "-");
				}
			}
		}
		++sp_site_no;
	}

	return make_pair(false, strand);
}


string sp_site_conversion (const string& sequence, unsigned start, unsigned end) {
	string sp_site;

	for ( unsigned i = start; i >= end; i-- ) {
		char base = sequence[i];
		base = toupper(base);

		switch(base) {
			case 'A':
				sp_site += "T";
				break;
			case 'C':
				sp_site += "G";
				break;
			case 'G':
				sp_site += "C";
				break;
			case 'T':
				sp_site += "A";
				break;
			default: //Not specific Base, e.g. N,...
				sp_site += "N"; //in any case not an accepted splice site
		}
	}

	return sp_site;
}


vector< pair <string, string> > parse_sp_sites(string raw_sp_sites) {
	const string SEP_SP_SITES = "#";
	const string SEP_DON_ACC_SITE = "_";

	vector< pair<string, string> > sp_sites;

	unsigned all_pos = raw_sp_sites.find("all");
	if (all_pos != raw_sp_sites.npos) { //"all" found
		return sp_sites; //empty sp_sites returned, no filtering
	}

	//filtering, so splice sites must be parsed further
	vector<string> paired_sp_sites = tokenize(SEP_SP_SITES, raw_sp_sites);

	for (unsigned i = 0; i < paired_sp_sites.size(); i++) {
		if (!paired_sp_sites[i].empty()) {
			pair<string, string> seperated_sp_site;
			vector<string> don_acc_site = tokenize(SEP_DON_ACC_SITE, paired_sp_sites[i]);

			if (don_acc_site.size() == 2) { //only accepted when one donor and one acceptor site
				to_lower(don_acc_site[0]); //lower case, unified format
				to_lower(don_acc_site[1]);

				seperated_sp_site = make_pair(don_acc_site[0], don_acc_site[1]);
				sp_sites += seperated_sp_site;
			}
			//else Incomplete or unaccepted splice site. No donor or acceptor side (_AG, GT_, or no separator)
		}
	}

	return sp_sites;
}

