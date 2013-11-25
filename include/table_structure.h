#ifndef _TABLE_STRUCTURE
#define _TABLE_STRUCTURE

#include <mysql++.h>
#include <ssqls.h>
#include <string>
#include <cstring>
#include <stdlib.h>
using namespace std;

/*
 * tables structure for comparative gene prediction, all in one table, one database for all species
 */
sql_create_6(genomes,
	     1,6,
	     int,seqid,
	     string,dnaseq,
	     string,seqname,
	     int,start,
	     int,end,
	     string,species)
sql_create_2(speciesnames,
	     1,2,
	     int,speciesid,
	     string,speciesname)
sql_create_3(seqnames,
	     1,3,
	     int,seqnr,
	     int,speciesid,
	     string,seqname)
sql_create_14(hints,
	      1,14,
	      int,hintid,
	      int,speciesid,
	      int,seqnr,
	      string,source,
	      int,start,
	      int,end,
	      double,score,
	      int,type,
	      string,strand,
	      string,frame,
	      int,priority,
	      string,grp,
	      int,mult,
	      string,esource)
// The database schema is an excerpt from ENSEMBL:
// http://www.ensembl.org/info/docs/api/core/core_schema.html
// The following is calling a very complex macro which will create
// the table structure row in a STL container.
sql_create_2(dna,
	     1, 2,
	     int,seq_region_id,
	     std::string, sequence)  
sql_create_4(seq_region,
	     1,4,
	     int,seq_region_id,
	     std::string,name,
	     std::string,coord_system_id,
	     int,length)
sql_create_6(assembly,
	     1, 6,
	     int, asm_seq_region_id,
	     int, cmp_seq_region_id,
	     int, asm_start,
	     int, asm_end,
	     int, cmp_start,
	     int, cmp_end)
// omit species name, assume every species has its own database,
// 'name' is the shared key in 'seq_region' table which refers to 'chrName' in augustus
// it could also be scaffold.'query_region' is defined in database 'cmpproject', so is 'augustus-gff'


#endif //_TABLE_STRUCTURE
