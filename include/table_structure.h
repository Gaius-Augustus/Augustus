#ifndef _TABLE_STRUCTURE
#define _TABLE_STRUCTURE

#include <mysql++.h>
#include <ssqls.h>
#include<string>
#include<cstring>
#include<stdlib.h>
using namespace std;
// The following is calling a very complex macro which will create
// table structure row in a STL container.
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
//ommit species name,assume every species has its own database,'name' is the shared key in 'seq_region' table which refers to 'chrName' in augustus,it could also be scaffold.'query_region' is defiend in database 'cmpproject',so is 'augustus-gff'


#endif //_TABLE_STRUCTURE
