#	Makefile for homGeneMapping
#
#	Created: 16-October-2015
#	Last modified: 16-October-2015 by Katharina Hoff (katharina.hoff@gmail.com)

all:
	(cd src; make)

clean:
	(rm -rf test; cd src; make clean)

test:
	rm -rf test
	mkdir -p test
	# test with hints
	(cd examples; ../../../bin/homGeneMapping --noDupes --gtfs=gtffilenames.tbl                                                                  --halfile=aln.hal --tmpdir=../test/tmp               --outdir=../test/outdir               --printHomologs=../test/homologs.txt               > ../test/homGeneMapping.out)
	# test without hints
	cat examples/gtffilenames.tbl | sed "s/gtf\t.*/gtf/g" > test/gtffilenames_without_hints.tbl
	(cd examples; ../../../bin/homGeneMapping --noDupes --gtfs=../test/gtffilenames_without_hints.tbl                                            --halfile=aln.hal --tmpdir=../test/tmp_without_hints --outdir=../test/outdir_without_hints --printHomologs=../test/homologs_without_hints.txt > ../test/homGeneMapping_without_hints.out)
	# test with hints provided by SQLite db
	../../bin/load2sqlitedb --species=hg19 --dbaccess=test/homGeneMapping_hints.db ../../examples/cgp/human.fa
	../../bin/load2sqlitedb --species=mm9  --dbaccess=test/homGeneMapping_hints.db ../../examples/cgp/mouse.fa
	../../bin/load2sqlitedb --noIdx --species=hg19 --dbaccess=test/homGeneMapping_hints.db examples/gtfs/human.hints.gff
	../../bin/load2sqlitedb --noIdx --species=mm9  --dbaccess=test/homGeneMapping_hints.db examples/gtfs/mouse.hints.gff
	../../bin/load2sqlitedb --makeIdx --dbaccess=test/homGeneMapping_hints.db
	(cd examples; ../../../bin/homGeneMapping --noDupes --gtfs=../test/gtffilenames_without_hints.tbl --dbaccess=../test/homGeneMapping_hints.db --halfile=aln.hal --tmpdir=../test/tmp_sqlite_hints --outdir=../test/outdir_sqlite_hints --printHomologs=../test/homologs_sqlite_hints.txt > ../test/homGeneMapping_sqlite_hints.out)