# wheat training from 11.06.2013

# genome:
# genome assembly of the IWGSC (International Wheat Genome Consortium)
# all contigs longer than 15kb (~5% of the complete assembly) were taken and masked using RepeatScout and RepeatMasker (~> 64% masked)
# training was then performed on the masked genome

# training data:
# 1286040 Triticum aestivum ESTs from NCBI + 6137 FL-cDNAs (Komugi-TriFLDB)

# estimated accuracy:
+-----------------------------------------+
|                SN           SP          |
+------------+----------------------------+
| gene level |  0.27        0.228         |
| exon level |  0.711       0.597         |
| nuc. level |  0.783       0.62          |
+------------+----------------------------+
