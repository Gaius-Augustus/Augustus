Project: similarity-score algorithm for block-profile and protein sequence
with intron information
Author:  Lars Gabriel

Table of Contents:
1. INTRODUCTION
2. INSTALLATION
3. RUNNING pp_simScore
4. FILE FORMATS
5. OPTIONS
6. EXAMPLES

1. INTRODUCTION

The program pp_simScore computes the similarity score and optimal alignments
of a block-profile and a protein sequence.
The algorithm can optionally take intron positions into account.

2. INSTALLATION
AUGUSTUS has to be installed (see README.md)

3. RUNNING pp_simScore
The program can be run with following command line from the directory src/

    ./pp_simScore [options] --fasta protein_sequence_file.fa --prfl protein_profile_file.prfl


4. FILE FORMATS

    Potein sequence:

        The file protein_sequence_file.fa has to be in FASTA format.
        It may contain an optional [Intron] section. This section denotes the
        intron positions in the protein sequence, which are specified as list of
        (j, f), where j is the index of the amino acid after witch the intron
        immediately occurs. The indices range from 0 to m - 1 if the protein
        sequence has a length of m.
        >protein sequence header
         XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
         XXXXXXX protein sequence XXXXXXXXXXX
         XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

         [Introns]
         # index of the position after which an intron occurs | residual nucleotides before the intron
         2 0
         5 1
         30 2
         104 1


    Block-Profile:

    The block profile file has to have following structure:

         [dist]
         min    max
         [block]
         B
         [intron profile]
         w
         inter-block_profile_list
         intra-block_profile_list

    This structure can be repeated in this file. The file has to end either
    in a [dist] section or a [dist] and than [intron profile] section.
    The [intron profile] sections are optional.

    [dist] explanation:
    min, max denote the distance interval of an inter-block section

    [block] explanation:
    B denotes a (20 x t) matrix for a block of t of the block-profile

    [intron profile] explanation:
    an intron profile describes the positions and frequencies of introns in and
    before the associated block
    w: number of protein family members used to build the intron profile
    inter-block_profile_list: list of (h, v),
    where   h denotes the number of introns which occurred within a family member,
            v the number of family members which have this number of introns
    intra-block_profile_list: list (s, f, v),
    where   s denotes the index of the position in the block after which an intron occurs,
            f denotes the number of nucleotides which are left before the intron (0,1,2)
            v the number of family members which have an intron at that position

5. OPTIONS

    --gap_inter:                gap costs for an alignment column that is a gap
                                in an inter-block section
                                default setting: -5

    --gap_intra:                gap costs for an alignment column that is a gap in a block
                                default setting: -50

    --gap_intron:               gap costs for an gap in intron positions
                                default setting: -5

    --epsilon_intron:           pseudocount parameter epsilon1, the pseudocount
                                is added to a relative intron frequency v/w with
                                (v+epsilon1)/(w+epsilon1+epsilon2)
                                default setting: 0.0000001

    --epsilon_noIntron:         pseudocount parameter epsilon2, the pseudocount
                                is added to a relative intron frequency v/w with
                                (v+epsilon1)/(w+epsilon1+epsilon2)
                                default setting: 0.1

    --intron_weight_intra:      value that is added to an intron score for a
                                match of intron positions in a block
                                default setting: 5

    --intron_weight_inter:      value that is added to an intron score for a
                                match of intron positions in an inter-block
                                default setting: 5

    --alignment:                number of optimal alignments that are computed
                                default setting: 1

    --help:                     print USAGE

    --out:                      denotes the output format, the following output options,
                                between " ", are implemented:

                                       "score" : output is the similarity score
                                       "matrix" : output are similarity matrix and similarity score
                                       "alignment" : output are the computed alignments to the console as

                                Alignment representation of P as symbols of
                                {AminoAcid, gap symbol or number of amino acids in inter-block}
                                Alignment representation of argmax of B as symbols of
                                {argmax AminoAcid for aligned block column, gap symbol or inter-block length}
                                Frequency of amino acid of P in aligned block column of B,
                                if alignment type is a match

                                        "matrix+alignment": output are similarity matrix,
                                                            similarity score and the computed
                                                            alignments in the format described above
                                        "db" : output are the computed alignment
                                               as list of alignment frames,
                                               an element of the list consists of:

                                - starting position of the first amino acid of the protein sequence
                                that is included in the alignment frame
                                - block number in which the alignment frame is located
                                - index of the first block column that is included in the alignment frame
                                - length of the frame (number of alignment columns)
                                - alignment type: 'm', 's'. 'p' or '-'

                                        "bp" : output is a list of translations from the index of a block
                                               to the number of the block in the .prfl file
                                        "consents" : output is the average of the argmax
                                                     of the block columns for the complete profile
                                        "interblock" : output is a list of all inter-block distance intervals

                                default setting: "score"

6. EXAMPLES

    Data:   intron-block-profile: examples/sim-score/EOG09150290.prfl
            protein sequence with intron positions: examples/sim-score/EDW03868.1.fa

    Example of the command line for the output of an alignment and the
    similarity score of the protein sequence and the intron-block-profile:

    ./pp_simScore --fasta EDW03868.1.fa --prfl EOG09150290.prfl --out alignment

    Example of the command line for the output of an alignment and the
    similarity score of the protein sequence and the intron-block-profile,
    without taking the intron positions into account:

    ./pp_simScore --fasta EDW03868.1.fa --prfl EOG09150290.prfl \
    --intron_weight_intra 0 --intron_weight_inter 0 --out alignment
