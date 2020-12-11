# About AUGUSTUS

AUGUSTUS is a gene prediction program written or maintained by Mario Stanke, Oliver Keller, Stefanie König, Lizzy Gerischer, Katharina Hoff, Giovanna Migliorelli, Lars Gabriel, Anica Hoppe, Tonatiuh Peña Centeno, Henry Mehlan, Daniel Honsel and Steffen Herbold. It can be used as an ab initio program, which means it bases its prediction purely on the sequence. AUGUSTUS may also incorporate hints on the gene structure coming from extrinsic sources such as EST, MS/MS, protein alignments and syntenic genomic alignments. Since version 3.0 AUGUSTUS can also predict the genes simultaneously in several aligned genomes (see [README-cgp.md](README-cgp.md)). AUGUSTUS has been trained for predicting genes, among other species, in: 

- Homo sapiens (human), 
- Drosophila melanogaster (fruit fly), 
- Arabidopsis thaliana (plant),
- Brugia malayi (nematode),
- Aedes aegypti (mosquito),
- Coprinus cinereus (fungus),
- Tribolium castaneum (beetle)
- Schistosoma mansoni (worm)
- Tetrahymena thermophila (ciliate)
- Galdieria sulphuraria (red algae)
- Zea mays (maize)
- Toxoplasma gondii (parasitic protozoa)
- Caenorhabditis elegans (worm)
- Aspergillus fumigatus
- Aspergillus nidulans
- Aspergillus oryzae
- Aspergillus terreus
- Botrytis cinerea
- Callorhinchus milii
- Candida albicans
- Candida guilliermondii
- Candida tropicalis
- Chaetomium globosum
- Coccidioides immitis
- Cryptococcus neoformans gattii
- Cryptococcus neoformans neoformans
- Danio rerio
- Debaryomyces hansenii
- Encephalitozoon cuniculi
- Eremothecium gossypii
- Fusarium graminearum
- Gallus gallus
- Histoplasma capsulatum
- Kluyveromyces lactis
- Laccaria bicolor
- Lodderomyces elongisporus
- Magnaporthe grisea
- Neurospora crassa
- Nicotiana attenuata (coyote tobacco)
- Petromyzon marinus (sea lamprey)
- Phanerochaete chrysosporium
- Pichia stipitis
- Rhizopus oryzae
- Saccharomyces cerevisiae
- Schizosaccharomyces pombe
- Staphylococcus aureus
- Ustilago maydis
- Yarrowia lipolytica
- Nasonia vitripennis (wasp)
- Solanum lycopersicum (tomato)
- Chlamydomonas reinhardtii (green algae)
- Amphimedon queenslandica (sponge)
- Acyrthosiphon pisum (pea aphid)
- Leishmania tarentolae (protozoa, intronless)
- Trichinella spiralis
- Theobroma cacao (cacao)
- Escherichia coli
- Thermoanaerobacter tengcongensis (a bacterium)
- Triticum aestivum (wheat)
- Ancylostoma ceylanicum
- Volvox carteri
- Mnemiopsis leidyi (Ctenophora)
- Nematostella vectensis (Cnidaria)
- Ciona intestinalis (Chordata)
- Strongylocentrotus purpuratus (Echinodermata)
- Pisaster ochraceus (starfish)
- Chiloscyllium punctatum (bamboo shark)
- Scyliorhinus torazame (cat shark)
- Rhincodon typus (whale shark)
  
The training annotation files of the following species are a courtesy of Jason Stajich (see also http://fungal.genome.duke.edu/): Aspergillus fumigatus, Aspergillus nidulans, Aspergillus oryzae, Aspergillus terreus, Botrytis cinerea, Candida albicans, Candida guilliermondii, Candida tropicalis, Chaetomium globosum, Coccidioides immitis, Coprinus cinereus, Cryptococcus neoformans gattii, Cryptococcus neoformans neoformans, Debaryomyces hansenii, Encephalitozoon cuniculi, Eremothecium gossypii, Fusarium graminearum, Histoplasma capsulatum, Kluyveromyces lactis, Laccaria bicolor, Lodderomyces elongisporus, Magnaporthe grisea, Neurospora crassa, Phanerochaete chrysosporium, Pichia stipitis, Rhizopus oryzae, Saccharomyces cerevisiae, Schizosaccharomyces pombe, Ustilago maydis, Yarrowia lipolytica.

The training for 'sealamprey' (Petromyzon marinus) was performed by Falk Hildebrand and Shigehiro Kuraku, based on the genome assembly (PMAR3.0) provided by the Genome Sequencing Center at Washington University School of Medicine (WUGSC) in St. Louis. The training is described in: "Molecular Evolution in the Lamprey Genomes and Its Relevance to the Timing of Whole Genome Duplications" T. Manousaki, H. Qiu, M. Noro, F. Hildebrand, A. Meyer and S. Kuraku in 'Jawless Fishes of the World, Volume 1', Cambridge Scholars Publishing http://www.cambridgescholars.com/jawless-fishes-of-the-world

The training for elephant shark (Callorhinchus milii) was performed by Tereza Manousaki and Shigehiro Kuraku, based on the genome assembly (made up of 1.4x whole genome shotgun reads) available at http://esharkgenome.imcb.a-star.edu.sg/resources.html. The training for 'japaneselamprey' (Lethenteron japonicum or Lethenteron camtschaticum) was performed by Shigehiro Kuraku.

The training for Pneumocystis jirovecii was performed by Marco Pagni, Philippe Hauser et al as described in Hauser PM, Burdet FX, Cissé OH, Keller L, Taffé P, Sanglard D, Pagni M., Comparative Genomics Suggests that the Fungal Pathogen Pneumocystis Is an Obligate Parasite Scavenging Amino Acids from Its Host's Lungs. PLoS One. 2010, Dec 20;5(12):e15152. PubMed PMID: 21188143; PubMed Central PMCID: PMC3004796. 

The parameter training for cacao was done by Don Gilbert (gilbertd at indiana.edu).

Parameters for Ancylostoma ceylanicum were trained by Erich Schwarz (Cornell University).

Parameters for camponotus_floridanus (ant) were contributed by Shishir K Gupta.

Parameters for Apis dorsata were contributed by Francisco Camara Ferreira.

Parameters for Mnemiopsis leidyi, Nematostella vectensis, Ciona intestinalis and Strongylocentrotus purpuratus were contributed by Joseph Ryan (University of Florida).

Parameters for Chiloscyllium punctatum (brownbanded bamboo shark), Scyliorhinus torazame (cloudy catshark) and Rhincodon typus (whale shark) were contributed by Shigehiro Kuraku; see Hara, Yamaguchi, et al., 2018. Nature Ecology and Evolution, doi: 10.1038/s41559-018-0673-5.
