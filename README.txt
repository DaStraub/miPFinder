MiPFinder: Cross-species genome-wide identification of evolutionary conserved microProteins.

Here is miPFinder described, a program to identify microProteins, novel regulators of protein function. 
MiPFinder starts with a set of protein sequences and considers information about protein size, sequence similarity and domain composition to create a list of miP candidates. 

This version (v1) was published.

The algorithm is available under the GPLv3 licence at https://github.com/DaStraub/miPFinder

GENERAL INFORMATION

The miPFinder program takes a single command line in the windows command prompt (e.g. “python miPFinder.py -f proteins.fasta -p ProteinGeneList.tsv -a annotation.tsv”). 
The minimum input requirement is a simple fasta file with all protein sequences (“-f”), however a file with protein annotations (“-a”) will aid the miP selection tremendously. 
For the addition of protein-protein interaction domain information, a Pfam domain database (“-d”) and a file specifying interaction domains (“-i”) is necessary. 
Moreover, a file specifying the protein-gene-relationship (“-p”) will allow for cis-miP detection, for filtering putative ancestors for their longest splice variant, and for the removal of redundant miP candidate splice variants. Parameters for the maximal miP and minimum ancestor length can be adjusted (“-M” and “-A”, respectively, standard setting: 140 and 250) as well as all cutoff values.
miPFinder is built with Python v2.7.9 running on Microsoft Windows 7 and using hmmer v3.1b1, blast+ v2.2.29, clustalw v2.1, but any python2, hmmer3, blast2, clustalw2 and Microsoft Windows version might be sufficient. Path to the dependencies (hmmer, blast, clustalw2) must be specified, if the accessory programs are not set as environment variables, using command line arguments “-H”,”-B”,”-C”, respectively. MipFinder will check the availability of specified input files and correct function of all dependencies before each run.
For further information about input requirements and data sources see below.

DETAILED INPUT

- USAGE
python miPFinder_v1.py -f FASTA.fasta [DEPENDENCIES, DATABASES, OTHER VARIABLES]

- REQUIRED

-- Sequence file ('-f', '--fasta')
FASTA file of all proteins, recommended: http://www.phytozome.net, http://www.ensembl.org/, ftp://ftp.ncbi.nlm.nih.gov/refseq ->protein.faa files

- DEPENDENCIES

-- BLAST path ('-B', '--blast', default='')
Path/to/blast-folder/, available at ftp://ftp.ncbi.nih.gov/blast/executables/blast+/'
Expects "makeblastdb.exe" and "blastp" commands to work when targeting that folder.

-- CLUSTALW2 path ('-C', '--clustalw', default='')
Path/to/clustalw2-folder/, available at http://www.clustal.org/clustal2/'
Expects "clustalw2.exe".

-- HMMER path ('-H', '--hmm', default='')
Path/to/hmmsearch-folder/, available at http://hmmer.org/'
Expects "hmmsearch.exe", "hmmbuild.exe" and "hmmscan.exe".

- DATABASES

-- PFAM database ('-d', '--PfamA')
Pfam database in .hmm format, e.g. Pfam-A_v28.hmm, available at ftp://ftp.ebi.ac.uk/pub/databases/Pfam/'.

-- iPFAM database ('-i', '--iPfam')
iPfam database, merged homodomain and heterodomain interaction file, available at http://www.ipfam.org/.
Expected format is 3 columns for homodomain interactions and 5 columns for heterodomain interactions.
Homodomain interaction information (c = column): Pfam ID (c1), either "both" or "interchain" (c3) therefore "intrachain" interaction domains will be skipped.
Heterodomain interaction information (c = column): Pfam ID (c1), Pfam ID (c3), either "both" or "interchain" (c5) therefore "intrachain" interaction domains will be skipped.

-- Protein-gene relation ('-p', '--ProteinGeneList')
List of identifiers of proteins and genes. 
Relevant information (c = column): Protein ID (c1), Gene ID (c2); tab separated

-- Annotation file ('-a', '--annotation')
Annotation file for FASTA file. Protein IDs have to match IDs in FASTA file.
Relevant information (c = column): Protein ID (c1), Description (c2); tab separated

-- Domain annotation file ('-m', '--hmmscanDB')
Domain annotation file in hmmscan's "--domtblout" format (http://hmmer.org/).
This file could be the entry point for a custom domain annotation file.
Relevant information (c = column): Protein ID (c4), Pfam accession (c2), E-value (c7), c-Evalue (c12), ali coord from (c18), ali coord to (c19), Pfam description (c23)

-- IDs of intrest ('-l', '--IDlist', dest='IDlist', default=known miPs)
Semicolon separated list of IDs that will be searched against miP candidate protein IDs and reported in the final result table.
known miPs = 'AT5G39860.1;AT1G26945.1;AT5G15160.1;AT3G28857.1;AT1G74500.1;AT3G47710.1;AT4G15248.1;AT3G21890.1;AT3G28917.1;AT1G74660.1;AT1G18835.1;AT1G14760.2;AT3G58850.1;AT2G42870.1;AT1G01380.1;AT2G30424.1;AT5G53200.1;AT2G30420.1;AT2G30432.1;AT2G46410.1;AT4G01060.1;AT3G52770.1'

- OTHER VARIABLES

-- Output file prefix ('-s', '--species', default='NA').
-- Maximum length of a microProtein in aminoacids ('-M', '--maxMIPlength', default=140).
-- Minimum length of a microProtein ancestor in aminoacids ('-A', '--minANCESTORlength', default=250).
-- E-value cutoff for Blast search ('-L', '--blastCUTOFF', default=1e-3).
-- Pfam domains: minimum overlap of a microProtein with an annotated domain ('-O', '--overlapCUTOFF', default=0.6).
-- E-value cutoff for HMMscan and HMMsearch ('-E', '--evalueCUTOFF', default=0.1).
-- c-Evalue cutoff for HMMscan and HMMsearch ('-V', '--cvalueCUTOFF', default=0.05).

OUTPUT

- FILES

-- "<prefix>_miPlist.txt"
Main output, contains miP candidates and its additional information.
For full functionality, paste the file contents into a MS Office Excel sheet using Windows and keep it in the same folder as the "alignment" folder.

-- "<prefix>_hmmscandb.txt"
Contains the domain table produced by HMMER. 
Reusing this file in subsequent miPFinder runs of the same FASTA file will save a lot of time.

-- "SkippedGroups.txt"
Contains skipped miP groups that contain >1000 small protein sequences. These groups are unlikely miP candidates.

- FOLDERS

-- "alignment"
Contains alignment files ("*.txt", "*_block.txt") that are linked in "<prefix>_miPlist.txt" and score ("*_individual.txt") files.
This folder must stay in the same directory as "<prefix>_miPlist.txt" for full functionality.

-- "blast"
Contains BLAST database, BLAST results and FASTA files of small and big proteins.
Reusing this folder in subsequent miPFinder runs of the same FASTA file will save time.

-- "tempfiles"
Contains CLUSTALW and HMMER files that are not used in subsequent miPFinder runs or provide any information in "<prefix>_miPlist.txt".
Can be helpful for trouble shooting in case of program failure.

CITATION

Daniel Straub, Stephan Wenkel: Cross-species genome-wide identification of evolutionary conserved microProteins. bioRxiv 061655; doi: http://dx.doi.org/10.1101/061655 