# getKASP_pipeline

These scripts make a simple pipeline to design KASP (Kompetitive Allele Specific PCR) primers.

Polymarker (http://polymarker.tgac.ac.uk/) is a great software to design KASP primers, but sometimes I have some specific requirements that Polymarker cannot meet, and it is difficult for me to modify its scripts because I do not know Ruby. That is why I wrote these simple scripts to just meet my requirements.

# Pseudo code
1. Read the polymarker input and get:
	- snp position
	- sequences and make them a fasta file for blast later
2. Blast using the fasta file and output blast results
3. Process blast output file to extract flanking sequences (250 bp each side)
4. Multiple sequence alignment of the homeologs
5. Use the msa file to design primers using primer3

# Dependencies

getKASP_pipeline needs following 3 software to find differences among homeologs and design primer.
1. Muscle: Multiple sequence alignment program (http://www.drive5.com/muscle/)
2. Primer3: program for designing PCR primers (http://primer3.sourceforge.net/)
3. blast+ package from NCBI (https://blast.ncbi.nlm.nih.gov/Blast.cgi)

Please make sure "muscle", "primer3_core" and "blastn" are in the software PATH. Otherwise, please modify specific scripts and give the software path.

# How it works
1. Find all the different sites that can differ all other sequences from the user provided alignment file;
2. Use these sites and the SNP site as SEQUENCE_FORCE_RIGHT_END in primer3 to design all possible left and right primers in the target sequence.

# Usage

I divided the pipeline into 6 steps:
- Script "parse_polymarker_input.py": parse the polymarker input and prepare a fasta file for blast
- Blast using system command "blastn"
- Script "getflanking.py": Parse the blast output file and output the homelog contigs and flanking ranges
- Split the range file for each marker with system command "awk"
- Get flanking sequences for each file with command "blastdbcmd"
- Get KASP primers using script "getkasp3.py"

You can run this step by step or run the whole pipeline with script "run_getkasp.py". I suggest run the 6 steps in the script "run_getkasp.py" stey by step to get familiar how it works first.

Example: `run_getkasp.py for_polymarker.csv`

Change the software paths, blast contig names and locations etc accordingly.

The "bin" folder has all the scripts for each step and softare primer3 and muscle in case your system does not have them.

# Acknowledgements
I borrowed ideas from the polymarker scripts (https://github.com/TGAC/bioruby-polyploid-tools), a great tool for Genome Specific KASPar design in polyploid species. Thanks to the author of Polymarker.

I also borrowed some codes from biopython (https://github.com/biopython/biopython/blob/master/Bio/Emboss/Primer3.py). Thanks to them too.

Thanks to the open source software **Primer3** (http://primer3.sourceforge.net/), **Muscle** (http://www.drive5.com/muscle/), and blast+ package from NCBI (https://blast.ncbi.nlm.nih.gov/Blast.cgi).
