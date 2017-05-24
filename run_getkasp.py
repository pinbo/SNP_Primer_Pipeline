#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  run_getkasp.py
#  
#  Copyright 2017 Junli Zhang <zhjl86@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

# Just to run all the steps at once for designing KASP primers
# Need to use in MacPro

# change the reference location accordingly.

# example: run_getkasp.py for_polymarker.csv 3

#########################
from glob import glob


def main(args):
	polymarker_input = args[1]
	genome_number =  args[2]
	script_path = os.path.dirname(os.path.realpath(__file__)) + "/bin/" # scripts folder
	#reference = "/Library/WebServer/Documents/blast/db/nucleotide/IWGSC_CSS_ABD-TGAC_v1.fa" # blast contig file
	reference = "/Library/WebServer/Documents/blast/db/nucleotide/161010_Chinese_Spring_v1.0_pseudomolecules.fasta"
	
	# step 1:
	cmd1 = script_path + "parse_polymarker_input.py " + polymarker_input
	print "Step 1: Parse polymarker input command:\n", cmd1
	call(cmd1, shell=True)
	
	#step 2: blast
	cmd2 = 'blastn -task blastn -db ' + reference + ' -query for_blast.fa -outfmt "6 std qseq sseq slen" -num_threads 3 -out blast_out.txt'
	print "Step 2: Blast command:\n", cmd2
	call(cmd2, shell=True)
	
	# Step 3: parse the blast output file and output the homelog contigs and flanking ranges
	cmd3 = script_path + "getflanking.py " + polymarker_input + " blast_out.txt temp_range.txt " + genome_number
	print "Step 3: Get the flanking range command:\n", cmd3
	call(cmd3, shell=True)
	
	# step 4: split file for each marker
	cmd4 = "awk  '{ print $2,$3,$4 > \"temp_marker_\"$1.txt }' temp_range.txt"
	print "Step 4: Flanking range for each marker command:\n", cmd4
	call(cmd4, shell=True)
	
	# step 5: get flanking sequences for each file
	cmd5 = "find . -iname \"temp_marker*\" | xargs -n1 basename | xargs -I {} sh -c 'blastdbcmd -entry_batch {} -db " + reference + " > flanking_{}.fa'"
	print "Step 5: Get flanking sequences for each marker command:\n", cmd5
	call(cmd5, shell=True)
	
	# step 6: get kasp
	cmd6 = script_path + "getkasp3.py 1" # add blast option
	print "Step 6: Get KASP primers for each marker command:\n", cmd6
	call(cmd6, shell=True)
	
	# step 7: parse the for_polymarker input and create the input for SNP2CAPS
	cmd7 = script_path + "parse_polymarker_input_for_CAPS.py " + polymarker_input
	print "Step 7: Create input file for SNP2CAPS:\n", cmd7
	call(cmd7, shell=True)
	
	# step 8: run SNP2CAPS script to find all potential Restriciton enzymes
	cmd8 = script_path + "SNP2CAPS.pl for_SNP2CAPS.fa " + script_path + "REgcg.txt EcoRV,BamHI > CAPS_output.txt"
	print "Step 8: run SNP2CAPS script to find all potential Restriciton enzymes:\n", cmd8
	call(cmd8, shell=True)
	
	# step 9: get CAPS markers
	cmd9 = script_path + "getCAPS.py 1" # add blast option
	print "Step 9: Get CAPS and dCAPS primers for each marker command:\n", cmd9
	call(cmd9, shell=True)
	
	# step 10: concatenate output files
	caps_files = glob("CAPS_output/selected_CAPS_primers*")
	kasp_files = glob("KASP_output/selected_KASP_primers*")
	alignment_files = glob("alignment_raw_*")
	with open("All_alignment_raw.fa", "w") as outfile:
		for f in alignment_files:
			with open(f) as infile:
				outfile.write(f + "\n")
				outfile.write(infile.read())
				outfile.write("\n\n")
	outfile.close()
	print "all output files are: ", caps_files
	cmd10 = "cat CAPS_output/selected_CAPS_primers* > Potential_CAPS_primers.tsv"
	cmd11 = "cat KASP_output/selected_KASP_primers* > Potential_KASP_primers.tsv"
	cmd12 = "cat alignment_raw_* > All_alignment_raw.fa"
	print "Concatenate all output files to single files\n", cmd10, "\n", cmd11, "\n", cmd12
	call(cmd10, shell=True)
	call(cmd11, shell=True)
	#call(cmd12, shell=True)
	
	print "\n\n\n KASP primers have been designed successfully!\n Check files beginning with 'select_primer' and CAPS_output.txt"
	return 0

if __name__ == '__main__':
	import sys, os
	from subprocess import call
	sys.exit(main(sys.argv))
