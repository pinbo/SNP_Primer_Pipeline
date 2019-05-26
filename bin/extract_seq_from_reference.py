#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  get_flanking_for_variaiton_sites.py
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
## Update on 2019-05-26 because BLAST+ update to version 2.9.0+
## I do not need add the range to the output, the default output already has it.

from subprocess import call
import sys



def main(args):
	reference_list = ["/Library/WebServer/Documents/blast/db/nucleotide/IWGSC_v2_ChrU.fa", 
	"/Library/WebServer/Documents/blast/db/nucleotide/IWGSC_CSS_AB-TGAC_UCW_v1.fa",
	"/Library/WebServer/Documents/blast/db/nucleotide/161010_Chinese_Spring_v1.0_pseudomolecules.fasta",
	"/Volumes/DATA2/databases/ncbi/nucleotide/Triticum_aestivum.TGACv1.dna.toplevel.fa",
	"/Users/galaxy/blastdb/IWGSC_v1.1_HC_20170706_cds.fasta"]
	infile = args[1]
	reference = reference_list[int(args[2]) - 1] # reference 1 or 2
	# step 1: get the flanking sequences
	flanking_file = "temp_flanking_seq.fa"
	cmd = "blastdbcmd -entry_batch " + infile + " -db " + reference + " > " + flanking_file
	call(cmd, shell=True)
	return 0

if __name__ == '__main__':
	import sys
	sys.exit(main(sys.argv))
