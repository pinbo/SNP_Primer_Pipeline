#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  check_hits_and_extract_seq
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

# check whether some sequences have more than 10 hits, which might be a repetitive elements and should be avoided"
# no input, will find temp_marker* files in the folder

### Imported
from subprocess import call
import sys, os
from glob import glob

# check hits number
def check_hit(infile):
	num_lines = sum(1 for line in open(infile) if line.strip())
	if num_lines > 10:
		print infile, "has more than 10 hits. Please check!!!\n"
		os.rename(infile, "problem_" + infile)
		return 0
	return 1

def extract_seq(range_file, reference):
	flanking_file = "flanking_" + range_file + ".fa"
	cmd = "blastdbcmd -entry_batch " + range_file + " -db " + reference + " > " + flanking_file
	print cmd
	call(cmd, shell=True)

def main():
	# get all the raw sequences
	raw = glob("temp_marker*") # All file names start from "flanking"
	raw.sort()
	reference = "/Library/WebServer/Documents/blast/db/nucleotide/161010_Chinese_Spring_v1.0_pseudomolecules.fasta"
	for ff in raw:
		if check_hit(ff):
			extract_seq(ff, reference)
	return 0

if __name__ == '__main__':
	main()
