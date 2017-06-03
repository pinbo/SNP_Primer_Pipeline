#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parse_polymarker_input.py
#  
#  Copyright 2017 Junli Zhang <junli@debian>
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
"""
Parse polymarker input and output a fasta file for blast.

./parse_polymarker_input.py for_polymarker.csv 
"""
import sys

iupac = {"[A/G]": "R", "[G/A]": "R", "[C/T]": "Y", "[T/C]": "Y", "[G/C]": "S", "[C/G]": "S", "[A/T]": "W", "[T/A]": "W", "[G/T]": "K", "[T/G]": "K", "[A/C]": "M", "[C/A]": "M"}


def main():
	polymarker_input = sys.argv[1]
	outfile = "for_blast.fa"
	out = open(outfile, "w")

	# get snp position
	for line in open(polymarker_input):
		line = line.strip()
		if not line:
			continue
		snpname, chrom, seq = line.replace(" ","").split(",")
		snpname = snpname.replace("_", "-") # in case there is already "_" in the snp name
		seq = seq.strip() # in case there is space in the input file
		pos = seq.find("[")
		snp = iupac[seq[pos:pos+5]]
		seq2 = seq[:pos] + snp + seq[pos+5:]
		out.write(">" + snpname + "_" + chrom + "_" + snp + "\n" + seq2 + "\n")
	out.close()

	return 0

if __name__ == '__main__':
	main()

