#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  getkasp
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

# example
# ./getflanking.py for_polymarker.csv blast_out.tsv temp_range.txt 3

### Imported
import sys

# i'll assume it's tab-delimited...
# and est is the query.
polymarker_input = sys.argv[1]
blast_file = sys.argv[2] # this is a special blast file with two columns in the last: qseq and sseq (query sequence and subject sequenc)
outfile = sys.argv[3]
genome_number =  int(sys.argv[4])
if genome_number not in [1, 2, 3]:
	sys.exit("Genome number need to be either 1, 2, or 3")

genomes = "ABD"
genomes = genomes[:genome_number] + "n" # final genomes + chrUn

# get snp position
snp_pos = {}
for line in open(polymarker_input):
	snp, chrom, seq = line.strip().replace(" ","").split(",") # in case there are spaces
	snp = snp.replace("_", "-") # in case there is already "_" in the snp name
	seq = seq.strip() # in case there are spaces in the input file
	snp_pos[snp] = seq.find("[") + 1

# function to find all the indexes (1-based) of a character
# for looking for gaps
def find(s, ch):
    return [i + 1 for i, ltr in enumerate(s) if ltr == ch]

# function to 
# take 500 bp up/down-stream.
# increased this value from 250 to 500 for dCAPS marker design later.
xstream = 500

flanking = {} # flanking information
snpinfo = {}

# two list below to store the matched subjects in order
# so I know which chromosme subject is the best hit
snp_list = [] # snp name
range_list = [] # range list for each subject

# blast fields
# IWB50236_7A_R	IWGSC_CSS_7DS_scaff_3919748	98.718	78	1	0	24	101	4891	4968	1.55e-30	138	CTCATCAAATGATTCAAAAATATCGATRCTTGGCTGGTGTATCGTGCAGACGACAGTTCGTCCGGTATCAACAGCATT	CTCATCAAATGATTCAAAAATATCGATGCTTGGCTGGTGTATCGTGCAGACGACAGTTCGTCCGGTATCAACAGCATT	5924
# Fields: 
# 1: query id, subject id, % identity, alignment length, mismatches, gap opens, 
# 7: q. start, q. end, s. start, s. end, evalue, bit score
# 13: q. sequence, s. sequence, s. length
snp_size_list = [] # max alignment length for each snp
for line in open(blast_file):
	if line.startswith('#'):
		continue
	fields = line.split("\t")
	query, subject = fields[:2]
	snp, qchrom = query.split("_")[0:2] # snp name, query chromosome name
	qchrom = qchrom[0:2] # no arm for pseudomolecule blast
	#schrom = subject.split("_")[2] # subject chromosome name with arm
	schrom = subject[-2:] # chr6A as in the pseudomolecule. No chromosome arm
	#print qchrom, schrom
	if schrom[1] not in genomes:
		continue
	#pct_identity = float(fields[2]) # big gaps cause low identity
	pct_identity = 100 - (float(fields[4]) + float(fields[5])) / float(fields[3]) * 100 # to avoid big gaps
	align_length = int(fields[3])
	if snp not in snp_size_list:
		snp_size_list.append(snp)
		min_align = max(50, align_length * 0.9) # to filter out those not very good alignment, since I will blast anyway later.
	# only get min-identity 90% and at least 50 bp alignment
	#print "snp, min_align", snp, min_align
	if pct_identity > 88 and align_length > min_align:
		qstart, qstop, sstart, sstop = [int(x) for x in fields[6:10]]
		qseq, sseq = fields[12:14]
		slen = int(fields[14]) # subject length
		if snp_pos[snp] < qstart or snp_pos[snp] > qstop: # if snp is not in the alignment
			continue
		# gap positions for qseq and sseq
		# gap will mess up the SNP position, so need to count
		qgap = find(qseq, "-")
		sgap = find(sseq, "-")
		
		temp = snp_pos[snp] - qstart # distance from qstart to the snp if there is no gap
		nqgap = nsgap = 0
		for n in qgap:
			if n < temp:
				nqgap += 1
				temp += 1

		for n in sgap:
			if n < temp: # distance from the begening to the snp
				nsgap += 1
		
		pos =  sstart + (temp - nsgap) # snp position in the subject
		strand = "plus"
		if sstart > sstop:
			strand = "minus"
			pos =  sstart - (temp - nsgap) # snp position in the subject
		
		up = max(1, pos - xstream)
		down = min(slen, pos + xstream)
		
		pos2 = pos - up + 1 # snp position in the extracted flanking sequences
		if sstart > sstop: # if minus strand
			pos2 = down - pos + 1 # snp position in the extracted flanking sequences
	
		if qchrom == schrom:
			snpinfo[query] = query + "_" + str(pos2)
		#flanking[query + "-" + subject] = "\t".join([subject, str(up) + "-" + str(down), strand])
		#flanking[query + "-" + subject + "-" + str(sstart)] = "\t".join([subject, str(up) + "-" + str(down), strand])
		snp_list.append(query)
		range_list.append("\t".join([subject, str(up) + "-" + str(down), strand]))

# find out which has too many hits
max_hit = 6
from collections import Counter
ct = Counter(snp_list) # count of each snp hits
for i in ct:
	print i, "has hits", ct[i]

# output
out = open(outfile, "w")

for i in range(len(snp_list)):
	snp = snp_list[i]
	if ct[snp] > max_hit:
		#print snp, ct[snp]
		continue
	rg = range_list[i]
	out.write(snpinfo[snp] + "\t" + rg + "\n")

#for k, v in flanking.items():
#	k2 = k.split("-")[0] # key for snpinfo
#	out.write(snpinfo[k2] + "\t" + v + "\n")
	
out.close()
