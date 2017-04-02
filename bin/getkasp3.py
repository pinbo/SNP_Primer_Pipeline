#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  getkasp
#
#  Copyright 2016 Junli Zhang <zhjl86@gmail.com>
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

# getkasp3.py use files with names starting with "flanking*" from step 5 of the pipeline
# getkasp3.py now uses all the variation sites that can differ all homeologs

# change the input wildcard and the paths of primer3 and muscle accordingly.
# NOTES: the output primer pair includes the common primer and the primer with SNP A or T, you need to change the 3' nucleartide to get the primer for the other SNP.

### Imported
from subprocess import call
import getopt, sys, os, re

#########################

from glob import glob

# get all the raw sequences
raw = glob("flanking_temp_marker*") # All file names start from "flanking"
raw.sort()

iupac = {"R": "AG", "Y": "TC", "S": "GC", "W": "AT", "K": "TG", "M": "AC"}

def kasp(seqfile):
	#flanking_temp_marker_IWB1855_7A_R_251.fa
	snpname, chrom, allele, pos =re.split("_|\.", seqfile)[3:7]
	getkasp_path = os.path.dirname(os.path.realpath(__file__))
	out = "selected_primers_" + snpname + ".txt"
	target = "2AS" # target sequence name
	product_min = 50
	product_max = 250
	alt_allele = iupac[allele][0] # choose A or T
	
	# software path
	muscle_path = "muscle"
	primer3_path = "primer3_core"

	########################
	# STEP 0: create alignment file and primer3output file
	RawAlignFile = "alignment_raw_" + snpname + ".fa"
	alignmentcmd = muscle_path + " -in " + seqfile + " -out " + RawAlignFile + " -quiet"
	print "Alignment command: ", alignmentcmd
	call(alignmentcmd, shell=True)
	
	# Primer3 input and output
	#outfile = open(out, 'w') # output file
	snp_site = int(pos) - 1 # 0-based
	########################
	# read alignment file
	fasta = {} # dictionary for alignment

	with open(RawAlignFile) as file_one:
		for line in file_one:
			if line.startswith(">"):
				sequence_name = line.rstrip().lstrip(">")
			else:
				fasta.setdefault(sequence_name, "")
				fasta[sequence_name] += line.rstrip()


	# get the variaiton site among sequences
	ids = [] # all other sequence names
	for kk in fasta.keys():
		if chrom in kk:
			target = kk
		else:
			ids.append(kk)
			
	print "The target: ", target
	print "The other groups: ", ids

	alignlen = len(fasta[target])
	print "Alignment length: ", alignlen
	
	seq_template = fasta[target].replace("-","") # remove all gaps

	variation = [] # variation sites that can differ ALL
	variation2 = [] # variation sites that can differ at least 2 homeologs
	
	# calculate gap number on the left and right
	gap_left_target = len(fasta[target]) - len(fasta[target].lstrip('-'))
	gap_left = max([len(v) - len(v.lstrip('-')) for k, v in fasta.items()])
	gap_right = min([len(v.rstrip('-')) for k, v in fasta.items()])
	print "gap_left_target, gap_left and gap_right: ", gap_left_target, gap_left, gap_right
	
	ngap = 0 # gap number
	diffarray = {} # a list of 0 or 1: the same as or different from the site in each sequences of ids
	#for i in range(alignlen):
	for i in range(gap_left, gap_right):
		j = i - gap_left_target # index number in target sequence without leading gaps
		nd = 0 # number of difference
		da = [0] * len(ids) # differ array
		m = 0 # counter of next loop
		nl1 = len(fasta[target][gap_left:i]) - len(fasta[target][gap_left:i].rstrip("-")) # number of gaps on the left of current site for the target
		nr1 = len(fasta[target][i:gap_right]) - len(fasta[target][i:gap_right].lstrip("-")) + 1
		if j - ngap < snp_site:
			pos_template = j - ngap + nl1 # position in the target template (no gaps)
		else:
			pos_template = j - ngap - nr1
		if (pos_template < snp_site and nl1 > 10) or (pos_template >= snp_site and nr1 > 10):
			if fasta[target][i] == "-":
				ngap += 1
			continue # go to next loop
		for k in ids:
			# below nl1 to nr2 is for handling gaps
			#nl1 = i - len(fasta[target][:i].rstrip("-")) # number of gaps on the left of current site for the target
			nl2 = len(fasta[k][gap_left:i]) - len(fasta[k][gap_left:i].rstrip("-")) # number of gaps on the left of current site for the homeolog
			#nr1 = len(fasta[target][i:]) - len(fasta[target][i:].lstrip("-")) + 1
			nr2 = len(fasta[k][i:gap_right]) - len(fasta[k][i:gap_right].lstrip("-")) + 1
			if j - ngap < snp_site: # variation site should be left end
				b1 = fasta[target][i:gap_right].replace("-","")[nl1] # target non-gap base
				b2 = fasta[k][i:gap_right].replace("-","")[nl2] # homeolog non-gap bas
			else:
				b1 = fasta[target][gap_left:i].replace("-","")[-nr1] # target non-gap base
				b2 = fasta[k][gap_left:i].replace("-","")[-nr2] # homeolog non-gap bas
			if b1 != b2:
				nd+=1
				da[m] = 1 # m sequence has variation from target
			m += 1
		
		#diffarray[j - ngap] = da # variation infor for each site
		if pos_template not in diffarray:
			diffarray[pos_template] = da
		else:
			diffarray[pos_template] = [sum(x) for x in zip(da,diffarray[pos_template])]
		if nd == len(ids): # different from all other sequences
			if pos_template not in variation: # in case multiple gaps
				variation.append(pos_template)
		if nd > 0: # different from at least 1 other sequence
			if pos_template not in variation2: # in case multiple gaps
				variation2.append(pos_template)
		# if site i is a gap, gap # + 1
		if fasta[target][i] == "-":
			ngap += 1
	
	print "Sites that can differ all\n", variation
	print "\nSites that can differ at least 1\n", variation2
	#############
	# loop to write primer3 input for each variation site
	# primer3 inputfile
	primer3input = "primer3.input." + snpname
	p3input = open(primer3input, 'w')

	#seq_template = fasta[target].replace("-","") # remove all gaps

	# because A and T give lower Tm, so use them as template
	if alt_allele in "ATat":
		seq_template = seq_template[:snp_site] +  alt_allele + seq_template[snp_site + 1:]

	for i in variation2:
		if i < snp_site:
			left_end = i
			right_end = snp_site
		else:
			left_end = snp_site
			right_end = i
		if right_end - left_end > product_max - 35: # suppose both primers are 18 bp
			continue
		settings = "PRIMER_TASK=generic" + "\n" + \
		"SEQUENCE_ID=" + target + "-" + str(i) + "\n" + \
		"SEQUENCE_TEMPLATE=" + seq_template + "\n" + \
		"PRIMER_PRODUCT_SIZE_RANGE=50-250" + "\n" + \
		"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + getkasp_path + "/primer3_config/"  + "\n" + \
		"PRIMER_MAX_SIZE=25" + "\n" + \
		"PRIMER_PAIR_MAX_DIFF_TM=6.0" + "\n" + \
		"PRIMER_FIRST_BASE_INDEX=1" + "\n" + \
		"PRIMER_LIBERAL_BASE=1" + "\n" + \
		"PRIMER_NUM_RETURN=5"  + "\n" + \
		"PRIMER_EXPLAIN_FLAG=1"  + "\n" + \
		"SEQUENCE_FORCE_LEFT_END=" + str(left_end + 1) + "\n" + \
		"SEQUENCE_FORCE_RIGHT_END=" + str(right_end + 1) + "\n" + \
		"="
		p3input.write(settings + "\n")

	p3input.close()

	# primer3 output file
	primer3output = "primer3.output." + snpname
	p3cmd = primer3_path + " -default_version=2 -output=" + primer3output + " " + primer3input
	print "Primer3 command 1st time: ", p3cmd
	call(p3cmd, shell=True)

	##########################
	# function to get reverse complement
	def ReverseComplement(seq):
		# too lazy to construct the dictionary manually, use a dict comprehension
		seq1 = 'ATCGTAGCatcgtagc'
		seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
		return "".join([seq_dict[base] for base in reversed(seq)])

	# classes
	class Primers(object):
		"""A primer set designed by Primer3"""
		def __init__(self):
			self.start = 0
			self.end = 0
			self.length = 0
			self.tm = 0.0
			self.gc = 0.0
			self.anys = 0.0
			self.three = 0.0
			self.hairpin = 0.0
			self.end_stability = 0.0
			self.seq = ""
			self.difthreeall = "NO" # whether 3' site can differ all
		def formatprimer(self):
			formatout = "\t".join(str(x) for x in [self.start, self.end, self.length, self.tm, self.gc, self.anys, self.three, self.end_stability, self.hairpin, self.seq, ReverseComplement(self.seq), self.difthreeall])
			return(formatout)

	class PrimerPair(object):
		"""A pair of primers designed by Primer3"""
		def __init__(self):
			self.left = Primers()
			self.right = Primers()
			self.compl_any = "NA"
			self.compl_end = "NA"
			self.penalty = "NA"
			self.product_size = 0



	### parse primer 3 primer check output
	primerpairs = {} # sequence ID is the key

	with open(primer3output) as infile:
		for line in infile:
			line = line.strip()
			if "SEQUENCE_ID" in line:
				seqid = line.split("=")[1]
				primerpairs[seqid] = PrimerPair()
			elif "PRIMER_PAIR_0_PENALTY" in line:
				primerpairs[seqid].penalty = line.split("=")[1]
			elif "PRIMER_PAIR_0_COMPL_ANY" in line:
				primerpairs[seqid].compl_any = line.split("=")[1]
			elif "PRIMER_PAIR_0_COMPL_END" in line:
				primerpairs[seqid].compl_end = line.split("=")[1]
			elif "PRIMER_PAIR_0_PRODUCT_SIZE" in line:
				primerpairs[seqid].product_size = int(line.split("=")[1])
			elif "PRIMER_LEFT_0_SEQUENCE" in line:
				primerpairs[seqid].left.seq = line.split("=")[1]
			elif "PRIMER_LEFT_0=" in line:
				primerpairs[seqid].left.start = int(line.split("=")[1].split(",")[0])
				primerpairs[seqid].left.length = int(line.split("=")[1].split(",")[1])
				primerpairs[seqid].left.end = primerpairs[seqid].left.start + primerpairs[seqid].left.length - 1
			elif "PRIMER_LEFT_0_TM" in line:
				primerpairs[seqid].left.tm = float(line.split("=")[1])
			elif "PRIMER_LEFT_0_GC_PERCENT" in line:
				primerpairs[seqid].left.gc = float(line.split("=")[1])
			elif "PRIMER_LEFT_0_SELF_ANY_TH" in line:
				primerpairs[seqid].left.anys = float(line.split("=")[1])
			elif "PRIMER_LEFT_0_SELF_END_TH" in line:
				primerpairs[seqid].left.three = float(line.split("=")[1])
			elif "PRIMER_LEFT_0_HAIRPIN_TH" in line:
				primerpairs[seqid].left.hairpin = float(line.split("=")[1])
			elif "PRIMER_LEFT_0_END_STABILITY" in line:
				primerpairs[seqid].left.end_stability = float(line.split("=")[1])
			elif "PRIMER_RIGHT_0_SEQUENCE" in line:
				primerpairs[seqid].right.seq = line.split("=")[1]
			elif "PRIMER_RIGHT_0=" in line:
				primerpairs[seqid].right.start = int(line.split("=")[1].split(",")[0])
				primerpairs[seqid].right.length = int(line.split("=")[1].split(",")[1])
				primerpairs[seqid].right.end = primerpairs[seqid].right.start - primerpairs[seqid].right.length + 1
			elif "PRIMER_RIGHT_0_TM" in line:
				primerpairs[seqid].right.tm = float(line.split("=")[1])
			elif "PRIMER_RIGHT_0_GC_PERCENT" in line:
				primerpairs[seqid].right.gc = float(line.split("=")[1])
			elif "PRIMER_RIGHT_0_SELF_ANY_TH" in line:
				primerpairs[seqid].right.anys = float(line.split("=")[1])
			elif "PRIMER_RIGHT_0_SELF_END_TH" in line:
				primerpairs[seqid].right.three = float(line.split("=")[1])
			elif "PRIMER_RIGHT_0_HAIRPIN_TH" in line:
				primerpairs[seqid].right.hairpin = float(line.split("=")[1])
			elif "PRIMER_RIGHT_0_END_STABILITY" in line:
				primerpairs[seqid].right.end_stability = float(line.split("=")[1])


	# write to file
	outfile = open(out, 'w')
	outfile.write("index\tproduct_size\ttype\tstart\tend\tlength\tTm\tGCcontent\tany\t3'\tend_stability\thairpin\tprimer_seq\tReverseComplement\t3'differall\tpenalty\tcompl_any\tcompl_end\n")

	#for pp in primerpairs:
	for i, pp in primerpairs.items():
		varsite = int(i.split("-")[-1]) # variation site
		if pp.product_size != 0:
			pl = pp.left
			pr = pp.right
			# check whether 3' can differ all
			if varsite in variation:
				pl.difthreeall = "YES"
				pr.difthreeall = "YES"
			if varsite < snp_site:
				pc = pl # pc is the common primer
				# rr: range to check; only check 10 bases from 3' end
				rr = range(pc.end - 10, pc.end) # pc.end is 1 based, so change to 0 based.
			else:
				pc = pr
				rr = range(pc.end -1, pc.end + 9)
			# sum of all the variation in each site
			aa = [sum(x) for x in zip(*(diffarray[k] for k in rr))]
			if min(aa) > 0: # if common primer can differ all
				outfile.write("\t".join([i, str(pp.product_size), "LEFT", pl.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end]) + "\n")
				outfile.write("\t".join([i, str(pp.product_size), "RIGHT", pr.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end]) + "\n")

	outfile.close()

# loop for all snp sequence files

for ff in raw:
	kasp(ff)
	
