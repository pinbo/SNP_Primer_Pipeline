#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  getcaps
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

# getcaps3.py use files with names starting with "flanking*" from step 5 of the pipeline
# getcaps3.py now uses all the variation sites that can differ all homeologs

# change the input wildcard and the paths of primer3 and muscle accordingly.
# NOTES: the output primer pair includes the common primer and the primer with SNP A or T, you need to change the 3' nucleartide to get the primer for the other SNP.

# 2018-03-29: from getCAPS.py, but modified to accept user input multiple sequences.
# two alleles now support indels or long haplotypes

### Imported
from subprocess import call
import getopt, sys, os, re
from glob import glob
#########################
## arguments
# enzyme_check_only.py -a GCCATACATGTAAAATTTAAGGAAAGACTGGGAAAGTTCTCTAAATGTTTATTTTCTAGT -b GCCATACATGTAAAATTTAAGGAAAGACTGGaAAAGTTCTCTAAATGTTTATTTTCTAGT -p 1900
usage = '''
Please make sure both sequence A and sequence B are given.

enzyme_check_only.py -a seqA -b seqB -p max_enzyme_price
-h: this message
-a: sequence A
-b: sequence B
-p: max_enzyme_price: max enzyme price ($ per 1000 U) (default is 200)

'''

# parameter or file names that need to be changed
max_price = 200
seqA = ""
seqB = ""


try:
	opts, args = getopt.getopt(sys.argv[1:], "a:b:p:h", ["help"])
except getopt.GetoptError as err:
	# print help information and exit:
	print str(err)  # will print something like "option -a not recognized"
	print usage
	sys.exit(2)
for o, a in opts:
	if o == "-a":
		seqA = a
	elif o in ("-h", "--help"):
		print usage
		sys.exit()
	elif o in ("-p"):
		max_price = int(a)
	elif o in ("-b"):
		seqB = a
	else:
		assert False, "unhandled option"
print "Options done"

if not seqA or not seqB:
	print usage
	sys.exit(1)


# code for my reference
iupac = {"R": "AG", "Y": "TC", "S": "GC", "W": "AT", "K": "TG", "M": "AC"}

# classes
class Restriction_Enzyme(object):
	def __init__(self, name, seq):
		self.name = name
		self.seq = seq.lower()
		self.length = len(seq)
		self.template_seq = ""
		self.primer_end_pos = [] # a list of end positions
		self.caps = "No"
		self.dcaps = "No"
		self.allpos = [] # all the match positions in the template
		self.change_pos = None
		self.price = int(name.split(',')[-1])
		self.primer_direction = "" # to use the end positions as left or right primer

# sequence class with allele sequence and position
class caps_seq(object):
	"""sequence with allele sequence and positions"""
	def __init__(self):
		self.refseq = "" # reference template
		self.ref = Primers() # ref allele: SNP_A
		self.alt = "" # alt allele: SNP_B
		self.compl_end = "NA"
		self.penalty = "NA"
		self.product_size = 0

# calculate GC content of a sequence
def Calc_GC(seq):
	t = 0.0 # float
	for a in seq:
		if a=='C' or a=='G':
			t += 1
	return t / len(seq) * 100

def parse_RE_file(RE_file):
	REs = {}
	with open(RE_file) as file_one:
		for line in file_one:
			enzyme, seq = line.rstrip().split("\t")
			REs[enzyme] = Restriction_Enzyme(enzyme, seq)
	return REs

def ReverseComplement(seq):
	s1 = "BDHKMNRSVWYATGCbdhkmnrsvwyatgc"
	s2 = "VHDMKNYSBWRTACGvhdmknysbwrtacg"
	seq_dict = {s1[i]:s2[i] for i in range(len(s1))}
	return "".join([seq_dict[base] for base in reversed(seq)])

def string_dif(s1, s2): # two strings with equal length
	return [i for i in xrange(len(s1)) if s1[i] != s2[i]]

def dif_region(s1, s2): # two strings do not need to have the same length
	s1r = s1[::-1] # reverse the string
	s2r = s2[::-1] # reverse the string
	L1 = [i for i in xrange(min([len(s1),len(s2)])) if s1[i] != s2[i]] # forward
	L2 = [i for i in xrange(min([len(s1),len(s2)])) if s1r[i] != s2r[i]] # reverse
	return [L1[0], len(s1) - L2[0] - 1] # differ positions for when counting forward and reversely on the template (the two numbers are both in forward direction)

# return the range of differences of two strings
def dif_region2(s1, s2): # two strings do not need to have the same length
	s1r = s1[::-1] # reverse the string
	s2r = s2[::-1] # reverse the string
	L1 = [i for i in xrange(min([len(s1),len(s2)])) if s1[i] != s2[i]] # forward
	L2 = [i for i in xrange(min([len(s1),len(s2)])) if s1r[i] != s2r[i]] # reverse
	return [L1[0], -L2[0]]

def seq2pattern(seq):
	iupac = {
		"B": "[CGT]",
		"D": "[AGT]",
		"H": "[ACT]",
		"K": "[GT]",
		"M": "[AC]",
		"N": "[ACGT]",
		"R": "[AG]",
		"S": "[CG]",
		"V": "[ACG]",
		"W": "[AT]",
		"Y": "[CT]"
	}
	seq = seq.upper()
	seq2 = ""
	for i in seq:
		if i in "ATGC":
			seq2 += i
		else:
			seq2 += iupac[i]
	#seq2 = '(?=(' +  seq2 + '))'
	return seq2.lower()

def check_pattern(enzyme, wild_seq, mut_seq): # check whether enzyme can match wild_seq after 1 base change but not mut_seq
	wild_seq = wild_seq.lower()
	mut_seq = mut_seq.lower()
	pos_L, pos_R = dif_region(wild_seq, mut_seq) # differnce region borders in the wild_seq
	#print "pos_L, pos_R ", pos_L, pos_R
	enzyme_name = enzyme.name
	enzyme_seq = enzyme.seq.strip("n") # some enzyme has sequence beginning and ending with n, such as TspRI,72, nncastgnn
	for i in range(len(enzyme_seq)):
		ss = seq2pattern(enzyme_seq[0:i]) + "[atgc]" + seq2pattern(enzyme_seq[i+1:]) # regular expression
		if find_substring(ss, wild_seq) != find_substring(ss, mut_seq): # even they are the same length, if the start position is different, it is still okay
		#print "Enzyme, Enzyme seq, pattern ", enzyme_name, enzyme_seq, ss
			for m in re.finditer(ss, wild_seq): # iterate all the matching places
				change_pos = m.start() + i # which was changed
				# if the left first different nt between wt and mut is in the enzyme recognization site and the change position is more than 1 nt from it.
				if pos_L in range(m.start(), m.end()) and pos_L - change_pos > 1:
					enzyme.primer_direction = "left" # use as left primer end positions
					enzyme.primer_end_pos += range(change_pos + 1, pos_L)
					change_pos = m.start() + i # which was changed
					print "One nt can be changed to fit enzyme", enzyme.name
					enzyme.dcaps = "Yes"
					enzyme.template_seq = wild_seq[:change_pos] + enzyme_seq[i].upper() + wild_seq[change_pos+1:]
					enzyme.change_pos = change_pos
					print "change position and primer end postions are ", change_pos, enzyme.primer_end_pos
					# break the loop if one change postion is found, because I do not think we need many
					break
				
				if pos_R in range(m.start(), m.end()) and change_pos - pos_R > 1:
					enzyme.primer_direction = "right" # use as right primer end positions
					enzyme.primer_end_pos += range(pos_R + 1, change_pos)
					change_pos = m.start() + i # which was changed
					print "One nt can be changed to fit enzyme", enzyme.name
					enzyme.dcaps = "Yes"
					enzyme.template_seq = wild_seq[:change_pos] + enzyme_seq[i].upper() + wild_seq[change_pos+1:]
					enzyme.change_pos = change_pos
					print "change position and primer end postions are ", change_pos, enzyme.primer_end_pos
					break
		# break the loop if dcaps found
		if enzyme.dcaps == "Yes":
			break
		#enzyme.primer_end_pos = list(set(enzyme.primer_end_pos)) # not necessary since I already made the break if one change positon were found
	return enzyme

def find_substring(substring, string): # find all the starting index of a substring
	return [m.start() for m in re.finditer(substring, string)]

# test whether an enzyme can be modifed to fit a dCAPS
def test_enzyme(enzyme, wild_seq, mut_seq): # enzyme is an Restriction_Enzyme object
	enzyme_seq = enzyme.seq
	enzyme_seq_RC = ReverseComplement(enzyme_seq) # 1
	enzyme_name = enzyme.name
	wild_seq = wild_seq.lower()
	mut_seq = mut_seq.lower()
	wild_allpos = find_substring(seq2pattern(enzyme_seq), wild_seq)
	mut_allpos = find_substring(seq2pattern(enzyme_seq), mut_seq)
	wild_allpos += find_substring(seq2pattern(enzyme_seq_RC), wild_seq) # also check reverse complement sequences of enzyme
	mut_allpos += find_substring(seq2pattern(enzyme_seq_RC), mut_seq)
	enzyme.allpos = list(set(wild_allpos))
	if len(wild_allpos) != len(mut_allpos): # snp cause digestion difference
		enzyme.caps = "Yes"
		enzyme.template_seq = wild_seq
		print "CAPS found with enzyme ", enzyme_name
		return enzyme
	# else no caps found, check dcaps
	#snp_pos = string_dif(wild_seq, mut_seq)[0] # snp position in the template
	#print "snp_pos is", snp_pos
	enzyme = check_pattern(enzyme, wild_seq, mut_seq)
	if enzyme.dcaps != "Yes":
		enzyme = check_pattern(enzyme, mut_seq, wild_seq)
	if enzyme.dcaps != "Yes" and enzyme_seq_RC != enzyme_seq:
		enzyme.seq = enzyme_seq_RC
		enzyme = check_pattern(enzyme, wild_seq, mut_seq)
	if enzyme.dcaps != "Yes" and enzyme_seq_RC != enzyme_seq:
		enzyme = check_pattern(enzyme, mut_seq, wild_seq)
	return enzyme

# function to count mismtaches
def mismatchn (s1, s2):
	return sum(c1!=c2 for c1,c2 in zip(s1,s2))

# function to extract sequences from a fasta file 
def get_fasta(infile):
	fasta = {} # dictionary for alignment
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line: # skip blank lines
				if line.startswith(">"):
					sequence_name = line.lstrip("> ").split()[0] # left strip > or space, so " > abc edf" will be "abc edf", then split by space to get "abc"
					fasta[sequence_name] = ""
				else:
					fasta[sequence_name] += line.replace(" ", "") # remove spaces in case
	return fasta

# in case multiple hit in the same chromosome in the psudomolecule
# the input fasta file are in order from the blast file
def get_fasta2(infile, target_chrom):
	fasta = {} # dictionary for alignment
	target = ""
	non_target_list = []
	n = 0 # add a number in the chromosome name
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line.startswith(">"):
				sequence_name = line.split()[0].lstrip(">")
				sequence_name += "-" + str(n)
				n += 1
				if not target and target_chrom in sequence_name:
					target = sequence_name
				else:
					non_target_list.append(sequence_name)
			else:
				fasta.setdefault(sequence_name, "")
				fasta[sequence_name] += line.rstrip()
			
	return fasta, target, non_target_list

# parse the flanking sequence input

def caps(seqA, seqB, max_price): # two alleles now support indels or long haplotypes, SNP_A is the template allele (or reference allel)
	snpname = "SNP"
	getcaps_path = os.path.dirname(os.path.realpath(__file__))
	out = "available_enzymes.txt"
	print "Output selected CAPS file name is: ", out

	################## Get CAPS information
	# step 1: read the enzyme file
	RE_file = getcaps_path + "/NEB_parsed_REs.txt" # this one removed some duplicated cuttings
	REs = parse_RE_file(RE_file) # get the list of restriction enzymes
	# step 2: get the list of enzymes that can be used for caps or dcaps
	wild_seq = seqA
	mut_seq = seqB
	pos_L, pos_R = dif_region2(wild_seq, mut_seq) # differnce region borders in the wild_seq, use this instead of snp_pos to decide use as left end or right end 
	var = "[" + wild_seq[pos_L:pos_R] + "/" + mut_seq[pos_L:pos_R] + "]" # variations, for examploe "[AT/TAG]"
	caps_list = []
	dcaps_list = []
	for k in REs:
		enzyme = REs[k]
		if enzyme.price > max_price:
			continue
		enzyme = test_enzyme(enzyme, wild_seq, mut_seq)
		if enzyme.caps == "Yes":
			caps_list.append(enzyme)
		elif enzyme.dcaps == "Yes":
			dcaps_list.append(enzyme)
	print "caps_list is ", [x.name for x in caps_list]
	print "dcaps_list is ", [x.name for x in dcaps_list]

	outfile = open(out, 'w')
	#outfile.write("CAPS cut information\n") # change to 1 based
	outfile.write("Sequence A:\t" + seqA + "\n")
	outfile.write("Sequence B:\t" + seqB + "\n\n")
	outfile.write("Enzyme\tEnzyme_seq\tChange_pos\tOther_cut_pos\tTemplate\n")
	for enzyme in dcaps_list + caps_list:
		seq = enzyme.template_seq
		seq = seq[:pos_L] + var + seq[pos_R:]
		outfile.write(enzyme.name + "\t" + enzyme.seq + "\t" + str(enzyme.change_pos) + "\t" + ", ".join([str(x + 1) for x in enzyme.allpos]) + "\t" + seq + "\n")
	# close outfile
	outfile.close()
	return 0


# Find CAPS/dCAPS enzymes

caps(seqA, seqB, max_price)

