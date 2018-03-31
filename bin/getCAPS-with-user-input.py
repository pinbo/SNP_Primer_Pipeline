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
# getCAPS-with-user-input seqfile target-name alleleA alleleB position price
usage = '''
Please make sure all the parameters are given.

getCAPS-with-user-input -i seqfile -t target-seq-name --ref Ref_allele --alt Alt_allele -p SNP/indel-position -m max_enzyme_price 
						--minTm min_Tm --maxTmm max_Tm --minSize min_primer_size --maxSize max_primer_size
-h or --help: this message
-i: iseqfile: has all the homeologs for designing genomic specific primers
-t: target-seq-name: the name of the target sequence
-p: SNP/indel-position: the postion of the mutation (1 based)
-m: max_enzyme_price: max enzyme price ($ per 1000 U) (200 for example)
-r or --ref: Ref_allele: the SNP or indel allele in the seqfile template (treated as reference sequence)
-a or --alt: Alt_allele: the altanative allele of the SNP or indel
--mintm <primer min Tm, default 58>
--maxtm <primer max Tm, default 62>
--minsize <primer min size, default 18>
--maxsize <primer max size, default 30>

'''

# parameter or file names that need to be changed

blast = 0 # 0 or 1, whether to blast
seqfile = ""
target = ""
SNP_A = ""
SNP_B = ""
snp_pos = 0
minTm = 58
maxTm = 62
minSize = 18
maxSize = 30
max_price = 200

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:p:m:t:r:a:h", ["help", "mintm=", "maxtm=", "minsize=", "maxsize=", "ref=", "alt="])
except getopt.GetoptError as err:
	# print help information and exit:
	print str(err)  # will print something like "option -a not recognized"
	print usage
	sys.exit(2)
for o, a in opts:
	if o == "-i":
		seqfile = a
	elif o in ("-h", "--help"):
		print usage
		sys.exit()
	elif o in ("-p"):
		snp_pos = int(a)
	elif o in ("-m"):
		max_price = int(a)
	elif o in ("-t"):
		target = a
	elif o in ("-r", "--ref"):
		SNP_A = a
	elif o in ("-a", "--alt"):
		SNP_B = a
	elif o in ("--mintm"):
		minTm = int(a)
	elif o in ("--mintm"):
		minTm = int(a)
	elif o in ("--maxtm"):
		maxTm = int(a)
	elif o in ("--minsize"):
		minSize = int(a)
	elif o in ("--maxsize"):
		maxSize = int(a)
	else:
		assert False, "unhandled option"
print "Options done"

if not target or not seqfile or not SNP_A or not SNP_B or not snp_pos:
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

class Primers(object):
	"""A primer set designed by Primer3"""
	def __init__(self):
		self.name = ""
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
		self.difnum = 0
		self.direction = ""
	def formatprimer(self):
		formatout = "\t".join(str(x) for x in [self.start, self.end, self.difnum, self.difthreeall, self.length, self.tm, self.gc, self.anys, self.three, self.end_stability, self.hairpin, self.seq, ReverseComplement(self.seq)])
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

#from sys import platform
def get_software_path(base_path):
	if sys.platform.startswith('linux'): # linux
		primer3_path = base_path + "/primer3_core"
		muscle_path = base_path + "/muscle"
	elif sys.platform == "win32" or sys.platform == "cygwin": # Windows...
		primer3_path = base_path + "/primer3_core.exe"
		muscle_path = base_path + "/muscle.exe"
	elif sys.platform == "darwin": # MacOSX
		primer3_path = base_path + "/primer3_core_darwin64"
		muscle_path = base_path + "/muscle3.8.31_i86darwin64"
	return primer3_path, muscle_path

# simple Tm calculator
def Tm(seq):
	t=0
	for a in seq:
		if a=='A' or a=='T':
			t=t+2
		if a=='C' or a=='G':
			t=t+4
	return t

# calculate GC content of a sequence
def Calc_GC(seq):
	t = 0.0 # float
	for a in seq:
		if a=='C' or a=='G':
			t += 1
	return t / len(seq) * 100

# function to find the segments with largest Tm
def FindLongestSubstring(s1, s2):
	longestStart = 0
	longestEnd = 0
	largestTm = 0
	start = 0
	gaps = [i for i, c in enumerate(s1) if c=='-' or s2[i]=='-']
	gaps.append(len(s1))
	for gap in gaps:
		end = gap
		tm = Tm(s1[start:end])
		if  tm > largestTm:
			longestStart = start
			longestEnd = end
			largestTm = tm
		start = gap + 1
	nL = len(s1[:longestStart].replace("-","")) # number of bases on the left of the longest segments
	nR = len(s1[longestEnd:].replace("-",""))  # number of bases on the right of the longest segments
	return [longestStart, longestEnd, nL, nR]

# get the list of sequences in the homeolog groups for comparison with current primer
def get_homeo_seq(fasta, target, ids, align_left, align_right):
	s1 = fasta[target] # primer sequence in the template with gaps
	seq2comp = [] # final sequence to compare for each other homeolog
	for k in ids:
		s2 = fasta[k]
		targetSeq = s1[align_left:(align_right + 1)]
		homeoSeq = s2[align_left:(align_right + 1)]
		#print "Targetseq ", targetSeq
		#print "homeoSeq  ", homeoSeq
		# Get the sequences for comparison
		indexL, indexR, nL, nR = FindLongestSubstring(targetSeq, homeoSeq)
		indexL += align_left
		indexR += align_left
		seqL = s2[:indexL].replace("-","")
		seqR = s2[indexR:].replace("-","")
		#print "indexL, indexR, nL, nR ", indexL, indexR, nL, nR
		#print "s2[indexL:indexR] ", s2[indexL:indexR]
		if len(seqL) < nL: # in case it does not have enough bases
			seqL = "-" * (nL - len(seqL)) + seqL
		if len(seqR) < nR:
			seqL = seqR + "-" * (nR - len(seqR))
		seqk = seqL[::-1][:nL][::-1] + s2[indexL:indexR] + seqR[:nR]
		seq2comp.append(seqk)
		#print k, "\t", seqk
	return seq2comp

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
		#print "find_substring(ss, wild_seq), ", i, find_substring(ss, wild_seq)
		#print "find_substring(ss, mutt_seq), ", i, find_substring(ss, mut_seq) 
		#if len(re.findall(ss, wild_seq)) != len(re.findall(ss, mut_seq)): # differnt length should be caused by the differnce in the SNP/indel
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

# function to blast and parse the output of each primer in the wheat genome
# depends on function: mismtachn
def primer_blast(primer_for_blast, outfile_blast):
	forblast = open("for_blast_primer.fa", 'w') # for blast against the gnome
	for k, v in primer_for_blast.items(): # k is the sequence and v is the number
		forblast.write(">" + v + "\n" + k + "\n")
	forblast.close()
	blast_hit = {} # matched chromosomes for primers: less than 2 mismatches in the first 4 bps from 3'
	### for blast
	reference = "/Library/WebServer/Documents/blast/db/nucleotide/161010_Chinese_Spring_v1.0_pseudomolecules.fasta"
	cmd2 = 'blastn -task blastn -db ' + reference + ' -query for_blast_primer.fa -outfmt "6 std qseq sseq qlen slen" -num_threads 3 -word_size 7 -out ' + outfile_blast
	print "Step 2: Blast command:\n", cmd2
	call(cmd2, shell=True)
	# process blast file
	# blast fields
	# IWB50236_7A_R	IWGSC_CSS_7DS_scaff_3919748	98.718	78	1	0	24	101	4891	4968	1.55e-30	138	CTCATCAAATGATTCAAAAATATCGATRCTTGGCTGGTGTATCGTGCAGACGACAGTTCGTCCGGTATCAACAGCATT	CTCATCAAATGATTCAAAAATATCGATGCTTGGCTGGTGTATCGTGCAGACGACAGTTCGTCCGGTATCAACAGCATT	101 5924
	# Fields: 
	# 1: query id, subject id, % identity, alignment length, mismatches, gap opens, 
	# 7: q. start, q. end, s. start, s. end, evalue, bit score
	# 13: q. sequence, s. sequence, q. length s. length
	for line in open(outfile_blast):
		if line.startswith('#'):
			continue
		fields = line.split("\t")
		query, subject, pct_identity, align_length= fields[:4]
		qstart, qstop, sstart, sstop = [int(x) for x in fields[6:10]]
		qseq, sseq = fields[12:14]
		qlen = int(fields[14])
		n1 = qlen - qstop
		if n1 < 2 and mismatchn(qseq[(n1 - 4):], sseq[(n1 - 4):]) + n1 < 2: # if less than 2 mismtaches in the first 4 bases from the 3' end of the primer
			blast_hit[query] = blast_hit.setdefault(query, "") + ";" + subject + ":" + str(sstart)
	return blast_hit

# function to extract sequences from a fasta file 
def get_fasta(infile):
	fasta = {} # dictionary for alignment
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line.startswith(">"):
				sequence_name = line.split()[0].lstrip(">")
			else:
				fasta.setdefault(sequence_name, "")
				fasta[sequence_name] += line.rstrip()
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

# parse primer3 output
def parse_primer3output(primer3output, primerpair_to_return):
	regex = "012345"[:primerpair_to_return]
	regex = "([" + regex + "])"
	primerpairs= {}
	with open(primer3output) as infile:
		for line in infile:
			line = line.strip()
			if "SEQUENCE_ID" in line:
				seqid = line.split("=")[1]
				for i in range(0, primerpair_to_return):
					primerpairs[seqid + "-" + str(i)] = PrimerPair()
				continue
			elif re.search("PRIMER_PAIR_" + regex + "_PENALTY", line):
				mm = re.search("PRIMER_PAIR_" + regex + "_PENALTY", line)
				primerpairs[seqid + "-" + mm.group(1)].penalty = line.split("=")[1]
				continue
			elif re.search( "PRIMER_PAIR_" + regex + "_COMPL_ANY", line):
				mm = re.search( "PRIMER_PAIR_" + regex + "_COMPL_ANY", line)
				primerpairs[seqid + "-" + mm.group(1)].compl_any = line.split("=")[1]
				continue
			elif re.search( "PRIMER_PAIR_" + regex + "_COMPL_END", line):
				mm = re.search( "PRIMER_PAIR_" + regex + "_COMPL_END", line)
				primerpairs[seqid + "-" + mm.group(1)].compl_end = line.split("=")[1]
			elif re.search( "PRIMER_PAIR_" + regex + "_PRODUCT_SIZE", line):
				mm = re.search( "PRIMER_PAIR_" + regex + "_PRODUCT_SIZE", line)
				primerpairs[seqid + "-" + mm.group(1)].product_size = int(line.split("=")[1])
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "_SEQUENCE", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_SEQUENCE", line)
				primerpairs[seqid + "-" + mm.group(1)].left.seq = line.split("=")[1]
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "=", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "=", line)
				primerpairs[seqid + "-" + mm.group(1)].left.start = int(line.split("=")[1].split(",")[0])
				primerpairs[seqid + "-" + mm.group(1)].left.length = int(line.split("=")[1].split(",")[1])
				primerpairs[seqid + "-" + mm.group(1)].left.end = primerpairs[seqid + "-" + mm.group(1)].left.start + primerpairs[seqid + "-" + mm.group(1)].left.length - 1
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "_TM", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_TM", line)
				primerpairs[seqid + "-" + mm.group(1)].left.tm = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "_GC_PERCENT", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_GC_PERCENT", line)
				primerpairs[seqid + "-" + mm.group(1)].left.gc = float(line.split("=")[1])
			elif re.search( "PRIMER_LEFT_" + regex + "_SELF_ANY_TH", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_SELF_ANY_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].left.anys = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "_SELF_END_TH", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_SELF_END_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].left.three = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "_HAIRPIN_TH", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_HAIRPIN_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].left.hairpin = float(line.split("=")[1])
			elif re.search( "PRIMER_LEFT_" + regex + "_END_STABILITY", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_END_STABILITY", line)
				primerpairs[seqid + "-" + mm.group(1)].left.end_stability = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_SEQUENCE", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_SEQUENCE", line)
				primerpairs[seqid + "-" + mm.group(1)].right.seq = line.split("=")[1]
			elif re.search( "PRIMER_RIGHT_" + regex + "=", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "=", line)
				primerpairs[seqid + "-" + mm.group(1)].right.start = int(line.split("=")[1].split(",")[0])
				primerpairs[seqid + "-" + mm.group(1)].right.length = int(line.split("=")[1].split(",")[1])
				primerpairs[seqid + "-" + mm.group(1)].right.end = primerpairs[seqid + "-" + mm.group(1)].right.start - primerpairs[seqid + "-" + mm.group(1)].right.length + 1
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_TM", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_TM", line)
				primerpairs[seqid + "-" + mm.group(1)].right.tm = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_GC_PERCENT", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_GC_PERCENT", line)
				primerpairs[seqid + "-" + mm.group(1)].right.gc = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_SELF_ANY_TH", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_SELF_ANY_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].right.anys = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_SELF_END_TH", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_SELF_END_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].right.three = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_HAIRPIN_TH", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_HAIRPIN_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].right.hairpin = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_END_STABILITY", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_END_STABILITY", line)
				primerpairs[seqid + "-" + mm.group(1)].right.end_stability = float(line.split("=")[1])
				continue
	return primerpairs

# function to find primer sequence variation site and highligh them in primer sequences
def format_primer_seq(primer, variation): # input is a primer object and variation list
	if primer.start < primer.end:
		start = primer.start
		end = primer.end
		seq = primer.seq
		#primer_range = range(primer.start - 1, primer.end)
	else:
		start = primer.end
		end = primer.start
		seq = ReverseComplement(primer.seq)
		#primer_range = range(primer.end - 1, primer.start)
	
	primer_range = range(start - 1, end)
	var_sites = set(variation).intersection(primer_range)
	var_sites_relative = [i - start + 1 for i in var_sites]
	#seq = primer.seq.lower()
	#seq = primer.seq
	for i in var_sites_relative:
		seq = seq[:i] + seq[i].upper() + seq[i+1:]
	if primer.start < primer.end:
		primer.seq = seq
	else:
		primer.seq = ReverseComplement(seq)
	primer.difnum = len(var_sites)
	return primer

# parse the flanking sequence input

def caps(seqfile, target, SNP_A, SNP_B, snp_pos, max_price): # two alleles now support indels or long haplotypes, SNP_A is the template allele (or reference allel)
	#flanking_temp_marker_IWB1855_7A_R_251.fa
	#snpname, chrom, allele, pos =re.split("_|\.", seqfile)[3:7]
	#chrom = chrom[0:2] # no arm
	#snp_pos = int(pos) - 1 # 0-based
	#print "snpname, chrom, allele, pos ", snpname, chrom, allele, pos
	snp_pos = snp_pos - 1 # 0-based
	snpname = re.sub(".fa|.fasta","", seqfile)
	getcaps_path = os.path.dirname(os.path.realpath(__file__))
	directory = "CAPS_output"
	if not os.path.exists(directory):
		os.makedirs(directory)
	out = directory + "/selected_CAPS_primers_" + snpname + ".txt"
	print "Output selected CAPS file name is: ", out
	
	# software path
	primer3_path, muscle_path = get_software_path(getcaps_path)
	
	# get target and ids and rename fasta seq names
	#fasta_raw, target, ids = get_fasta2(seqfile, chrom)
	fasta_raw = get_fasta(seqfile)
	seq_names = fasta_raw.keys()
	ids = [i for i in seq_names if i != target] # other non-target sequence names
	print "ids is ", ids
	
	seq_template = fasta_raw[target]
	################## Get CAPS information
	# step 1: read the enzyme file
	RE_file = getcaps_path + "/NEB_parsed_REs.txt" # this one removed some duplicated cuttings
	REs = parse_RE_file(RE_file) # get the list of restriction enzymes
	# step 2: get the list of enzymes that can be used for caps or dcaps
	#SNP_A, SNP_B = iupac[allele] # SNP 2 alleles
	print "SNP_A, SNP_B ", SNP_A, SNP_B # two alleles now support indels or long haplotypes
	wild_seq = seq_template[:snp_pos] + SNP_A + seq_template[snp_pos+len(SNP_A):]
	mut_seq = seq_template[:snp_pos] + SNP_B + seq_template[snp_pos+len(SNP_A):]
	pos_L, pos_R = dif_region(wild_seq, mut_seq) # differnce region borders in the wild_seq, use this instead of snp_pos to decide use as left end or right end 
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
	variation = [] # variation sites that can differ ALL
	variation2 = [] # variation sites that can differ at least 2 homeologs
	
	# STEP 0: create alignment file and primer3output file
	RawAlignFile = "alignment_raw_" + seqfile
	alignmentcmd = muscle_path + " -in " + seqfile + " -out " + RawAlignFile + " -quiet"
	print "Alignment command: ", alignmentcmd
	call(alignmentcmd, shell=True)
	
	#########################
	if not len(ids): # if there is no homeologs found, such as when ploidy is 1
		# loop to write primer3 input for each variation site
		# primer3 inputfile
		primer3input = directory + "/primer3.input." + snpname
		p3input = open(primer3input, 'w')
		# for dcaps
		product_min = 70
		product_max = 250
		n = 0 # number to see how many records were written to primer3 input file
		for enzyme in dcaps_list:
			for primer_end_pos in enzyme.primer_end_pos: # a list of potential end positions
				if enzyme.primer_direction == "right":
					left_end = -1000000
					right_end = primer_end_pos + 1
				else:
					left_end = primer_end_pos + 1
					right_end = -1000000
				#if right_end - left_end < product_min - 35 or right_end - left_end > product_max - 35: # suppose both primers are 18 bp
				#	continue
				settings = "PRIMER_TASK=generic" + "\n" + \
				"SEQUENCE_ID=" + snpname + "-dCAPS-" + enzyme.name + "-" + enzyme.seq + "-" + str(primer_end_pos+1) + "\n" + \
				"SEQUENCE_TEMPLATE=" + enzyme.template_seq + "\n" + \
				"PRIMER_PRODUCT_SIZE_RANGE=150-200 200-250 70-150" + "\n" + \
				"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + getcaps_path + "/primer3_config/"  + "\n" + \
				"PRIMER_MAX_SIZE=" + maxSize + "\n" + \
				"PRIMER_MIN_SIZE=" + minSize  + "\n" + \
				"PRIMER_MAX_TM=" + maxTm + "\n" + \
				"PRIMER_MIN_TM=" + minTm + "\n" + \
				"PRIMER_PAIR_MAX_DIFF_TM=6.0" + "\n" + \
				"PRIMER_FIRST_BASE_INDEX=1" + "\n" + \
				"PRIMER_LIBERAL_BASE=1" + "\n" + \
				"PRIMER_NUM_RETURN=5"  + "\n" + \
				"PRIMER_EXPLAIN_FLAG=1"  + "\n" + \
				"SEQUENCE_FORCE_LEFT_END=" + str(left_end) + "\n" + \
				"SEQUENCE_FORCE_RIGHT_END=" + str(right_end) + "\n" + \
				"="
				n += 1
				p3input.write(settings + "\n")
		# for caps
		product_min = 300
		product_max = 900
		left_end = -1000000
		right_end = -1000000
		for enzyme in caps_list:
			settings = "PRIMER_TASK=generic" + "\n" + \
			"SEQUENCE_ID=" + snpname + "-CAPS-" + enzyme.name + "-" + enzyme.seq  + "\n" + \
			"SEQUENCE_TEMPLATE=" + enzyme.template_seq + "\n" + \
			"PRIMER_PRODUCT_SIZE_RANGE=" + str(product_min) + "-" + str(product_max) + "\n" + \
			"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + getcaps_path + "/primer3_config/"  + "\n" + \
			"PRIMER_MAX_SIZE=" + maxSize + "\n" + \
			"PRIMER_MIN_SIZE=" + minSize  + "\n" + \
			"PRIMER_MAX_TM=" + maxTm + "\n" + \
			"PRIMER_MIN_TM=" + minTm + "\n" + \
			"PRIMER_PAIR_MAX_DIFF_TM=6.0" + "\n" + \
			"PRIMER_FIRST_BASE_INDEX=1" + "\n" + \
			"PRIMER_LIBERAL_BASE=1" + "\n" + \
			"PRIMER_NUM_RETURN=5"  + "\n" + \
			"PRIMER_EXPLAIN_FLAG=1"  + "\n" + \
			"SEQUENCE_FORCE_LEFT_END=" + str(left_end) + "\n" + \
			"SEQUENCE_FORCE_RIGHT_END=" + str(right_end) + "\n" + \
			"SEQUENCE_TARGET=" + str(snp_pos - 20) + ",40" + "\n" + \
			"="
			n += 1
			p3input.write(settings + "\n")

		p3input.close()
		
		if n == 0:
			print "No primer3 input were found"
			outfile = open(out, 'w')
			outfile.write("\nCAPS cut information for SNP " + snpname + "\n") # change to 1 based
			outfile.write("Enzyme\tEnzyme_seq\tChange_pos\tOther_cut_pos\n")
			for enzyme in caps_list:
				outfile.write(enzyme.name + "\t" + enzyme.seq + "\t" + ", ".join([str(x + 1) for x in enzyme.allpos]) + "\n")
			outfile.write("\ndCAPS cut information for SNP " + snpname + "\n") # change to 1 based
			for enzyme in dcaps_list:
				outfile.write(enzyme.name + "\t" + enzyme.seq + "\t" + str(enzyme.change_pos) + "\t" + ", ".join([str(x + 1) for x in enzyme.allpos]) + "\n")
			outfile.close()
			return 0
		# primer3 output file
		primer3output = directory + "/primer3.output." + snpname
		p3cmd = primer3_path + " -default_version=2 -output=" + primer3output + " " + primer3input
		print "Primer3 command 1st time: ", p3cmd
		call(p3cmd, shell=True)
		primerpairs = parse_primer3output(primer3output, 5)
	else: # if there are homeolog sequences	
	########################
		
		########################
		# read alignment file
		fasta = get_fasta(RawAlignFile)

		## get the variaiton site among sequences
		#ids = [] # all other sequence names
		#for kk in fasta.keys():
			#key_chr = kk.split("_")[2] # sequence chromosome name
			#if chrom in key_chr or key_chr in chrom: # 3B contig names do not have chromosome arm
				#target = kk
			#else:
				#ids.append(kk)
				
		print "The target: ", target
		print "The other groups: ", ids

		alignlen = len(fasta[target])
		print "Alignment length: ", alignlen
		
		# get the target ID template base coordinate in the alignment
		t2a = {} # template to alignment
		a2t = {}
		ngap = 0 # gaps
		for i in range(alignlen):
			if fasta[target][i] == "-":
				ngap += 1
				continue
			t2a[i - ngap] = i
			a2t[i] = i - ngap

		print "last key of t2a", i - ngap
		print "last key of a2t", i
		
		seq_template = fasta[target].replace("-","") # remove all gaps
		variation = [] # variation sites that can differ ALL
		variation2 = [] # variation sites that can differ at least 2 homeologs
		
		# calculate gap number on the left and right
		gap_left_target = len(fasta[target]) - len(fasta[target].lstrip('-'))
		gap_left = max([len(v) - len(v.lstrip('-')) for k, v in fasta.items()])
		gap_right = min([len(v.rstrip('-')) for k, v in fasta.items()])
		print "gap_left_target, gap_left and gap_right: ", gap_left_target, gap_left, gap_right
		
		diffarray = {} # a list of 0 or 1: the same as or different from the site in each sequences of ids
		for i in range(gap_left, gap_right): # exclude 20 bases on each side
			b1 = fasta[target][i]
			if b1 == "-":  # target non-gap base
				ngap += 1
				continue
			pos_template = a2t[i] # position in the target template (no gaps)
			if pos_template < 20 or pos_template > len(seq_template) - 20:
				continue
			nd = 0 # number of difference
			da = [0] * len(ids) # differ array
			m = 0 # counter of next loop
			if pos_template < snp_pos:
				align_left = t2a[pos_template - 19] # 20 bp left of current position
				align_right = i
			else:
				align_left = i # 20 bp left of current position
				align_right = t2a[pos_template + 19]
			seq2comp = get_homeo_seq(fasta, target, ids, align_left, align_right) # list of sequences for comparison
			for k in seq2comp:
				if pos_template < snp_pos:
					b2 = k[-1] # homeolog non-gap bas
				else:
					b2 = k[0]
				if b1 != b2:
					nd += 1
					da[m] = 1 # m sequence has variation from target
				m += 1

			# for each site pos_template
			diffarray[pos_template] = da
			if nd == len(ids): # different from all other sequences
				if pos_template not in variation: # in case multiple gaps
					variation.append(pos_template)
			if nd > 0: # different from at least 1 other sequence
				if pos_template not in variation2: # in case multiple gaps
					variation2.append(pos_template)
		
		print "Sites that can differ all\n", variation
		print "\nSites that can differ at least 1\n", variation2
		#print "\nKeys of diffarray: ", diffarray.keys()
		
		########################################################
		# loop to write primer3 input for each variation site
		# primer3 inputfile
		primer3input = directory + "/primer3.input." + snpname
		p3input = open(primer3input, 'w')
		# for dcaps
		product_min = 70
		product_max = 350
		n = 0 # number to see how many records were written to primer3 input file
		for enzyme in dcaps_list:
			for primer_end_pos in enzyme.primer_end_pos: # a list of potential end positions
				for i in variation:
					#if primer_end_pos > snp_pos:
					if enzyme.primer_direction == "right":
						left_end = i + 1
						right_end = primer_end_pos + 1
					else:
						left_end = primer_end_pos + 1
						right_end = i + 1
					if right_end - left_end < product_min - 35 or right_end - left_end > product_max - 35: # suppose both primers are 18 bp
						continue
					settings = "PRIMER_TASK=generic" + "\n" + \
					"SEQUENCE_ID=" + snpname + "-dCAPS-" + enzyme.name + "-" + enzyme.seq + "-" + str(i+1) + "-" + str(primer_end_pos+1) + "\n" + \
					"SEQUENCE_TEMPLATE=" + enzyme.template_seq + "\n" + \
					"PRIMER_PRODUCT_SIZE_RANGE=" + str(product_min) + "-" + str(product_max) + "\n" + \
					"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + getcaps_path + "/primer3_config/"  + "\n" + \
					"PRIMER_MAX_SIZE=" + maxSize + "\n" + \
					"PRIMER_MIN_SIZE=" + minSize  + "\n" + \
					"PRIMER_MAX_TM=" + maxTm + "\n" + \
					"PRIMER_MIN_TM=" + minTm + "\n" + \
					"PRIMER_PAIR_MAX_DIFF_TM=6.0" + "\n" + \
					"PRIMER_FIRST_BASE_INDEX=1" + "\n" + \
					"PRIMER_LIBERAL_BASE=1" + "\n" + \
					"PRIMER_NUM_RETURN=5"  + "\n" + \
					"PRIMER_EXPLAIN_FLAG=1"  + "\n" + \
					"SEQUENCE_FORCE_LEFT_END=" + str(left_end) + "\n" + \
					"SEQUENCE_FORCE_RIGHT_END=" + str(right_end) + "\n" + \
					"="
					n += 1
					p3input.write(settings + "\n")
		# for caps
		product_min = 300
		product_max = 900
		for enzyme in caps_list:
			for i in variation: # sites that can differ all
				if i < snp_pos:
					left_end = i + 1
					right_end = -1000000
				else:
					left_end = -1000000
					right_end = i + 1
				settings = "PRIMER_TASK=generic" + "\n" + \
				"SEQUENCE_ID=" + snpname + "-CAPS-" + enzyme.name + "-" + enzyme.seq + "-" + str(i+1) + "\n" + \
				"SEQUENCE_TEMPLATE=" + enzyme.template_seq + "\n" + \
				"PRIMER_PRODUCT_SIZE_RANGE=" + str(product_min) + "-" + str(product_max) + "\n" + \
				"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + getcaps_path + "/primer3_config/"  + "\n" + \
				"PRIMER_MAX_SIZE=" + maxSize + "\n" + \
				"PRIMER_MIN_SIZE=" + minSize  + "\n" + \
				"PRIMER_MAX_TM=" + maxTm + "\n" + \
				"PRIMER_MIN_TM=" + minTm + "\n" + \
				"PRIMER_PAIR_MAX_DIFF_TM=6.0" + "\n" + \
				"PRIMER_FIRST_BASE_INDEX=1" + "\n" + \
				"PRIMER_LIBERAL_BASE=1" + "\n" + \
				"PRIMER_NUM_RETURN=5"  + "\n" + \
				"PRIMER_EXPLAIN_FLAG=1"  + "\n" + \
				"SEQUENCE_FORCE_LEFT_END=" + str(left_end) + "\n" + \
				"SEQUENCE_FORCE_RIGHT_END=" + str(right_end) + "\n" + \
				"SEQUENCE_TARGET=" + str(snp_pos - 20) + ",40" + "\n" + \
				"="
				n += 1
				p3input.write(settings + "\n")	

		p3input.close()
		
		if n == 0:
			print "No primer3 input were found"
			outfile = open(out, 'w')
			outfile.write("Sites that can differ all for " + snpname + "\n")
			outfile.write(", ".join([str(x + 1) for x in variation])) # change to 1 based
			outfile.write("\nCAPS cut information for SNP " + snpname + "\n") # change to 1 based
			outfile.write("Enzyme\tEnzyme_seq\tChange_pos\tOther_cut_pos\n")
			for enzyme in caps_list:
				outfile.write(enzyme.name + "\t" + enzyme.seq + "\t" + ", ".join([str(x + 1) for x in enzyme.allpos]) + "\n")
			outfile.write("\ndCAPS cut information for SNP " + snpname + "\n") # change to 1 based
			for enzyme in dcaps_list:
				outfile.write(enzyme.name + "\t" + enzyme.seq + "\t" + str(enzyme.change_pos) + "\t" + ", ".join([str(x + 1) for x in enzyme.allpos]) + "\n")
			outfile.close()
			return 0
		# primer3 output file
		primer3output = directory + "/primer3.output." + snpname
		p3cmd = primer3_path + " -default_version=2 -output=" + primer3output + " " + primer3input
		print "Primer3 command 1st time: ", p3cmd
		call(p3cmd, shell=True)

		##########################
		### parse primer 3 primer check output
		primerpairs = parse_primer3output(primer3output, 1)
		#primerpairs = {} # sequence ID is the key

	#######################################
	# write to file
	outfile = open(out, 'w')
	outfile.write("index\tproduct_size\ttype\tstart\tend\tdiff_number\t3'differall\tlength\tTm\tGCcontent\tany\t3'\tend_stability\thairpin\tprimer_seq\tReverseComplement\tpenalty\tcompl_any\tcompl_end\tPrimerID\tmatched_chromosomes\n")
	# Get primer list for blast
	primer_for_blast = {}
	final_primers = {}
	nL = 0 # left primer count
	nR = 0 # right primer count
	for i, pp in primerpairs.items():
		if pp.product_size != 0:
			pl = pp.left
			pr = pp.right
			if pl.seq not in primer_for_blast:
				nL += 1
				pl.name = "L" + str(nL)
				primer_for_blast[pl.seq] = pl.name # use seq as keys
			else:
				pl.name = primer_for_blast[pl.seq]
			if pr.seq not in primer_for_blast:
				nR += 1
				pr.name = "R" + str(nR)
				primer_for_blast[pr.seq] = pr.name # because a lot of same sequences
			else:
				pr.name = primer_for_blast[pr.seq]
			pp.left = pl
			pp.right = pr
			final_primers[i] = pp
				
	# blast primers
	blast_hit = {}
	outfile_blast = directory + "/primer_blast_out_" + snpname + ".txt"
	if blast and len(primer_for_blast) > 0:
		blast_hit = primer_blast(primer_for_blast, outfile_blast) # chromosome hit for each primer
	# write output file
	for i, pp in final_primers.items():
		#varsite = int(i.split("-")[-1]) # variation site
		#pl = pp.left
		#pr = pp.right
		pl = format_primer_seq(pp.left, variation)
		pr = format_primer_seq(pp.right, variation)
		# check whether 3' can differ all: not necessary for here, because I only used sites that can differ all.
		#if varsite in variation:
		pl.difthreeall = "YES"
		pr.difthreeall = "YES"
		outfile.write("\t".join([i, str(pp.product_size), "LEFT", pl.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, pl.name, blast_hit.setdefault(pl.name, "")]) + "\n")
		outfile.write("\t".join([i, str(pp.product_size), "RIGHT", pr.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, pr.name, blast_hit.setdefault(pr.name, "")]) + "\n")
	
	outfile.write("\n\nSites that can differ all for " + snpname + "\n")
	outfile.write(", ".join([str(x + 1) for x in variation])) # change to 1 based
	# write CAPS cut information
	outfile.write("\n\nCAPS cut information for snp " + snpname + "\n") # change to 1 based
	outfile.write("Enzyme\tEnzyme_seq\tChange_pos\tOther_cut_pos\n")
	for enzyme in dcaps_list + caps_list:
		outfile.write(enzyme.name + "\t" + enzyme.seq + "\t" + str(enzyme.change_pos) + "\t" + ", ".join([str(x + 1) for x in enzyme.allpos]) + "\n")
	# close outfile
	outfile.write("\n\n\n")
	outfile.close()

# design primers

caps(seqfile, target, SNP_A, SNP_B, snp_pos, max_price)
	
