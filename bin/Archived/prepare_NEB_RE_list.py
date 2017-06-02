#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  prepare_NEB_RE_list.py
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

# the file is sorted by price
neb_re = "/home/junli/Documents/Git/getKASP_pipeline/bin/NEB_enzymes.txt"

REs = {}
prices = {}
with open(neb_re) as file_one:
	next(file_one)
	for line in file_one:
		info = line.rstrip().split("\t")
		enzyme = info[1]
		price = info[4]
		seq = info[6]
		REs[seq] = REs.setdefault(seq, "") + "," + enzyme
		if seq not in prices:
			prices[seq] = price
out = open("NEB_parsed_REs.txt", "w")
for k, v in REs.items():
	v = v.strip(',')
	out.write( v + "," + prices[k] + "\t" + k + "\n")

out.close()



