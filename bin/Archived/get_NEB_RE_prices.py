#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  get_NEB_RE_prices.py
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

from bs4 import BeautifulSoup
import urllib2

# get links for all NEB restriction enzymes
def get_link():
	enzyme_links = []
	url = "https://www.neb.com/tools-and-resources/selection-charts/frequencies-of-restriction-sites"
	content = urllib2.urlopen(url).read()
	soup = BeautifulSoup(content, 'html.parser')
	table = soup.find('table')
	table_body = table.find('tbody')
	for link in table_body.find_all('a'):	
		if "-HF" in link.contents[0]: # enzyme name
			continue
		enzyme_links.append("https://www.neb.com" + link.get('href'))
	return enzyme_links

# get the price for a specific NEB enzyme
def get_price(url):
	#url = "https://www.neb.com/products/r0117-aatii"
	content = urllib2.urlopen(url).read()
	soup = BeautifulSoup(content, 'html.parser')
	enzyme = soup.title.string.strip().split(" ")[0]
	# price
	table = soup.find('table', attrs={'class':"items add-to-cart-list"})
	table_body = table.find('tbody')
	rows = table_body.find_all('tr')
	if not rows: # some enzymes are discontinued
		return enzyme + "\tNA\tNA"
	cols = rows[0].find_all('td')
	enzyme_size =  cols[1].get_text().split(" ")[0] # quantity
	enzyme_price =  cols[3].get_text().strip("$") # price
	return "\t".join([enzyme, enzyme_size, enzyme_price])


def main():
	enzyme_links = get_link()
	print "Got all enzyme links"
	out = open("NEB_RE_prices.txt", 'w')
	out.write("Enzyme\tSize\tPrice\n")
	for url in enzyme_links[82:]:
		price = get_price(url)
		out.write(price.encode("utf-8") + "\n")
	print "Got all prices!"
	out.close()
	return 0

if __name__ == '__main__':
	import sys
	sys.exit(main())
