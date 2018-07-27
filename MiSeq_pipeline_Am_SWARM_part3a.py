#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "9, July 18th, 2018"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from Bio import SeqIO
seqlist = []
OTUlist = []
duplicatelist = []



def removeoutgroup(outputpath, outgrouptree, tree):
	print("Remove OTUs based on outgroup tree\n")
	if not os.path.exists(outputpath + 'outgroup_removal/'): 
		print('ERROR! No outgroup_removal folder. Run part 2!')	
		main()
	os.system('python3 Miseq_scripts/8_remove_outgroup_from_tree.py outputs/outgroup_removal/' +outgrouptree + ' outputs/outgroup_removal/' +tree + ' outputs/OTUs/SWARM_postout.txt outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta')
	
def main():
	c = sys.argv[1]

	pathA = os.getcwd()
	path = pathA + "/Rawdata/"


#	print("The trees needs to be in the outgroup_removal folder and in newick format\n")
#	c = input('what is the name of your tree (for example:  RAxML_labelledTree_masked_<project>.tre and the same name with _outgroup at the end for the outgroup tree :RAxML_labelledTree_masked_<project>_outgroup.tre ) \n Do Not add the full path or the script will break \n')
	try:
		tree = c
		outgrouptree = c.split('.tre')[0]+'_outgroup.tre'
	except ValueError:
		c = ""	
	if c == "":
		print ('Your input is empty.  Try again. ')

	outputpath = pathA + '/outputs/'
	if not os.path.exists(outputpath): 
		print('ERROR! No output folder. run part 2!')
		main()
	removeoutgroup(outputpath, outgrouptree, tree)
main()