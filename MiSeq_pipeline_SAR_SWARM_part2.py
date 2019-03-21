# python3 MiSeq_pipeline_SAR_SWARM_part2.py list_sample.txt (data_name)
# then reply to prompts (read guide before)
#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "9, June 1st, 2018"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from Bio import SeqIO
from sys import argv
seqlist = []
OTUlist = []
duplicatelist = []



def makesinglefastafile(file, path, outputpath, listsample):
	if not os.path.exists(outputpath + 'dereplicated/'): 
		os.makedirs(outputpath + 'dereplicated/') 	

	os.system('vsearch --derep_fulllength '+path+'/'+ file + ' --sizeout --fasta_width 0 --output outputs/dereplicated/'+file.split('.assembled.fas')[0]+'_dereplicated.fas')
	os.system('python3 Miseq_scripts/1_pool_rename.py ' + outputpath + 'dereplicated/'+file.split('.assembled.fas')[0]+'_dereplicated.fas ' + listsample)			# pool all the reads together in a huge file

def PickOTUSwarm(dSWARM , path,outputpath, listsample, readcutoff): #, dataname):
	if not os.path.exists(outputpath + 'OTUs/'): 
		os.makedirs(outputpath + 'OTUs/') 	
	#pick OTUs using SWARM
	print ("Pick OTUs")
	os.system('vsearch --derep_fulllength outputs/dereplicated/readpooled.fas --sizein --sizeout --fasta_width 0 --output outputs/OTUs/dereplicated_seqfile.fas --uc outputs/OTUs/dereplicated_seqfile.map.txt')
	os.system('python3 Miseq_scripts/2b_check_primer.py outputs/OTUs/dereplicated_seqfile.fas')
	os.system('vsearch --derep_fulllength outputs/OTUs/dereplicated_seqfile_primer.fas --sizein --sizeout --fasta_width 0 --output outputs/OTUs/dereplicated_seqprimer.fas --uc outputs/OTUs/dereplicated_seqprimer.map.txt')
	os.system('swarm -t 2 -s outputs/OTUs/statSWARM -d '+  str(dSWARM) +' -z outputs/OTUs/dereplicated_seqprimer.fas > outputs/OTUs/derepseqfile_output.swarm')
	print("Merge SWARM and dereplicate list")
	os.system('python3 Miseq_scripts/3_postSwarm_v2.py outputs/OTUs/derepseqfile_output.swarm outputs/OTUs/dereplicated_seqfile.map.txt outputs/OTUs/dereplicated_seqprimer.map.txt outputs/OTUs/dereplicated_seqprimer.fas')
	print ("Add read numbers")
	os.system('python3 Miseq_scripts/4_Add_numbers_v2.py outputs/OTUs/SWARM_postout.fas outputs/OTUs/SWARM_postout.txt '+listsample)# +' '+ dataname) #  '+ str(runref))
	print("Prepare files for Chimeras check")
	if not os.path.exists(outputpath + 'chimeras/'): 
		os.makedirs(outputpath + 'chimeras/') 	
	os.system('python3 Miseq_scripts/5a_Pre_Uchime_v.py outputs/OTUs/SWARM_postout_nosingleton.fas outputs/OTUs/OTUtable_temp.txt '+str(readcutoff))
	print("Chimera check using uchime_denovo implemented in vsearch")
	os.system('vsearch --uchime3_denovo outputs/chimeras/Seq_reads.fas --nonchimera outputs/chimeras/Seq_reads_nochimera_nosingleton.fas --uchimeout outputs/chimeras/chimeratable.txt')
	os.system('python3 Miseq_scripts/5b_Post_Uchime_v.py outputs/chimeras/Seq_reads.fas')
	os.system('python3 Miseq_scripts/5c_Water_remove_contaminant.py outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed.fas')
	
def RunBlast(AssTaxo, outputpath, idmin, qcov, readcutoff):
	if AssTaxo == 0:
		print("No taxonomic assignemt")
		print("Pipeline over")
		
	if AssTaxo == 1:
		if not os.path.exists(outputpath + 'taxonomic_assignment/'): 
			os.makedirs(outputpath + 'taxonomic_assignment/') 
		print("Run BLAST")
		os.system('python3 Miseq_scripts/6_BLASTn_Vsearch.py outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta ' +  str(idmin) + " "+ str(qcov) + ' SAR '+str(readcutoff))
#		os.system('python3 Miseq_scripts/6_BLASTn_V3_differential.py outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta ' + str(idmin) + " "+ str(qcov) + ' '+str(readcutoff) + ' ' + str(diffcutoff))

def makealignment(AssTaxo, outputpath):
	print("make alignment for outgroup removal. Take a while\n")
	if not os.path.exists(outputpath + 'outgroup_removal/'): 
		os.makedirs(outputpath + 'outgroup_removal/') 
	if AssTaxo == 1:
		os.system('mafft --addfragments outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_vsearch.fasta --thread 10 --reorder --mapout SAR_db/SSU_SAR_EUK_v14.3_mafft.fasta > outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_BLASTed_TA.fasta')
		os.system('python3 Miseq_scripts/7_remove_column.py outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_BLASTed_TA.fasta fasta SAR_db/SSU_SAR_EUK_v14.3_mask_75.txt')
		print("Now, you have to make a tree using this alignment (OTUseq_TA.fasta in outgroup_removal folder) and a constraint tree (in SAR_db folder)")
		print("commandline for tree builidng")
		print("raxmlHPC-PTHREADS-AVX2 -f v -s <alignment> -m GTRGAMMAI -t <constraint tree> -n <name of the output> -T 10")
		os.system('raxmlHPC-PTHREADS-AVX2 -f v -s outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_BLASTed_TA_masked.fas -m GTRGAMMAI -t SAR_db/SSU_SAR_EUK_v14.3_RAxML_constraint_rooted.tre -n test_MiSeq2018.tre -T 10')
	if AssTaxo == 0:
		os.system('mafft --addfragments outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta --thread 10 --reorder --mapout SAR_db/SSU_SAR_EUK_v14.3_mafft.fasta > outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_TA.fasta')
		os.system('python3 Miseq_scripts/7_remove_column.py outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_TA.fasta fasta SAR_db/SSU_SAR_EUK_v14.3_mask_75.txt')
		print("Now, you have to make a tree using this alignment (OTUseq_TA.fasta in outgroup_removal folder) and a constraint tree (in SAR_db folder)")
		print("commandline for tree builidng")
		print("raxmlHPC-PTHREADS-AVX2 -f v -s <alignment> -m GTRGAMMAI -t <constraint tree> -n <name of the output> -T 10")
		os.system('raxmlHPC-PTHREADS-AVX2 -f v -s outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_TA_masked.fas -m GTRGAMMAI -t SAR_db/SSU_SAR_EUK_v14.3_RAxML_constraint_rooted.tre -n test_MiSeq2018.tre -T 10')
	
def main():
	samplefile = sys.argv[1]
#	dname = sys.argv[2]
	listsamp = []
	pathA = os.getcwd()
	path = pathA + "/outputs/convertPEARfiles/"
	try:
		listsample = samplefile
	except ValueError:
		samplefile = ""	
	if samplefile == "":
		print ('Your input samplefile is empty.  Try again. ')
	for samp in open(listsample,'r'):
		if samp.split('\t')[0] not in listsamp:
			listsamp.append(samp.split('\t')[0])
# 	try:
# 		dataname = dname
# 	except ValueError:
# 		dname = ""	
# 	if dname == "":
# 		print ('Your input dataname is empty.  Try again. ')
	i = input('What percentage would you like to cluster your OTUs with SWARM (hit return for default of 1) ')
	try:
		num = int(i) + 1
	except TypeError:
		print ('Your input must be a number.  Try again. ')
		main()
	except ValueError:
		i = ""	
	if i == "":
		dSWARM = 1
	else:
		dSWARM = int(i)
	print ("you want to use SWARM at ", str(dSWARM))
	r = input('what is the minimum number of read for each OTU? (hit return for default of 100) :')
	try:
		num = float(r) + 1
	except TypeError:
		print('Your input must be a number. Try again.')
		main()
	except ValueError:
		r = ""
	if r == "":
		readcutoff = 100
	else:
		readcutoff = int(r)
	print(readcutoff)

	BLAST = input('Do you need to assign taxonomy using the BLAST tool? (yes or no)')
	if BLAST[0] == 'y':
		AssTaxo = 1
		x = input('what is the minimum identity cut off? (hit return for default of 90%) :')
		try:
			num = float(x) + 1
		except TypeError:
			print('Your input must be a number. Try again.')
			main()
		except ValueError:
			x = ""
		if x == "":
			idmin = 90.00
		else:
			idmin = float(x)
		print(idmin)

		y = input('what is the minimum query coverage cut off? (hit return for default of 70%) :')
		try:
			num = float(y) + 1
		except TypeError:
			print('Your input must be a number. Try again.')
			main()
		except ValueError:
			y = ""
		if y == "":
			qcov = 70.00
		else:
			qcov = float(y)
		print(qcov)

		z = input('Hit return when you are ready to continue. ')


	elif BLAST[0] == 'n':
		AssTaxo = 0; idmin = 90.00; qcov = 70.00; taxa = 'na'; readcutoff = 0
		z = input('Hit return when you are ready to continue. ')

	else:
		print ('Please answer yes or no. ')
		main()
	
	outputpath = pathA + '/outputs/'
	if not os.path.exists(outputpath): 
		os.makedirs(outputpath) 	
	temppath = pathA + '/temp/'
	if not os.path.exists(temppath): 
		os.makedirs(temppath) 	
	filnum = 0
	for file in os.listdir(path):
		filnum += 1
		makesinglefastafile(file, path, outputpath,listsample)
	if len(listsamp) != filnum:
		print(int(len(listsamp)), "<>",int(filnum),  "ISSUE with sample list! PLEASE CHECK !")
	else:
		PickOTUSwarm(dSWARM , path, outputpath, listsample, readcutoff)#, dataname)
		RunBlast(AssTaxo, outputpath, idmin, qcov, readcutoff)
		makealignment(AssTaxo, outputpath)
main()
