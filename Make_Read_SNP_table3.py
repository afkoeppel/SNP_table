#!/usr/bin/python

#Imports the os, sys, ,re, subprocess and string modules.
import sys
import re
from subprocess import *
import os
import string
import pybedtools
import pysam

#Arguments:  
# 1) working directory (absolute path)
path = (sys.argv[1])

#Sets up fixed variables (change as needed)
mastersamfilename = "SNP_readlist2_sorted.sam"
headerfilename = path + "/header.txt"
snpfile = "/m/cphg-nextgen1/cphg-nextgen1/projects/afk2s/Farber_rnaSeq/C3H_snps1.vcf"
genomefile = "/m/cphg-nextgen1/cphg-nextgen1/projects/afk2s/Farber_rnaSeq/mm9.fa"
counter = 0
totalcount = int(sys.argv[2])
writefilename = "Read_SNP_table.txt"
writefile = open(writefilename,'w')

#samfile = open(mastersamfilename, 'r')
#samBigList = samfile.readlines()
samfile = pysam.Samfile(mastersamfilename, "rs" )

#Read in SNP file with pybedtools
snpbed = pybedtools.BedTool(snpfile)

while counter <= totalcount:
	asam = samBigList[counter]
	print("cat " + headerfilename + " > read.sam")
	sts = call("cat " + headerfilename + " > read.sam", shell=True)
	littlesam=open("read.sam",'a')
	littlesam.write(asam)
	littlesam.close()
	print("samtools view -bS read.sam -o read.bam")
	sts = call("samtools view -bS read.sam -o read.bam", shell=True)
	readbed = pybedtools.BedTool('read.bam')
	readbed.intersect(snpbed, bed=True).saveas('readhits.bed')
#	print("bedtools intersect -abam -wo -a read.bam -b " + snpfile + " > readhits.bed")
#	sts = call("bedtools intersect -abam -wo -a read.bam -b " + snpfile + " > readhits.bed", shell=True)
	hitbedfile=open("readhits.bed",'r')
	bedList = hitbedfile.readlines()
	numsnps = len(bedList)
	print(numsnps)
	if numsnps > 1:
		snpList=[]
		genCallList=[]
		for abed in bedList:
			abed = abed.rstrip()
			lineSplit = string.split(abed,'\t')
			readid = lineSplit[3]
			snpList.append(lineSplit[2])	
		#samtools mpileup -d ../mm9.fa -l readhits.bed SNP_readlist.bam > SNP_pile.txt 
		print("samtools mpileup -f " + genomefile + " read.bam > SNP_pile.txt")
		sts = call("samtools mpileup -f " + genomefile + " read.bam > SNP_pile.txt", shell=True)
		pilefile = open("SNP_pile.txt",'r')
		pileList = pilefile.readlines()
		for apile in pileList:
			apile = apile.rstrip()
			splitPile = string.split(apile,'\t')
			position = splitPile[1]
			baselist = splitPile[4]
			if position in snpList:
				nummatch = len(re.findall('[,.]',baselist))
				mismatch = len(re.findall('[AaCcTtGg]',baselist))
				numcarrot = len(re.findall('[><]',baselist))
				numenns = len(re.findall('[Nn]',baselist))
				if nummatch == 1 and mismatch ==0 and numcarrot==0 and numenns ==0:
					genCallList.append('ref')
				elif nummatch == 0 and mismatch ==1 and numcarrot==0 and numenns ==0:
					genCallList.append('alt')
				elif nummatch == 0 and mismatch ==0 and numcarrot==1 and numenns ==0:
					genCallList.append('carrot')
				elif nummatch == 0 and mismatch ==0 and numcarrot==0 and numenns ==1:
					genCallList.append('N')
				else:
					genCallList.append('WTF?')
		x=0
		numcall = len(genCallList)
		while x < numcall:
			gencall1 = genCallList[x]
			snp1 = snpList[x]
			y = x + 1
			while y < numcall:
				snp2 = snpList[y]
				gencall2 = genCallList[y]
				if gencall1 =='alt' or gencall1=='ref':
					if gencall2=='alt' or gencall2=='ref':
						writefile.write(readid + '\t' + snp1 + '\t' + snp2 + '\t' + gencall1 + '\t' + gencall2 + '\n')
				y+=1
			x+=1 
		pilefile.close()	
	hitbedfile.close()
	counter += 1
#Closes write file:
writefile.close()
