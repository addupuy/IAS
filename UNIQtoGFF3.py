#!/usr/bin/python

import fileinput, glob, os, sys

twoSideCut = int(sys.argv[1])
oneSideCut = int(sys.argv[2])
readcut = int(sys.argv[3])
mypath = os.getcwd()
os.chdir(mypath)
os.system("mkdir GFF3")
filterSettings = open(mypath + "/GFF3/filter_settings.txt","wa")
filterString = "Sites with reads on both sides: IRL or IRR >= " + str(twoSideCut) + "\nSites with reads on only one side: IRL or IRR >= " + str(oneSideCut) + "\nMinimum read number = " +str(readcut) + "\n"
filterSettings.write(filterString)
filterSettings.close()
for file in glob.glob("*.uniq"):
	os.chdir(mypath + "/GFF3")
	output = open(file[:-4]+"gff3","wa")
	os.chdir('..')
	myinput = open(file, "rU").readlines()
	count = 0
	for line in myinput:
		mylist = line.split('\t')
		sampleID = mylist[0]
		if mylist[1][0:3] == "chr":
			chrom = mylist[1]
		else:
			chrom = "chr" + mylist[1]
		IRL = float(mylist[8])
		IRR = float(mylist[9].replace('\n',''))
		IRLread = int(mylist[3])
		IRRread = int(mylist[4])
		score = max(IRL,IRR)
		if mylist[5] == "+":
			IDstring = "ID=" + str(sampleID) +"*"+ str(count) +";color=#154360"
		else:
			IDstring = "ID=" + str(sampleID) +"*"+ str(count) +";color=#D35400"
		if mylist[5] == "+":
			start = mylist[2]
			stop = int(start) + 500
		elif mylist[5] == '-':
			if int(mylist[2]) < 500:
				start = 1
			else:
				start = (int(mylist[2]) - 500)
			stop = mylist[2]
		if IRL > 0 and IRR > 0:
			if IRL >= twoSideCut or IRR >= twoSideCut:
				mystring = chrom + '\tT2/Onc3\tinsertion site\t' + str(start) + '\t' + str(stop) + '\t' + str(score) + '\t' + str(mylist[5]) + '\t.\t' + IDstring + '\t' + str(mylist[3]) + '\t' + str(mylist[4]) + '\t' + str(mylist[8]) + '\t' + str(mylist[9])
				output.write(mystring)
		elif IRL >= oneSideCut or IRR >= oneSideCut:
			if IRLread >= readcut or IRRread >= readcut:
				mystring = chrom + '\tT2/Onc3\tinsertion site\t' + str(start) + '\t' + str(stop) + '\t' + str(score) + '\t' + str(mylist[5]) + '\t.\t' + IDstring + '\t' + str(mylist[3]) + '\t' + str(mylist[4]) + '\t' + str(mylist[8]) + '\t' + str(mylist[9])
				output.write(mystring)
		count += 1
output.close()