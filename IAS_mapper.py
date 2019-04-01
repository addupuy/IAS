#!/usr/bin/python

import sys,os,linecache,itertools
barcode = sys.argv[1]

def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

def cutadapt_IRL1(inputR1, inputR2):
	os.system("cutadapt --discard-untrimmed -g GTATGTAAACTTCCGACTTCAACTG -o "+trimirlf1+" -p "+trimirlr1+" "+inputR1+" "+inputR2)

def cutadapt_IRL2(inputR1, inputR2):
	os.system("cutadapt -G ^GTAATACGACTCACTATAGGGCTCCGCTTAAGGGAC -o "+trimirlf2+" -p "+trimirlr2+" "+inputR1+" "+inputR2)
	
def cutadapt_IRL3(inputR1, inputR2):
	os.system("cutadapt -a GTCCCTTAAGCGGAGCCCTATAGTGAGTCGTATTAC -A CAGTTGAAGTCGGAAGTTTACATAC -o "+trimirlf3+" -p "+trimirlr3+" "+inputR1+" "+inputR2)

def cutadapt_IRR1(inputR1, inputR2):
	os.system("cutadapt --discard-untrimmed -g GTATGTAAACTTCCGACTTCAACTG -o "+trimirrf1+" -p "+trimirrr1+" "+inputR1+" "+inputR2)

def cutadapt_IRR2(inputR1, inputR2):
	os.system("cutadapt -G ^GTAATACGACTCACTATAGGGCTCCGCTTAAGGGAC -o "+trimirrf2+" -p "+trimirrr2+" "+inputR1+" "+inputR2)

def cutadapt_IRR3(inputR1, inputR2):
	os.system("cutadapt -a GTCCCTTAAGCGGAGCCCTATAGTGAGTCGTATTAC -A CAGTTGAAGTCGGAAGTTTACATAC -o "+trimirrf3+" -p "+trimirrr3+" "+inputR1+" "+inputR2)
	
def mapreads(genome, inputR1, inputR2, mysammy):
	os.system("hisat2 --no-hd --no-sq --quiet -x "+genome+" -q -1 "+inputR1+" -2 "+inputR2+" -S "+mysammy)

def cigar_caller(mystring):
	mynumb = ""
	total = 0
	for char in mystring:
		if char.isdigit() == True:
			mynumb += char
		elif char == "M":
			total += int(mynumb)
			mynumb = ""
		elif char == "I":
			total -= int(mynumb)
			mynumb = ""
		elif char == "D":
			total += int(mynumb)*2
			mynumb = ""
	return total
		
def samparse(irlsam,irrsam):
	results = {}
	ligation = {}
	read1_flags = [69,73,89,99,81,83,89]
	read2_flags = [161,163,169,145,147]
	read1_flags_pos = [69,73,89,99]
	read1_flags_neg = [81,83,89]
	read2_flags_pos = [161,163,169]
	read2_flags_neg = [145,147]
	IRLmax,IRRmax = 0,0
	with open(irlsam, "rU") as SAM:
		count = 0
		linecount = file_len(irlsam)
		for line, line2 in itertools.izip_longest(SAM,SAM,fillvalue=''):
			if linecount - count > 2:
				count += 2
				mapcheck = False
				mylist = line.split("\t")
				mylist2 = line2.split("\t")
				flag1,flag2 = int(mylist[1]), int(mylist2[1])
				if flag1 in read1_flags:
					myflag = flag1
					chrom, address = mylist[2], int(mylist[3])
					length = cigar_caller(mylist[5])
					myseq = mylist[9]
				elif flag2 in read1_flags:
					myflag = flag2
					chrom, address = mylist2[2], int(mylist2[3])
					length = cigar_caller(mylist2[5])
					myseq = mylist2[9]
				else:
					length = 0
				#Runs if first line contains the mapping data for the transposon read
				#Executes if transposon IRL read mapped to + strand
				if length >= 12:
					if myflag in read1_flags_pos:
						if myseq[0:2] == "TA":
							mytag = chrom +":"+ str(address) +"-" #flips strand for IRL
							if mytag in results:
								results[mytag][0] += 1
								if results[mytag][0] > IRLmax:
									IRLmax = results[mytag][0]
								mapcheck = True
							else:
								results[mytag] = [1,0]
								mapcheck = True
							if mytag in ligation:
								if length not in ligation[mytag][0]:
									ligation[mytag][0].append(length)
							else:
								ligation[mytag] = [[],[]]
								ligation[mytag][0].append(length)
					elif myflag in read1_flags_neg:
						if myseq[-2:] == "TA":
							mytag = chrom +":"+ str(address+(length-2)) +"+" #flips strand for IRL
							if mytag in results:
								results[mytag][0] += 1
								if results[mytag][0] > IRLmax:
									IRLmax = results[mytag][0]
								mapcheck = True
							else:
								results[mytag] = [1,0]
								mapcheck = True
							if mytag in ligation:
								if length not in ligation[mytag][0]:
									ligation[mytag][0].append(length)
							else:
								ligation[mytag] = [[],[]]
								ligation[mytag][0].append(length)
					if mapcheck == False: #paired read only if the tnp read cannot be mapped precisely
						if flag1 in read2_flags:
							myflag = flag1
							chrom, address = mylist[2], int(mylist[3])
							length = cigar_caller(mylist[5])
							myseq = mylist[9]
						elif flag2 in read2_flags:
							myflag = flag2
							chrom, address = mylist2[2], int(mylist2[3])
							length = cigar_caller(mylist2[5])
							myseq = mylist2[9]
						if length >= 10:
							if myflag in read2_flags_pos:
								if myseq[0:2] == "TA":
									mytag = chrom +":"+ str(address) +"-"
									if mytag in results:
										results[mytag][0] += 1
										if results[mytag][0] > IRLmax:
											IRLmax = results[mytag][0]
										mapcheck = True
									else:
										results[mytag] = [1,0]
										mapcheck = True
									if mytag in ligation:
										if length not in ligation[mytag][0]:
											ligation[mytag][0].append(length)
									else:
										ligation[mytag] = [[],[]]
										ligation[mytag][0].append(length)
							elif myflag in read2_flags_neg:
								if myseq[-2:] == "TA":
									mytag = chrom +":"+ str(address+(length-2)) +"+" #flips strand since it's the IRL read
									if mytag in results:
										results[mytag][0] += 1
										if results[mytag][0] > IRLmax:
											IRLmax = results[mytag][0]
										mapcheck = True
									else:
										results[mytag] = [1,0]
										mapcheck = True
									if mytag in ligation:
										if length not in ligation[mytag][0]:
											ligation[mytag][0].append(length)
									else:
										ligation[mytag] = [[],[]]
										ligation[mytag][0].append(length)
	with open(irrsam, "rU") as SAM:
		count = 0
		linecount = file_len(irrsam)
		for line, line2 in itertools.izip_longest(SAM,SAM,fillvalue=''):
			if linecount - count > 2:
				count += 2
				mapcheck = False
				mylist = line.split("\t")
				mylist2 = line2.split("\t")
				flag1, flag2 = int(mylist[1]), int(mylist2[1])
				if flag1 in read1_flags:
					myflag = flag1
					chrom, address = mylist[2], int(mylist[3])
					length = cigar_caller(mylist[5])
					myseq = mylist[9]
				elif flag2 in read1_flags:
					myflag = flag2
					chrom, address = mylist2[2], int(mylist2[3])
					length = cigar_caller(mylist2[5])
					myseq = mylist2[9]
				else:
					length = 0
				#Runs if first line contains the mapping data for the transposon read
				#Executes if transposon IRL read mapped to + strand
				if length >= 12:
					if myflag in read1_flags_pos:
						if myseq[0:2] == "TA":
							mytag = chrom +":"+ str(address) +"+" #same strand for IRR
							if mytag in results:
								results[mytag][1] += 1
								if results[mytag][1] > IRRmax:
									IRRmax = results[mytag][1]
								mapcheck = True
							else:
								results[mytag] = [0,1]
								mapcheck = True
							if mytag in ligation:
								if length not in ligation[mytag][1]:
									ligation[mytag][1].append(length)
							else:
								ligation[mytag] = [[],[]]
								ligation[mytag][1].append(length)
					elif myflag in read1_flags_neg:
						if myseq[-2:] == "TA":
							mytag = chrom +":"+ str(address+(length-2)) +"-" #flips strand for IRR
							if mytag in results:
								results[mytag][1] += 1
								if results[mytag][1] > IRRmax:
									IRRmax = results[mytag][1]
								mapcheck = True
							else:
								results[mytag] = [0,1]
								mapcheck = True
							if mytag in ligation:
								if length not in ligation[mytag][1]:
									ligation[mytag][1].append(length)
							else:
								ligation[mytag] = [[],[]]
								ligation[mytag][1].append(length)
					if mapcheck == False: #paired read only analyzed if the transposon read cannot be mapped precisely
						if flag1 in read2_flags:
							myflag = flag1
							chrom, address = mylist[2], int(mylist[3])
							length = cigar_caller(mylist[5])
							myseq = mylist[9]
						elif flag2 in read2_flags:
							myflag = flag2
							chrom, address = mylist2[2], int(mylist2[3])
							length = cigar_caller(mylist2[5])
							myseq = mylist2[9]
						if length >= 10:
							if myflag in read2_flags_pos:
								if myseq[0:2] == "TA":
									mytag = chrom +":"+ str(address) +"-"
									if mytag in results:
										results[mytag][1] += 1
										if results[mytag][1] > IRRmax:
											IRRmax = results[mytag][1]
										mapcheck = True
									else:
										results[mytag] = [0,1]
										mapcheck = True
									if mytag in ligation:
										if length not in ligation[mytag][1]:
											ligation[mytag][1].append(length)
									else:
										ligation[mytag] = [[],[]]
										ligation[mytag][1].append(length)
							elif myflag in read2_flags_neg:
								if myseq[-2:] == "TA":
									mytag = chrom +":"+ str(address+(length-2)) +"+" #flips strand since it's the IRL read
									if mytag in results:
										results[mytag][1] += 1
										if results[mytag][1] > IRRmax:
											IRRmax = results[mytag][1]
										mapcheck = True
									else:
										results[mytag] = [0,1]
										mapcheck = True
									if mytag in ligation:
										if length not in ligation[mytag][1]:
											ligation[mytag][1].append(length)
									else:
										ligation[mytag] = [[],[]]
										ligation[mytag][1].append(length)
	myfile = irlsam[:-8]
	output = open(myfile + ".uniq","wa")
	for key in results:
		mylist = key.split(":")
		chrom,address = mylist[0],mylist[1][:-1]
		ligationlist = ligation[key]
		IRLlp = len(ligationlist[0])
		IRRlp = len(ligationlist[1])
		strand = mylist[1][-1:]
		myvals = results[key]
		IRL,IRR = myvals[0],myvals[1]
		try:
			IRLnorm = round((float(IRL)/IRLmax*100),3)
		except:
			IRLnorm = 0
		try:
			IRRnorm = round((float(IRR)/IRRmax*100),3)
		except:
			IRRnorm = 0
		mystring = myfile +"\t"+ chrom +"\t"+ str(address) +"\t"+  str(IRL) +"\t"+ str(IRR) +"\t"+ strand +"\t"+ str(IRLlp) +"\t"+ str(IRRlp) +"\t"+ str(IRLnorm) +"\t"+ str(IRRnorm) +"\n"
		output.write(mystring)
	output.close()

mybarcode = open(barcode,"rU").readlines()
for line in mybarcode:
	infolist = line.replace("\n","").split("\t")
	mysample = infolist[0]
	irlF,irlR,irrF,irrR = infolist[1],infolist[2],infolist[3],infolist[4]
	genome = infolist[5]
	trimirlf1 = irlF[:-9] + "_5trim.fastq.gz"
	trimirlr1 = irlR[:-9] + "_5trim.fastq.gz"
	trimirlf2 = irlF[:-9] + "_53trim.fastq.gz"
	trimirlr2 = irlR[:-9] + "_53trim.fastq.gz"
	trimirlf3 = irlF[:-9] + "_trim.fastq.gz"
	trimirlr3 = irlR[:-9] + "_trim.fastq.gz"
	irlsam = mysample + "_IRL.sam"
	trimirrf1 = irrF[:-9] + "_5trim.fastq.gz"
	trimirrr1 = irrR[:-9] + "_5trim.fastq.gz"
	trimirrf2 = irrF[:-9] + "_53trim.fastq.gz"
	trimirrr2 = irrR[:-9] + "_53trim.fastq.gz"
	trimirrf3 = irrF[:-9] + "_trim.fastq.gz"
	trimirrr3 = irrR[:-9] + "_trim.fastq.gz"
	irrsam = mysample + "_IRR.sam"

	#Process IRL reads ---> trim adaptor ---> map reads
	cutadapt_IRL1(irlF, irlR)
	cutadapt_IRL2(trimirlf1, trimirlr1)
	os.system("rm " + trimirlf1)
	os.system("rm " + trimirlr1)
	cutadapt_IRL3(trimirlf2, trimirlr2)
	os.system("rm " + trimirlf2)
	os.system("rm " + trimirlr2)
	mapreads(genome, trimirlf3, trimirlr3, irlsam)


	#Process IRR reads ---> trim adaptor ---> map reads
	cutadapt_IRR1(irrF, irrR)
	cutadapt_IRR2(trimirrf1, trimirrr1)
	os.system("rm " + trimirrf1)
	os.system("rm " + trimirrr1)
	cutadapt_IRR3(trimirrf2, trimirrr2)
	os.system("rm " + trimirrf2)
	os.system("rm " + trimirrr2)
	mapreads(genome, trimirrf3, trimirrr3, irrsam)

	#Combine mapped IRL and IRR reads ---> generate .UNIQ file
	samparse(irlsam,irrsam)
	os.system("rm " + irlsam)
	os.system("rm " + irrsam)