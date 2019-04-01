#!/usr/bin/python
from __future__ import division
from scipy.stats import chisquare
from scipy.stats import skew
from scipy.stats import kurtosis
import cPickle as pickle
import statsmodels.sandbox.stats.multicomp as sm
import os,sys,glob,time

myfile = sys.argv[1]
promoter = int(sys.argv[2])
geneset = sys.argv[3]  #can be "ensembl", "ucsc", or "refseq"
genome = "GRCh38"

#counts number of TAs in reference genome using FASTA files in given directory
def countTA(genome,pipepath):
	os.chdir(pipepath +"/"+ genome)
	myta = open("TA_count.txt","rU").readlines()
	TAcount = int(myta[0].replace("\n",""))
	return TAcount

#takes GFF3 file as input, builds and returns a dictionary of insertion sites for each sample in the file
#structure of dictionary is:  sampleDB[sampleID][chr][list of sites]
#a second dictionary is also built and returned that contains the # of insertions for each sample
def buildsampleDB(myfile):
	sampleDB,samplecounts = {},{}
	platelist = []
	myinput = open(myfile,"rU").readlines()
	for line in myinput:
		mylist = line.replace("\n","").split("\t")
		chrom = mylist[0]
		address = int(mylist[3])
		mysample = mylist[8]
		plate = mysample.split("*")[0]
		if plate not in platelist:
			platelist.append(plate)
		strand = mylist[6]
		if strand == "-":
			address = address * -1
		if mysample not in sampleDB:
			sampleDB[mysample] = {}
			sampleDB[mysample][chrom] = [address]
		elif chrom not in sampleDB[mysample]:
			sampleDB[mysample][chrom] = [address]
		else:
			sampleDB[mysample][chrom].append(address)
		if mysample not in samplecounts:
			samplecounts[mysample] = 1
		else:
			samplecounts[mysample] += 1
	samplenumber = len(platelist)
	return sampleDB,samplecounts,samplenumber

#this simply builds and returns a list of chromosomes found in the reference genome specified by the user
def buildchromList(genome,pipepath):
	os.chdir(pipepath +"/"+ genome)
	chromlist = []
	myinput = open("chromosomes.txt","rU").readlines()
	for line in myinput:
		mychrom = line.replace("\n","")
		chromlist.append(mychrom)
	return chromlist

#this function simply counts and returns the number genes that have at least one annotated insertion event
#the function takes an input of the results dictionary built in the pipeline below
def countgenehits(resultsDB):
	genehits = 0
	for key in resultsDB:
		hitinfo = resultsDB[key]
		hits = int(len(hitinfo)-1)
		if hits > 1:
			genehits += 1
	return genehits

#this function takes four inputs and returns 6 different values that are used to assess the quality of
#each gene hit. These values are used to detect strand bias in the insertion sites and to identify
#false positive hits
def hiteval(gene,hits,resultsDB,geneonlyDB):
	if gene in geneonlyDB:
		genebodyhits = len(geneonlyDB[gene]) - 1
	else:
		genebodyhits = 0
	promoterhits = hits - genebodyhits
	if gene in geneonlyDB:
		genebodyplus = int(geneonlyDB[gene][0][4])
		genebodyneg = int(geneonlyDB[gene][0][5])
	else:
		genebodyplus = 0
		genebodyneg = 0
	promoterplus = int(resultsDB[gene][0][4]) - genebodyplus
	promoterneg = int(resultsDB[gene][0][5]) - genebodyneg
	return genebodyhits,promoterhits,genebodyplus,genebodyneg,promoterplus,promoterneg

#this function takes the ordered list of sites within each gene and calculates the standard deviation of the
#distance between the sites. This is meant to detect clustering within genes that is indicitive of a gene
#truncation event. 	
def cluster(sitelist):
	distance = []
	for x in range(0,len(sitelist)-1):
		space = sitelist[0] - sitelist[x]
		distance.append(int(space))
	geneskew = skew(distance)
	genekurt = kurtosis(distance)
	return geneskew,genekurt

mypath = os.getcwd() #sets path variable containing the GFF3 file
timestamp = str(time.clock())
resultsdir = "gCIS_results_" + timestamp.split(".")[1]
os.mkdir(mypath +"/"+ resultsdir)
pipepath = os.path.dirname(os.path.realpath(__file__)) #sets path variable for directory containing pipeline
genepath = mypath +"/"+ genome +"/"+ geneset
parameters = open(mypath +"/"+ resultsdir +"/"+ resultsdir + "_parameter_settings.txt","wa")
parameters.write("input filename: " + myfile +"\npromoter size (bp): " + str(promoter) +"\ngenome: "+genome +"\ngene annotation: " +geneset+"\n")
parameters.close()

#these are the only allowed promoter sizes since the gene dictionaries are pre-built using these values
promoters = [0,5000,10000,15000,20000,25000,30000,35000,40000] 
promoterref = {0:3,5000:4,10000:5,15000:6,20000:7,25000:8,30000:9,35000:10,40000:11}
promindex = promoterref[promoter]
if promoter not in promoters:
	sys.exit("allowed promoter sizes are: 0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, and 40000")
if geneset not in ["ensembl","ucsc","refseq"]:
	sys.exit("allowed gene annotation types: ensembl, ucsc, refseq")

sampleDB,samplecounts,samplenumber = buildsampleDB(myfile)  #builds sample dictionary
TAcount = countTA(genome,pipepath) #assigns variable containing # of TA sites in reference genome
chromlist = buildchromList(genome,pipepath) #generates list of chromosomes in the reference genome

#generate output file for all genes with at least one insertion event
allhits = open(mypath +"/"+ resultsdir +"/"+ myfile[:-5] + "_all_genes.txt","wa")
allhits.write("gene_symbol\tchr\tstart (bp)\tstop (bp)\tstrand\t# of hits\t# of samples\tp-value\tFDR\t# in gene body [+,-]\t# in promoter [+,-]\tskewness\tkurtosis\tfunctional prediction\tsamples with insertion\n")

#generate output file for filtered results (i.e. adjusted p-value <= 0.05)
trimhits = open(mypath +"/"+ resultsdir +"/"+ myfile[:-5] + "_filtered_genes.txt","wa")
trimhits.write("gene_symbol\tchr\tstart (bp)\tstop (bp)\tstrand\t# of hits\t# of samples\tp-value\tFDR\t# in gene body [+,-]\t# in promoter [+,-]\tskewness\tkurtosis\tfunctional prediction\tsamples with insertion\n")
os.chdir(pipepath +"/"+ genome)

#this begins a section that builds the results and stats dictionaries containing the data for the run
resultsDB = {}
statsDB = {}
for entry in chromlist:
	mygenes = pickle.load(open(genepath +"/"+ entry,"rb"))   #[start,stop,strand,0k,5k,10k,15k,20k,25k,30k,35k,40k]
	chrom = entry.replace("_gene_TAs.p","")
	for gene in mygenes:
		geneinfo = mygenes[gene][0:3]
		geneTAs = mygenes[gene][promindex]
		geneprob = geneTAs / TAcount
		geneinfo.append(chrom) #structure of list = [start,stop,strand,chrom]
		geneinfo.extend([0,0,'','']) #structure of new list = [start,stop,strand,chrom,+stand count,-strand count,skewness,kurtosis]
		strand = geneinfo[2]
		if strand == "+":
			start = int(mygenes[gene][0]) - promoter
			stop = int(mygenes[gene][1])
		elif strand == "-":
			start = int(mygenes[gene][0])
			stop = int(mygenes[gene][1]) + promoter
		genetotalprob = 0
		sitelist = [start,stop]
		for sample in sampleDB:
			if chrom in sampleDB[sample]:
				myhits = sampleDB[sample][chrom]
				for hit in myhits:
					if abs(int(hit)) >= start and abs(int(hit)) <= stop:
						sitelist.append(abs(int(hit)))
						if gene not in resultsDB:
							resultsDB[gene] = [geneinfo,sample]
							if int(hit) < 0:
								resultsDB[gene][0][5] += 1
							else:
								resultsDB[gene][0][4] += 1
						elif sample not in resultsDB[gene]:
							resultsDB[gene].append(sample)
							if int(hit) < 0:
								resultsDB[gene][0][5] += 1
							else:
								resultsDB[gene][0][4] += 1
			sampleprob = geneprob * samplecounts[sample]
			genetotalprob += sampleprob
		statsDB[gene] = [genetotalprob]
		sitelist.sort()
		if len(sitelist) > 3:
			hitskew,hitkurt = cluster(sitelist)
			resultsDB[gene][0][6] = hitskew
			resultsDB[gene][0][7] = hitkurt
		elif len(sitelist) == 1:
			resultsDB[gene][0][6] = "---"
			resultsDB[gene][0][7] = "---"
			
#build dictionary containing hits in gene body only (i.e no promoter)
#also adds expected number of hits in gene body only to statsDB to be used to identify false positives later
geneonlyDB = {}
for entry in chromlist:
	#mygenes[gene] = [start,stop,strand,0k,5k,10k,15k,20k,25k,30k,35k,40k]
	mygenes = pickle.load(open(genepath +"/"+ entry,"rb")) 
	chrom = entry.replace("_gene_TAs.p","")
	for gene in mygenes:
		geneinfo = mygenes[gene][0:3]
		geneTAs = mygenes[gene][3]
		geneprob = geneTAs / TAcount
		geneinfo.append(chrom) #structure of list = [start,stop,strand,chrom]
		geneinfo.extend([0,0]) #structure of new list = [start,stop,strand,chrom,+stand count,-strand count]
		start = int(mygenes[gene][0])
		stop = int(mygenes[gene][1])
		genetotalprob = 0
		for sample in sampleDB:
			if chrom in sampleDB[sample]:
				myhits = sampleDB[sample][chrom]
				for hit in myhits:
					if abs(int(hit)) >= start and abs(int(hit)) <= stop:
						if gene not in geneonlyDB:
							geneonlyDB[gene] = [geneinfo,sample]
							if int(hit) < 0:
								geneonlyDB[gene][0][5] += 1
							else:
								geneonlyDB[gene][0][4] += 1
						elif sample not in geneonlyDB[gene]:
							geneonlyDB[gene].append(sample)
							if int(hit) < 0:
								geneonlyDB[gene][0][5] += 1
							else:
								geneonlyDB[gene][0][4] += 1
			sampleprob = geneprob * samplecounts[sample]
			genetotalprob += sampleprob
		statsDB[gene].append(genetotalprob)

genecount = countgenehits(resultsDB) #needed for p-value correction below
printDB = {}
peas,hitgenes = [],[]
for key in resultsDB:
	#sets variables for printing results to output files
	hitinfo = resultsDB[key]
	hits = int(len(hitinfo)-1)
	samplesize = len(sampleDB)
	samplelist = ",".join(map(str,hitinfo[1:]))
	samplehits = len(hitinfo)-1
	hitskew = hitinfo[0][6]
	hitkurt = hitinfo[0][7]
	samplenohits = samplesize - samplehits
	genebodyhits,promoterhits,genebodyplus,genebodyneg,promoterplus,promoterneg = hiteval(key,hits,resultsDB,geneonlyDB)
	
	#generates strings for output
	genebodystring = str(genebodyhits) +" ["+ str(genebodyplus) +","+ str(genebodyneg) +"]"
	promoterstring = str(promoterhits) +" ["+ str(promoterplus) +","+ str(promoterneg) +"]"
	
	#performs chi-square test for each gene
	expecthits = statsDB[key][0]
	if expecthits == 0:
		expecthits = 0.0001
	expectnohits = samplesize - expecthits
	mychitest = chisquare([samplehits,samplenohits],[expecthits,expectnohits])
	
	#determines if each gene is significant without promoter hits
	#this is needed to identify false positive genes with promoter hits without orientation bias
	if key in geneonlyDB:
		geneonlyhits = len(geneonlyDB[key])-1
	else:
		geneonlyhits = 0.0001
	geneonlynohits = samplesize - geneonlyhits
	geneonlyexpect = statsDB[key][1]
	geneonlyexpectrest = samplesize - geneonlyexpect
	geneonlychisq = chisquare([geneonlyhits,geneonlynohits],[geneonlyexpect,geneonlyexpectrest])	
	
	#adjusts p-value for all genes with > 1 hit that have a "0" p-value
	if samplehits > 1:
		rawp = mychitest[1]
		if rawp == 0:
			rawp = 10**-300
	else:
		rawp = "----"
	
	#sets variables to estimate prediction of functional impact of insertion on each gene
	plusperc = float(((genebodyplus + promoterplus) / hits) * 100)
	negperc = float(((genebodyneg + promoterneg) / hits) * 100)
	genebodyplus = float((genebodyplus / hits)*100)
	genebodyneg = float((genebodyneg / hits)*100)

	#prepares final string to write to outpute file
	mystring = str(hitinfo[0][3]) +"\t"+ str(hitinfo[0][0]) +"\t"+ str(hitinfo[0][1]) +"\t"+ str(hitinfo[0][2]) +"\t"+ str(samplehits) +"\tsample#\tp-value\tFDR\t"+ genebodystring +"\t"+ promoterstring +"\t"+ str(hitskew) +"\t"+ str(hitkurt) +"\tprediction\t"+ samplelist +"\n"	
	#printDB guide = mystring, rawp, FDR, hits, strand, plusperc, negperc, genebodyplus, genebodyneg, mychitest[1], geneonlychisq[1], prediction
	printDB[key] = [mystring, rawp, "", hits, str(hitinfo[0][2]), plusperc, negperc, genebodyplus, genebodyneg, mychitest[1], geneonlychisq[1],""]
	if rawp != "----":
		peas.append(rawp)
		hitgenes.append(key)

# calculate FDR for p-values
multitest = sm.multipletests(peas, alpha=0.05, method="fdr_bh", is_sorted=False, returnsorted=False)
fdrs = multitest[1]
for x in range(0,len(hitgenes)):
	thisgene = hitgenes[x]
	pval,myfdr = peas[x],fdrs[x]
	if pval == printDB[thisgene][1]:
		printDB[thisgene][2] = myfdr
	else:
		print "error: p-values don't match"

for key in printDB:
	pval = printDB[key][1]
	if printDB[key][2] != "":
		myfdr = printDB[key][2]
	else:
		myfdr = "----"
	prediction = "----"
	if printDB[key][11] == "" and printDB[key][3] > 4:
		if printDB[key][4] == "+":
			if printDB[key][5] >= 75 and myfdr <= 10**-5:
				prediction = "over-expression"
			elif printDB[key][5] < 75 and printDB[key][5] > 25:
				if (printDB[key][7] + printDB[key][8]) >= 90 and myfdr <= 10**-5:
					prediction = "disruption"
				else:
					prediction = "false positive"
		elif printDB[key][4] == "-":
			if printDB[key][6] >= 75 and myfdr <= 10**-5:
				prediction = "over-expression"
			elif printDB[key][6] < 75 and printDB[key][6] > 25:
				if (printDB[key][7] + printDB[key][8]) >= 90 and myfdr <= 10**-5:
					prediction = "disruption"
				else:
					prediction = "false positive"
	elif printDB[key][11] == "" and printDB[key][3] < 3:
		printDB[key][11] = "----"
	mystring = printDB[key][0]
	newstring1 = mystring.replace("p-value", str(pval))
	newstring2 = newstring1.replace("FDR", str(myfdr))
	newstring3 = newstring2.replace("prediction", prediction)
	tempstring = newstring3.replace("\n","").split("\t")[13]
	hitlist = tempstring.split(",")
	thesamplelist = []
	for item in hitlist:
		sample = item.split("*")[0]
		if sample not in thesamplelist:
			thesamplelist.append(sample)
	finalstring = newstring3.replace("sample#", str(len(thesamplelist)))
	if (len(thesamplelist) / samplenumber) >= 0.05:
		if prediction != "----" and prediction != "false positive" and myfdr != "----" and myfdr <= 10**-5:
			trimhits.write(key +"\t"+ finalstring)
			allhits.write(key +"\t"+ finalstring)
		else:
			allhits.write(key +"\t"+ finalstring)
	else:
		allhits.write(key +"\t"+ finalstring)
trimhits.close()
allhits.close()