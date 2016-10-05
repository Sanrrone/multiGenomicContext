from __future__ import with_statement 

# ==============================================================================
# 						multiGenomicContext
#
# Author: Sandro Valenzuela (sandrolvalenzuead@gmail.com) 
#
# Please type "python multiGenomicContext.py -h" for usage help
#
# ==============================================================================

__author__ = 'Sandro Valenzuela (sandrolvalenzuead@gmail.com)'
__version__ = '1.0'
__date__ = '10 August 2016'

import sys, os, re, subprocess, csv, glob
from operator import itemgetter
from collections import deque
from optparse import OptionParser
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

def printPlotStep(outfilename,totalgenes,totalgenomes):
	print outfilename
	plotstep=open("plotstep.R", 'w')
	plotstep.write("""
rm(list=ls())
library(ggplot2)
library(genoPlotR)
args<-commandArgs()
outfilename<-args[6]
totalgenes<-as.numeric(args[7])
totalgenomes<-as.numeric(args[8])

temp = list.files(pattern="*.DNASEGcsv")
if (length(temp)>1) {
  gbknames<- lapply(as.list(temp),function(x){strsplit(x = gsub(pattern = "[.]",replacement = " ",x = x),split = c(" "))[[1]][1]})
  
  df<-lapply(temp, read.csv, header = FALSE)
  df<-lapply(df,function(x){colnames(x)<-c("name", "start",  "end" ,"strand"  ,"col" ,"lty" ,"lwd" ,"pch" ,"cex", "gene_type");x})
  df<-lapply(df,function(x){dna_seg(x)})
  names(df)<-gbknames
  
}else{
  df<-read.csv(temp,header = F)
  colnames(df)<-c("name", "start",  "end" ,"strand"  ,"col" ,"lty" ,"lwd" ,"pch" ,"cex", "gene_type")
  df<-list(dna_seg(df))
}

uniqnames<-unique(do.call(rbind.data.frame, df)["name"])
uniqnames<-sort(uniqnames[,1])
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
colors<-sample(color, length(uniqnames))
df2color<-data.frame(as.matrix(uniqnames),as.matrix(colors))
df2color<-t(df2color)
colnames(df2color)<-df2color[1,]
df2color<- df2color[-1,]
df<-lapply(df,function(x){x["col"]<-df2color[x$name];x})

uniqnames<-gsub(pattern = "_",x = as.matrix(uniqnames),replacement = " ")
uniqnames<-gsub(pattern = "[.]",x = as.matrix(uniqnames),replacement = ",")
if(totalgenes>totalgenomes){
	pdf(file=outfilename, width = totalgenes, height = totalgenes)
}else{
	pdf(file=outfilename, width = totalgenes, height = totalgenomes)
}

par(mar=c(2,2,2,0))
plot(c(0,1000), c(0,1000), type="n", axes=FALSE, xlab="", ylab="")

legend("center", legend = c(as.matrix(uniqnames)),ncol = 1,xpd = NA, cex = 0.8,
       bty="n",fill=c(as.matrix(colors)),border = c("white"),title = "Genes")

plot_gene_map(dna_segs = df,dna_seg_label_cex = 0.9)

dev.off()""")

	plotstep.close()

	subprocess.call(["Rscript", "plotstep.R", str(outfilename), str(totalgenes), str(totalgenomes)])
	
	filenames = glob.glob('*.DNASEGcsv')
	for filename in filenames:
		os.remove(filename)
	
	os.remove("plotstep.R")

	return None

def foundGenomicContext(gene,faafile,upstream,downstream,GCX): #function to search genomic context
	#the faa files are with genes in order and formatted >gene|contig|position
	gene_list = [] #to save up and downstream genes
	faa_sequences = SeqIO.parse(open(faafile),'fasta')
	for proteins in faa_sequences:
		name = str(proteins.id)
		gene_list.append(name)

	#get the position of our gene
	gene_position=gene_list.index(gene)

	#calculate the total num of genes to print
	downstream=downstream+upstream
	#backup gene_position to change the color
	ourgene_position=gene_position
	
	#save the contigname of out gene
	contigname=str(gene_list[gene_position]).split("|")[2]

	#check if we are close to the begining of the list
	if (gene_position-upstream)<0:
		gene_position=0
		upstream=upstream-gene_position
	else:
		gene_position=gene_position-upstream


	outname=str(faafile).replace(".faa","")
	dna_segs=open(str(outname+".DNASEGcsv"),"w")

	while gene_position<len(gene_list) and downstream>=0:

		#only prints genes in the same contig of our gene
		if str(gene_list[gene_position]).split("|")[2] == contigname:
			genid=str(gene_list[gene_position]).split("|")[0]
			name=str(gene_list[gene_position]).split("|")[1]
			contig=str(gene_list[gene_position]).split("|")[2]
			pos1=str(gene_list[gene_position]).split("|")[3].split(":")[0]
			pos1=str(pos1).replace(">","").replace("<","")
			pos2=str(gene_list[gene_position]).split("|")[3].split(":")[1]
			pos2=str(pos2).replace(">","").replace("<","")
			strand=str(gene_list[gene_position]).split("|")[3]
			strand=str(strand).split(":")[2].replace("+","1").replace("-","-1")

			if gene_position == ourgene_position:
				color="red"
				GCX.write("%s,%s,%s,%s,%s,%s,%s\n" % (faafile,str(genid+"_match"),contig,name,pos1,pos2,strand))

			else:
				color="gray"
				GCX.write("%s,%s,%s,%s,%s,%s,%s\n" % (faafile,genid,contig,name,pos1,pos2,strand))


			#print name,pos1,pos2,strand,color

			dna_segs.write("%s,%s,%s,%s,%s,1,1,8,1,arrows\n" % (name, pos1, pos2, strand, color))

		gene_position=gene_position+1
		downstream=downstream-1

	dna_segs.close()

	return None


def which(program):#function to check if some program exists 
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def Makefaa(gbk):
	wd=os.getcwd()
	gbkname=gbk.replace("/"," ").split()[len(gbk.replace("/"," ").split())-1]
	protdict={}
	contiggene={}
	location={}
	location2={}
	location3={}
	faa= open(str(wd+"/"+gbkname+ ".faa"), 'w')
	recs = [rec for rec in SeqIO.parse(gbk, "genbank")]
	for rec in recs:
		contigname=rec[0:].id
		feats = [feat for feat in rec.features if feat.type == "CDS"]

		for feat in feats:
			if "product" in feat.qualifiers:
				cdsname=str(feat.qualifiers["locus_tag"]).replace("'","").replace("[","").replace("]","")
				product=str(feat.qualifiers["product"]).replace("'","").replace("[","").replace("]","").replace(" ","_").replace(",",".")
				locat1=str(feat.location).replace("[","").replace("]","").replace("(",":").replace(")","").split(":")[0]
				locat2=str(feat.location).replace("[","").replace("]","").replace("(",":").replace(")","").split(":")[1]
				locat3=str(feat.location).replace("[","").replace("]","").replace("(",":").replace(")","").split(":")[2]

				if "translation" in feat.qualifiers:	
					translation=str(feat.qualifiers["translation"]).replace("'","").replace("[","").replace("]","")
					band=0
					for seq in protdict.items():
						if seq[1] == translation:
							band=1

					if band==0:
						faa.write(">%s|%s|%s|%s\n" % (cdsname, product, contigname, str(locat1)+":"+str(locat2)+":"+locat3))
						faa.write("%s\n" % (translation))
	
	faa.close
	return str(gbkname + ".faa")

def main():

	parser = OptionParser(usage = "Usage: python multiGenomicContext.py -f protein.fasta -l MYgbklist.txt")
	parser.add_option("-f","--proteinFasta",dest="fastaProtein",help="default:none. your protein in fasta format to search on the gbk")
	parser.add_option("-l","--gbklist",dest="gbkList",help="List of the gbk (remember also have the files), also you can give the complete path in the list")
	parser.add_option("-u","--upstreamGenes",dest="Upstream",help="default:5 number of genes to search upstream on the gbks",default=4)
	parser.add_option("-d","--downstreamGenes",dest="Downstream",help="default:5 number of genes to search downstream on the gbks",default=4)
	parser.add_option("-e","--evalue",dest="evalue",help="default:1e-5 e-value for blastp search",default=1e-5)
	parser.add_option("-i","--identity",dest="Identity",help="default:85 range 1-100 % of identity on the blastp alignment to consider the gene exists on the genome",default=85)
	parser.add_option("-a","--alignmentLength",dest="alignL",help="default:75 range 1-100 % of aligment length to consider the gene exists on the genome",default=75)

	(options,args) = parser.parse_args()

	Inputprotein = options.fastaProtein
	gbkList= options.gbkList
	Upstream = int(options.Upstream)
	Downstream = int(options.Downstream)
	Evalue=str(options.evalue)
	Identity=int(options.Identity)
	alignL=int(options.alignL)


	#check variables
	if not Inputprotein:
		print "No input provided, use -h for help"
		sys.exit()

	if not gbkList:
		print "No gbk list provided, use -h for help"
		sys.exit()

	#searching for blastp
	blastpBIN=which("blastp")
	if blastpBIN == None:
		print "No blastp found, install it before continue"
		sys.exit()

	fasta_sequences = SeqIO.parse(open(Inputprotein),'fasta')
	gbks = open(gbkList,'r')

#################################################################################
	#get proteins from gbks
	print "Making .faa from gbk files"
	faafiles = [] #create list to save .faa
	for gbk in gbks:
		gbk=gbk.rstrip()#delete \n character
		name=Makefaa(gbk) #makefaa return the name of .faa (and create the file)
		faafiles.append(name)

	gbks.close()
#################################################################################

	#walk through the fastas
	for fasta in fasta_sequences:
		name, sequence = str(fasta.id), str(fasta.seq)
		print "Find",name,"in faa files"
		#making an individual fasta with the protein
		tmp=open('tmp.faa','w')
		tmp.write(">%s\n%s\n" % (name,sequence))
		tmp.close()

		GCX=open(str(name+".csv"),'w')
		GCX.write("source,genId,contig,name,start,end,strand\n")
		for faa in faafiles:
			subprocess.call([blastpBIN, "-query", "tmp.faa", "-subject", str(faa), "-out", "tmp.out", "-evalue", Evalue, "-outfmt", "10", "-max_target_seqs", "1", "-max_hsps", "1"])
			#now we check if the results pass the filter to consider the gene "exists" in the genome
			if os.path.getsize("tmp.out")>0:
				tmp=open("tmp.out","r")
				uniquerow=next(csv.reader(tmp))
				tmp.close()
				os.remove("tmp.out")
				#uniquerow[1] is the name of protein that match with our query (header of the fasta to be specific)
				#uniquerow[2] is identity
				#uniquerow[3] is alignment coverage (length)
				if uniquerow[2]>=Identity and (float(uniquerow[3])/len(sequence))>=(alignL/100.0):
					#if we are here, so, the protein exist in the gbk, the next step is find the genes up and down stream of the gbk
					#call the function
					foundGenomicContext(uniquerow[1],faa,Upstream,Downstream,GCX)

			else:
				print "No match found on gbk",str(">"+name),"for",faa
				if os.path.isfile("tmp.out"):
					os.remove("tmp.out")

		os.remove("tmp.faa")
		GCX.close()
		#call plot step
		#sys.exit()
		printPlotStep(str(name+".pdf"), Upstream+Downstream+1, len(faafiles))

	print "Clean files"
	for faa in faafiles:
		os.remove(faa)

	print "Done"



if __name__ == '__main__':
	main()
	sys.exit()