#!/usr/bin/python
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

import sys, os, re, subprocess, csv, glob
from operator import itemgetter
from collections import deque
from optparse import OptionParser
from Bio import SeqIO

def printPlotStep(outfilename,globalA,cleanProcess):
	print outfilename
	plotstep=open("plotstep.R", 'w')
	plotstep.write("""
rm(list=ls())
library(ggplot2)
library(genoPlotR)
args<-commandArgs()
outfilename<-args[6]
globalA<-ifelse(toupper(args[7])=="TRUE",TRUE,FALSE)

temp = list.files(pattern="*.DNASEGcsv")
nfiles<-length(temp)

#parse names
gbknames<- lapply(as.list(temp),function(x){strsplit(x = x,split = "[.]DNASEGcsv")[[1]][1]})
gbknames<- lapply(gbknames,function(x){gsub(pattern = "[.]gbff",replacement = "",x = x)})

#read regions of interest
df<-lapply(temp, read.csv, header = F, stringsAsFactors = F)


df<-lapply(df,function(x){
  x<-unique(x)
  colnames(x)<-c("name", "start",  "end" ,"strand"  ,"col" ,"lty" ,"lwd" ,"pch" ,"cex", "gene_type","locus_tag","contig")
  name<-strsplit(x = x$name,split = "_|-")
  x$name<-unlist(lapply(name,function(y){
    halfy<-round(length(y)/2,0)
    if(halfy+1>=3){
      y<-paste0(paste(y[1:halfy],collapse = " "),"\n",
                paste(y[(length(y)-(halfy-1)):length(y)],collapse = " ")
      )
    }else{
      y<-paste(y,collapse = " ")
    }
    y
  }))
  x
})

df<-lapply(df,function(x){
  x<-x[order(x$contig),]
  prevContig<-x[1,"contig"]
  newx<-list()
  
  for(c in unique(x$contig)){
    tmp<-subset(x, contig == c)
    tmp<-tmp[order(tmp$start),]
    if(tmp[1,"contig"] == prevContig){
      endpos<-tmp[nrow(tmp),"end"] + 3001
    }else{
      if(nrow(tmp)>=2){
        prevEnd<-0
        nextStart<-0
        for(i in 1:(nrow(tmp)-1)){
          distGene<- nextStart - prevEnd
          posdif<-tmp[i,"end"] - tmp[i,"start"]
          tmp[i,"start"] <- endpos + distGene
          tmp[i,"end"] <- tmp[i,"start"] + posdif
          endpos<-tmp[i,"end"]
          prevEnd<-tmp[i,"end"]
          nextStart<-tmp[i+1,"start"]
        }
      }else{
        posdif<-tmp[1,"end"] - tmp[1,"start"]
        tmp[1,"start"] <- endpos
        tmp[1,"end"] <- tmp[1,"start"] + posdif
        endpos<-tmp[1,"end"]
      }
      
      prevContig<- c
    }
    newx[[c]]<-tmp
  }
  x<-bind_rows(newx)
})

annot<-lapply(df,function(x){

  annotation(x1=x$start+10,x2=x$end-5,text=x$name,rot=replicate(nrow(x),35))
  #annotation(x1=x$start,x2=x$end,text=x$locus_tag,rot=replicate(nrow(x),35))
})

xlims<-lapply(df,function(x){
  x<-x[order(x$start),]
  lims<-c(ifelse(min(x[1,"start"])-100<=0,0,min(x[,"start"])-100))
  for(i in 1:(nrow(x)-1)){
    if(x[i+1,"start"]-x[i,"end"] >= 3000){
      lims<-c(lims,x[i,"end"]+100,
              x[i+1,"start"]-100)
    }
  }
  lims<-c(lims,max(x[,"end"])+50)
})

#set unique colors for genes
uniqnames<-unique(do.call(rbind.data.frame, df)["name"])
uniqnames<-sort(uniqnames[,1])
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
colors<-sample(color, length(uniqnames))
df2color<-data.frame(as.matrix(uniqnames),as.matrix(colors))
df2color<-t(df2color)
colnames(df2color)<-df2color[1,]
df2color<- df2color[-1,]
df<-lapply(df,function(x){
  x["fill"]<-df2color[x$name]
  x["col"]<-"black"
  #x["gene_type"]<-"side_blocks"
  #tmp<-x["name"]
  #x["name"]<-x["locus_tag"]
  #x["locus_tag"]<-tmp
  x
})


df<-lapply(df,function(x){dna_seg(x)})

if(nfiles>1){
  names(df)<-gbknames
}


wformula=as.integer(log(max(sapply(annot,nrow)))*log(max(sapply(annot,nrow)))*2)+max(sapply(xlims, length))
hformula=as.integer(log(nfiles)*nfiles)+1
pdf(file=outfilename, width = wformula, height = hformula)

par(mar=c(0,3,2,3))
plot(c(0,1000), c(0,1000), type="n", axes=FALSE, xlab="", ylab="")

legend("center", legend = gsub("_"," ",c(as.matrix(uniqnames))),ncol = as.integer(length(uniqnames)/20)+1,xpd = NA, 
       cex = 0.8, bty="n",fill=c(as.matrix(colors)),border = c("white"),title = "Genes")

if(globalA){
  #read mauve comparison
  mauvebb<-read_mauve_backbone("tmpbb.mauve")
  plot_gene_map(dna_segs = mauvebb$dna_segs,dna_seg_label_cex = 0.8,
                comparisons = mauvebb$comparisons,
                dna_seg_scale = T)
}else{
    plot_gene_map(dna_segs = df,dna_seg_label_cex = 0.5,annotation_height = round(2+log(hformula),0),
                annotations = annot, xlims = xlims,
                scale = F, dna_seg_scale = T,plot_new=T)
}


dev.off()


""")

	plotstep.close()
	RBIN=which("Rscript")
	subprocess.call([RBIN, "plotstep.R", str(outfilename), str(globalA)])
	if not cleanProcess:
		return None


	filenames = glob.glob('*.DNASEGcsv')
	for filename in filenames:
		os.remove(filename)
	
	os.remove("plotstep.R")

	return None

def foundGenomicContext(gene,faafile,upstream,downstream,GCX,dna_segs): #function to search genomic context
	#the faa files are with genes in order and formatted >gene|locustag|contig|position
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
				#GCX.write("%s,%s,%s,%s,%s,%s,%s\n" % (faafile,str(genid+" (input)"),contig,name,pos1,pos2,strand))
				GCX.write("%s,%s,%s,%s,%s,%s,%s\n" % (faafile,genid,contig,name,pos1,pos2,strand))

			else:
				color="gray"
				GCX.write("%s,%s,%s,%s,%s,%s,%s\n" % (faafile,genid,contig,name,pos1,pos2,strand))


			#print name,pos1,pos2,strand,color
			if gene_position == ourgene_position:
				#dna_segs.write("%s,%s,%s,%s,%s,1,1,8,1,arrows,%s,%s\n" % (str(name+" (input)"), pos1, pos2, strand, color, genid, contig))
				dna_segs.write("%s,%s,%s,%s,%s,1,1,8,1,arrows,%s,%s\n" % (name, pos1, pos2, strand, color, genid, contig))

			else:
				dna_segs.write("%s,%s,%s,%s,%s,1,1,8,1,arrows,%s,%s\n" % (name, pos1, pos2, strand, color, genid, contig))

		gene_position=gene_position+1
		downstream=downstream-1


	return None


def which(program):#function to check if some program exists 
  
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
				if "gene" in feat.qualifiers:
					product=str(feat.qualifiers["gene"]).replace("'","").replace("[","").replace("]","")

				else:
					product=str(feat.qualifiers["product"]).replace("'","").replace("[","").replace("]","").replace(" ","_").replace(",",".")

				locationStart = list()
				locationEnd = list()
				locationStrand = list()
				i=0
				if "join" in str(feat.location):
					for gensection in str(feat.location).replace("[","").replace("]","").replace("{","").replace("}","").replace("join","").replace(" ","").replace(">","").replace("<","").split(","):
						locationStart.append(gensection.split(":")[0])
						locationEnd.append(gensection.split(":")[1].split("(")[0])
						locationStrand.append(gensection.split(":")[1].split("(")[1].replace(")",""))

						i=i+1

				else:
					featlocation=str(feat.location).replace(">","").replace("<","").replace("[","").replace("]","").replace("(",":").replace(")","")
					locationStart.append(featlocation.split(":")[0])
					locationEnd.append(featlocation.split(":")[1])
					locationStrand.append(featlocation.split(":")[2])
					i=i+1

				if "translation" in feat.qualifiers:
					translation=str(feat.qualifiers["translation"]).replace("'","").replace("[","").replace("]","")
					band=0
					for seq in protdict.items():
						if seq[1] == translation:
							band=1

					if band==0:
						for i in range(0,i):
							#print str(cdsname+" || "+product+" || "+contigname+" || "+locationStart[i]+" || "+locationEnd[i]+" || "+locationStrand[i])
							faa.write(">%s|%s|%s|%s\n" % (cdsname, product, contigname, str(int(locationStart[i])+1)+":"+locationEnd[i]+":"+locationStrand[i]))
							faa.write("%s\n" % (translation))
	
	faa.close()
	return str(gbkname + ".faa")

def main():

	parser = OptionParser(usage = "Usage: python multiGenomicContext.py -f protein.fasta -g mygbff.gbff [-w if your protein fasta have 2>= proteins]")
	#parser.add_option("-g","--global", dest="globalA",help="default:False plot gbk-gbk aligment instead only a region. use only with -l parameter", default=False, action='store_true')
	parser.add_option("-f","--proteinFasta",dest="proteinFasta",help="default:none. your protein in fasta format to search on the gbk/gbff")
	parser.add_option("-g","--gbklist",dest="gbkList",help="Comma separated gbk, for example: mygbk1.gbk,mygbk2.gbk,mygbk3.gbk")
	parser.add_option("-u","--upstreamGenes",dest="Upstream",help="default:5 number of genes to search upstream on the gbks",default=4)
	parser.add_option("-d","--downstreamGenes",dest="Downstream",help="default:5 number of genes to search downstream on the gbks",default=4)
	parser.add_option("-e","--evalue",dest="evalue",help="default:1e-5 e-value for blastp search",default=1e-5)
	parser.add_option("-i","--identity",dest="Identity",help="default:85 range 1-100 % of identity on the blastp alignment to consider the gene exists on the genome",default=85)
	parser.add_option("-a","--alignmentLength",dest="alignL",help="default:75 range 1-100 % of aligment length to consider the gene exists on the genome",default=85)
	parser.add_option("-b","--blastpBIN", dest="blastpBIN",help="default:/usr/bin/blastp blastp binary path", default="/usr/bin/blastp")
	parser.add_option("-m","--progressiveMauveBIN", dest="progressiveMauveBIN",help="default:/usr/bin/progressiveMauve mauve binary path", default="/usr/bin/progressiveMauve")
	parser.add_option("-c","--cleanProcessOff", dest="cleanProcess",help="default: True plot this kind of files is complex, so if you turn this flag False, you will have the R file to manipulate the plots", default=True, action='store_false')
	parser.add_option("-w","--wholeGenomicInput", dest="wholeGenomicInput",help="default: False by default the script plot one chart/csv per input sequence, with this parameter all proteins input are in the same chart (per gbk)", default=False, action='store_true')

	(options,args) = parser.parse_args()

	#globalA = options.globalA
	Inputprotein = options.proteinFasta
	gbkList= options.gbkList
	Upstream = int(options.Upstream)
	Downstream = int(options.Downstream)
	Evalue=str(options.evalue)
	Identity=int(options.Identity)
	alignL=int(options.alignL)
	blastpBIN=options.blastpBIN
	mauveBIN=options.progressiveMauveBIN
	cleanProcess=options.cleanProcess
	wholeGenomicInput=options.wholeGenomicInput
    
	globalA=False
	#check variables
	if not Inputprotein:
		if globalA is not True:
			print "No input provided (-f), use -h for help"
			sys.exit()
	else:
		if not os.path.isfile(Inputprotein):
			print str("* "+Inputprotein+"doesn't exist, check the file directory")
		if globalA:
			print "-g/--global is only vaild with -l/--gbklist option"
			sys.exit()

	if gbkList is None:
		print "No gbk list provided (-l), use -h for help"
		sys.exit()
	else:
		gbkList=str(gbkList).split(",")

	#searching for blastp
	if which(blastpBIN) is None:
		print "No blastp found, install it before continue or use --blastpBIN for custom binary path"
		sys.exit()

	#searching for mauve
	if globalA and which(mauveBIN) is None:
		print "No mauveAligner found, install it before continue or use progressiveMauveBIN for custom progressiveMauve binary path"
		sys.exit()

	if Upstream+Downstream >= 433 and globalA is not True:
		#max number of colors for R script
		print "too much genes for plot, use the option --global for entire sequences"
		sys.exit()

	RBIN=which("Rscript")
	if RBIN == None:
		print "No Rscript binary found, install it before continue"
		sys.exit()

	Inputprotein=os.path.abspath(Inputprotein)
	inputProteins = SeqIO.parse(open(Inputprotein),'fasta')

#################################################################################
	#get proteins from gbks
	print "Making .faa from gbk files"
	for i in range(0,len(gbkList)):
		if os.path.isfile(gbkList[i]):
			gbkList[i]=os.path.abspath(gbkList[i])
		else:
			print str("* "+gbkList[i]+"doesn't exist")
			sys.exit()

	gbkfaafiles = [] #create list to save .faa
	for gbk in gbkList:
		gbk=gbk.rstrip() #delete \n character
		name=Makefaa(gbk) #makefaa return the name of .faa (and create the file)
		gbkfaafiles.append(name)

#################################################################################

	#genome-genome aligment
	if globalA:
		command=str(mauveBIN+" --output=tmp.mauve --backbone-output=tmpbb.mauve --seed-family --muscle-args='-refine' "+str(" ".join(gbkList)))
		subprocess.call(command, shell=True)
		printPlotStep(str(name+".pdf"), globalA,cleanProcess)

	else:
		if wholeGenomicInput:
			#walking through the fastas and genes
			for faa in gbkfaafiles:
				GCX=open(str(faa+".csv"),'w')
				GCX.write("source,genId,contig,name,start,end,strand\n")

				outname=str(faa).replace(".faa","")
				dna_segs=open(str(outname+".DNASEGcsv"),"w")
				print "working on "+faa
				for fasta in SeqIO.parse(open(Inputprotein),'fasta'):
					name, qsequence = str(fasta.id), str(fasta.seq)
					#making an individual fasta with the protein
					tmp=open('tmp.faa','w')
					tmp.write(">%s\n%s\n" % (name,qsequence))
					tmp.close()


					command=str(blastpBIN+" -query tmp.faa -subject "+str(faa)+" -out tmp.out -evalue "+Evalue+" -outfmt 10")
					subprocess.call(command, shell=True)
					os.remove("tmp.faa")
					#now we check if the results pass the filter to consider the gene "exists" in the genome
					if os.path.getsize("tmp.out")>0:
						tmp=open("tmp.out","r")
						for uniquerow in csv.reader(tmp,delimiter=','):
							#uniquerow[0] is our query protein
							#uniquerow[1] is the name of protein that match with our query (header of the fasta to be specific)
							#uniquerow[2] is identity
							#uniquerow[3] is alignment coverage (length)
							if uniquerow[2]>=Identity and (float(uniquerow[3])/len(qsequence))>=(alignL/100.0):
								#if we are here, so, the protein exist in the gbk, the next step is find the genes up and down stream of the gbk
								#call the function
								#print(uniquerow)
								foundGenomicContext(uniquerow[1],faa,0,0,GCX,dna_segs)
						tmp.close()
						os.remove("tmp.out")

					else:
						print "No match found in gbk",str(">"+name),"for",faa
						if os.path.isfile("tmp.out"):
							os.remove("tmp.out")

				GCX.close()
				dna_segs.close()
				#call plot step
			printPlotStep(str("Gcontext"+".pdf"), globalA, cleanProcess)
		else:
			#walking through the fastas and genes
			for fasta in inputProteins:
				name, sequence = str(fasta.id), str(fasta.seq)
				print "Find",name,"in faa files"
				#making an individual fasta with the protein
				tmp=open('tmp.faa','w')
				tmp.write(">%s\n%s\n" % (name,sequence))
				tmp.close()

				GCX=open(str(name+".csv"),'w')
				GCX.write("source,genId,contig,name,start,end,strand\n")
				for faa in gbkfaafiles:
					command=str(blastpBIN+" -query tmp.faa -subject "+str(faa)+" -out tmp.out -evalue "+Evalue+" -outfmt 10 -max_target_seqs 1 -max_hsps 1")
					subprocess.call(command, shell=True)
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
							outname=str(faa).replace(".faa","")
							dna_segs=open(str(outname+".DNASEGcsv"),"w")
							foundGenomicContext(uniquerow[1],faa,Upstream,Downstream,GCX,dna_segs)
							dna_segs.close()

					else:
						print "No match found on gbk",str(">"+name),"for",faa
						if os.path.isfile("tmp.out"):
							os.remove("tmp.out")

				os.remove("tmp.faa")
				GCX.close()
				#call plot step
				printPlotStep(str(name+".pdf"), globalA, cleanProcess)

	print "Clean files"
	for faa in gbkfaafiles:
		os.remove(faa)

	if os.path.exists("*.mauve"):
		os.remove("*.mauve")

	print "Done"



if __name__ == '__main__':
	main()
	sys.exit()