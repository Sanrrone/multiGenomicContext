# multiGenomicContext
---------------------
multiGenomicContext is a python + R script that plot the genomic context of a protein on the genome that you want (or many genomes). You only needs two things: a fasta file and a list of your gbk to find the genomic context (and the gbk's too).

# Output
![Banner](https://github.com/Sanrrone/multiGenomicContext/blob/master/example/sample.png)


#Requisites
* Python >= 2.7 with the module [Biopython](http://biopython.org/wiki/Download)
* [blastp binary](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* R (tested in 3.3.1), with modules:
	* ggplot2: ```install.packages("ggplot2")```
	* genoPlotR: ```install.packages("genoPlotR", repos="http://R-Forge.R-project.org")```

# Usage

multiGenomicContext have a minimal use:
	
	python multiGenomicContext.py -f protein.fasta -l gbklist.txt
	
and a complete use:

	python multiGenomicContext.py -f protein.fasta -l gbklist.txt -u 4 -d 4 -e 1e-5 -i 85 -a 75 -s 15
	
where the options are:

* -f: The protein sequence in fasta format (also can be a multifasta of proteins).
* -l: The list of gbk (if you have one, also put that name in a file)
* -u: Number of genes to put the genomic context in upstream search (default: 4)
* -d Number of genes to put the genomic context in downstream search (default: 4)
* -e E-value for blastp search (default: 1e-5)
* -i Identity % of the alignment on blastp results to consider the gene exists on the genome (default: 85)
* -a Alignment length (%) between gene and the match for blastp search to consider the gene "exists" on the genome (default 75)
* -s Number of character of labels genes (default 15)

# How can I do a list?
The simple way is write a txt name by name. Or do in a terminal:
		
	ls -1 *.gbk > myGbkList.txt

# Notes
* multiGenomicContext search genes on the gbk because the cds are ordered, but this true only for one chromosome assembly, for gbk files where two or more contigs exists, it's show genomic context for the same contig of the gene.
* We strongly recommend use one software to annotate all gbk's to conserve the same names outcome.
