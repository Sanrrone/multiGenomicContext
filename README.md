# multiGenomicContext
----------------------
multiGenomicContext is a python + R script to plot the genomic context of a protein on the genome that you want (or many genomes). You only needs two things: a fasta file and the gbk to find the genomic context.

# Output
![Banner](https://github.com/Sanrrone/multiGenomicContext/blob/master/example/sample.png)


# Requisites
* Python >= 2.7 with the module [Biopython](http://biopython.org/wiki/Download)
* [blastp binary](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* R (tested in 3.3.1), with modules:
	* ggplot2: ```install.packages("ggplot2")```
	* genoPlotR: ```install.packages("genoPlotR", repos="http://R-Forge.R-project.org")```
* optional: [Mauve](http://darlinglab.org/mauve/download.html) for gbk-gbk whole comparison.

# Usage

multiGenomicContext have a minimal use:
	
	python multiGenomicContext.py -f protein.fasta -l mygbk1.gbk,mygbk2.gbk
	
and a complete use:

	python multiGenomicContext.py -f protein.fasta -l mygbk1.gbk,mygbk2.gbk -u 4 -d 4 -e 1e-5 -i 85 -a 75 
	
where the options are:

* -f: The protein sequence in fasta format (also can be a multifasta of proteins).
* -l: The list of gbk (if you have one, also put that name in a file)
* -u: Number of genes to put the genomic context in upstream search (default: 4)
* -d Number of genes to put the genomic context in downstream search (default: 4)
* -e E-value for blastp search (default: 1e-5)
* -i Identity % of the alignment on blastp results to consider the gene exists on the genome (default: 85)
* -a Alignment length (%) between gene and the match for blastp search to consider the gene "exists" on the genome (default 75)

There are more options available for a more customizable way.

* -g: Plot gbk-gbk aligment (use -m for exact binary progressiveMauve path), this option will ignore "-f" and only should be used with "-l"
* -b: blastp binary path
* -m: progressive mauve align. This option is to plot the entire gbk regions (not genes), is used to see synteny between genomes. (use -g)
* -c: for clean files off, by default the script will remove all tmp files on the way, if you used "-c" you will conserve all files including the .R for the plot.

# Notes
* multiGenomicContext search genes on the gbk because the cds are ordered, but this true only for one chromosome assembly, for gbk files where two or more contigs exists, it's show genomic context for the same contig of the gene.
* We strongly recommend use one software to annotate all gbk's to conserve the same names outcome.
* -g option use prograssiveMauve that only accept "gbk" extension, not "gbff".
* Sometimes the names are too large for the area plot, maybe you should fix it with illustrator :).

# External useful tools
check for these tools to extract some useful information from your data:


* [fetchMyLineage](https://github.com/Sanrrone/fetchMyLineage): Return the complete lineage of your organism just providing the genus and species names.

* [extractSeq](https://github.com/Sanrrone/extractSeq): Extract and size defined sequence from and specific contig, from and specific genome.

* [plotMyGBK](https://github.com/Sanrrone/plotMyGBK): Plot your GBK in a circular graph with COG categories.

* [pasteTaxID](https://github.com/Sanrrone/pasteTaxID): fetch the taxonomic IDs to your fastas.

* [GGisy](https://github.com/Sanrrone/GGisy): Plot synteny of two sequence (you can use two genomes), and see the identity of the matched regions.

* [getS2](https://github.com/Sanrrone/getS2): obtain the order parameter to each residue of your simulation.
