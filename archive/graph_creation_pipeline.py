import os, subprocess

genome_file = "/cs/work/home/aekuosma/genomes/Homo_sapiens_hg38/hg38.fa"
annotation_file = "/cs/work/home/aekuosma/genomes/Homo_sapiens_hg38/annotation/Homo_sapiens_genes_hg38_Ensembl_ID.fixed.gtf"

genes = list()

inhandle = open("added_gene.txt")

while inhandle:
	line = inhandle.readline().strip("\n")
	if line == "":
		break
	genes.append(line)


for gene in genes:
	cmd = "mkdir %s" % (gene)
	p = subprocess.Popen(cmd, shell = True, stdout = None, stderr = None)
	p.wait()

	cmd = "grep '%s' %s > %s/%s.tmp" % (gene, annotation_file, gene, gene)
	p = subprocess.Popen(cmd, shell = True, stdout = None, stderr = None)
	p.wait()

	cmd = "python /cs/work/home/aekuosma/software/colinearchaining/create_input_files.py %s/%s.tmp %s %s/%s" % (gene, gene, genome_file, gene, gene)
	p = subprocess.Popen(cmd, shell = True, stdout = None, stderr = None)
	p.wait()

