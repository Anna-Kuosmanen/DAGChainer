## Input a gene annotation (gtf) file, output a list of the size of the sequence graph made from each gene
## Prints format "gene name\tno of vertices"

import sys, os, subprocess

if len(sys.argv) < 2:
	print "Give an annotation file."
	sys.exit()

gtf_in = open(sys.argv[1])

genes = list()

# Collect all gene IDs
while gtf_in:
	line = gtf_in.readline().strip("\n")
	if line == "":
		break
	line_parts = line.split('\t')
	
	description_parts = line_parts[8].split(';')
	description_dict = dict()

	for index, part in enumerate(description_parts):
		if index == len(description_parts)-1:
			continue
		key_value_pair = part.strip().split(" ")
		description_dict[key_value_pair[0].strip()] = key_value_pair[-1].strip()

	gene_id = description_dict["gene_id"].strip('"')
	if gene_id not in genes:
		genes.append(gene_id)

gtf_in.close()


for gene in genes:
	cmd = "grep '%s' %s > temp_gene.tmp" % (gene, sys.argv[1])
	p = subprocess.Popen(cmd, shell = True, stdout = None, stderr = None)
	p.wait()

	exons = list()

	gtf_in = open("temp_gene.tmp", "r")

	## Save as tuples of (chrom, start, end) to prune duplicates
	while True:
		line = gtf_in.readline()
		if line == "":
			break
		line_parts = line.split("\t")

		if line_parts[2] == "exon":
			exon = (line_parts[0], line_parts[3], line_parts[4])
			if exon not in exons:
				exons.append(exon)

	gtf_in.close()

	length = 0;

	for exon in exons:
		length += (int(exon[2])-int(exon[1]) +1)

	print "%s\t%s\t%s" % (gene, exons[0][0], length)





