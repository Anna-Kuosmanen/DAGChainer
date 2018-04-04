import os, subprocess

genes = list()

inhandle = open("added_gene.txt")

while inhandle:
	line = inhandle.readline().strip("\n")
	if line == "":
		break
	genes.append(line)


for gene in genes:
	print gene
	cmd = "/cs/work/home/aekuosma/software/colinearchaining/extract_patterns %s/%s.graph.fixed %s/%s.nodes.fixed %s 100 1000 > %s/patterns.txt" % (gene, gene, gene, gene, gene, gene)
	p = subprocess.Popen(cmd, shell = True, stdout = None, stderr = None)
	p.wait()

