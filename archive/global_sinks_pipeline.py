import os, subprocess

genes = list()

inhandle = open("chosen_genes.txt")

while inhandle:
	line = inhandle.readline().strip("\n")
	if line == "":
		break
	genes.append(line)


for gene in genes:
	cmd = "python /cs/work/home/aekuosma/software/colinearchaining/add_global_sinks.py %s/%s" % (gene, gene)
	p = subprocess.Popen(cmd, shell = True, stdout = None, stderr = None)
	p.wait()

