## Sample genes from the vertex count file (gene\tchrom\tvertexcount)
## Give as parameter the vertex count file, number of genes to sample,
## and the range of the graph size

import sys
from random import uniform

if len(sys.argv) < 5:
	print "Please give as parameters the vertex count file, number of genes to sample and the range of the graph size (min, max)."
	sys.exit()

genes_in = open(sys.argv[1])
number_to_sample = int(sys.argv[2])
sizemin = int(sys.argv[3])
sizemax = int(sys.argv[4])

candidate_genes = list()

while genes_in:

	line = genes_in.readline().strip("\n")
	if line == "":
		break
	parts = line.split("\t")

	# Drop the genes for alternative loci (they have longer names)
	if len(parts[1]) > 5:
		continue

	if (int(parts[2]) >= sizemin and int(parts[2]) <= sizemax):
		candidate_genes.append(parts[0])

genes_in.close()
chosen_genes = list()

while len(chosen_genes) < number_to_sample and len(candidate_genes) > 0:
	next_index = int(uniform(0,len(candidate_genes)))
	next_chosen = candidate_genes[next_index]
	chosen_genes.append(next_chosen)
	del candidate_genes[next_index]


for gene in chosen_genes:
	print gene

