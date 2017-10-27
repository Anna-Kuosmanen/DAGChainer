## Input: annotation (for one gene) and FASTA and output file prefix
## Output: graph file and nodes file for creating the sequence graph

import os, sys
from Bio import SeqIO

if len(sys.argv) < 4:
	print "Give annotation file and FASTA file and output file prefix"
	sys.exit(1)

output_prefix = sys.argv[3]

# List of lists, each corresponding to the exons in one transcript
transcripts = list()

exons = list()

gtf_in = open(sys.argv[1])

while True:
	line = gtf_in.readline()
	if line == "":
		break
	line_parts = line.split("\t")
	if line_parts[2] == "exon":
 		exons.append(line)

gtf_in.close()

while(len(exons) > 0):

	# The first exon decides the transcript name to look for
	exon = exons[0]

	# List of exons for this transcript
	tuples = list()

	exon_parts = exon.split('\t')

	chrom = exon_parts[0]

	description_parts = exon_parts[8].split(';')
	description_dict = dict()
	for index, part in enumerate(description_parts):
		if index == len(description_parts)-1:
			continue
		key_value_pair = part.strip().split(" ")
		description_dict[key_value_pair[0].strip()] = key_value_pair[-1].strip()

	tuples.append((exon_parts[0],exon_parts[3],exon_parts[4]))

	del exons[0]

	to_delete = list()

	# Then find all the other exons with this transcript name
	for index in range(len(exons)):
		exon = exons[index]

		exon_parts = exon.split('\t')
		description_parts = exon_parts[8].split(';')

		description_dict_2 = dict()
		for i, part in enumerate(description_parts):
			if i == len(description_parts)-1:
				continue
			key_value_pair = part.strip().split(" ")
			description_dict_2[key_value_pair[0].strip()] = key_value_pair[-1].strip()


		if(description_dict["transcript_id"] == description_dict_2["transcript_id"]):
			to_delete.append(index)
			tuples.append((exon_parts[0],exon_parts[3],exon_parts[4]))

	for index in range(len(to_delete)):
		del exons[int(to_delete[index]-index)]

	transcripts.append(tuples)

# Combine exon lists and sort by coordinates
combo_list = list()

for transcript in transcripts:
	for exon in transcript:
		if exon not in combo_list:
			combo_list.append(exon)

combo_list.sort(key=lambda tup: tup[1])

# Create a dict that maps the exon tuple (chr, start, end) to the node number
tuples_to_nodes = dict()
# Create another dict that maps the node numbers to exon tuples
nodes_to_tuples = dict()

for index, item in enumerate(combo_list):
	tuples_to_nodes[item] = index
	nodes_to_tuples[index] = item


# Go over the transcript lists to find arcs and write to graph file
graph_out = open("%s.graph" % (output_prefix), "w")

# List of lists, one list per node
arc_list = list()

for i in range(len(combo_list)):
	temp_list = list()
	arc_list.append(temp_list)


for transcript in transcripts:

	for index, exon in enumerate(transcript):
		if index == len(transcript)-1:
			continue
		first_node = tuples_to_nodes[exon]
		second_node = tuples_to_nodes[transcript[index+1]]
		if second_node not in arc_list[first_node]:
			arc_list[first_node].append(second_node)


# First write the number of nodes
graph_out.write("%s\n" % (len(combo_list)))

for item in arc_list:
	int_list = sorted([int(x) for x in item])
	graph_out.write("%s\n" % (" ".join([str(x) for x in int_list])))




# Extract the sequences and write to nodes file
for record in SeqIO.parse(sys.argv[2], "fasta"):
	if record.id == chrom:
		sequence = str(record.seq)

nodes_out = open("%s.nodes" % (output_prefix), "w")

for index, exon in enumerate(combo_list):
	seq = sequence[int(exon[1]):int(exon[2])+1].upper()
	nodes_out.write("%s\n" % (index))
	nodes_out.write("%s\n" % (seq))




graph_out.close()
nodes_out.close()
