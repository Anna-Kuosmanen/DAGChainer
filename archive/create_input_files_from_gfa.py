## TODO Reverse complementing, how to read the links? Is first +/- the inneighbors and second +/- the outneighbor?

import os, sys

if len(sys.argv) < 3:
	print "Give GFA file and output file prefix"
	sys.exit(1)


inhandle = open(sys.argv[1])
graph_out = open("%s.graph" % (sys.argv[2]), "w")
nodes_out = open("%s.nodes" % (sys.argv[2]), "w")

vertex_lines = list()
arc_lines = list()


while inhandle:
	line = inhandle.readline().strip("\n")
	if line == "":
		break
	parts = line.split("\t")
	# Vertex
	if parts[0] == "S":
		vertex_lines.append(line)
	# Edge
	if parts[0] == "L":
		arc_lines.append(line)


for vertex in vertex_lines:
	parts = vertex.split("\t")
	index = parts[1]
	seq = parts[2]
	nodes_out.write("%s\n" % (index))
	nodes_out.write("%s\n" % (seq))

# List of lists, one list per node
arc_list = list()

for i in range(len(vertex_lines)):
	temp_list = list()
	arc_list.append(temp_list)


for arc in arc_lines:
	print arc
	parts = arc.split("\t")
	start = int(parts[1])
	end = int(parts[3])
	arc_list[start].append(end)



# First write the number of nodes
graph_out.write("%s\n" % (len(vertex_lines)))

for item in arc_list:
	int_list = sorted(item)
	graph_out.write("%s\n" % (" ".join([str(x) for x in int_list])))



inhandle.close()
graph_out.close()
nodes_out.close()

