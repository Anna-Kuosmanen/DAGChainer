## Helper script to add global sink to every graph, to prevent them from being disjoint.

import os, sys

if len(sys.argv) < 2:
	print "Give the prefix of the graph and nodes files."
	sys.exit(1)

prefix = sys.argv[1]

graph_in = open("%s.graph" % (prefix))
nodes_in = open("%s.nodes" % (prefix))

graph_out = open("%s.graph.fixed" % (prefix), "w")
nodes_out = open("%s.nodes.fixed" % (prefix), "w")

# For graph file, every row that has no outgoing edges add an edge to the new vertex
line = graph_in.readline().strip("\n")

maxvertex = int(line)

graph_out.write("%s\n" % (str(maxvertex+1)))

for i in range(maxvertex):
	line = graph_in.readline().strip("\n")
	if line == "":
		graph_out.write("%s\n" % (str(maxvertex)))
	else:
		graph_out.write("%s\n" % (line))

graph_out.write("\n")

graph_in.close()
graph_out.close()

# For nodes file add a dummy vertex
while nodes_in:
	line = nodes_in.readline().strip("\n")
	if line == "":
		break
	nodes_out.write("%s\n" % (line))

nodes_out.write("%s\nN\n" % (str(maxvertex)))

nodes_in.close()
nodes_out.close()
