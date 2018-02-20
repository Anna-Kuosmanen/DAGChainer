#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include <string.h>
#include "../util/utils.h"

using namespace lemon;
using namespace std;

int split(string filename, string output_folder);

int main(int argc, char* argv[])
{
  if (argc != 3){
    cerr << "Usage: " << argv[0] << " GRAPH_FILENAME OUTPUT_FOLDER/" << endl;
    return 1;
  }
  if(!file_exists(argv[1])){
    cerr << "ERROR: input file not found\n";
    return 1;
  }
  if(!directory_exists(argv[2])){
    cerr << "ERROR: output directory not found\n";
    return 1;
  }
  return split(argv[1], argv[2]);
}

int split(string filename, string output_folder)
{

  ListDigraph g;

  ListDigraph::NodeMap<int> node_labels(g);
  ListDigraph::ArcMap<int> arc_labels(g);
  ListDigraph::ArcMap<int> arc_weights(g);

  digraphReader(g, filename)
    .nodeMap("label", node_labels)
    .arcMap("label", arc_labels)
    .arcMap("weight", arc_weights)
    .run();
  
  Undirector<ListDigraph> undirected(g);
	ListDigraph::NodeMap<int> components(g);
	stronglyConnectedComponents(undirected, components);

	int num_subgraphs = 0;
	for(ListDigraph::NodeIt n(g); n != INVALID; ++n){
		if(components[n] > num_subgraphs) num_subgraphs = components[n];
	}
	num_subgraphs++;
	ListDigraph::NodeMap<ListDigraph::Node> map(g);

	for(int i = 0; i < num_subgraphs; i++){

		ListDigraph temp;
    ListDigraph::NodeMap<int> temp_node_labels(g);
    ListDigraph::ArcMap<int> temp_arc_labels(g);
    ListDigraph::ArcMap<int> temp_arc_weights(g);

		for(ListDigraph::NodeIt n(g); n != INVALID; ++n){
			if(components[n] == i){
				map[n] = temp.addNode();
        temp_node_labels[map[n]] = node_labels[n];
			}
		}
		for(ListDigraph::NodeIt n(g); n != INVALID; ++n){
			if(components[n] == i){
				for(ListDigraph::OutArcIt o(g, n); o != INVALID; ++o){
					ListDigraph::Arc temp_arc = temp.addArc(map[g.source(o)], map[g.target(o)]);
          temp_arc_labels[temp_arc] = arc_labels[o];
          temp_arc_weights[temp_arc] = arc_weights[o];
				}
			}
		}
    string output_filename = output_folder + filename + "_split_" + to_string(i);
      DigraphWriter<ListDigraph>(temp, output_filename)
        .nodeMap("label", temp_node_labels)
        .arcMap("label", temp_arc_labels)
        .arcMap("weight", temp_arc_weights)
        .run();
	}
  return 0;
}