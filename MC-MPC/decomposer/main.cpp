#include <lemon/list_graph.h>
#include "decomposition.h"
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include "MPC.h"
#include "../util/utils.h"
#include <string.h>

using namespace lemon;
using namespace std;

int decompose(string filename, string output_folder);

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

  return decompose(argv[1], argv[2]);
}

int decompose(string filename, string output_folder)
{
  ListDigraph graph;

  ListDigraph::NodeMap<int> node_labels(graph);
  ListDigraph::ArcMap<int> arc_labels(graph);
  ListDigraph::ArcMap<int> arc_weights(graph);

  digraphReader(graph, filename)
    .nodeMap("label", node_labels)
    .arcMap("label", arc_labels)
    .arcMap("weight", arc_weights)
    .run();
  
  ListDigraph::Node s = add_source(graph);
  ListDigraph::Node t = add_sink(graph);

  ListDigraph::ArcMap<int> minflow(graph);
  find_minflow(graph, minflow, s, t);

  ListDigraph::ArcMap<int> decomposition(graph);
  decompose_graph(graph, minflow, s, t, decomposition);

  //sink and source nodes added to the graph should not be included in the output
  graph.erase(s);
  graph.erase(t);

  //find the highest # of a decomposition
  int num_decompositions = 0;
  for(ListDigraph::ArcIt ai(graph); ai != INVALID; ++ai){
    if(num_decompositions < decomposition[ai]) {
      num_decompositions = decomposition[ai];
    }
  } 

  //this index is just for output file names. there is not necessarily a decomposition for every normal index number
  int decomposition_index = 0;

  for(int i = 0; i <= num_decompositions; i++)
  {
    ListDigraph temp;
    ListDigraph::NodeMap<int> temp_node_labels(temp);
    ListDigraph::ArcMap<int> temp_arc_labels(temp);
    ListDigraph::ArcMap<int> temp_arc_weights(temp);

    //mapping from original graph to a decomposed part
    ListDigraph::Node null_node = graph.addNode();
    ListDigraph::NodeMap<ListDigraph::Node> mapping(graph, null_node);

    bool decomposition_found = false;

    for(ListDigraph::ArcIt a(graph); a != INVALID; ++a){

      if(decomposition[a] == i){
        decomposition_found = true;
        ListDigraph::Node source = graph.source(a);
        ListDigraph::Node target = graph.target(a);

        if(mapping[source] == null_node){
          mapping[source] = temp.addNode();
          temp_node_labels[mapping[source]] = node_labels[source];
        }
        if(mapping[target] == null_node){
          mapping[target] = temp.addNode();
          temp_node_labels[mapping[target]] = node_labels[target];
        }

        ListDigraph::Arc temp_arc = temp.addArc(mapping[source], mapping[target]);
        temp_arc_labels[temp_arc] = arc_labels[a];
        temp_arc_weights[temp_arc] = arc_weights[a];
      }
    }
    if(decomposition_found == true){
      string output_filename = output_folder + "decomp_" + to_string(decomposition_index);
      DigraphWriter<ListDigraph>(temp, output_filename)
        .nodeMap("label", temp_node_labels)
        .arcMap("label", temp_arc_labels)
        .arcMap("weight", temp_arc_weights)
        .run();
      decomposition_index++;
    }
  }
  return 0;
}