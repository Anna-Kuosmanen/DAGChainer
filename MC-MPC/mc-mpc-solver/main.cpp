#include <lemon/list_graph.h>
#include <lemon/concepts/digraph.h>
#include <lemon/cost_scaling.h>
#include <lemon/network_simplex.h>
#include <string.h>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <lemon/connectivity.h>
#include <lemon/circulation.h>
#include "../util/utils.h"
#include "../decomposer/decomposition.h"
#include "../decomposer/MPC.h"
#include <time.h>

using namespace lemon;
using namespace lemon::concepts;
using namespace std;

string solve_minflow(string filename);
string solve_with_decomposition(string filename);
string solve_with_decomposition2(string filename);

int main(int argc, char* argv[])
{
  if (argc != 3 && argc != 4){
    cerr << "Usage: " << argv[0] << " GRAPH_FILENAME OUTPUT_FILENAME (--decomp)" << endl;
    return 1;
  }

  if(!file_exists(argv[1])){
    cerr << "ERROR: input file not found\n";
    return 1;
  }
  string result;

  if(argc == 4){
    if(strcmp(argv[3],"--decomp") == 0){
      result = solve_with_decomposition(argv[1]);
    } else {
      cerr << "last argument not recognized\n";
      return 1;
    }
  }else{
    result = solve_minflow(argv[1]);
  }
  ofstream output_file;
  output_file.open(argv[2]);
  output_file << result;
  output_file.close();
  return 0;
}

void transform_to_minflow_graph(ListDigraph& g, ListDigraph::ArcMap<int>& demands)
{
	//splits each node into connected nodes, where the first has all the incoming and the second all the outgoing edges

	for (ListDigraph::NodeIt n(g); n != INVALID; ++n)
	{
		ListDigraph::Node new_node = g.split(n, false);
		ListDigraph::Arc new_arc = g.addArc(n, new_node);
		demands[new_arc] = 1;
	}
}

string solve_with_decomposition(string filename)
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
  clock_t start = clock();
  cout << "start ";
  ListDigraph::ArcMap<int> minflow(graph);
  find_minflow(graph, minflow, s, t);
  clock_t after_minflow = clock();
  cout << "minflow ";

  ListDigraph::ArcMap<int> decomposition(graph);
  decompose_graph(graph, minflow, s, t, decomposition);
  clock_t after_decomp = clock();
  cout << "decomp " << endl;
  cout << "minflow: " << to_string(((float) (after_minflow-start)/CLOCKS_PER_SEC)) << " decomp: " << to_string(((float) (after_decomp-after_minflow)/CLOCKS_PER_SEC)) << endl;


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

  cout << "num decompositions: " << to_string(num_decompositions) << endl;
  
  clock_t copy_time = 0;
  clock_t solve_time = 0;

  //this index is just for output file names. there is not necessarily a decomposition for every normal index number
  int decomposition_index = 0;

  string result = "";

  for(int i = 0; i <= num_decompositions; i++)
  {
    clock_t start2 = clock();
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
    copy_time += (clock()-start2);
    clock_t start3 = clock();
    //solve in decomposed parts
    if(decomposition_found){
      ListDigraph::ArcMap<bool> original_arc(temp, false);
      int arc_weight_sum = 0;
      for(ListDigraph::ArcIt ai(temp); ai != INVALID; ++ai){
        original_arc[ai] = true;
        arc_weight_sum += temp_arc_weights[ai];
      }
      ListDigraph::ArcMap<int> demands(temp);
	    transform_to_minflow_graph(temp, demands);

	    ListDigraph::Node s = add_source(temp);
      ListDigraph::Node t = add_sink(temp);

      ListDigraph::Arc sink_to_source = temp.addArc(t, s);
      temp_arc_weights[sink_to_source] = arc_weight_sum;
  
	    //initialize cost scaling algorithm and give it the demands and the weights of the edges
      /**
	    CostScaling<ListDigraph> costScaling(temp);
	    costScaling.lowerMap(demands);
	    costScaling.costMap(temp_arc_weights);

	    costScaling.run();

      ListDigraph::ArcMap<int> minflow(temp);
      costScaling.flowMap(minflow);
      **/

      NetworkSimplex<ListDigraph> solver(temp);
      solver.lowerMap(demands);
      solver.costMap(arc_weights);
      solver.run();
      ListDigraph::ArcMap<int> minflow(temp);
      solver.flowMap(minflow);
      
      for(ListDigraph::ArcIt ai(temp); ai != INVALID; ++ai){
      //we have added many extra arcs which are not needed in the output
        if(original_arc[ai]){
          result += (to_string(temp_arc_labels[ai]) + " " + to_string(minflow[ai]) + " " +   to_string(temp_arc_weights[ai]) + "\n");
        }
      }
    }
    solve_time += (clock()-start3);
  }
  cout << "copy: " << to_string(((float) copy_time)/CLOCKS_PER_SEC) << " solve: " << to_string(((float) solve_time)/CLOCKS_PER_SEC) << endl;
  return result; 

}

string solve_minflow(string filename)
{

  clock_t start = clock();

  ListDigraph g;
  ListDigraph::NodeMap<int> node_labels(g);
  ListDigraph::ArcMap<int> arc_labels(g);
  ListDigraph::ArcMap<int> arc_weights(g);
  digraphReader(g, filename)
    .nodeMap("label", node_labels)
    .arcMap("label", arc_labels)
    .arcMap("weight", arc_weights)
    .run();
  
  ListDigraph::ArcMap<bool> original_arc(g, false);
  int arc_weight_sum = 0;
  for(ListDigraph::ArcIt ai(g); ai != INVALID; ++ai){
    original_arc[ai] = true;
    arc_weight_sum += arc_weights[ai];
  }

  clock_t start2 = clock();
	// transform the graph so that every node is split in two, and demand for the edge between the duplicates is 1

	ListDigraph::ArcMap<int> demands(g);
	transform_to_minflow_graph(g, demands);

	ListDigraph::Node s = add_source(g);
  ListDigraph::Node t = add_sink(g);

  ListDigraph::Arc sink_to_source = g.addArc(t, s);
  arc_weights[sink_to_source] = arc_weight_sum;
  
	//initialize cost scaling algorithm and give it the demands and the weights of the edges
  /**
	CostScaling<ListDigraph> costScaling(g);
	costScaling.lowerMap(demands);
	costScaling.costMap(arc_weights);

	costScaling.run();

  ListDigraph::ArcMap<int> minflow(g);
  costScaling.flowMap(minflow);
  **/
  
  NetworkSimplex<ListDigraph> solver(g);
  solver.lowerMap(demands);
  solver.costMap(arc_weights);
  solver.run();
  ListDigraph::ArcMap<int> minflow(g);
  solver.flowMap(minflow);
 
  string result = "";
  for(ListDigraph::ArcIt ai(g); ai != INVALID; ++ai){
    //we have added many extra arcs which are not needed in the output
    if(original_arc[ai]){
      result += (to_string(arc_labels[ai]) + " " + to_string(minflow[ai]) + " " + to_string(arc_weights[ai]) + "\n");
    }
  }
  clock_t start3 = clock();
  cout << "normal: reading: " << to_string(((float) start2-start) / CLOCKS_PER_SEC) << " solving: " << to_string(((float) start3-start2) / CLOCKS_PER_SEC) << endl;

  return result;
}