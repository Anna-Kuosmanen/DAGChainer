#include <lemon/list_graph.h>
#include <lemon/lgf_writer.h>
#include <lemon/dfs.h>
#include "../util/utils.h"
#include <string.h>
#include <time.h>

using namespace lemon;
using namespace std;

#define MAX_WEIGHT 1000

void createKPathGraph(string output_filename, int k, int n, int m);
void createKPathGraphAlt(string output_filename, int k, int n, int m);
void createRandomDAG(string output_filename, int n, float prob);

int main(int argc, char* argv[])
{
  if (argc != 6 && argc != 5 && argc != 4){
    cerr << "Usage: " << argv[0] << " OUTPUT_FILENAME NUM_PATHS PATH_LENGTH NUM_EXTRA_ARCS (--alt)"
     << endl << "or" << endl << "OUTPUT_FILENAME NUM_NODES ARC_PROBABILITY";
    return 1;
  }

  if(argc == 4){
    int num_nodes;
    float arc_probability;
  try {
    num_nodes = stoi(argv[2]);
    arc_probability = stof(argv[3]);
    createRandomDAG(argv[1], num_nodes, arc_probability);
    return 0;

  } catch (const std::exception& e) {
    cerr << "ERROR: error reading numerical parameters\n";
    return 1;
    }
  }
  int num_paths, path_length;
  int extra_arcs;
  try {
    num_paths = stoi(argv[2]);
    path_length = stoi(argv[3]);
    extra_arcs = stoi(argv[4]);
  } catch (const std::exception& e) {
    cerr << "ERROR: error reading numerical parameters\n";
    return 1;
  }
  if(argc == 5){
    createKPathGraph(argv[1], num_paths, path_length, extra_arcs);
  }
  else if(argc == 6){
    createKPathGraphAlt(argv[1], num_paths, path_length, extra_arcs);
  }

  return 0;
}

void createRandomDAG(string output_filename, int n, float prob)
{
   ListDigraph g;
   ListDigraph::NodeMap<int> node_labels(g);
   
   ListDigraph::Node* nodes = new ListDigraph::Node[n];


   for(int i = 0; i < n; ++i){
      ListDigraph::Node temp = g.addNode();
      node_labels[temp] = i;
      nodes[i] = temp;
   }
   srand(time(NULL));

   for(int i = 0; i < n; i++){
     for(int j = i+1; j < n; j++){
       if(rand()%100 < prob*100){
          g.addArc(nodes[i], nodes[j]);
       }
     }
   }

   ListDigraph::ArcMap<int> arc_labels(g);
   ListDigraph::ArcMap<int> arc_weights(g);

   int index = 0;
   for(ListDigraph::ArcIt ai(g); ai != INVALID; ++ai){
     arc_labels[ai] = index;
     arc_weights[ai] = rand() % MAX_WEIGHT;
     index++;
  }

  DigraphWriter<ListDigraph>(g, output_filename)
    .nodeMap("label", node_labels)
    .arcMap("label", arc_labels)
    .arcMap("weight", arc_weights)
    .run();

}

/**
Normal k-path graphs where extra arcs are always "forward", which guarantees acyclicity
**/
void createKPathGraph(string output_filename, int k, int n, int m)
{
	srand(time(NULL));

  ListDigraph g;

	ListDigraph::Node* nodes[k];
	for (int i = 0; i < k; ++i){
		nodes[i] = (ListDigraph::Node*) calloc(n, sizeof(ListDigraph::Node));
	}
	for (int i = 0; i < k; ++i){
		for (int j = 0; j < n; ++j){
			nodes[i][j] = g.addNode();
			if(j != 0){
        ListDigraph::Arc temp_arc = g.addArc(nodes[i][j-1], nodes[i][j]);
      }
		}
	}
  cout << "m:" << to_string(m) << endl;
  for (int i = 0; i < m; ++i){
 		int k1 = rand()%k;
 		int k2 = rand()%k;
 		int n1 = rand()%(n-1);
 		int n2 = (rand()%(n-n1-1))+n1+1;
 
 		if(findArc(g, nodes[k1][n1], nodes[k2][n2]) == INVALID){
 			ListDigraph::Arc temp_arc = g.addArc(nodes[k1][n1], nodes[k2][n2]);
 		}
 	}
   
  //Shuffle nodes

  int num_nodes = countNodes(g);
  ListDigraph::Node* node_arr = new ListDigraph::Node[num_nodes];
  int node_index = 0;
  for(ListDigraph::NodeIt n(g); n != INVALID; ++n){
    node_arr[node_index] = n;
    node_index++;
  }
  random_shuffle(node_arr, node_arr + num_nodes);


  //Shuffle arcs
  int num_arcs = countArcs(g);
  ListDigraph::Arc* arcs = new ListDigraph::Arc[num_arcs];
  int arc_index = 0;
  for(ListDigraph::ArcIt a(g); a != INVALID; ++a) { 
    arcs[arc_index] = a;
    arc_index++;
  }
  random_shuffle(arcs, arcs + num_arcs);

  //construct shuffled graph
  ListDigraph shuffled;
  ListDigraph::NodeMap<int> shuffled_node_labels(shuffled);
  ListDigraph::ArcMap<int> shuffled_arc_labels(shuffled);
  ListDigraph::ArcMap<int> arc_weights(shuffled);

  ListDigraph::NodeMap<ListDigraph::Node> mapping(g);
  int index = 0;
  for(int i = 0; i< num_nodes; i++) {
    mapping[node_arr[i]] = shuffled.addNode();
    shuffled_node_labels[mapping[node_arr[i]]] = index;
    index++;
  }
  index = 0;
  for(int i = 0; i< num_arcs; i++) {
    ListDigraph::Node n1 = mapping[g.source(arcs[i])];
    ListDigraph::Node n2 = mapping[g.target(arcs[i])];
    ListDigraph::Arc temp = shuffled.addArc(n1, n2);
    shuffled_arc_labels[temp] = index;
    arc_weights[temp] = rand()%MAX_WEIGHT;
    index++;
  }

  DigraphWriter<ListDigraph>(shuffled, output_filename)
    .nodeMap("label", shuffled_node_labels)
    .arcMap("label", shuffled_arc_labels)
    .arcMap("weight", arc_weights)
    .run();
}


/**
Alternative model of k-path-graphs. In this model arcs are added between any two random nodes, if the new arc would not add a cycle to the graph. 
**/
void createKPathGraphAlt(string output_filename, int k, int n, int m)
{
	srand(time(NULL));

  ListDigraph g;

	ListDigraph::Node* nodes[k];
	for (int i = 0; i < k; ++i){
		nodes[i] = (ListDigraph::Node*) calloc(n, sizeof(ListDigraph::Node));
	}
	for (int i = 0; i < k; ++i){
		for (int j = 0; j < n; ++j){
			nodes[i][j] = g.addNode();
			if(j != 0){
        ListDigraph::Arc temp_arc = g.addArc(nodes[i][j-1], nodes[i][j]);
      }
		}
	}
  cout << "m:" << to_string(m) << endl;
  Dfs<ListDigraph> dfs(g);

	for (int i = 0; i < m; ++i){
		int k1 = rand()%k;
		int k2 = rand()%k;
		int n1 = rand()%n;
		int n2 = rand()%n;

    ListDigraph::Node from_node = nodes[k1][n1];
    ListDigraph::Node to_node = nodes[k2][n2];

		if(findArc(g, from_node, to_node) == INVALID && !dfs.run(to_node, from_node)){
			ListDigraph::Arc temp_arc = g.addArc(nodes[k1][n1], nodes[k2][n2]);
		}
	}
   
  //Shuffle nodes

  int num_nodes = countNodes(g);
  ListDigraph::Node* node_arr = new ListDigraph::Node[num_nodes];
  int node_index = 0;
  for(ListDigraph::NodeIt n(g); n != INVALID; ++n){
    node_arr[node_index] = n;
    node_index++;
  }
  random_shuffle(node_arr, node_arr + num_nodes);


  //Shuffle arcs
  int num_arcs = countArcs(g);
  ListDigraph::Arc* arcs = new ListDigraph::Arc[num_arcs];
  int arc_index = 0;
  for(ListDigraph::ArcIt a(g); a != INVALID; ++a) { 
    arcs[arc_index] = a;
    arc_index++;
  }
  random_shuffle(arcs, arcs + num_arcs);

  //construct shuffled graph
  ListDigraph shuffled;
  ListDigraph::NodeMap<int> shuffled_node_labels(shuffled);
  ListDigraph::ArcMap<int> shuffled_arc_labels(shuffled);
  ListDigraph::ArcMap<int> arc_weights(shuffled);

  ListDigraph::NodeMap<ListDigraph::Node> mapping(g);
  int index = 0;
  for(int i = 0; i< num_nodes; i++) {
    mapping[node_arr[i]] = shuffled.addNode();
    shuffled_node_labels[mapping[node_arr[i]]] = index;
    index++;
  }
  index = 0;
  for(int i = 0; i< num_arcs; i++) {
    ListDigraph::Node n1 = mapping[g.source(arcs[i])];
    ListDigraph::Node n2 = mapping[g.target(arcs[i])];
    ListDigraph::Arc temp = shuffled.addArc(n1, n2);
    shuffled_arc_labels[temp] = index;
    arc_weights[temp] = rand()%MAX_WEIGHT;
    index++;
  }

  DigraphWriter<ListDigraph>(shuffled, output_filename)
    .nodeMap("label", shuffled_node_labels)
    .arcMap("label", shuffled_arc_labels)
    .arcMap("weight", arc_weights)
    .run();
}