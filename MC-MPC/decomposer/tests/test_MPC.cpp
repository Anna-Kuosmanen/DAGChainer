#include <lemon/list_graph.h>
#include <lemon/cost_scaling.h>
#include "../MPC.h"
#include <stdlib.h>
#include <time.h>
#include "../catch.hpp"
#include "../../util/utils.h"

using namespace std;
using namespace lemon;


TEST_CASE("simple test case succeeds"){
  ListDigraph g;
  ListDigraph::Node s = g.addNode();
  ListDigraph::Node t = g.addNode();
  ListDigraph::Node a = g.addNode();
  ListDigraph::Node b = g.addNode();
  
  ListDigraph::Arc sa = g.addArc(s, a);
  ListDigraph::Arc sb = g.addArc(s, b);
  ListDigraph::Arc ab = g.addArc(a, b);
  ListDigraph::Arc bt = g.addArc(b, t);
  ListDigraph::Arc at = g.addArc(a, t);
  
  
  ListDigraph::ArcMap<int> flow(g);
  
  find_minflow(g, flow, s, t);
  
  REQUIRE(flow[sa] == 1);
  REQUIRE(flow[ab] == 1);
  REQUIRE(flow[bt] == 1);
  REQUIRE(flow[sb] == 0);
  REQUIRE(flow[at] == 0);
        
}
TEST_CASE("Another simple test case succeeds"){
  ListDigraph g;
  ListDigraph::Node s = g.addNode();
  ListDigraph::Node t = g.addNode();
  ListDigraph::Node a = g.addNode();
  ListDigraph::Node b = g.addNode();
  ListDigraph::Node c = g.addNode();
  ListDigraph::Node d = g.addNode();
  ListDigraph::Node e = g.addNode();
  
  ListDigraph::Arc sa = g.addArc(s, a);
  ListDigraph::Arc ac = g.addArc(a, c);
  ListDigraph::Arc cb = g.addArc(c, b);
  ListDigraph::Arc bt = g.addArc(b, t);
  ListDigraph::Arc sd = g.addArc(s, d);
  ListDigraph::Arc dc = g.addArc(d, c);
  ListDigraph::Arc ce = g.addArc(c, e);
  ListDigraph::Arc et = g.addArc(e, t);
  ListDigraph::Arc ab = g.addArc(a, b);
  ListDigraph::Arc de = g.addArc(d, e);
  

  ListDigraph::ArcMap<int> flow(g);
  
  find_minflow(g, flow, s, t);
  
  REQUIRE(flow[bt]+flow[et] == 2);
        
}

//Simple function for generating acyclic graphs
void createRandomGraph(ListDigraph& g, int num_nodes, float edge_prob){

	srand(time(NULL));
	ListDigraph::NodeMap<int> labels(g);

	for(int i=0; i<num_nodes; i++){
		ListDigraph::Node new_node = g.addNode();
		labels[new_node] = i;
	}
	for(ListDigraph::NodeIt n(g); n != INVALID; ++n){
		for(ListDigraph::NodeIt v(g); v != INVALID; ++v){

			//no edges from bigger nodes to smaller to ensure acyclicity,
			//and no edges from node to itself
			//+ an attempt to create longer graphs
			if(labels[n] >= labels[v] || labels[n] < labels[v]-20) continue;

			if(rand()%100 <= edge_prob*100){
				g.addArc(n, v);
			}
		}
	}
}

void check_flow_on_nodes(ListDigraph& g, ListDigraph::ArcMap<int>& flow, ListDigraph::Node s, ListDigraph::Node t){
  for(ListDigraph::NodeIt n(g); n != INVALID; ++n){
    if(n == s || n == t) continue;
    int incoming_flow = 0;
    for(ListDigraph::InArcIt i(g,n); i != INVALID; ++i){
      incoming_flow += flow[i];
    }
    REQUIRE(incoming_flow > 0);
    for(ListDigraph::OutArcIt o(g,n); o != INVALID; ++o){
      incoming_flow -= flow[o];
    }
    REQUIRE(incoming_flow == 0);
  }
}



TEST_CASE("feasible minflow is generated for a random graph"){
  srand(time(NULL));
  ListDigraph g;
  createRandomGraph(g, 100, 0.9);
 
  ListDigraph::Node s, t;
  
  s = add_source(g);
  t = add_sink(g);

  ListDigraph::ArcMap<int> flow(g, 0);
  
  find_minflow(g, flow, s, t);
    
  check_flow_on_nodes(g, flow, s, t);
  
}


TEST_CASE("feasible minflow is generated for a random graph 2"){
  srand(time(NULL));
  ListDigraph g;
  createRandomGraph(g, 1000, 0.9);
 
  ListDigraph::Node s, t;
  
  s = add_source(g);
  t = add_sink(g);

  ListDigraph::ArcMap<int> flow(g, 0);
  
  find_minflow(g, flow, s, t);
    
  check_flow_on_nodes(g, flow, s, t);
  
}

TEST_CASE("feasible minflow is generated for a random graph 3"){
  srand(time(NULL));
  ListDigraph g;
  createRandomGraph(g, 1000, 0.6);
 
  ListDigraph::Node s, t;
  
  s = add_source(g);
  t = add_sink(g);

  ListDigraph::ArcMap<int> flow(g, 0);
  
  find_minflow(g, flow, s, t);
    
  check_flow_on_nodes(g, flow, s, t);
  
}