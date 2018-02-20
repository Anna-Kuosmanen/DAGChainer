#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include <stack>
#include <stdlib.h>
#include "../util/utils.h"


using namespace lemon;
using namespace std;

void find_feasible_flow(ListDigraph& g, ListDigraph::ArcMap<int>& flow, ListDigraph::Node s, ListDigraph::Node t){
  int num_nodes = countNodes(g);
  ListDigraph::NodeMap<int> top_order(g);
  topologicalSort(g, top_order);
  //ListDigraph::Node top_nodes[num_nodes];
  ListDigraph::Node* top_nodes = new ListDigraph::Node[num_nodes];

  ListDigraph::NodeMap<bool> covered(g);
  covered[s] = true;
  covered[t] = true;

  //top nodes will have nodes in descending topological order
  for(ListDigraph::NodeIt n(g); n != INVALID; ++n){
    top_nodes[num_nodes-top_order[n]-1] = n;
  }

  //this node map tells us the maximum number of uncovered nodes reachable from this node via single path
  ListDigraph::NodeMap<int> uncovered(g, 0);

  while(true){

    //update uncovered
    for(int i = 0; i < num_nodes; i++){
      ListDigraph::Node temp = top_nodes[i];
      uncovered[temp] = 0;
      for (ListDigraph::OutArcIt o(g, temp); o != INVALID; ++o)
      {
       uncovered[temp] = max(uncovered[g.target(o)], uncovered[temp]);
      }
      if(!covered[temp]){
        uncovered[temp]++;
      }
    }

    ListDigraph::Node current = s;
    while(current != t){
      int max_uncovered = 0;
      ListDigraph::Arc nextArc;
      for(ListDigraph::OutArcIt o(g, current); o != INVALID; ++o){
        if(uncovered[g.target(o)] > max_uncovered){
          nextArc = o;
          max_uncovered = uncovered[g.target(o)];
        }
      }
      //if there is no uncovered nodes at the source, every node is covered and we are ready
      //else we just choose an arbitrary edge (there are no reachable uncovered nodes)
      if(max_uncovered == 0){
        if(current == s) return;
        ListDigraph::OutArcIt o(g, current);
        nextArc = o;
      }

      flow[nextArc]++;
      current = g.target(nextArc);
      covered[current] = true;
    }
  }


}

//returns true if an augmenting path is found, false otherwise. This uses no demands on arcs
bool find_augmenting_path(ListDigraph& g, ListDigraph::ArcMap<int>& flow, ListDigraph::Node s, ListDigraph::Node t)
{

  ListDigraph::NodeMap<bool> visited(g, false);
  visited[s] = true;

  //the boolean value tells if it is a forward arc
  stack<pair<ListDigraph::Arc, bool> > path;
  ListDigraph::Node extraNode = g.addNode();
  ListDigraph::Arc extraArc = g.addArc(extraNode, s);
  visited[extraNode] = true;
  path.push(make_pair(extraArc, true));

  ListDigraph::NodeMap<int> flow_on_node(g);
  for(ListDigraph::NodeIt n(g); n != INVALID; ++n){
    for(ListDigraph::OutArcIt o(g, n); o != INVALID; ++o){
      flow_on_node[n] += flow[o];
    }
  }

  ListDigraph::ArcMap<pair<bool, bool> > visitedArc(g, make_pair(false, false));

  loop: while(!path.empty()){

    ListDigraph::Arc currentArc = path.top().first;
    ListDigraph::Node currentNode;
    if(path.top().second){
      currentNode = g.target(currentArc);
      visitedArc[currentArc].first = true;
    }
    else{
      currentNode = g.source(currentArc);
      visitedArc[currentArc].second = true;
    }

    visited[currentNode] = true;

    if(currentNode == t){
      
      while(!path.empty()){
        if(path.top().second){
         flow[path.top().first]--;
        }else{
          flow[path.top().first]++;
        }
        path.pop();
      }

        g.erase(extraArc);
        g.erase(extraNode);
      return true;
    }

    //we can use forward arcs only if there is extra flow on that node
    if(flow_on_node[currentNode] > 0 || !path.top().second){

      for(ListDigraph::OutArcIt o(g, currentNode); o != INVALID; ++o){
        if(flow[o] > 0 && !visitedArc[o].first && o != currentArc){
          path.push(make_pair(o, true));
          flow_on_node[g.target(o)]--;
          goto loop;
        }
      }
    }
    //and the backward arcs
    for(ListDigraph::InArcIt i(g, currentNode); i != INVALID; ++i){
      if(!visitedArc[i].second && i != currentArc){
        path.push(make_pair(i, false));
        flow_on_node[g.source(i)]++;
        goto loop;
      }
    }

    path.pop();
  }

  g.erase(extraArc);
  g.erase(extraNode);
  return false;
}

void find_minflow(ListDigraph& g, ListDigraph::ArcMap<int>& flow, ListDigraph::Node s, ListDigraph::Node t)
{
  find_feasible_flow(g, flow, s, t);
  //drawGraphToFileWithArcMap(g, flow, "feasible.dot");
  int index = 0;
  while(find_augmenting_path(g, flow, s, t)){
    index++;
  }
  //cout << "INDEX: " << to_string(index) << endl;
  //drawGraphToFileWithArcMap(g, flow, "minflow.dot");
}

