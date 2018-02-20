#include <lemon/list_graph.h>
#include <time.h>
#include <stack>
#include <climits>

using namespace lemon;
using namespace std;


bool can_reach_another_node(int path_index, ListDigraph::NodeMap<int*>& reachable, vector<ListDigraph::Node>* paths, int* ants, int num_paths)
{
	for(int i = 0; i< num_paths; i++){
		if(i == path_index) continue;
		if(reachable[paths[path_index][ants[path_index]]][i] <= ants[i]) return true;
	}
	return false;
}

bool is_MAC(ListDigraph::NodeMap<int*>& reachable, vector<ListDigraph::Node>* paths, int* ants, int num_paths)
{
	for(int i = 0; i< num_paths; i++){
		if(can_reach_another_node(i, reachable, paths, ants, num_paths)) return false;
	}
	return true;
}

void decompose_graph(ListDigraph& g, ListDigraph::ArcMap<int>& minFlow, ListDigraph::Node s, ListDigraph::Node t, ListDigraph::ArcMap<int>& decomposition){

	int num_paths = 0;
	for(ListDigraph::OutArcIt o(g, s); o != INVALID; ++o){
		num_paths += minFlow[o];
	}

	//Extract paths from the graph according to the minFlow

	vector<ListDigraph::Node>* paths = new vector<ListDigraph::Node>[num_paths];
	ListDigraph::NodeMap<int*> reachable(g);

	//paths are picked up one path at time
	for (int i = 0; i < num_paths; ++i)
	{
		ListDigraph::Node node = s;
		int path_index = 0;
		while(node != t){
			if(reachable[node] == NULL){
				reachable[node] = (int*) calloc(num_paths, sizeof(int));
				for (int i = 0; i < num_paths; ++i)
				{
					reachable[node][i] = INT_MAX;
				}
			}
			reachable[node][i] = path_index;

			for (ListDigraph::OutArcIt o(g, node); o != INVALID; ++o)
			{
				if(minFlow[o] > 0){
					minFlow[o] -= 1;
					paths[i].push_back(node);
					node = g.target(o);
					break;
				}
			}
			path_index++;
		}

		if(reachable[t] == NULL){
			reachable[t] = (int*) calloc(num_paths, sizeof(int));
		}
		reachable[t][i] = path_index;
		paths[i].push_back(t);
	}

	// Create reachability tables for nodes using a DFS. We know the reachability table for the sink node 
	// (it can't reach any other nodes), so we can recursively deduce other reachability tables from that

	stack<ListDigraph::Node> node_stack;
	node_stack.push(s);
	ListDigraph::NodeMap<bool> visited(g);

	loop: while(!node_stack.empty()){
		ListDigraph::Node node = node_stack.top();

		// If any of the child nodes of a node is still unvisited, the child is added to the stack and the loop is restarted
		for(ListDigraph::OutArcIt o(g, node); o != INVALID; ++o){
			ListDigraph::Node child = g.target(o);
			if(!visited[child]){
				node_stack.push(child);
				goto loop;
			}
		}
		// All the child nodes are visited, so we can calculate the reachability tables from the reachability tables of the child nodes
		for(ListDigraph::OutArcIt o(g, node); o != INVALID; ++o){
			ListDigraph::Node child = g.target(o);
			for (int i = 0; i < num_paths; ++i)
			{
				reachable[node][i] = min(reachable[node][i], reachable[child][i]);
			}
		}
		// The node can be marked as visited because it's reachability table is computed
		visited[node] = true;
		node_stack.pop();
	}

	// Then we create the decomposition by running a set of k "ants" in the graph, each ant following a different path
	// As we pass through a node, we color that node with the current color
	// When the configuration of ants forms a MAC, we increase the color by one, causing the following nodes to be colored
	// with a different color. 

	int color = 1;

	ListDigraph::NodeMap<int> node_color(g, 0);
	int ants[num_paths];
	for (int i = 0; i < num_paths; ++i){
		ants[i] = 0;
		node_color[paths[i][0]] = color;
	}

	while(true){

		if(is_MAC(reachable, paths, ants, num_paths)){
			
			// When a MAC is found, all the upcoming nodes will belong to the next part
			color++;

			for(int i = 0; i < num_paths; ++i){
				//update the colors of current MAC nodes to the new decomposition color
				node_color[paths[i][ants[i]]] = color;
			}
			// Move every ant forwards by 1
			for (int i = 0; i < num_paths; ++i)
			{
				ants[i]++;
				//after moving an ant, color all the incoming arcs that belong to same part with current color. if there are arcs that span from earlier parts, mark them as discarded
				ListDigraph::Node current = paths[i][ants[i]];
				node_color[current] = color;
				for(ListDigraph::InArcIt ai(g, current); ai != INVALID; ++ai){
					if(node_color[g.source(ai)] == color || node_color[g.source(ai)] == 0){
						decomposition[ai] = color;
					} else {
						decomposition[ai] = -1;
					}
				}
			}
		}

		else{
			// If current ant configuration is not a MAC, we have to move ants forwards.
			// By moving each ant forwards until it can't reach any other ant we create a configuration of ants,
			// which is not necessarily a MAC, but we won't miss any MACs.
			for (int i = 0; i < num_paths; ++i)
			{
				while(can_reach_another_node(i, reachable, paths, ants, num_paths) && paths[i][ants[i]] != t)
				{
					ants[i]++;
					//after moving an ant, color all the incoming arcs that belong to same part with current color. if there are arcs that span from other parts, mark them as discarded
					ListDigraph::Node current = paths[i][ants[i]];
					node_color[current] = color;
					for(ListDigraph::InArcIt ai(g, current); ai != INVALID; ++ai){
						if(node_color[g.source(ai)] == color || node_color[g.source(ai)] == 0){
							decomposition[ai] = color;
						} else {
							decomposition[ai] = -1;
						}
					}
				}
			}
		}

		// terminate only when all ants are in the sink node
		bool terminate = true;
		for (int i = 0; i < num_paths; ++i)
		{
			if(paths[i][ants[i]] != t){
				terminate = false;
			}
		}
		if(terminate) break;
	}
}