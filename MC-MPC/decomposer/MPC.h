#ifndef MPC_H_
#define MPC_H_

#include <lemon/list_graph.h>

using namespace lemon;

void find_feasible_flow(ListDigraph& g, ListDigraph::ArcMap<int>& flow, ListDigraph::Node s, ListDigraph::Node t);

void find_minflow(ListDigraph& g, ListDigraph::ArcMap<int>& flow, ListDigraph::Node s, ListDigraph::Node t);

#endif /* MPC_H_ */
