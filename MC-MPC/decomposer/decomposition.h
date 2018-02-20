#ifndef DECOMPOSITION_H_
#define DECOMPOSITION_H_

#include <lemon/list_graph.h>

using namespace lemon;
using namespace std;

void decompose_graph(ListDigraph& g, ListDigraph::ArcMap<int>& minFlow, ListDigraph::Node s, ListDigraph::Node t, ListDigraph::ArcMap<int>& decomposition);

#endif /* DECOMPOSITION_H_ */
