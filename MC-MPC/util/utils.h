#ifndef UTILS_H_
#define UTILS_H_

#include <lemon/list_graph.h>
#include <string.h>

using namespace lemon;
using namespace std;

ListDigraph::Node add_source(ListDigraph& g);
ListDigraph::Node add_sink(ListDigraph& g);

bool file_exists(std::string filename);
bool directory_exists(std::string filename);

void drawGraphToFile(ListDigraph& g, std::string filename);
void drawGraphToFileWithArcMap(ListDigraph& g, ListDigraph::ArcMap<int>& map, string filename);


#endif /* UTILS_H_ */
