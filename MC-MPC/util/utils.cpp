#include <lemon/list_graph.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <fstream>

using namespace std;
using namespace lemon;

//adds and returns a new source node, which has an outbound arc to all previous source nodes of the graph
ListDigraph::Node add_source(ListDigraph& g){
  ListDigraph::Node s = g.addNode();
  for(ListDigraph::NodeIt n(g); n != INVALID; ++n){
    
        if(countInArcs(g, n) == 0 && n != s){
            g.addArc(s, n);
        }
    }
    
  return s;
}

//adds a sink node so that every sink node in the graph has an outbound arc to the new sink, which is returned
ListDigraph::Node add_sink(ListDigraph& g){

  ListDigraph::Node t = g.addNode();  
  for(ListDigraph::NodeIt n(g); n != INVALID; ++n){
    
        if(countOutArcs(g, n) == 0 && n != t){
            g.addArc(n, t);
        }
    }
    
  return t;
}

bool directory_exists(string directory_name)
{
  //http://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c
  struct stat info;
  if( stat( directory_name.c_str(), &info ) != 0 )
    return false;
  else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows 
    return true;
  else
    return false;
}

bool file_exists(string filename)
{
  std::ifstream infile(filename.c_str());
  if(infile.good()){
    return true;
  }
  return false;
}


//functions for visualizing graphs in dot format
void drawGraphToFile(ListDigraph& g, string filename){
	ofstream myfile;
	myfile.open(filename.c_str());
	myfile << "digraph g {\n";
	for (ListDigraph::ArcIt a(g); a!= INVALID; ++a)
	{
		myfile << g.id(g.source(a)) << " -> " << g.id(g.target(a)) <<  "\n";
	}
	myfile << "}\n";
	myfile.close();
}
/**
void drawGraphToFileWithArcMap(ListDigraph& g, ListDigraph::ArcMap<int>& map, string filename){
	ofstream myfile;
	myfile.open(filename);
	myfile << "digraph g {\n";
	for (ListDigraph::ArcIt a(g); a!= INVALID; ++a)
	{
		myfile << g.id(g.source(a)) << " -> " << g.id(g.target(a)) << " [label=\"" << map[a] << "\"] \n";
	}
	myfile << "}\n";
	myfile.close();
}**/

void drawGraphToFileWithArcMap(ListDigraph& g, ListDigraph::ArcMap<int>& map, string filename){
	ofstream myfile;
	myfile.open(filename.c_str());
	myfile << "digraph g {\n";
	for (ListDigraph::ArcIt a(g); a!= INVALID; ++a)
	{

    if(map[a] == 1){
      myfile << g.id(g.source(a)) << " -> " << g.id(g.target(a)) << " [label=\"" << map[a] << "\", color=blue] \n";
    } else {
		  myfile << g.id(g.source(a)) << " -> " << g.id(g.target(a)) << " [label=\"" << map[a] << "\"] \n";
    }
	}
	myfile << "}\n";
	myfile.close();
}
