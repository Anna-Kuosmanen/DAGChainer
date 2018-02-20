#include <lemon/list_graph.h>
#include "../decomposition.h"
#include "../MPC.h"
#include <stdlib.h>
#include <time.h>
#include "../catch.hpp"
#include <lemon/bfs.h>

using namespace std;
using namespace lemon;


TEST_CASE("decomposing an easy test case")
{
  ListDigraph g;
  ListDigraph::Node s = g.addNode();
  ListDigraph::Node t = g.addNode();
  ListDigraph::Node a = g.addNode();
  ListDigraph::Node b = g.addNode();

  ListDigraph::Arc sa = g.addArc(s, a);
  ListDigraph::Arc bt = g.addArc(b, t);
  ListDigraph::Arc sb = g.addArc(s, b);
  ListDigraph::Arc at = g.addArc(a, t);

  ListDigraph::ArcMap<int> minflow(g);
  find_minflow(g, minflow, s, t);

  ListDigraph::ArcMap<int> decomposition(g);

  decompose_graph(g, minflow, s, t, decomposition);

  REQUIRE(decomposition[sa] == decomposition[sb]);
  REQUIRE(decomposition[at] == decomposition[bt]);
  REQUIRE(decomposition[sa] != decomposition[bt]);
}

TEST_CASE("decomposing an easy test case 2")
{
  ListDigraph g;
  ListDigraph::Node s = g.addNode();
  ListDigraph::Node t = g.addNode();

  //node g's variable name is 'i'

  ListDigraph::Node a = g.addNode();
  ListDigraph::Node b = g.addNode();
  ListDigraph::Node c = g.addNode();
  ListDigraph::Node d = g.addNode();
  ListDigraph::Node e = g.addNode();
  ListDigraph::Node f = g.addNode();
  ListDigraph::Node i = g.addNode();
  ListDigraph::Node h = g.addNode();

  ListDigraph::Arc sa = g.addArc(s, a);
  ListDigraph::Arc sb = g.addArc(s, b);
  ListDigraph::Arc sc = g.addArc(s, c);
  ListDigraph::Arc ad = g.addArc(a, d);
  ListDigraph::Arc bd = g.addArc(b, d);
  ListDigraph::Arc ce = g.addArc(c, e);
  ListDigraph::Arc df = g.addArc(d, f);
  ListDigraph::Arc dg = g.addArc(d, i);
  ListDigraph::Arc eh = g.addArc(e, h);
  ListDigraph::Arc eg = g.addArc(e, i);
  ListDigraph::Arc ft = g.addArc(f, t);
  ListDigraph::Arc gt = g.addArc(i, t);
  ListDigraph::Arc ht = g.addArc(h, t);

  ListDigraph::ArcMap<int> minflow(g);
  find_minflow(g, minflow, s, t);

  ListDigraph::ArcMap<int> decomposition(g);

  decompose_graph(g, minflow, s, t, decomposition);

  REQUIRE(decomposition[sa] == decomposition[sb]);
  REQUIRE(decomposition[sb] == decomposition[sc]);

  REQUIRE(decomposition[ad] == decomposition[eh]);
  REQUIRE(decomposition[dg] == decomposition[eh]);
  REQUIRE(decomposition[eh] == decomposition[df]);

  REQUIRE(decomposition[ft] == decomposition[gt]);
  REQUIRE(decomposition[gt] == decomposition[ht]);

  REQUIRE(decomposition[sa] != decomposition[dg]);
  REQUIRE(decomposition[ad] != decomposition[ft]);
  REQUIRE(decomposition[bd] != decomposition[sc]);
  REQUIRE(decomposition[ht] != decomposition[eh]);
}