# Generator

This is a generator for randomized k-path graphs. K-path graphs are guaranteed to have a path cover of at most k paths. The parameters given to the generator are

- k: number of paths
- n: length of one path. This makes whole graph to have k*n nodes
- m: number of extra arcs. Extra arcs are added between random nodes, as long as the graph retains it acyclicity.

## Usage

Generator takes four arguments, output graph file name, k, n and m. Graphs are output in [LEMON graph format](http://lemon.cs.elte.hu/pub/tutorial/a00018.html). Nodes have labels and arcs have labels and weights. Arc weights are randomized between 0 and `MAX_WEIGHT` (defined in code).

```
make
./generator output_graph 5 100 100
```
