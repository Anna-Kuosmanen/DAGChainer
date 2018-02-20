# Minimum-cost minimum path cover solver

This program takes a graph as a input and outputs the minimum cotst path cover in output file. The result is given in the following format: One line for each arc of the graph, and each line contains arc's label and number of paths going through that arc, separated by space. Output for `example_graph`in the directory:

```
3 0
1 1
5 1
2 1
4 1
6 1
7 1
```

## Usage

Solver takes two arguments, input graph file name and output file name. Input graphs are in [LEMON graph format](http://lemon.cs.elte.hu/pub/tutorial/a00018.html). Only @nodes and @arcs sections are required. Nodes should have labels and arcs should have labels and weights. See the example `example_graph`.

```
make
mc-mpc example_graph output.txt
```
