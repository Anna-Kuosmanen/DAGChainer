# Split

This is a tool for splitting the input graph to its connected components, so that every connected component is in its own file. Splitting preserves node and arc labeling and arc weights. This is designed to be used with directed graphs.

## Usage

Split takes two arguments, input graph file name and output folder. Input (and output) graphs are in [LEMON graph format](http://lemon.cs.elte.hu/pub/tutorial/a00018.html). Only @nodes and @arcs sections are required. Nodes should have labels and arcs should have labels and weights. See the example `example_graph`. Decomposition parts are outputted in individual files named `<filename> + _split_ + <index>`.

```
mkdir output
make
./split example_graph output/
```
