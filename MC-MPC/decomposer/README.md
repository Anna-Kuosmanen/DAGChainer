# Decomposer

This is a decomposer for directed acyclic graphs, optimized for graphs with minimum path covers of small size.

## Usage

Decomposer takes two arguments, input graph file name and output folder. Input (and output) graphs are in [LEMON graph format](http://lemon.cs.elte.hu/pub/tutorial/a00018.html). Only @nodes and @arcs sections are required. Nodes should have labels and arcs should have labels and weights. See the example `example_graph`. Decomposition parts are outputted in individual files named `decomp_ + <index>`.

```
mkdir output
make
./decompose example_graph output/
```

## Running tests

Tests use [catch](https://github.com/philsquared/Catch).  

```
wget https://github.com/philsquared/Catch/releases/download/v1.8.2/catch.hpp
make test
./decomposition_tests
```