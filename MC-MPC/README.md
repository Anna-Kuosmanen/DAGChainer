# Faster algorithms to minimum path cover by graph decomposition

This is the repository for my master's thesis. It consists of several smaller programs which are used together

`decomposer` decomposes directed acyclic graphs into smaller components using maximum anti-chains as separators. Cutting a graph along maximum anti-chains allows us to solve some problems, like minimum path cover, in decomposed parts instead of in the full graph. 

`split` is a tool for splitting directed graphs' connected components into different files

`mc-mpc-solver` solves minimum-cost minimum path cover in a graph by solving minimum flow

`generator` is a generator tool for randomized k-path graphs

`util` contains common helper functions used by other programs