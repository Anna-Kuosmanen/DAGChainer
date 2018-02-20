#!/bin/bash
### compile all programs

cd decomposer && make && cd ..
cd generator && make && cd ..
cd mc-mpc-solver && make && cd ..
cd split && make && cd ..

### move all binaries to bin/

mkdir -p bin
cp decomposer/decompose bin/
cp generator/generator bin/
cp split/split bin/
cp mc-mpc-solver/mc-mpc bin/