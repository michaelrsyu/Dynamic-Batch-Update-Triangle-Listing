# Dynamic-Batch-Update-Triangle-Listing

## Input Format
We assume the input graph text file is named *graph.txt* and follows the SNAP edge-list format.

## Pre-process
We pre-process the text file into a *graph-sort.bin* that contains the graph in binary format. We do so by running the following...

```bash
$ ./run txt-to-bin <directory containing graph.txt> 
```

## Usage

After *graph-sort.bin* is generated, the remaining can be run with the following commands...

```bash
$ ./run <batch-size> <alg> <thread-count> 
```
### Arugement Description:
0) run
1) Float [0.0-100.0] for size of batch as a percentage of graph size
2) Number corresponding to the algorithm - DPTL+[4], DPTL[2], Naive[0]
3) Path to directory containing graph-sort.bin 
4) Thread count 
