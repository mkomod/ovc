# Integrating multi-omics

Code to reproduce the analysis of our report "Integrating Multi-omics through Sparse Canonical Correlation Analysis". 

## Usage

```
$ cd data
$ gunzip *
$ cd ../R
$ R
> source("main.R")
```

Figures and tables will be added to the `./figures` and `./tables` directories respectively.


## Notes

Part of the analysis are quite computationally intensive. By default we have disabled the re-computation of certain values, however, they can be re-run by editing main and setting `RUN_ALL <- TRUE`. We advise running the analysis on a computer cluster.

