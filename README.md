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

Part of the analysis are quite computationally intensive. They can be re-run by editing main and setting `RUN_ALL <- TRUE`. However, we advise running the analysis on a computer cluster.

