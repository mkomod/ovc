l1_range <- seq(0.1, 3, length.out=30)
l2_range <- seq(0.1, 2, length.out=20)
l1l2 <- expand.grid(l1_range, l2_range)

cl <- parallel::makeCluster(getOption("cl.cores", CORES))
parallel::clusterExport(cl, 
	    varlist=c("l1l2", "TCGA.rna.std", "TCGA.rad.std", "PERMUTATIONS"))

res <- parallel::parSapplyLB(cl, 1:nrow(l1l2), function(i) {
    l1 <- l1l2[i, 1]
    l2 <- l1l2[i, 2]
    p <- rcca::cca_permutation_validation(TCGA.rna.std, TCGA.rad.std, l1, l2,
		PERMUTATIONS, 1000, 1e-6, F)
    cat("(", l1, ",", l2, "): ", sprintf("%.3f", p), "\n", sep="")
    return(p)
})

rcca.pv <- cbind(l1l2, res)
