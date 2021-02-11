sink(file="../tables/rad_models.tex")
cat.models(c("m.1", "m.2", "m.3"), labs=1:3)
sink()

sink(file="../tables/rna_models.tex")
cat.models(c("m.4", "m.5", "m.6"))
sink()

sink(file="../tables/pca_rad_models.tex")
cat.models(c("p.1", "p.2", "p.3"), labs=c(1, "1, 2", "1, 2, 3"))
sink()

sink(file="../tables/pca_rna_models.tex")
cat.models(c("p.4", "p.5", "p.6"))
sink()

sink(file="../tables/pca_rna_rad_models.tex")
cat.models(c("p.7", "p.8", "p.9"))
sink()

sink(file="../tables/other_models.tex")
cat.models(c("m.hl", "m.cust"))
sink()

sink(file="../tables/pca_pv.tex")
N.pcs <- 8
cat("&", paste(1:N.pcs, collapse=" & "), "\\\\\n \\hline \n",
    "Radiomics &", paste(sprintf("%.3f", 
    (cumsum(TCGA.rad.eigen$values) / sum(TCGA.rad.eigen$values))[1:N.pcs]), 
	     collapse=" & "),"\\\\\n",
    "mRNA &", paste(sprintf("%.3f", 
    (cumsum(TCGA.rna.eigen$values) / sum(TCGA.rna.eigen$values))[1:N.pcs]),
	     collapse=" & "),"\n")
sink()


sink(file="../tables/networks.tex")
    cat.networks(mods, mod.names, ncol(TCGA.rna.std))
sink()
