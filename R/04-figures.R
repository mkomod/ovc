# Correlation Structure ---
rgb.palette <- colorRampPalette(c("red", "black", "blue"), space = "rgb")
lattice::levelplot(TCGA.rna.cor, 
		   main="Correlation matrix of mRNA expression data", 
		   col.regions=rgb.palette(150), cuts=100,  pretty=T)
lattice::levelplot(TCGA.rad.cor,
		   main="Correlation matrix of mRNA expression data", 
		   col.regions=rgb.palette(150), cuts=100)
lattice::levelplot(TCGA.rna.rad.cor,
		   xlab="", ylab="", 
		   col.regions=rgb.palette(150), cuts=100, xaxt="n", pretty=T)


# CCA permutation validation ---
pdf(file="../figures/rcca_tuning_hyperparameters.pdf", width=9, height=7.5)
par(family="Times", cex.lab=1.3, cex.axis=1.3)
filled.contour(z=(rcca.pv$z), x=rcca.pv$x, y=rcca.pv$y, 
	       color.palette=function(n) hcl.colors(n),
	       key.title = title(main = "p-value", cex.main=1, font.main=3),
	       xlab=expression(lambda[1]), ylab=expression(lambda[2]))
dev.off()


# CCA Weights ---
pdf(file="../figures/rcca_u_1-3_v_1-3.pdf", width=9, height=6.5)
w1 <- TCGA.cca$w1
w2 <- TCGA.cca$w2
par(mfrow=c(2, 3), family="Times", mar=c(5, 3, 3, 1))
for (i in 1:3) {
    pch.col <- fifelse(w1[ , i] == 0, rgb(0, 0, 1, 0.04), rgb(0, 0, 1, 0.8))
    plot(w1[ , i], pch=20, ylim=c(-0.1, 0.1), col=pch.col,
	 main=bquote(w[1]^(.(i))), ylab="", las=1,
	 cex.lab=1.4, cex.axis=1.4, cex.main=1.7)
    mtext(paste0("(", letters[i], ")"), side=3, adj=-0.15, line=0.9, cex=1.3)
}
for (i in 1:3) {
    pch.col <- fifelse(w2[ , i] == 0, rgb(0, 0, 1, 0.3), rgb(0, 0, 1, 0.8))
    plot(w2[ , i], pch=20, ylim=c(-0.1, 0.1), col=pch.col,
	 main=bquote(w[2]^(.(i))), ylab="", las=1,
	 cex.lab=1.4, cex.axis=1.4, cex.main=1.7)
    mtext(paste0("(", letters[i + 3], ")"), side=3, adj=-0.15, line=0.9, cex=1.3)
}
dev.off()


# PCA Projections ---
pdf(file="../figures/pca_projections.pdf", width=9, height=5)
{
    par(mfrow=c(1, 2), family="Times", mar=c(5, 4.2, 3, 1))
    plot(TCGA.rna.pc.proj[ , 1], TCGA.rna.pc.proj[ , 2],
	 main="mRNA", font.main=1,
	 xlab=bquote("Principal component 1"),
	 ylab=bquote("Principal component 2"), yaxt="n", pch=20, cex=0.8)
    mtext(paste0("(", letters[1], ")"), side=3, adj=-0.15, line=0.9, cex=1.3)
    axis(2, las=1)
    plot(TCGA.rad.pc.proj[ , 1], TCGA.rad.pc.proj[ , 2],
	 main="Radiomics", font.main=1,
	 xlab=bquote("Principal component 1"),
	 ylab=bquote("Principal component 2"), yaxt="n", pch=20, cex=0.8)
    mtext(paste0("(", letters[2], ")"), side=3, adj=-0.15, line=0.9, cex=1.3)
    axis(2, las=1)
}
dev.off()


# PCA cross projections ---
pdf(file="../figures/pca_cross_projections.pdf", width=9, height=3.5)
par(mfrow=c(1, 3), family="Times", mar=c(5, 4, 3, 1), cex.lab=1.5, cex.main=1.5)
for (i in 1:3) {
    plot(TCGA.rna.pc.proj[ , i], TCGA.rad.pc.proj[ , i],
	 main=paste("Principal component", i), font.main=1,
	 xlab="mRNA",
	 ylab="Radiomics", yaxt="n", pch=20)
    mtext(paste0("(", letters[i], ")"), side=3, adj=-0.15, line=0.9, cex=1.3)
    axis(2, las=1)
}
dev.off()


# CCA Projections ---
pdf(file="../figures/cca_projections.pdf", width=9, height=3.5)
par(mfrow=c(1, 3), family="Times", mar=c(5, 4.2, 3, 1), cex.lab=1.5, cex.main=1.5)
for (i in 1:3) {
    plot(TCGA.rna.proj[ , i], TCGA.rad.proj[ , i],
	 main=paste("Canonical Variates", i), font.main=1,
	 ylab="Radiomics", 
	 xlab="mRNA", yaxt="n", pch=20, cex=0.8)
    mtext(paste0("(", letters[i], ")"), side=3, adj=-0.13, line=0.8, cex=1.3)
    axis(2, las=1)
}
dev.off()


# Networks ---
pdf(file="../figures/networks.pdf", width=9, height=6.5)
par(mfrow=c(2, 3), family="Times", mar=c(5, 3, 3, 1))
for(i in seq_along(1:6)) {
    SmCCNet::plotMultiOmicsNetwork(A_bar, cor(cbind(TCGA.rna.std, TCGA.rad.std)), 
	mods, i, P1=ncol(TCGA.rna.std), NetLayout="lgl",
	FeatureLabel=c(colnames(TCGA.rna.std), colnames(TCGA.rad.std)))
    # readline(sprintf("[%d / %d] Press Enter for next: ", i, length(mods)))
    mtext(paste0("(", letters[i], ")"), side=3, adj=-0.13, line=0.8, cex=1.3)
}
dev.off()
