# Cross validation
rCCA_CV <- function(X1, X2, l1_range, l2_range, folds = 5, 
		    niter = 1000, threshold = 1.0e-6, verbose= TRUE, seed=1) 
{
    set.seed(seed)

    training_loss <- array(dim=c(length(l1_range), length(l2_range), folds))
    testing_loss <- array(dim=c(length(l1_range), length(l2_range), folds))

    perm <- sample(1:nrow(X1), nrow(X1), replace=F)
    X1 <- X1[perm, ]                   # shuffle the rows
    X2 <- X2[perm, ]

    fold.size <- ceiling(nrow(X1) / folds)
    for (fold in 1:folds) {
	sub <- (fold.size * (fold - 1) + 1):(min(nrow(X1), (fold.size * fold)))
	X1.train <- X1[ setdiff(1:nrow(X1), sub) , ]
	X2.train <- X2[ setdiff(1:nrow(X2), sub) , ]
	X1.test <- X1[sub , ]
	X2.test <- X2[sub , ]
	B.train <- t(X1.train) %*% X2.train
	B.test <- t(X1.test) %*% X2.test
	
	for (l1.index in seq_along(l1_range)) {
	    l1 <- l1_range[l1.index]
	    for (l2.index in seq_along(l2_range)) {
		l2 <- l2_range[l2.index]
		res <- rCCA(X1.train, X2.train, l1, l2, K=1, niter=niter, 
			    threshold=threshold, verbose=FALSE)
		w1 <- res$w1
		w2 <- res$w2
		training_loss[l1.index, l2.index, fold] <- 
		    t(w1)%*% B.train %*% w2
		
		loss <- cor(X1.test %*% w1, X2.test %*% w2)
		if (is.na(loss)) loss <- 0
		testing_loss[l1.index, l2.index, fold] <- loss
		    
	    }
	}
	print(training_loss[ , , fold])
	print(testing_loss[ , , fold])
    }

    training_loss_mean <- apply(training_loss, c(1,2), mean) 
    testing_loss_mean <- apply(testing_loss, c(1,2), mean) 
    return(list(training_loss = training_loss_mean, 
		testing_loss = testing_loss_mean))
}


# Permutation based validaiton
rCCA_permute <- function(X1, X2, l1_range, l2_range, permutations=50, seed=1, 
    niter = 1000, threshold = 1.0e-6, verbose= TRUE, progress=TRUE) 
{
    set.seed(seed)
    p_values  <- matrix(nrow=length(l1_range), ncol=length(l2_range))
    for (l1.index in seq_along(l1_range)) {
	for (l2.index in seq_along(l2_range)) {
	    l1 <- l1_range[l1.index]
	    l2 <- l2_range[l2.index]
	    res <- rCCA(X1, X2, l1, l2, K=1, niter=niter, threshold=threshold, 
			verbose=FALSE)
	    d <- cor(X1 %*% res$w1, X2 %*% res$w2)[1, 1]

	    cat("(", l1, ", ", l2, "):", sep="")
	    d.ps <- sapply(1:permutations, function(p) {
		perm <- sample(1:nrow(X1), nrow(X1), replace=F)
		X1.p <- X1[perm, ]
		res.p <- rCCA(X1.p, X2, l1, l2, K=1, niter=niter, 
				 threshold=threshold, verbose=FALSE)
		d.p <- cor(X1.p %*% res.p$w1, X2 %*% res.p$w2)
		if (is.na(d.p)) d.p <- 0
		return(d.p)
	    })
	    p_value <- mean(d.ps >= d)
	    p_values[l1.index, l2.index] <- p_value
	    cat(" p-value: ", p_value, "\n")
	}
    }

    return(list(p_values = p_values, l1_range = l1_range, l2_range = l2_range))
}


# Compute robust weights for w1 and w2 using column subsampling of X1 and X2
robustCCA <- function(X1, X2, l1, l2, niter=1000, verbose=F, subsample=0.7,
		    n_samples=500, seed=1, standardize=T, progress=T) 
{
    set.seed(seed)

    X1.colnames <- colnames(X1)
    X2.colnames <- colnames(X2)
    X1.subsamp_size <- ceiling(subsample * ncol(X1))
    X2.subsamp_size <- ceiling(subsample * ncol(X2))
    
    W1 <- matrix(nrow=n_samples, ncol=ncol(X1))
    W2 <- matrix(nrow=n_samples, ncol=ncol(X2))

    for (i in 1:n_samples) {
	X1.subsamp_cols <- sort(sample(1:ncol(X1), X1.subsamp_size, replace=F))
	X2.subsamp_cols <- sort(sample(1:ncol(X2), X2.subsamp_size, replace=F))
	X1.sub <- X1[ , X1.subsamp_cols]
	X2.sub <- X2[ , X2.subsamp_cols]
	if (standardize) {
	    X1.sub <- scale(X1.sub , T, T)
	    X2.sub <- scale(X2.sub , T, T)
	}

	res <- rCCA(X1.sub, X2.sub, l1, l2, K=1, niter=niter, verbose=verbose)
	w1 <- res$w1
	w2 <- res$w2
	W1[i, X1.subsamp_cols] <- w1
	W2[i, X2.subsamp_cols] <- w2
	if (progress) cat(".")
    }

    return(list(W1=W1, W2=W2)) 
}


# Construct a similarity matrix for SmCCNet (see: Shi et al.)
similarityMatrix <- function(W1, W2) 
{
    W <- cbind(W1, W2)
    A.bar <- matrix(0, nrow=ncol(W), ncol=ncol(W))
    A.count <- matrix(0, nrow=ncol(W), ncol=ncol(W))
    for (i in 1:nrow(W)) {
	A <- outer(abs(W[i, ]), abs(W[i, ]))
	A.count <- A.count + !is.na(A)
	A[is.na(A)] <- 0
	A.bar <- A.bar + A
    }
    diag(A.bar) <- 0
    A.bar <- A.bar / A.count
    A.bar <- A.bar / max(A.bar)
    Matrix::Matrix(A.bar, sparse=T)
}


constructSimilarityMatrix <- function(X1, X2, l1, l2) 
{
    res <- robustCCA(X1, X2, l1, l2)
    Abar <- similarityMatrix(res$W1, res$W2)
    return(list(Abar=Abar, W1=res$W1, W2=res$W2))
}


ccaVariableLimits <- function(X1, X2, X1.n_covariates, X2.n_covariates, 
	X1.names, X2.names, l1, l2) 
{
    sapply(X1.n_covariates, function(n_X1) {
	X1.sub <- X1[ , X1.names[1:n_X1]]
	sapply(X2.n_covariates, function(n_X2) {
	    X2.sub <- X2[ , X2.names[1:n_X2]]
	    res <- rcca::rCCA(X1.sub, X2.sub, l1, l2)
	    return(res)
	})
    })
}


ccaCorrGrid <- function(X1, X2, X1.n_covariates, X2.n_covariates, 
	X1.names, X2.names, l1, l2) 
{
    sapply(X1.n_covariates, function(n_X1) {
	X1.sub <- X1[ , X1.names[1:n_X1]]
	sapply(X2.n_covariates, function(n_X2) {
	    X2.sub <- X2[ , X2.names[1:n_X2]]
	    res <- rcca::rCCA(X1.sub, X2.sub, l1, l2)
	    return(res$loss)
	})
    })
}


grid.as.matrix <- function(g) 
{
    x <- g[ , 1]
    y <- g[ , 2]
    z <- matrix(g[ , 3], nrow=length(unique(x)))
    return(list(x=unique(x), y=unique(y), z=z))
}


plot.persp <- function(x, y, z, ...) 
{
    jet.colors <- colorRampPalette( c("blue", "green") )
    n_col <- 100
    nrz <- nrow(z)
    ncz <- ncol(z)
    color <- jet.colors(n_col)
    z_facet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    facetcol <- cut(z_facet, n_col)
    persp(x, y, z, col = color[facetcol], ...)
}


ccaSensitivity <- function(X1, X2, l1l2)
{
    W1 <- matrix(ncol=0, nrow=ncol(X1))
    W2 <- matrix(ncol=0, nrow=ncol(X2))
    for (i in 1:nrow(l1l2)) {
	l1 <- l1l2[i, 1]
	l2 <- l1l2[i, 2]
	res <- rCCA(X1, X2, l1, l2, K=1)
	w1 <- res$w1
	w2 <- res$w2
	W1 <- cbind(W1, w1)
	W2 <- cbind(W2, w2)
    }
    return(list(W1=W1, W2=W2))
}


plot.ccaSensitivity <- function(cca.sen, w, cca.actual=NULL, l1l2=NULL, l1=NULL,
	l2=NULL, n.sd=2)
{
    if (!is.null(l1)) {
	fixed <- which(l1l2[, 1] == l1)
    } else if (!is.null(l2)) {
	fixed <- which(l1l2[, 2] == l2)
    } else {
	fixed <- 1:ncol(cca.sen[[1]])
    }

    W <- cca.sen[[w]][ , fixed]
    W.mean <- apply(W, 1, mean)
    W.sd <- apply(W, 1, sd)
    W.u <- W.mean + n.sd * W.sd
    W.l <- W.mean - n.sd * W.sd
    pch.color <- fifelse(abs(W.mean) <= 5e-3, 
		     rgb(0, 0, 0, 0.01), rgb(0, 0, 0, 0.8))
    arr.color <- fifelse(abs(W.mean) <= 5e-3, 
		     rgb(0, 0, 0.8, 0.01), rgb(0, 0, 0.8, 0.8))
    plot(W.mean, ylim=c(min(W.l), max(W.u)), pch=20, cex=0.7, col=pch.color)
    arrows(1:nrow(W), W.l, 1:nrow(W), W.u, code=3, angle=90, length=0.02,
    col=arr.color) 
    if (!is.null(cca.actual)) {
	W.actual <- cca.actual[[w]][ , 1]
	pch.color <- fifelse(abs(W.actual) <= 5e-3, 
			 rgb(1, 0, 0, min(0.3, 10/nrow(W))), rgb(1, 0, 0, 0.8))
	points(W.actual, pch=20, col=pch.color)
    }
}


concordance.gonen <- function(model)
{
    b <- matrix(coef(model), ncol=1)
    X <- model.matrix(model)
    n <- nrow(X)
    h <- 0.5 * sd(X %*% b) * n^(-1/3)
    v.b <- model$var
    conc <- 0
    conc.s <- 0
    d.conc.s <- rep(0, ncol(b))
    for (i in 1:(n-1)) {
	x.i <- X[i, ]
	for (j in (i+1):n) {
	    x.j <- X[j, ]
	    x.ij <- x.i - x.j;               	x.ji <- -x.ij
	    b.ij <- as.double(t(b) %*% x.ij);   b.ji <- -b.ij
	    u.ij <- pnorm(- b.ij / h );         u.ji <- pnorm(- b.ji / h )
	    conc <- conc +
		(b.ji < 0) / (1 + exp(b.ji)) + (b.ij < 0) / (1 + exp(b.ij))
	    conc.s <- conc.s + u.ji / (1 + exp(b.ji)) + u.ij / (1 + exp(b.ij))
	    d.conc.s <- d.conc.s +
		-x.ji/h * dnorm(-b.ji/h) * (1 + exp(b.ji))^-1 +
		u.ji * -x.ji * exp(b.ji) * (1 + exp(b.ji))^-2 +
		-x.ij/h * dnorm(-b.ij/h) * (1 + exp(b.ij))^-1 +
		u.ij * -x.ij * exp(b.ij) * (1 + exp(b.ij))^-2
	}
    }
    conc <- 2 * conc / (n * (n-1))
    conc.s <- 2 * conc.s / (n * (n-1))
    d.conc.s <- 2 * d.conc.s / (n * (n-1))

    var.s <- 0
    for (i in 1:(n-1)) {
	x.i <- X[i, ]
	u.j <- 0; u.k <- 0
	for (j in (i+1):n) {
	    x.j <- X[j, ]
	    x.ij <- x.i - x.j;               	x.ji <- -x.ij
	    b.ij <- as.double(t(b) %*% x.ij);   b.ji <- -b.ij
	    u.ij <- pnorm(- b.ij / h );         u.ji <- pnorm(- b.ji / h )
	    u.j <- u.j + u.ji / (1 + exp(b.ji)) + u.ij / (1 + exp(b.ij)) - conc.s
	}
	for (k in 1:i) {
	    x.k <- X[k, ]
	    x.ik <- x.i - x.k;               	x.ki <- -x.ik
	    b.ik <- as.double(t(b) %*% x.ik);   b.ki <- -b.ik
	    u.ik <- pnorm(- b.ik / h );         u.ki <- pnorm(- b.ki / h )
	    u.k <- u.k + u.ki / (1 + exp(b.ki)) + u.ik / (1 + exp(b.ik)) - conc.s
	}
	var.s <- var.s + u.j * u.k
    }
    var.s <- (2 / (n * (n-1)))^2 * var.s + t(d.conc.s) %*% v.b %*% d.conc.s
    var.s <- as.double(var.s)
    return(list(c=conc, c.s=conc.s, se=sqrt(var.s)))
}


concordance.harrell <- function(model)
{
    b <- matrix(coef(model), ncol=1)
    X <- model.matrix(model)
    m.y <- as.numeric(model$y)
    y <- m.y[1:(length(m.y)/2)]
    d <- m.y[(length(m.y)/2 + 1):(length(m.y))]
    n <- nrow(X)
    num <- 0
    den <- 0
    
    for (i in 2:nrow(X)) {
	x.i <- X[i, ]
	y.i <- y[i]
	d.i <- d[i]
	for (j in 1:(i-1)) {
	    x.j <- X[j, ]
	    y.j <- y[j]
	    d.j <- d[j]
	    x.ij <- t(b) %*% x.i - t(b) %*% x.j
	    num <- num +
		(y.i < y.j) * (x.ij > 0) * d.i + (y.j < y.i) * (x.ij < 0) * d.j
	    den <- den +
		(y.i < y.j) * d.i + (y.j < y.i) * d.j
	}
    }
    return(num / den)
}


cat.models <- function(models) 
{
    for (i in 1:length(models)) {
	m <- get(models[i])
	con <- concordance(m)
	kon <- concordance.gonen(m)
	cat(sprintf("%.3f", con$con), 
	    "\\ (", sprintf("%.3f", sqrt(con$var)), ") & ",
	    sprintf("%.3f", kon$c), "\\ (", sprintf("%.3f", sqrt(kon$se)) , ")",
	    ifelse(i == length(models), "\n", "\\\\\n"), sep="")
    }
}
