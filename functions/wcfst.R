#########################         wc.calc        #########################       
#
# SOURCED FROM https://github.com/kjgilbert/QstFstComp/blob/master/R/QstFstComp.R
# returns the a, b, and c values required to compute Fst according to Weir and Cockerham 1984
#	modified from code previously implemented in hierfstat package 
#	(Jerome Goudet, http://cran.r-project.org/web/packages/hierfstat/index.html)
#
#	ndat		data frame with first column indicating population of origin and following representing loci
#	diploid 	Whether data are diploid

wc.calc <- function(ndat, diploid = TRUE){
  if (!diploid){
    dum <- ndat[, -1]
    nd <- max(dum, na.rm = TRUE)
    modu <- 1000
    if(nd < 10)		modu <- 10
    if(nd < 100)	modu <- 100
    dum <- dum * modu + dum
    ndat <- data.frame(ndat[, 1], dum)
  }
  pop <- ndat[, 1]
  ni <- length(pop)
  dat <- ndat
  loc.names <- names(dat)[-1]
  n <- t(ind.count(dat)) ## NEED in.count to give number of inds genotyped per locus and per pop
  nt <- apply(n, 1, sum, na.rm = TRUE)
  untyped.loc <- which(nt == 0)
  typed.loc <- which(nt != 0)
  if(length(untyped.loc) > 0){
    dat <- dat[, -(untyped.loc + 1)]
    n <- t(ind.count(dat))
    nt <- apply(n, 1, sum, na.rm = TRUE)
  }
  alploc <- nb.alleles(cbind(rep(1, ni), dat[, -1]))
  np <- dim(n)[2]
  npl <- apply(n, 1, tempfun <- function(x) sum(!is.na(x)))
  nl <- dim(n)[1]
  p <- pop.freq(dat, diploid)
  pb <- pop.freq(cbind(rep(1, length(pop)), dat[, -1]), diploid)
  n <- matrix(unlist(n), ncol = np)
  nal <- n[rep(1:nl, alploc), ]
  nc <- (nt - apply(n^2, 1, sum, na.rm = TRUE)/nt)/(npl - 1)
  ntal <- rep(nt, alploc)
  ncal <- rep(nc, alploc)
  p <- matrix(unlist(lapply(p, t)), ncol = np, byrow = TRUE)
  pb <- matrix(unlist(pb), ncol = 1)
  if(diploid){
    dum <- getal.b(dat[, -1])
    all.loc <- apply(dum, 2, tempfun1 <- function(y) as.numeric(dimnames(table(y))[[1]]))
    hetpl <- apply(dum, 2, fun <- function(z) {
      lapply(as.numeric(dimnames(table(z))[[1]]), who.is.het <- function(y) apply(z == 
                                                                                    y, 1, ind.is.het <- function(x) xor(x[1], x[2])))
    })
    mho <- lapply(hetpl, tempfun2 <- function(x) matrix(unlist(lapply(x, 
                                                                      tempfun3 <- function(y) tapply(y, pop, sum, na.rm = TRUE))), 
                                                        ncol = np))
    mho <- matrix(unlist(mho), ncol = np, byrow = TRUE)
    mhom <- (2 * nal * p - mho)/2
  }else{mhom <- nal * p}
  SSG <- apply(nal * p - mhom, 1, sum, na.rm = TRUE)
  dum <- nal * (p - 2 * p^2) + mhom
  SSi <- apply(dum, 1, sum, na.rm = TRUE)
  dum1 <- nal * (sweep(p, 1, pb))^2
  SSP <- 2 * apply(dum1, 1, sum, na.rm = TRUE)
  ntalb <- rep(npl, alploc)
  MSG <- SSG/ntal
  MSP <- SSP/(ntalb - 1)
  MSI <- SSi/(ntal - ntalb)
  sigw <- MSG
  sigb <- 0.5 * (MSI - MSG)
  siga <- 1/2/ncal * (MSP - MSI)
  
  abc.mat <- cbind(siga, sigb, sigw)
  # this returns a matrix of a, b, and c values in columns with one allele per row
  return(abc.mat) 
}