library(MASS)

d3m <- function(case,control,bsn,qn=101){
  
  wasser <- function(x, y, paranum=11, q=1){
    xq <- quantile(x, probs = seq(0,1,length=paranum),na.rm = T)
    yq <- quantile(y, probs = seq(0,1,length=paranum),na.rm = T)
    d <- sum((xq - yq)^q)
    d
  }
  
  case <- as.numeric(case)
  case <- case[!is.na(case)]
  control <- as.numeric(control)
  control <- control[!is.na(control)]
  d <- wasser(case,control,qn, 2)
  
  ncase <- length(case)
  ncontrol <- length(control)
  nsample <- ncase + ncontrol
  data <- c(case,control)
  
  shuffleID <- sapply(1:bsn,function(j)sample(nsample,nsample,replace=FALSE))
  bootd <- rep(0,bsn)
  for(i in 1:ncol(shuffleID)){
    set.seed(i)
    id <- shuffleID[,i]
    bootd[i] <- wasser(data[id][1:ncase],data[id][(ncase+1):nsample],101,2)
  }
  rm(shuffleID)
  
  pval <- 1 - ecdf(bootd)(d)
  if(pval==0){
    threshold <- quantile(bootd,.995)
    tx <- bootd[bootd >= threshold] - threshold
    fit <- fitdistr(x=tx, densfun="exponential")
    scale <- fit$estimate
    pval_exp <- pexp(d,scale,lower.tail=F)
    pval_threshold <- pexp(threshold,scale,lower.tail=F)
    r <- pval_threshold / .005
    pval <- as.numeric(pval_exp / r)
  }
  pval
}
