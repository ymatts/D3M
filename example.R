source("d3m.R")
source("plot_d3m.R")
data <- read.table("sample.txt")
case <- data[,1:160]
control <- data[,161:300]
pval <- rep(0,nrow(data))
for(i in 1:nrow(data)){
  pval[i] <- d3m(case = case[i,],control = control[i,],bsn = 2000,qn = 101)
  print(i)
}
plot.d3m(case[2,],control[2,])
