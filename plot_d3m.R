plot.d3m <- function(case,control){
  library(beanplot)
  par(mar=c(2,2,2,1))
  par(mfrow = c(1,2))
  case <- as.numeric(case)
  control <- as.numeric(control)
  beanplot(case,control,side="both",col = list("black",c("grey","white")),ylim=c(0,1),ll=0.05,names = c("case","control"),horizontal = T,log = "",main="Density of beta value")
  grid()
  qqplot(case,control,xlim=c(0,1),ylim=c(0,1),main="Q-Q plot")
  abline(a = 0,b = 1,col="red")
}
