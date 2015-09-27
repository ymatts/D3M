#' Subsetting Top Significant Distribution Pairs and Output Graphics.
#' @param obj D3M object resulted from d3m.
#' @param id identifier of each row. Default is sequential integer number starting from 1.
#' @param ntop the number of top distribution pairs to be investigated. Default is 10.
#' @param qval if TRUE, ordering the result by false discovery rate.
#' @param plot.it whether output the grphics based on Q-Q plot by JPEG format. Defarult is plot.it = FALSE.
#' @return top.cases top significant distributions in case group.
#' @return top.control top significant distributions in control group.
#' @return pval p-values related to top significant distribution pairs.
#' @return Q-Q plots of distribution pairs.
#' @details this function extract a subset of the most significant distribution pairs based on p-values. The graphical representation is based on Q-Q plots, which represents the shape difference being tested by d3m function.
#' @author Yusuke Matsui & Teppei Shimamura
#' @examples
#' library(D3M)
#' nrep <- 12
#' cases <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep)); cases <- do.call("rbind",cases)
#' control <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep)); control <- do.call("rbind",control)
#' obj <- d3m(cases,control,paranum = 101, q = 2, bsn = 1000)
#' topD3M(obj,ntop = 10,qval = TRUE, plot.it = FALSE)
#' @export
#' 

topD3M <- function(obj, id = seq_along(obj[[1]]), ntop = 10, qval = TRUE, plot.it = FALSE){
  
  if(qval){
      pvals <- p.adjust(obj[[1]],method = "BH")
  }else{
      pvals <- obj[[1]]
  }
  
  obj[[5]] <- id
  
  order.pvals <- order(pvals,decreasing = F)
  
  #topcases <- obj$cases[1:ntop,]
  
  topcases <- obj[[3]][order.pvals[1:ntop],]
  
  topcontrol <- obj[[4]][order.pvals[1:ntop],]
  
  TargetID <- obj[[5]][order.pvals[1:ntop]]
  
  pval <- pvals[order.pvals[1:ntop]]
  
  para <- par()
  
  if(plot.it){
      
      
      grDevices::jpeg("topD3M.jpg",800,600)
      
      #qq <- apply(topcases,1,stats::qqplot,topcontrol,plot.it=F)
      qx <- apply(topcases,1,function(x)stats::quantile(x,probs=seq(0,1,.05),type = 4))
      qy <- apply(topcontrol,1,function(x)stats::quantile(x,probs=seq(0,1,.05),type = 4))
      
      #graphics::plot(qq[[1]],xlim = c(0,1),ylim = c(0,1),type = "n",ann=F)
      graphics::plot(qx[,1],qy[,1],xlim = c(0,1),ylim = c(0,1),type = "n",ann=F)
      
      for(i in 1:ntop){
          
          #graphics::lines(qq[[i]],lwd = 1,col = grDevices::rgb(red = 0,0,1,alpha = 0.2))
          graphics::lines(qx[,i],qy[,i],lwd = 1,col = grDevices::rgb(red = 0,0,1,alpha = 0.2))
          
          #graphics::points(qq[[i]],pch = 2,col = grDevices::rgb(red = 0,0,1,alpha = 0.2))
          graphics::points(qx[,i],qy[,i],pch = 2,col = grDevices::rgb(red = 0,0,1,alpha = 0.2))
          
      }
      
      graphics::abline(a = 0,b = 1,col="red")
      
      graphics::title(main=paste("Q-Q plot of top ", ntop, " significant sites",sep = ""))
      
      graphics::mtext(text = "cases",side = 1,line = 3)
      
      graphics::mtext(text = "control",side = 2,line = 3)
  }
  
  on.exit(par(para))
  #on.exit(dev.off())
  
  return(list(top.cases = topcases, top.contorl = topcontrol,pval = pval,TargetID = TargetID))

}

