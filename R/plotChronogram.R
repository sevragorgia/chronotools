# Sergio Vargas, 2018
#
#
#

#' plotChronogram.ageHistogram
#' 
#' This function allows you to add a histogram to each node to display the posterior distribution of the node age
#' @param tree The tree to plot. Branch lengths are expected to be dates
#' @param age.range A matrix specifying the age range (max, min) of each node in the tree
#' @param age.brackets A matrix specifying the breaks of the age histogram at each node in the tree
#' @param alpha.vector A matrix specifying the density for each age bracket for each node. The density is used to calculate the age.brackt alpha.
#' @param base.alpha A value between 0.0 and 1.0 (defaults to 0.8) defining the highest alpha.
#' @param hist.scaler A user specified number controlling the appearance of the histograms at each node (arbitrarily set to 10).
#' @keywords chronogram, histogram, age distributions
#' plotChronogram.ageHistogram()

plotChronogram.ageHistogram<-function(tree,age.range,age.brackets,alpha.vector,...){
  
  args<-list(...)
  
  if(is.null(args$ftype)) args$ftype<-"i"
  
  fsize<-if(!is.null(args$fsize)) args$fsize else 1
  
  if(is.null(args$direction)) args$direction<-"leftwards"
  
  box.width<-if(!is.null(args$box.width)) args$box.width else 0.5
  
  base.alpha<-if(!is.null(args$base.alpha)) args$base.alpha else 0.8
  
  hist.scaler<-if(!is.null(args$hist.scaler)) args$hist.scaler else 10
  
  if(!is.null(args$cex)){
    cex<-args$cex
    args$cex<-NULL
  } else cex<-1.2
  
  if(!is.null(args$bar.col)){
    bar.col<-args$bar.col
    args$bar.col<-NULL
  } else bar.col<-"blue"
  
  par(mar=c(0,0,0,0))
  
  plot.new()		
  
  th<-max(nodeHeights(tree))
  
  h<-max(th,max(age.range))
  
  if(is.null(args$xlim)){
    m<-min(min(nodeHeights(tree)),min(age.range))
    d<-diff(c(m,h))
    pp<-par("pin")[1]
    
    sw<-fsize*(max(strwidth(tree$tip.label,units="inches")))+
      1.37*fsize*strwidth("W",units="inches")
    alp<-optimize(function(a,d,sw,pp) (a*1.04*d+sw-pp)^2,
                  d=d,sw=sw,pp=pp,
                  interval=c(0,1e6))$minimum
    args$xlim<-if(args$direction=="leftwards") c(h,m-sw/alp) else 
      c(m,h+sw/alp)
  }
  
  
  args$tree<-tree
  args$add<-TRUE
  
  do.call(plotTree,args=args)
  
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  
  for(i in 1:tree$Nnode+Ntip(tree)){
    
    alphas<-round((alpha.vector[[i-Ntip(tree)]]-min(alpha.vector[[i-Ntip(tree)]]))/(max(alpha.vector[[i-Ntip(tree)]])-min(alpha.vector[[i-Ntip(tree)]])),2)
    densities<-alpha.vector[[i-Ntip(tree)]]
    ages<-age.brackets[[i-Ntip(tree)]]
    
    for(j in 1:(length(age.brackets[[i-Ntip(tree)]])-1)){
      
      rect(ages[j], obj$yy[i], ages[j+1], obj$yy[i]+(densities[j]*hist.scaler), border=make.transparent(bar.col, base.alpha), col=make.transparent(bar.col, alphas[j]*base.alpha))
      
    }
  }
}


#' plotChronogram.ageBar
#' 
#' This function allows you to add a heat bar to each node to display the posterior density of the node age
#' @param tree The tree to plot. Branch lengths are expected to be dates
#' @param age.range A matrix specifying the age range (max, min) of each node in the tree
#' @param age.brackets A matrix specifying the breaks of the age histogram at each node in the tree
#' @param alpha.vector A matrix specifying the density for each age bracket for each node. The density is used to calculate the age.brackt alpha.
#' @param base.alpha A value between 0.0 and 1.0 (defaults to 0.8) defining the highest alpha.
#' @param box.width A user specified number controlling the appearance of the histograms at each node (arbitrarily set to 10).
#' @keywords chronogram, histogram, age distributions
#' plotChronogram.ageBar()

plotChronogram.ageBar<-function(tree,age.range,age.brackets,alpha.vector,...){
  
  args<-list(...)
  
  if(is.null(args$ftype)) args$ftype<-"i"
  
  fsize<-if(!is.null(args$fsize)) args$fsize else 1
  
  if(is.null(args$direction)) args$direction<-"leftwards"
  
  base.alpha<-if(!is.null(args$base.alpha)) args$base.alpha else 0.8
  
  box.width<-if(!is.null(args$box.width)) args$box.width else 0.2
  
  if(!is.null(args$cex)){
    cex<-args$cex
    args$cex<-NULL
  } else cex<-1.2
  
  if(!is.null(args$bar.col)){
    bar.col<-args$bar.col
    args$bar.col<-NULL
  } else bar.col<-"blue"
  
  par(mar=c(0,0,0,0))
  
  plot.new()		
  
  th<-max(nodeHeights(tree))
  
  h<-max(th,max(age.range))
  
  if(is.null(args$xlim)){
    m<-min(min(nodeHeights(tree)),min(age.range))
    d<-diff(c(m,h))
    pp<-par("pin")[1]
    
    sw<-fsize*(max(strwidth(tree$tip.label,units="inches")))+
      1.37*fsize*strwidth("W",units="inches")
    alp<-optimize(function(a,d,sw,pp) (a*1.04*d+sw-pp)^2,
                  d=d,sw=sw,pp=pp,
                  interval=c(0,1e6))$minimum
    args$xlim<-if(args$direction=="leftwards") c(h,m-sw/alp) else 
      c(m,h+sw/alp)
  }
  
  
  args$tree<-tree
  args$add<-TRUE
  
  do.call(plotTree,args=args)
  
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  
  for(i in 1:tree$Nnode+Ntip(tree)){
    
    alphas<-round((alpha.vector[[i-Ntip(tree)]]-min(alpha.vector[[i-Ntip(tree)]]))/(max(alpha.vector[[i-Ntip(tree)]])-min(alpha.vector[[i-Ntip(tree)]])),2)
    densities<-alpha.vector[[i-Ntip(tree)]]
    ages<-age.brackets[[i-Ntip(tree)]]
    
    for(j in 1:(length(age.brackets[[i-Ntip(tree)]])-1)){
      
      if(alphas[j] != 0) rect(ages[j], obj$yy[i]-box.width, ages[j+1], obj$yy[i]+box.width, border=make.transparent(bar.col, base.alpha), col=make.transparent(bar.col, alphas[j]*base.alpha))
      
    }
  }
}



#' plotChronogram.quantiles
#' 
#' This function allows you to add a heat bar to each node to display the posterior density of the node age
#' @param tree The tree to plot. Branch lengths are expected to be dates
#' @param age.range A matrix specifying the age range (max, min) of each node in the tree
#' @param age.quantiles A matrix specifying the quantiles to be plot at each node in the tree
#' @param median.ages A vector with the nodes median ages
#' @param base.alpha A value between 0.0 and 1.0 (defaults to 0.8) defining the alpha.
#' @keywords chronogram, histogram, age distributions
#' plotChronogram.quantiles()

plotChronogram.quantiles<-function(tree,age.range,age.quantiles,median.ages,...){
  
  args<-list(...)
  
  if(is.null(args$ftype)) args$ftype<-"i"
  
  fsize<-if(!is.null(args$fsize)) args$fsize else 1
  
  if(is.null(args$direction)) args$direction<-"leftwards"
  
  base.alpha<-if(!is.null(args$base.alpha)) args$base.alpha else 0.8
  
  if(!is.null(args$bar.width)){
    bar.width<-args$bar.width
    args$bar.width<-NULL
  } else bar.width<-11
  
  
  if(!is.null(args$cex)){
    cex<-args$cex
    args$cex<-NULL
  } else cex<-1.2
  
  if(!is.null(args$bar.col)){
    bar.col<-args$bar.col
    args$bar.col<-NULL
  } else bar.col<-"blue"
  
  par(mar=c(0,0,0,0))
  
  plot.new()		
  
  th<-max(nodeHeights(tree))
  
  h<-max(th,max(age.range))
  
  if(is.null(args$xlim)){
    m<-min(min(nodeHeights(tree)),min(age.range))
    d<-diff(c(m,h))
    pp<-par("pin")[1]
    
    sw<-fsize*(max(strwidth(tree$tip.label,units="inches")))+
      1.37*fsize*strwidth("W",units="inches")
    alp<-optimize(function(a,d,sw,pp) (a*1.04*d+sw-pp)^2,
                  d=d,sw=sw,pp=pp,
                  interval=c(0,1e6))$minimum
    args$xlim<-if(args$direction=="leftwards") c(h,m-sw/alp) else 
      c(m,h+sw/alp)
  }
  
  
  args$tree<-tree
  args$add<-TRUE
  
  do.call(plotTree,args=args)
  
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  
  for(i in 1:tree$Nnode+Ntip(tree)){
    
    #plot range
    lines(x=c(age.range[i-Ntip(tree),1],age.range[i-Ntip(tree),2]), y=rep(obj$yy[i],2),lwd=bar.width,lend=0, col=make.transparent(bar.col,base.alpha/2))
    lines(x=c(age.quantiles[i-Ntip(tree),1],age.quantiles[i-Ntip(tree),2]), y=rep(obj$yy[i],2),lwd=bar.width,lend=0, col=make.transparent(bar.col,base.alpha))
    
    #plot median
    points(median.ages[i-Ntip(tree)],obj$yy[i],pch=19, col=make.transparent(bar.col, base.alpha), cex=cex)
    
  }
}


###########################################################################################################################################

my_plotTree.errorbars<-function(tree,age.range,age.brackets,alpha.vector,mean.ages,base.alpha, hist.scale,...){
  
  args<-list(...)
  
  if(is.null(args$ftype)) args$ftype<-"i"
  
  fsize<-if(!is.null(args$fsize)) args$fsize else 1
  
  if(is.null(args$direction)) args$direction<-"leftwards"
  
  box.width<-if(!is.null(args$box.width)) args$box.width else 0.5
  
  if(!is.null(args$cex)){
    cex<-args$cex
    args$cex<-NULL
  } else cex<-1.2
  
  if(!is.null(args$bar.col)){
    bar.col<-args$bar.col
    args$bar.col<-NULL
  } else bar.col<-"blue"
  
  par(mar=c(0,0,0,0))
  
  plot.new()		
  
  th<-max(nodeHeights(tree))
  
  h<-max(th,max(age.range))
  
  if(is.null(args$xlim)){
    m<-min(min(nodeHeights(tree)),min(age.range))
    d<-diff(c(m,h))
    pp<-par("pin")[1]
    
    sw<-fsize*(max(strwidth(tree$tip.label,units="inches")))+
      1.37*fsize*strwidth("W",units="inches")
    alp<-optimize(function(a,d,sw,pp) (a*1.04*d+sw-pp)^2,
                  d=d,sw=sw,pp=pp,
                  interval=c(0,1e6))$minimum
    args$xlim<-if(args$direction=="leftwards") c(h,m-sw/alp) else 
      c(m,h+sw/alp)
  }
  
  
  args$tree<-tree
  args$add<-TRUE
  
  do.call(plotTree,args=args)
  
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  
  for(i in 1:tree$Nnode+Ntip(tree)){
    
    alphas<-round((alpha.vector[[i-Ntip(tree)]]-min(alpha.vector[[i-Ntip(tree)]]))/(max(alpha.vector[[i-Ntip(tree)]])-min(alpha.vector[[i-Ntip(tree)]])),2)
    densities<-alpha.vector[[i-Ntip(tree)]]
    ages<-age.brackets[[i-Ntip(tree)]]
    
    #color.gradient<-sapply(alphas, function(my.alpha) make.transparent(bar.col, my.alpha*base.alpha))
    #gradient.rect(age.range[i-Ntip(tree),1], obj$yy[i]-box.width, age.range[i-Ntip(tree),2], obj$yy[i]+box.width, col=rev(color.gradient),nslices=length(alpha.vector), border=NA)
    #points(mean.ages[i-Ntip(tree)],obj$yy[i],pch=19, col=make.transparent(bar.col, base.alpha), cex=cex)
    
    for(j in 1:(length(age.brackets[[i-Ntip(tree)]])-1)){
      
      #if(alphas[j] != 0) rect(ages[j], obj$yy[i]-box.width, ages[j+1], obj$yy[i]+box.width, border=make.transparent(bar.col, base.alpha), col=make.transparent(bar.col, alphas[j]*base.alpha))
      rect(ages[j], obj$yy[i], ages[j+1], obj$yy[i]+(densities[j]*hist.scale), border=make.transparent(bar.col, base.alpha), col=make.transparent(bar.col, alphas[j]*base.alpha))
      
      
    }
    
  }
  #lines(x=c(CI[i-Ntip(tree),1],CI[i-Ntip(tree),2]),
  #y=rep(obj$yy[i],2),lwd=bar.width,lend=0,
  #      col=make.transparent(bar.col,0.4))
  
  
  
  #points(obj$xx[1:tree$Nnode+Ntip(tree)],
  #       obj$yy[1:tree$Nnode+Ntip(tree)],pch=19,col=bar.col,
  #       cex=cex)
}


