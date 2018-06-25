
#displacement = matrix with x,y displacement from node coordinates
labelnodes<-function(text,node=NULL,...){

  displacement<-if(hasArg(displacement)) list(...)$displacement else matrix(rep(c(0,0),length(text)),length(text),2)

  if(hasArg(cex)) cex<-list(...)$cex
  else cex<-1

  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)

  h<-cex*strheight("A")
  w<-cex*strwidth(text)

  if(is.null(node)){
      cat("No node provided. Nothing to do...\n")
  }else{
    for(i in 1:length(text)){
      ii<-node[i]
      text(obj$xx[ii]-displacement[i,1],obj$yy[ii]+displacement[i,2],label=text[i],cex=cex)
    }
  }
}

# new params
#
# show.names
#
#
geo.legend<-function(leg=NULL,colors=NULL,alpha=0.2,...){
  if(hasArg(cex)) cex<-list(...)$cex
  else cex<-par()$cex
  if(hasArg(plot)) plot<-list(...)$plot
  else plot<-TRUE
  if(hasArg(show.lines)) show.lines<-list(...)$show.lines
  else show.lines<-TRUE

  if(hasArg(show.names)) show.names<-list(...)$show.names
  else show.names<-TRUE
  if(hasArg(show.dates)) show.dates<-list(...)$show.dates
  else show.dates<-TRUE


  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(colors)){
    colors<-setNames(c(
      rgb(255,242,127,255,maxColorValue=255),
      rgb(255,230,25,255,maxColorValue=255),
      rgb(253,154,82,255,maxColorValue=255),
      rgb(127,198,78,255,maxColorValue=255),
      rgb(52,178,201,255,maxColorValue=255),
      rgb(129,43,146,255,maxColorValue=255),
      rgb(240,64,40,255,maxColorValue=255),
      rgb(103,165,153,255,maxColorValue=255),
      rgb(203,140,55,255,maxColorValue=255),
      rgb(179,225,182,255,maxColorValue=255),
      rgb(0,146,112,255,maxColorValue=255),
      rgb(127,160,86,255,maxColorValue=255),
      rgb(247,67,112,255,maxColorValue=255)),
      c("Quaternary","Neogene","Paleogene",
        "Cretaceous","Jurassic","Triassic",
        "Permian","Carboniferous","Devonian",
        "Silurian","Ordovician","Cambrian",
        "Precambrian"))
  }
  if(is.null(leg)){
    leg<-rbind(c(2.588,0),
               c(23.03,2.588),
               c(66.0,23.03),
               c(145.0,66.0),
               c(201.3,145.0),
               c(252.17,201.3),
               c(298.9,252.17),
               c(358.9,298.9),
               c(419.2,358.9),
               c(443.8,419.2),
               c(485.4,443.8),
               c(541.0,485.4),
               c(4600,541.0))
    rownames(leg)<-c("Quaternary","Neogene","Paleogene",
                     "Cretaceous","Jurassic","Triassic",
                     "Permian","Carboniferous","Devonian",
                     "Silurian","Ordovician","Cambrian",
                     "Precambrian")
    t.max<-max(obj$xx)
    ii<-which(leg[,2]<=t.max)
    leg<-leg[ii,]
    leg[max(ii),1]<-t.max
  }
  colors<-sapply(colors,make.transparent,alpha=alpha)
  if(plot){
    y<-c(rep(0,2),rep(par()$usr[4],2))
    ylabel<--1/25*obj$Ntip
    for(i in 1:nrow(leg)){
      strh<-strheight(rownames(leg)[i])
      polygon(c(leg[i,1:2],leg[i,2:1]), y,
              col=colors[rownames(leg)[i]],border=NA)
      if(show.lines){
        lines(x=rep(leg[i,1],2),y=c(0,par()$usr[4]),
              lty="dotted",col="grey")
        if(show.dates){#sevra
          text(x=leg[i,1],y=-1,label=paste(leg[i,1]),col="black",cex=cex)
        }
      }
      if(show.names){
        polygon(x=c(leg[i,1],
                  mean(leg[i,])-0.8*cex*get.asp()*strh,
                  mean(leg[i,])-0.8*cex*get.asp()*strh,
                  mean(leg[i,])+0.8*cex*get.asp()*strh,
                  mean(leg[i,])+0.8*cex*get.asp()*strh,
                  leg[i,2]),y=c(0,ylabel,par()$usr[3],
                                par()$usr[3],ylabel,0),
                  col=colors[rownames(leg)[i]],border=NA)
        text(x=mean(leg[i,])+
               if(obj$direction=="leftwards") 0.12*strh else -0.12*strh,
               y=ylabel,labels=rownames(leg)[i],
               srt=90,adj=c(1,0.5),cex=cex)
        if(show.lines){
          lines(x=c(leg[i,1],mean(leg[i,])-0.8*cex*
                    get.asp()*strheight(rownames(leg)[i])),
                    y=c(0,ylabel),lty="dotted",col="grey")
          lines(x=c(leg[i,2],mean(leg[i,])+0.8*cex*
                    get.asp()*strheight(rownames(leg)[i])),
                    y=c(0,ylabel),lty="dotted",col="grey")
          lines(x=rep(mean(leg[i,])-0.8*cex*
                    get.asp()*strheight(rownames(leg)[i]),2),
                    y=c(ylabel,par()$usr[3]),lty="dotted",col="grey")
          lines(x=rep(mean(leg[i,])+0.8*cex*
                    get.asp()*strheight(rownames(leg)[i]),2),
                    y=c(ylabel,par()$usr[3]),lty="dotted",col="grey")
        }
      }
    }
  }
  invisible(list(leg=leg,colors=colors))
}


##################################################################
#
#
# Original code; forked 13.06.2018
#
#
#
#########################################################################################3





## some utility functions
## written by Liam J. Revell 2011, 2012, 2013, 2014, 2015, 2016, 2017

## function to rescale a tree according to an EB model
## written by Liam J. Revell 2017


## function to expand clades in a plot by a given factor
## written by Liam J. Revell 2017
expand.clade<-function(tree,node,factor=5){
	cw<-reorder(tree)
	tips<-setNames(rep(1,Ntip(tree)),cw$tip.label)
	get.tips<-function(node,tree){
    		dd<-getDescendants(tree,node)
    		tree$tip.label[dd[dd<=Ntip(tree)]]
	}
	desc<-unlist(lapply(node,get.tips,tree=cw))
	for(i in 2:Ntip(cw)){
		tips[i]<-tips[i-1]+
			if(names(tips)[i]%in%desc){
				1
			} else if(names(tips)[i-1]%in%desc){
				1
			} else 1/factor
	}
	obj<-list(tree=tree,tips=tips)
	class(obj)<-"expand.clade"
	obj
}

## S3 print method for the object class "expand.clade"
print.expand.clade<-function(x,...){
	cat("An object of class \"expand.clade\" consisting of:\n")
	cat(paste("(1) A phylogenetic tree (x$tree) with",Ntip(x$tree),
		"tips and\n   ",x$tree$Nnode,"internal nodes.\n"))
	cat("(2) A vector (x$tips) containing the desired tip-spacing.\n\n")
}

## S3 plot method for the object class "expand.clade"
plot.expand.clade<-function(x,...){
	args<-list(...)
	args$tree<-x$tree
	args$tips<-x$tips
	if(inherits(args$tree,"simmap")) do.call(plotSimmap,args)
	else do.call(plotTree,args)
}

## draw a box around a clade
## written by Liam J. Revell 2017

cladebox<-function(tree,node,color=NULL,...){
	if(is.null(color)) color<-make.transparent("yellow",0.2)
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	h<-max(nodeHeights(tree))
	parent<-tree$edge[which(tree$edge[,2]==node),1]
	x0<-max(c(obj$xx[node]+obj$xx[parent])/2,obj$xx[node]-0.05*h)
	x1<-obj$x.lim[2]
	dd<-getDescendants(tree,node)
	y0<-min(range(obj$yy[dd]))-0.5
	y1<-max(range(obj$yy[dd]))+0.5
	polygon(c(x0,x1,x1,x0),c(y0,y0,y1,y1),col=color,
		border=0)
}



## generic function to convert an object of class "simmap" to "phylo"
## written by Liam J. Revell 2016
as.phylo.simmap<-function(x,...){
	x$maps<-NULL
	x$mapped.edge<-NULL
	if(!is.null(x$node.states)) x$node.states<-NULL
	if(!is.null(x$states)) x$states<-NULL
	if(!is.null(x$Q)) x$Q<-NULL
	if(!is.null(x$logL)) x$logL<-NULL
	if(!is.null(attr(x,"map.order"))) attr(x,"map.order")<-NULL
	class(x)<-setdiff(class(x),"simmap")
	x
}

## generic function to convert an object of class "multiSimmap" to "multiPhylo"
## written by Liam J. Revell 2016
as.multiPhylo.multiSimmap<-function(x,...){
	obj<-lapply(x,as.phylo)
	class(obj)<-setdiff(class(x),"multiSimmap")
	obj
}

## generic function to convert object of class "phylo" to "multiPhylo"
## written by Liam J. Revell 2016
as.multiPhylo.phylo<-function(x,...){
	obj<-list(x)
	class(obj)<-"multiPhylo"
	obj
}

as.multiPhylo<-function(x,...){
	if (identical(class(x),"multiPhylo")) return(x)
	UseMethod("as.multiPhylo")
}

## function to label clades
## written by Liam J. Revell 2014, 2015
cladelabels<-function(tree=NULL,text,node,offset=NULL,wing.length=NULL,cex=1,
	orientation="vertical"){
	lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(is.null(tree)){
		wing.length<-1
		if(is.null(offset)) offset<-8
		tree<-list(edge=lastPP$edge,
			tip.label=1:lastPP$Ntip,
			Nnode=lastPP$Nnode)
		H<-matrix(lastPP$xx[tree$edge],nrow(tree$edge),2)
		tree$edge.length<-H[,2]-H[,1]
		class(tree)<-"phylo"
	}
	if(is.null(offset)) offset<-0.5
	xx<-mapply(labelSubTree,node,text,
		MoreArgs=list(tree=tree,pp=lastPP,offset=offset,wl=wing.length,cex=cex,
		orientation=orientation))
}

## internal function used by cladelabels
## written by Liam J. Revell 2014, 2015
labelSubTree<-function(tree,nn,label,pp,offset,wl,cex,orientation){
	if(is.null(wl)) wl<-1
	tree<-reorder(tree)
	tips<-getDescendants(tree,nn)
	tips<-tips[tips<=length(tree$tip.label)]
	ec<-0.7 ## expansion constant
	sw<-pp$cex*max(strwidth(tree$tip.label[tips]))
	sh<-pp$cex*max(strheight(tree$tip.label))
	cw<-mean(strwidth(LETTERS)*cex)
	h<-max(sapply(tips,function(x,tree)
		nodeHeights(tree)[which(tree$edge[,2]==x),2],
		tree=tree))+sw+offset*cw
	y<-range(pp$yy[tips])
	lines(c(h,h),y+ec*c(-sh,sh))
	lines(c(h-wl*cw,h),
		c(y[1]-ec*sh,y[1]-ec*sh))
	lines(c(h-wl*cw,h),
		c(y[2]+ec*sh,y[2]+ec*sh))
	text(h+cw,mean(y),
		label,srt=if(orientation=="horizontal") 0 else 90,
		adj=if(orientation=="horizontal") 0 else 0.5,cex=cex)
}

## get all the extant/extinct tip names
## written by Liam J. Revell 2012, 2015

getExtant<-function(tree,tol=1e-8){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	H<-nodeHeights(tree)
	tl<-max(H)
	x<-which(H[,2]>=(tl-tol))
	y<-tree$edge[x,2]
	y<-y[y<=length(tree$tip)]
	z<-tree$tip.label[y]
	return(z)
}

getExtinct<-function(tree,tol=1e-8) setdiff(tree$tip.label,getExtant(tree,tol))


## function to re-root a phylogeny along an edge
## written by Liam J. Revell 2011-2016

reroot<-function(tree,node.number,position=NULL,interactive=FALSE,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(interactive){
		plotTree(tree,...)
		cat("Click where you would like re-root the plotted tree\n")
		flush.console()
		obj<-get.treepos(message=FALSE)
		node.number<-obj$where
		position<-tree$edge.length[which(tree$edge[,2]==node.number)]-obj$pos
	}
	if(is.null(position)) position<-tree$edge.length[which(tree$edge[,2]==node.number)]
	tt<-splitTree(tree,list(node=node.number,bp=position))
	p<-tt[[1]]
	d<-tt[[2]]
	tip<-if(length(which(p$tip.label=="NA"))>0) "NA" else p$tip.label[which(p$tip.label%in%tree$node.label)]
	p<-root(p,outgroup=tip,resolve.root=T)
	bb<-which(p$tip.label==tip)
	p$tip.label[bb]<-"NA"
	ee<-p$edge.length[which(p$edge[,2]==bb)]
	p$edge.length[which(p$edge[,2]==bb)]<-0
	cc<-p$edge[which(p$edge[,2]==bb),1]
	dd<-setdiff(p$edge[which(p$edge[,1]==cc),2],bb)
	p$edge.length[which(p$edge[,2]==dd)]<-p$edge.length[which(p$edge[,2]==dd)]+ee
	obj<-paste.tree(p,d)
	if(interactive) plotTree(obj,...)
	obj
}

## function to add an arrow pointing to a tip or node in the tree
## written by Liam J. Revell 2014, 2017

add.arrow<-function(tree=NULL,tip,...){
	lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(!is.null(tree)){
		if(inherits(tree,"contMap")) tree<-tree$tree
		else if(inherits(tree,"densityMap")) tree<-tree$tree
	}
	if(is.numeric(tip)){
		ii<-tip
		if(!is.null(tree)&&ii<=length(tree$tip.label)) tip<-tree$tip.label[ii]
		else tip<-""
	} else if(is.character(tip)&&!is.null(tree)) ii<-which(tree$tip.label==tip)
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-1
	strw<-lastPP$cex*(strwidth(tip)+offset*mean(strwidth(c(LETTERS,letters))))
	if(hasArg(arrl)) arrl<-list(...)$arrl
	else {
		if(lastPP$type=="fan") arrl<-0.3*max(lastPP$xx)
		else if(lastPP$type=="phylogram") arrl<-0.15*max(lastPP$xx)
	}
	if(hasArg(hedl)) hedl<-list(...)$hedl
	else hedl<-arrl/3
	if(hasArg(angle)) angle<-list(...)$angle
	else angle<-45
	arra<-angle*pi/180
	asp<-if(lastPP$type=="fan") 1 else (par()$usr[4]-par()$usr[3])/(par()$usr[2]-par()$usr[1])
	if(hasArg(col)) col<-list(...)$col
	else col<-"black"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(lastPP$type=="fan") theta<-atan2(lastPP$yy[ii],lastPP$xx[ii])
	else if(lastPP$type=="phylogram") theta<-0
	segments(x0=lastPP$xx[ii]+cos(theta)*(strw+arrl),
		y0=lastPP$yy[ii]+sin(theta)*(strw+arrl),
		x1=lastPP$xx[ii]+cos(theta)*strw,
		y1=lastPP$yy[ii]+sin(theta)*strw,
		col=col,lwd=lwd,lend="round")
	segments(x0=lastPP$xx[ii]+cos(theta)*strw+cos(theta+arra/2)*hedl,
		y0=lastPP$yy[ii]+sin(theta)*strw+sin(theta+arra/2)*hedl*asp,
		x1=lastPP$xx[ii]+cos(theta)*strw,
		y1=lastPP$yy[ii]+sin(theta)*strw,
		col=col,lwd=lwd,lend="round")
	segments(x0=lastPP$xx[ii]+cos(theta)*strw+cos(theta-arra/2)*hedl,
		y0=lastPP$yy[ii]+sin(theta)*strw+sin(theta-arra/2)*hedl*asp,
		x1=lastPP$xx[ii]+cos(theta)*strw,
		y1=lastPP$yy[ii]+sin(theta)*strw,
		col=col,lwd=lwd,lend="round")
	invisible(list(x0=lastPP$xx[ii]+cos(theta)*(strw+arrl),
		y0=lastPP$yy[ii]+sin(theta)*(strw+arrl),
		x1=lastPP$xx[ii]+cos(theta)*strw,
		y1=lastPP$yy[ii]+sin(theta)*strw))
}

## function finds the height of a given node
## written by Liam Revell 2014, 2015, 2016
nodeheight<-function(tree,node,...){
	if(hasArg(root.edge)) root.edge<-list(...)$root.edge
	else root.edge<-FALSE
	if(root.edge) ROOT<-if(!is.null(tree$root.edge)) tree$root.edge else 0
	else ROOT<-0
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(node==(length(tree$tip.label)+1)) h<-0
	else {
		a<-setdiff(c(getAncestors(tree,node),node),length(tree$tip.label)+1)
		h<-sum(tree$edge.length[sapply(a,function(x,e) which(e==x),e=tree$edge[,2])])
	}
	h+ROOT
}

## function gets ancestor node numbers, to be used internally by
## written by Liam J. Revell 2014
getAncestors<-function(tree,node,type=c("all","parent")){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	type<-type[1]
	if(type=="all"){
		aa<-vector()
		rt<-length(tree$tip.label)+1
		currnode<-node
		while(currnode!=rt){
			currnode<-getAncestors(tree,currnode,"parent")
			aa<-c(aa,currnode)
		}
		return(aa)
	} else if(type=="parent"){
		aa<-tree$edge[which(tree$edge[,2]==node),1]
		return(aa)
	} else stop("do not recognize type")
}

# returns the heights of each node
# written by Liam J. Revell 2011, 2012, 2013, 2015, 2016
# modified by Klaus Schliep 2017
nodeHeights<-function(tree,...){
    if(hasArg(root.edge)) root.edge<-list(...)$root.edge
    else root.edge<-FALSE
    if(root.edge) ROOT<-if(!is.null(tree$root.edge)) tree$root.edge else 0
    else ROOT<-0
    nHeight <- function(tree){
        tree <- reorder(tree)
        edge <- tree$edge
        el <- tree$edge.length
        res <- numeric(max(tree$edge))
        for(i in seq_len(nrow(edge))) res[edge[i,2]] <- res[edge[i,1]] + el[i]
        res
    }
    nh <- nHeight(tree)
    return(matrix(nh[tree$edge], ncol=2L)+ROOT)
}

# function gets sister node numbers or names
# written by Liam J. Revell 2013, 2015
getSisters<-function(tree,node,mode=c("number","label")){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	mode<-mode[1]
	if(is.character(node)) node<-match(node,c(tree$tip.label,tree$node.label))
	sisters<-tree$edge[which(tree$edge[,1]==tree$edge[which(tree$edge[,2]==node),1]),2]
	sisters<-setdiff(sisters,node)
	if(mode=="number") return(sisters)
	else if(mode=="label"){
		result<-list()
		n<-length(tree$tip.label)
		if(is.null(tree$node.label)&&any(sisters>n)) result$nodes<-sisters[which(sisters>n)]
		else if(any(sisters>n)) result$nodes<-tree$node.label[sisters[which(sisters>n)]-n]
		if(any(sisters<=n)) result$tips<-tree$tip.label[sisters[which(sisters<=n)]]
		return(result)
	}
}

##required

## borrowed from mapplots
get.asp<-function(){
  pin<-par("pin")
  usr<-par("usr")
  asp<-(pin[2]/(usr[4]-usr[3]))/(pin[1]/(usr[2]-usr[1]))
  asp
}

# function reorders simmap tree
# written Liam Revell 2011, 2013, 2015
reorderSimmap<-function(tree,order="cladewise",index.only=FALSE,...){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  ii<-reorder.phylo(tree,order,index.only=TRUE,...)
  if(!index.only){
    if(inherits(ii,"phylo")) ii<-whichorder(ii$edge[,2],tree$edge[,2]) ## bug workaround
    tree$edge<-tree$edge[ii,]
    tree$edge.length<-tree$edge.length[ii]
    if(!is.null(tree$maps)){
      tree$maps<-tree$maps[ii]
      tree$mapped.edge<-tree$mapped.edge[ii,]
    }
    attr(tree,"order")<-order
    return(tree)
  } else return(ii)
}

## S3 reorder method for objects of class "simmap"
reorder.simmap<-function(x,...) reorderSimmap(x,...)


## function to make a color (e.g., "blue") transparent with alpha level alpha
make.transparent<-function(color,alpha){
  RGB<-col2rgb(color)[,1]/255
  rgb(RGB[1],RGB[2],RGB[3],alpha)
}


# gets descendant node numbers
# written by Liam Revell 2012, 2013, 2014
getDescendants<-function(tree,node,curr=NULL){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  if(length(curr)==0&&node<=length(tree$tip.label)) curr<-node
  w<-which(daughters>length(tree$tip.label))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}


## function to add a geological or other temporal legend to a plotted tree
## written by Liam J. Revell 2017
phytools.geo.legend<-function(leg=NULL,colors=NULL,alpha=0.2,...){
  if(hasArg(cex)) cex<-list(...)$cex
  else cex<-par()$cex
  if(hasArg(plot)) plot<-list(...)$plot
  else plot<-TRUE
  if(hasArg(show.lines)) show.lines<-list(...)$show.lines
  else show.lines<-TRUE
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(colors)){
    colors<-setNames(c(
      rgb(255,242,127,255,maxColorValue=255),
      rgb(255,230,25,255,maxColorValue=255),
      rgb(253,154,82,255,maxColorValue=255),
      rgb(127,198,78,255,maxColorValue=255),
      rgb(52,178,201,255,maxColorValue=255),
      rgb(129,43,146,255,maxColorValue=255),
      rgb(240,64,40,255,maxColorValue=255),
      rgb(103,165,153,255,maxColorValue=255),
      rgb(203,140,55,255,maxColorValue=255),
      rgb(179,225,182,255,maxColorValue=255),
      rgb(0,146,112,255,maxColorValue=255),
      rgb(127,160,86,255,maxColorValue=255),
      rgb(247,67,112,255,maxColorValue=255)),
      c("Quaternary","Neogene","Paleogene",
        "Cretaceous","Jurassic","Triassic",
        "Permian","Carboniferous","Devonian",
        "Silurian","Ordovician","Cambrian",
        "Precambrian"))
  }
  if(is.null(leg)){
    leg<-rbind(c(2.588,0),
               c(23.03,2.588),
               c(66.0,23.03),
               c(145.0,66.0),
               c(201.3,145.0),
               c(252.17,201.3),
               c(298.9,252.17),
               c(358.9,298.9),
               c(419.2,358.9),
               c(443.8,419.2),
               c(485.4,443.8),
               c(541.0,485.4),
               c(4600,541.0))
    rownames(leg)<-c("Quaternary","Neogene","Paleogene",
                     "Cretaceous","Jurassic","Triassic",
                     "Permian","Carboniferous","Devonian",
                     "Silurian","Ordovician","Cambrian",
                     "Precambrian")
    t.max<-max(obj$xx)
    ii<-which(leg[,2]<=t.max)
    leg<-leg[ii,]
    leg[max(ii),1]<-t.max
  }
  colors<-sapply(colors,make.transparent,alpha=alpha)
  if(plot){
    y<-c(rep(0,2),rep(par()$usr[4],2))
    ylabel<--1/25*obj$Ntip
    for(i in 1:nrow(leg)){
      strh<-strheight(rownames(leg)[i])
      polygon(c(leg[i,1:2],leg[i,2:1]),y,
              col=colors[rownames(leg)[i]],border=NA)
      if(show.lines){
        lines(x=rep(leg[i,1],2),y=c(0,par()$usr[4]),
              lty="dotted",col="grey")
        lines(x=c(leg[i,1],mean(leg[i,])-0.8*cex*
                    get.asp()*strheight(rownames(leg)[i])),
              y=c(0,ylabel),lty="dotted",col="grey")
        lines(x=c(leg[i,2],mean(leg[i,])+0.8*cex*
                    get.asp()*strheight(rownames(leg)[i])),
              y=c(0,ylabel),lty="dotted",col="grey")
        lines(x=rep(mean(leg[i,])-0.8*cex*
                      get.asp()*strheight(rownames(leg)[i]),2),
              y=c(ylabel,par()$usr[3]),lty="dotted",col="grey")
        lines(x=rep(mean(leg[i,])+0.8*cex*
                      get.asp()*strheight(rownames(leg)[i]),2),
              y=c(ylabel,par()$usr[3]),lty="dotted",col="grey")
      }
      polygon(x=c(leg[i,1],
                  mean(leg[i,])-0.8*cex*get.asp()*strh,
                  mean(leg[i,])-0.8*cex*get.asp()*strh,
                  mean(leg[i,])+0.8*cex*get.asp()*strh,
                  mean(leg[i,])+0.8*cex*get.asp()*strh,
                  leg[i,2]),y=c(0,ylabel,par()$usr[3],
                                par()$usr[3],ylabel,0),
              col=colors[rownames(leg)[i]],border=NA)
      text(x=mean(leg[i,])+
             if(obj$direction=="leftwards") 0.12*strh else -0.12*strh,
           y=ylabel,labels=rownames(leg)[i],
           srt=90,adj=c(1,0.5),cex=cex)
    }
  }
  invisible(list(leg=leg,colors=colors))
}


## function to return a node index interactively from a plotted tree
## written by Liam J. Revell 2017
getnode<-function(...){
  if(hasArg(env)) env<-list(...)$env
  else env<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(hasArg(show.pt)) show.pt<-list(...)$show.pt
  else show.pt<-FALSE
  xy<-unlist(locator(n=1))
  if(show.pt) points(xy[1],xy[2])
  d<-sqrt((xy[1]-env$xx)^2+(xy[2]-env$yy)^2)
  ii<-which(d==min(d))[1]
  ii
}

## function mostly to interactively label nodes by clicking
## written by Liam J. Revell 2017
phytool_labelnodes<-function(text,node=NULL,interactive=TRUE,
                     shape=c("circle","ellipse","rect"),...){
  shape<-shape[1]
  if(hasArg(circle.exp)) circle.exp<-list(...)$circle.exp
  else circle.exp<-1.3
  if(hasArg(rect.exp)) rect.exp<-list(...)$rect.exp
  else rect.exp<-1.6
  if(hasArg(cex)) cex<-list(...)$cex
  else cex<-1
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  h<-cex*strheight("A")
  w<-cex*strwidth(text)
  rad<-circle.exp*h*diff(par()$usr[1:2])/diff(par()$usr[3:4])
  if(is.null(node)){
    if(!interactive){
      cat("No nodes provided. Setting interactive mode to TRUE.\n")
      interactive<-TRUE
    }
    node<-vector(length=length(text))
  }
  for(i in 1:length(text)){
    if(interactive){
      cat(paste("Click on the node you would like to label ",
                text[i],".\n",sep=""))
      flush.console()
      ii<-getnode(env=obj)
      node[i]<-ii
    } else ii<-node[i]
    if(shape=="circle")
      draw.circle(obj$xx[ii],obj$yy[ii],rad,col="white")
    else if(shape=="ellipse")
      draw.ellipse(obj$xx[ii],obj$yy[ii],0.8*w[i],h,
                   col="white")
    else if(shape=="rect")
      rect(xleft=obj$xx[ii]-0.5*rect.exp*w[i],
           ybottom=obj$yy[ii]-0.5*rect.exp*h,
           xright=obj$xx[ii]+0.5*rect.exp*w[i],
           ytop=obj$yy[ii]+0.5*rect.exp*h,col="white",
           ljoin=1)
    text(obj$xx[ii],obj$yy[ii],label=text[i],cex=cex)
  }
  invisible(node)
}

# function splits tree at split
# written by Liam Revell 2011, 2014, 2015

splitTree<-function(tree,split){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  if(split$node>length(tree$tip.label)){
    # first extract the clade given by shift$node
    tr2<-extract.clade(tree,node=split$node)
    tr2$root.edge<-tree$edge.length[which(tree$edge[,2]==split$node)]-split$bp
    #now remove tips in tr2 from tree
    tr1<-drop.clade(tree,tr2$tip.label)
    nn<-if(!is.null(tree$node.label)) c(tree$node.label,"NA") else "NA"
    tr1$tip.label[which(tr1$tip.label%in%nn)]<-"NA"
    tr1$edge.length[match(which(tr1$tip.label=="NA"),tr1$edge[,2])]<-split$bp
  } else {
    # first extract the clade given by shift$node
    tr2<-list(edge=matrix(c(2L,1L),1,2),tip.label=tree$tip.label[split$node],edge.length=tree$edge.length[which(tree$edge[,2]==split$node)]-split$bp,Nnode=1L)
    class(tr2)<-"phylo"
    # now remove tip in tr2 from tree
    tr1<-tree
    tr1$edge.length[match(which(tr1$tip.label==tr2$tip.label[1]),tr1$edge[,2])]<-split$bp
    tr1$tip.label[which(tr1$tip.label==tr2$tip.label[1])]<-"NA"
  }
  trees<-list(tr1,tr2)
  class(trees)<-"multiPhylo"
  trees
}

# function drops entire clade
# written by Liam Revell 2011, 2015

drop.clade<-function(tree,tip){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  nn<-if(!is.null(tree$node.label)) c(tree$node.label,"NA") else "NA"
  tree<-drop.tip(tree,tip,trim.internal=FALSE)
  while(sum(tree$tip.label%in%nn)>1)
    tree<-drop.tip(tree,tree$tip.label[tree$tip.label%in%nn],
                   trim.internal=FALSE)
  tree
}

# function pastes subtree onto tip
# written by Liam Revell 2011, 2015
paste.tree<-function(tr1,tr2){
  if(!inherits(tr1,"phylo")||!inherits(tr2,"phylo")) stop("tr1 & tr2 should be objects of class \"phylo\".")
  if(length(tr2$tip)>1){
    temp<-tr2$root.edge; tr2$root.edge<-NULL
    tr1$edge.length[match(which(tr1$tip.label=="NA"),tr1$edge[,2])]<-tr1$edge.length[match(which(tr1$tip.label=="NA"),tr1$edge[,2])]+temp
  }
  tr.bound<-bind.tree(tr1,tr2,where=which(tr1$tip.label=="NA"))
  return(tr.bound)
}

# function rotates a node or multiple nodes
# written by Liam J. Revell 2013, 2015
rotateNodes<-function(tree,nodes,polytom=c(1,2),...){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  n<-length(tree$tip.label)
  if(nodes[1]=="all") nodes<-1:tree$Nnode+n
  for(i in 1:length(nodes)) tree<-rotate(tree,nodes[i],polytom)
  if(hasArg(reversible)) reversible<-list(...)$reversible
  else reversible<-TRUE
  if(reversible){
    ii<-which(tree$edge[,2]<=n)
    jj<-tree$edge[ii,2]
    tree$edge[ii,2]<-1:n
    tree$tip.label<-tree$tip.label[jj]
  }
  return(tree)
}



