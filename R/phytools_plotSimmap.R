## skimmed version to serve as chronotools companion and allow customization

plotPhylogram<-function(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,
                                   add,offset,direction,setEnv,xlim,ylim,placement,tips,split.vertical,lend,
                                   asp,plot){
  if(split.vertical&&!setEnv){
    cat("split.vertical requires setEnv=TRUE. Setting split.vertical to FALSE.\n")
    spit.vertical<-FALSE
  }
  # set offset fudge (empirically determined)
  offsetFudge<-1.37
  # reorder
  cw<-reorderSimmap(tree)
  pw<-reorderSimmap(tree,"postorder")
  # count nodes and tips
  n<-Ntip(cw)
  m<-cw$Nnode
  # Y coordinates for nodes
  Y<-matrix(NA,m+n,1)
  # first, assign y coordinates to all the tip nodes
  if(is.null(tips)) Y[cw$edge[cw$edge[,2]<=n,2]]<-1:n
  else Y[cw$edge[cw$edge[,2]<=n,2]]<-if(is.null(names(tips)))
    tips[sapply(1:Ntip(cw),function(x,y) which(y==x),y=cw$edge[cw$edge[,2]<=n,2])]
  else tips[gsub(" ","_",cw$tip.label)]
  # get Y coordinates of the nodes
  nodes<-unique(pw$edge[,1])
  for(i in 1:m){
    if(placement=="intermediate"){
      desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
      Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
    } else if(placement=="centered"){
      desc<-getDescendants(pw,nodes[i])
      desc<-desc[desc<=Ntip(pw)]
      Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
    } else if(placement=="weighted"){
      desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
      n1<-desc[which(Y[desc]==min(Y[desc]))]
      n2<-desc[which(Y[desc]==max(Y[desc]))]
      v1<-pw$edge.length[which(pw$edge[,2]==n1)]
      v2<-pw$edge.length[which(pw$edge[,2]==n2)]
      Y[nodes[i]]<-((1/v1)*Y[n1]+(1/v2)*Y[n2])/(1/v1+1/v2)
    } else if(placement=="inner"){
      desc<-getDescendants(pw,nodes[i])
      desc<-desc[desc<=Ntip(pw)]
      mm<-which(abs(Y[desc]-median(Y[1:Ntip(pw)]))==min(abs(Y[desc]-
                                                              median(Y[1:Ntip(pw)]))))
      if(length(mm>1)) mm<-mm[which(Y[desc][mm]==min(Y[desc][mm]))]
      Y[nodes[i]]<-Y[desc][mm]
    }
  }
  # compute node heights
  H<-nodeHeights(cw)
  # open plot
  par(mar=mar)
  if(is.null(offset)) offset<-0.2*lwd/3+0.2/3
  if(!add) plot.new()
  ###
  if(is.null(xlim)){
    pp<-par("pin")[1]
    sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
      offsetFudge*fsize*strwidth("W",units="inches")
    alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=H,sw=sw,pp=pp,
                  interval=c(0,1e6))$minimum
    xlim<-if(direction=="leftwards") c(min(H)-sw/alp,max(H)) else c(min(H),max(H)+sw/alp)
  }
  if(is.null(ylim)) ylim=range(Y)
  if(direction=="leftwards") H<-max(H)-H
  plot.window(xlim=xlim,ylim=ylim,asp=asp)
  if(plot){
    ####
    if(!split.vertical){#sevra: this plots the vertical lines in the phylogram
      for(i in 1:m) lines(H[which(cw$edge[,1]==nodes[i]),1],
                          Y[cw$edge[which(cw$edge[,1]==nodes[i]),2]],col=colors[match(nodes[i],cw$edge[,1])],lwd=lwd)
                          #sevra: original
                          #col=colors[names(cw$maps[[match(nodes[i],cw$edge[,1])]])[1]],lwd=lwd)
    }
    for(i in 1:nrow(cw$edge)){
      x<-H[i,1]
      for(j in 1:length(cw$maps[[i]])){
        if(direction=="leftwards")
          lines(c(x,x-cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),col=colors[i],lwd=lwd,lend=lend)
        #sevra: original
        #col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
        else lines(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),col=colors[i],lwd=lwd,lend=lend)
        #sevra: original
        #col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
        if(pts) points(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
                       pch=20,lwd=(lwd-1))
        x<-x+if(direction=="leftwards") -cw$maps[[i]][j] else cw$maps[[i]][j]
        j<-j+1
      }
    }
    if(node.numbers){
      symbols(if(direction=="leftwards") max(H) else 0,
              mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),
              rectangles=matrix(c(1.2*fsize*strwidth(as.character(Ntip(cw)+1)),
                                  1.4*fsize*strheight(as.character(Ntip(cw)+1))),1,2),inches=FALSE,
              bg="white",add=TRUE)
      text(if(direction=="leftwards") max(H) else 0,
           mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),Ntip(cw)+1,
           cex=fsize)
      for(i in 1:nrow(cw$edge)){
        x<-H[i,2]
        if(cw$edge[i,2]>Ntip(cw)){
          symbols(x,Y[cw$edge[i,2]],
                  rectangles=matrix(c(1.2*fsize*strwidth(as.character(cw$edge[i,2])),
                                      1.4*fsize*strheight(as.character(cw$edge[i,2]))),1,2),inches=FALSE,
                  bg="white",add=TRUE)
          text(x,Y[cw$edge[i,2]],cw$edge[i,2],cex=fsize)
        }
      }
    }
    if(direction=="leftwards") pos<-if(par()$usr[1]>par()$usr[2]) 4 else 2
    if(direction=="rightwards") pos<-if(par()$usr[1]>par()$usr[2]) 2 else 4
    for(i in 1:n) if(ftype) text(H[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],pos=pos,
                                 offset=offset,cex=fsize,font=ftype)
  }
  if(setEnv){
    PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
             show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
             font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
             x.lim=xlim,y.lim=ylim,
             direction=direction,tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
             edge=cw$edge,xx=sapply(1:(Ntip(cw)+cw$Nnode),
                                    function(x,y,z) y[match(x,z)],y=H,z=cw$edge),yy=Y[,1])
    assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
  }
  if(plot) if(split.vertical) splitEdgeColor(cw,colors,lwd)
}





########################################################################################
#
#
#
# Original source forked 13.06.2018
#
#
############################################

## functions plot stochastic character mapped trees
## written by Liam Revell 2011-2017

plotSimmap<-function(tree,colors=NULL,fsize=1.0,ftype="reg",lwd=2,
	pts=FALSE,node.numbers=FALSE,mar=NULL,add=FALSE,offset=NULL,direction="rightwards",
	type="phylogram",setEnv=TRUE,part=1.0,xlim=NULL,ylim=NULL,nodes="intermediate",
	tips=NULL,maxY=NULL,hold=TRUE,split.vertical=FALSE,lend=2,asp=NA,plot=TRUE){
	if(inherits(tree,"multiPhylo")){
		par(ask=TRUE)
		for(i in 1:length(tree)) plotSimmap(tree[[i]],colors=colors,fsize=fsize,ftype=ftype,
			lwd=lwd,pts=pts,node.numbers=node.numbers,mar,add,offset,direction,type,
			setEnv,part,xlim,ylim,nodes,tips,maxY,hold,split.vertical,lend,asp,plot)
	} else {
		# check tree
		if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\"")
		if(is.null(tree$maps)) stop("tree should contain mapped states on edges.")
		# check font
		ftype<-which(c("off","reg","b","i","bi")==ftype)-1
		if(!ftype) fsize=0
		# check colors
		if(is.null(colors)){
			st<-sort(unique(unlist(sapply(tree$maps,names))))
			colors<-palette()[1:length(st)]
			names(colors)<-st
			if(length(st)>1){
				cat("no colors provided. using the following legend:\n")
				print(colors)
			}
		}
		# swap out "_" character for spaces (assumes _ is a place holder)
		tree$tip.label<-gsub("_"," ",tree$tip.label)
		# get margin
		if(is.null(mar)) mar=rep(0.1,4)
		if(hold) null<-dev.hold()
		if(type=="phylogram"){
			if(direction%in%c("upwards","downwards"))
				updownPhylogram(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,add,offset,
					direction,setEnv,xlim,ylim,nodes,tips,split.vertical,lend,asp,plot)
			else plotPhylogram(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,add,offset,
					direction,setEnv,xlim,ylim,nodes,tips,split.vertical,lend,asp,plot)
		} else if(type=="fan"){
			plotFan(tree,colors,fsize,ftype,lwd,mar,add,part,setEnv,xlim,ylim,tips,
				maxY,lend,plot)
		} else if(type=="cladogram"){
			plotCladogram(tree,colors,fsize,ftype,lwd,mar,add,offset,direction,xlim,ylim,
				nodes,tips,lend,asp,plot)
		}
		if(hold) null<-dev.flush()
	}
}

# function to plot simmap tree in type "phylogram"
# written by Liam J. Revell 2011-2017
phylotools_plotPhylogram<-function(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,
	add,offset,direction,setEnv,xlim,ylim,placement,tips,split.vertical,lend,
	asp,plot){
	if(split.vertical&&!setEnv){
		cat("split.vertical requires setEnv=TRUE. Setting split.vertical to FALSE.\n")
		spit.vertical<-FALSE
	}
	# set offset fudge (empirically determined)
	offsetFudge<-1.37
	# reorder
	cw<-reorderSimmap(tree)
	pw<-reorderSimmap(tree,"postorder")
	# count nodes and tips
	n<-Ntip(cw)
 	m<-cw$Nnode
	# Y coordinates for nodes
	Y<-matrix(NA,m+n,1)
	# first, assign y coordinates to all the tip nodes
	if(is.null(tips)) Y[cw$edge[cw$edge[,2]<=n,2]]<-1:n
	else Y[cw$edge[cw$edge[,2]<=n,2]]<-if(is.null(names(tips)))
		tips[sapply(1:Ntip(cw),function(x,y) which(y==x),y=cw$edge[cw$edge[,2]<=n,2])]
		else tips[gsub(" ","_",cw$tip.label)]
	# get Y coordinates of the nodes
	nodes<-unique(pw$edge[,1])
	for(i in 1:m){
		if(placement=="intermediate"){
			desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
			Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
		} else if(placement=="centered"){
			desc<-getDescendants(pw,nodes[i])
			desc<-desc[desc<=Ntip(pw)]
			Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
		} else if(placement=="weighted"){
			desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
			n1<-desc[which(Y[desc]==min(Y[desc]))]
			n2<-desc[which(Y[desc]==max(Y[desc]))]
			v1<-pw$edge.length[which(pw$edge[,2]==n1)]
			v2<-pw$edge.length[which(pw$edge[,2]==n2)]
			Y[nodes[i]]<-((1/v1)*Y[n1]+(1/v2)*Y[n2])/(1/v1+1/v2)
		} else if(placement=="inner"){
			desc<-getDescendants(pw,nodes[i])
			desc<-desc[desc<=Ntip(pw)]
			mm<-which(abs(Y[desc]-median(Y[1:Ntip(pw)]))==min(abs(Y[desc]-
				median(Y[1:Ntip(pw)]))))
			if(length(mm>1)) mm<-mm[which(Y[desc][mm]==min(Y[desc][mm]))]
			Y[nodes[i]]<-Y[desc][mm]
		}
	}
	# compute node heights
	H<-nodeHeights(cw)
	# open plot
	par(mar=mar)
	if(is.null(offset)) offset<-0.2*lwd/3+0.2/3
	if(!add) plot.new()
	###
	if(is.null(xlim)){
		pp<-par("pin")[1]
		sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
			offsetFudge*fsize*strwidth("W",units="inches")
		alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=H,sw=sw,pp=pp,
			interval=c(0,1e6))$minimum
		xlim<-if(direction=="leftwards") c(min(H)-sw/alp,max(H)) else c(min(H),max(H)+sw/alp)
	}
	if(is.null(ylim)) ylim=range(Y)
	if(direction=="leftwards") H<-max(H)-H
	plot.window(xlim=xlim,ylim=ylim,asp=asp)
	if(plot){
		####
		if(!split.vertical){
			for(i in 1:m) lines(H[which(cw$edge[,1]==nodes[i]),1],
				Y[cw$edge[which(cw$edge[,1]==nodes[i]),2]],col=colors[names(cw$maps[[match(nodes[i],
				cw$edge[,1])]])[1]],lwd=lwd)
		}
		for(i in 1:nrow(cw$edge)){
			x<-H[i,1]
			for(j in 1:length(cw$maps[[i]])){
				if(direction=="leftwards")
					lines(c(x,x-cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
						col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
				else lines(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
						col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
				if(pts) points(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
					pch=20,lwd=(lwd-1))
				x<-x+if(direction=="leftwards") -cw$maps[[i]][j] else cw$maps[[i]][j]
				j<-j+1
			}
		}
		if(node.numbers){
			symbols(if(direction=="leftwards") max(H) else 0,
				mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),
				rectangles=matrix(c(1.2*fsize*strwidth(as.character(Ntip(cw)+1)),
				1.4*fsize*strheight(as.character(Ntip(cw)+1))),1,2),inches=FALSE,
				bg="white",add=TRUE)
			text(if(direction=="leftwards") max(H) else 0,
				mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),Ntip(cw)+1,
				cex=fsize)
			for(i in 1:nrow(cw$edge)){
				x<-H[i,2]
				if(cw$edge[i,2]>Ntip(cw)){
					symbols(x,Y[cw$edge[i,2]],
						rectangles=matrix(c(1.2*fsize*strwidth(as.character(cw$edge[i,2])),
						1.4*fsize*strheight(as.character(cw$edge[i,2]))),1,2),inches=FALSE,
						bg="white",add=TRUE)
					text(x,Y[cw$edge[i,2]],cw$edge[i,2],cex=fsize)
				}
			}
		}
		if(direction=="leftwards") pos<-if(par()$usr[1]>par()$usr[2]) 4 else 2
		if(direction=="rightwards") pos<-if(par()$usr[1]>par()$usr[2]) 2 else 4
		for(i in 1:n) if(ftype) text(H[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],pos=pos,
			offset=offset,cex=fsize,font=ftype)
	}
	if(setEnv){
		PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
			show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
			font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
			x.lim=xlim,y.lim=ylim,
			direction=direction,tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
			edge=cw$edge,xx=sapply(1:(Ntip(cw)+cw$Nnode),
			function(x,y,z) y[match(x,z)],y=H,z=cw$edge),yy=Y[,1])
		assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
	}
	if(plot) if(split.vertical) splitEdgeColor(cw,colors,lwd)
}

# function plots a tree; in the new version this is just a wrapper for plotSimmap
# written by Liam Revell 2012-2017
plotTree<-function(tree,...){
	if(hasArg(color)) color<-list(...)$color
	else color<-NULL
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-1.0
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-"reg"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(pts)) pts<-list(...)$pts
	else pts<-FALSE
	if(hasArg(node.numbers)) node.numbers<-list(...)$node.numbers
	else node.numbers<-FALSE
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-NULL
	if(hasArg(add)) add<-list(...)$add
	else add<-FALSE
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-NULL
	if(hasArg(type)) type<-list(...)$type
	else type<-"phylogram"
	if(hasArg(direction)) direction<-list(...)$direction
	else direction<-"rightwards"
	if(hasArg(setEnv)) setEnv<-list(...)$setEnv
	else setEnv<-TRUE
	if(hasArg(part)) part<-list(...)$part
	else part<-1.0
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-NULL
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-NULL
	if(hasArg(nodes)) nodes<-list(...)$nodes
	else nodes<-"intermediate"
	if(hasArg(tips)) tips<-list(...)$tips
	else tips<-NULL
	if(hasArg(maxY)) maxY<-list(...)$maxY
	else maxY<-NULL
	if(hasArg(hold)) hold<-list(...)$hold
	else hold<-TRUE
	if(hasArg(lend)) lend<-list(...)$lend
	else lend<-2
	if(hasArg(asp)) asp<-list(...)$asp
	else asp<-NA
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	if(inherits(tree,"multiPhylo")){
		par(ask=TRUE)
		if(!is.null(color)) names(color)<-"1"
		for(i in 1:length(tree)) plotTree(tree[[i]],color=color,fsize=fsize,ftype=ftype,
			lwd=lwd,pts=pts,node.numbers=node.numbers,mar=mar,add=add,offset=offset,
			direction=direction,type=type,setEnv=setEnv,part=part,xlim=xlim,ylim=ylim,
			nodes=nodes,tips=tips,maxY=maxY,hold=hold,lend=lend,asp=asp,plot=plot)
	} else {
		if(is.null(tree$edge.length)) tree<-compute.brlen(tree)
		tree$maps<-as.list(tree$edge.length)
		for(i in 1:length(tree$maps)) names(tree$maps[[i]])<-c("1")
		if(!is.null(color)) names(color)<-"1"
		plotSimmap(tree,colors=color,fsize=fsize,ftype=ftype,lwd=lwd,pts=pts,
			node.numbers=node.numbers,mar=mar,add=add,offset=offset,direction=direction,
			type=type,setEnv=setEnv,part=part,xlim=xlim,ylim=ylim,nodes=nodes,tips=tips,maxY=maxY,
			hold=hold,lend=lend,asp=asp,plot=plot)
	}
}

