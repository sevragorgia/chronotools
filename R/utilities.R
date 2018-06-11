
#
# Sergio Vargas
#
#



#extract node ages from a collection of trees
get_phylobayes_dates<-function(tree_set){
  
  
  phylo_dates<-data.frame(Node=c(1:tree_set[[1]]$Nnode))
  
  for(tree in 1:length(tree_set)){
    
    phylo_dates<-cbind(phylo_dates, tree_set[[tree]]$node.label)
    
  }
  
  phylo_dates<-apply(t(phylo_dates[,-1]), 2, function(x)as.numeric(as.character(x)))  
  
}




