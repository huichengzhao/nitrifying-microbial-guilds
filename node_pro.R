node_pro<-function(igraph,outdir){  
  # 节点度
  igraph.degree<-igraph::degree(igraph)
  # 节点度中心性
  igraph.cen.degree<-centralization.degree(igraph)$res
  # 节点介数中心性
  igraph.betweenness<-centralization.betweenness(igraph)$res
  # 节点中心性
  igraph.closeness<-centralization.closeness(igraph)$res
  
  igraph.node.pro <- cbind(igraph.degree,igraph.closeness,igraph.betweenness,igraph.cen.degree)
  colnames(igraph.node.pro)<-c("igraph.degree","igraph.closeness","igraph.betweenness","igraph.cen.degree")
  igraph.node.pro
}