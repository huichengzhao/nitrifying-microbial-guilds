site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 2.1 安装CRAN来源常用包
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list <- c("reshape2","ggplot2","devtools","bindrcpp", "VennDiagram",
                  "ggthemes","agricolae","dplyr","igraph", "psych","sqldf")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
package_list <- c("digest","AnnotationDbi", "impute", "GO.db", "preprocessCore","WGCNA","multtest")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 加载otu矩阵转化为igraph对象函数
source("matrix2igraph.R")
# 加载网络性质的计算函数
source("net_pro.R")
# 加载节点性质的计算函数
source("node_pro.R")

# 存储otu-sample矩阵的文件名,列为sample,行为otu或者gene
otu_sample_file <- "data/otutab_sort.txt"
# 存储otu分类信息的文件名
otu_tax_file<-"data/taxonomy.txt"
# 实验设计
# design_file<-"data/design.txt"
# 设定构建网络的r和p阈值
r.threshold=0.65
p.threshold=0.01
# 调整vertex大小
size=3
# 设定节点颜色所采用的信息列，如taxonomy.txt第二列Phylum
gcol=2
# 设定节点颜色所显示分类群的最大个数
maxtaxnum=10
# 设定节点标签所采用的信息列，如taxonomy.txt第二列Phylum
glab=2

#读取otu矩阵文件
otu <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 

# 统计矩阵中每行为0的元素个数
n0 <- apply(otu == 0, 1, sum)
# 统计0元素大于2个的行号
i0 <- which(n0 > 18)
# 删除0元素大于2个的行
otu=otu[-i0, ]


# 转置
otu<-t(otu)
#otu筛选,可选步骤，如果otu数量太多，建议筛选按丰度筛选
# otu <- otu[,colSums(otu)/sum(otu)>=(0.001/100)]

#读取otu分类信息和实验设计文件 
otu_tax <- read.table(otu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
# design <- read.table(design_file, header=T, row.names= 1, sep="\t")
# 物种注释+丰度
otu_tax<-as.data.frame(otu_tax[colnames(otu),])
otu_abundance <- colSums(otu)
otu_pro <- cbind(otu_tax,otu_abundance)

# matrix2igraph函数提供3个参数，分别是otu矩阵，相关性r的阈值，相关性p的阈值，
igraph<-matrix2igraph(otu,r.threshold,p.threshold)

# 将igraph weight属性赋值到igraph.weight,用于后边做图
igraph.weight <- E(igraph)$weight

# 做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
E(igraph)$weight <- NA
# igraph<-remove.edge.attribute(igraph,"weight")#把边值删除

# 按相关类型设置边颜色
# 如果构建网络时，weighted=NULL,此步骤不能统计
sum(igraph.weight>0)# number of postive correlation
sum(igraph.weight<0)# number of negative correlation

# set edge color，postive correlation 设定为red, negative correlation设定为blue
E.color <- igraph.weight
E.color <- ifelse(E.color>0, "grey",ifelse(E.color<0, "grey","grey"))
E(igraph)$color <- as.character(E.color)

# 添加OTU注释信息，如分类单元和丰度
# 另外可以设置vertices size, vertices color来表征更多维度的数据
# 根据otu的分类地位、多度等性质设置顶点颜色、大小、线粗细等
# set vertices size
igraph.otu <- as.data.frame(otu_pro[V(igraph)$name,]) # 筛选对应OTU属性
igraph.size <- log(as.numeric(igraph.otu$otu_abundance),1000)*size
V(igraph)$size <- igraph.size

# 设置节点颜色
node.col<-igraph.otu[,gcol]
node.col<-as.character(node.col)
fre.tax <- names(sort(table(node.col),decreasing =T)[1:maxtaxnum])
node.col[!(node.col %in% fre.tax)]<-"Rare_groups"
node.col<-as.factor(node.col)
otu.tax.levels<-levels(node.col)

library(randomcoloR)
levels(node.col) <- distinctColorPalette(length(otu.tax.levels))


## levels(node.col) <- rainbow(length(otu.tax.levels)) #   直接修改levles可以连值全部对应替换
V(igraph)$color <- as.character(node.col)

# 创建存放结果文件夹
# dir.create("network_results")
# pdf(file = "network_results/co-occurrence_network_tax.pdf",width = 13)
# # 按物种分类信息为节点配色
# set.seed(123)
# plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
#      edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# legend(1.2,1.4, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
# legend(1.2,0.8, legend=c("pos.cor","neg.cor"),col=c("red","blue"),lty=1,lwd=2, bty="n")
# 
# dev.off()

# 可以添加label项可显示标签，vertex.label.cex选项可以改变标签大小；
pdf(file = "network_results/co-occurrence_network_withtaxlab.pdf",width = 13)
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")

dev.off()

# plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_with_kk,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")


# legend(1.2,0.8, legend=c("pos.cor","neg.cor"),col=c("red","blue"),lty=1,lwd=2, bty="n")

# 根据自己需要修改节点标签
# otu_lab<-as.data.frame(otu_tax[colnames(otu),])
# otu_lab <- as.data.frame(otu_lab[V(igraph)$name,]) # 筛选对应OTU属性
# node.label<-as.character(otu_lab[,glab])
# V(igraph)$label<-node.label


# set.seed(123)
# plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# legend(1.2,1.4, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
# legend(1.2,0.8, legend=c("pos.cor","neg.cor"),col=c("red","blue"),lty=1,lwd=2, bty="n")

dev.off()

# 模块性 modularity
fc <- cluster_fast_greedy(igraph,weights =NULL)# cluster_walktrap   cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
modularity <- modularity(igraph,membership(fc))
# 按照模块为节点配色
comps <- membership(fc)
colbar <- distinctColorPalette(max(comps))
# colbar <- rainbow(max(comps))
V(igraph)$color <- colbar[comps] 

pdf(file = "network_results/co-occurrence_network_modules.pdf")
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.width=0.01,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1,1,legend=c(paste0("Modules",c(1:max(comps)))),col=colbar,pch=16,cex=1.2,bty="n")
dev.off()

# pdf(file = "network_results/co-occurrence_network_layout.pdf",width = 13)
# 
# set.seed(123)
# plot(igraph,main="Co-occurrence network layout_with_kk",layout=layout_with_kk,vertex.frame.color=NA,vertex.label=NA,
#      edge.lty=1,margin=c(0,0,0,0))
# 
# set.seed(123)
# plot(igraph,main="Co-occurrence network layout.fruchterman.reingold",layout=layout.fruchterman.reingold,vertex.frame.color=NA,vertex.label=NA,
#      edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# 
# dev.off()

# 最后添加删除color和label项可显示标签和点颜色边框
# pdf(file = "network_results/co-occurrence_network_vertex.frame.color.pdf",width = 13)
# 
# # 显示节点的外边框vertex.frame.color，默认为黑色
# set.seed(123)
# plot(igraph,main="Co-occurrence network",vertex.frame.color="blue",vertex.label.cex=0.5,
#      edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# 
# # 还可以改变标签离节点距离vertex.label.dist；标签颜色vertex.label.color；标签角度vertex.label.degree
# 
# dev.off()

# pdf(file = "network_results/co-occurrence_network_degree_distribution.pdf")
# 
# # random null model
# rand.g<- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = c("gnm"))
# 
# # degree_distribution
# plot(degree_distribution(igraph, cumulative = FALSE),col="red",cex=1.5,cex.axis=1,xlab="Degree",ylab="Relative frequency",main="The distribution of degree for co-occurrence network")
# points(degree_distribution(rand.g, cumulative = FALSE),col="blue",cex=1.5,type="p")
# legend("topright",c("Co-occurrence network","Erd??s–Rényi network"),col=c("red","blue"),pch=1,cex=1, bty="n")
# 
# dev.off()

# 展示大随机网络度分布图
# g <- erdos.renyi.game(1000, 10000,type = c("gnm"))
# degree_distribution(g)
# 
# plot(degree_distribution(g, cumulative = FALSE),xlab="Degree",main="The distribution of degree for co-occurrence network")

# net_pro:自定义函数，提供需要计算网络性质的igraph对象，结果返回计算好的网络性质
# netpro_result<-net_pro(igraph)
# write.csv(netpro_result,"network_results/igraph.network.pro.csv")
# 
# # node_pro:自定义函数，提供需要计算网络性质的igraph对象，结果返回计算好的节点性质
# nodepro_result<-node_pro(igraph)
# write.csv(nodepro_result,"network_results/igraph.node.pro.csv")
# 
# # igraph对象写出到文件
# write_graph(igraph, "network_results/igraph_edgelist.txt", format="edgelist")
# E(igraph)$weight<-igraph.weight 
# write_graph(igraph, "network_results/igraph_col.txt", format="ncol")#含节点名字
# igraph保存文件的读取
# read_graph("network_results/igraph_edgelist.txt", format ="edgelist")
# read_graph("network_results/igraph_col.txt", format ="ncol")
## 模块数量
max(comps)
## 每个OTU所属模块
comps

module = data.frame(nodename=names(comps),nodecommunity=as.vector(comps))
write.csv(module,"network_results/module.csv")

## 按照某个模块构建子图
modules.g<-induced_subgraph(igraph,comps==1)
E(modules.g)$weight<-NA

pdf(file = "network_results/co-occurrence_network_degree_distribution.pdf")
plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
dev.off()
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_with_kk,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_with_lgl,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_with_mds,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_with_sugiyama,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_with_graphopt,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_with_gem,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_with_fr,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_randomly,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_on_grid,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_nicely,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_as_tree,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_as_star,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
# plot(modules.g,main="Module1 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
#      layout=layout_in_circle,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

modules.g<-induced_subgraph(igraph,comps==2)
E(modules.g)$weight<-NA
pdf(file = "network_results/comps-2.pdf")
plot(modules.g,main="Module2 co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
dev.off()


modules.g<-induced_subgraph(igraph,comps==3)
E(modules.g)$weight<-NA
pdf(file = "network_results/comps-3.pdf")
plot(modules.g,main="Module3 co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
dev.off()

modules.g<-induced_subgraph(igraph,comps==4)
E(modules.g)$weight<-NA
pdf(file = "network_results/comps-4.pdf")
plot(modules.g,main="Module4 co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     layout=layout_with_kk,edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
dev.off()


modules.g<-induced_subgraph(igraph,comps==5)
E(modules.g)$weight<-NA
pdf(file = "network_results/comps-5.pdf")
plot(modules.g,main="Module5 co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
dev.off()

modules.g<-induced_subgraph(igraph,comps==6)
E(modules.g)$weight<-NA
pdf(file = "network_results/comps-6.pdf")
plot(modules.g,main="Module6 co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
dev.off()


modules.g<-induced_subgraph(igraph,comps==7)
E(modules.g)$weight<-NA
pdf(file = "network_results/comps-7.pdf")
plot(modules.g,main="Module7 co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
dev.off()


modules.g<-induced_subgraph(igraph,comps==8)
E(modules.g)$weight<-NA
pdf(file = "network_results/comps-8.pdf")
plot(modules.g,main="Module8 co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
dev.off()


modules.g<-induced_subgraph(igraph,comps==9)
E(modules.g)$weight<-NA
pdf(file = "network_results/comps-9.pdf")
plot(modules.g,main="Module9 co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
dev.off()


modules.g<-induced_subgraph(igraph,comps==9)
E(modules.g)$weight<-NA
pdf(file = "network_results/comps-9.pdf")
plot(modules.g,main="Module9 co-occurrence network",vertex.frame.color=NA,vertex.label.cex=0.6,
     edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,1.2, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=1,bty="n")
dev.off()

# 子图（subgraph)的创建
# 创建包含特定顶点的子图
# subg<-induced_subgraph(igraph, c("OTU_10","OTU_20","OTU_33"))
# # 创建包含特定边的子图
# g <- make_ring(10)
# g2 <- induced_subgraph(g, 1:7)
# g3 <- subgraph.edges(g, 1:5)

rand.g.netpro_result<-c()
for (i in 1:1000){
  #####random null model
  rand.g<- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = c("gnm"))
  suppressWarnings(suppressMessages(tem_netpro_result<-net_pro(rand.g)))
  rand.g.netpro_result<-cbind(rand.g.netpro_result,tem_netpro_result)
}
write.csv(rand.g.netpro_result,"network_results/rand.g.1000.result.csv")

# # 对随机矩阵结果求均值和sd值
result_summary<-cbind(rowMeans(rand.g.netpro_result),apply(rand.g.netpro_result,1,sd))
colnames(result_summary)<-c("Means","SD")
write.csv(result_summary,"network_results/rand.g.1000.result.summary.csv")