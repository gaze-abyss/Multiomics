####################update 2023/5/24
##########downstream analysis for 16s
work_dir <- '/path/test/'

library(igraph)
library(ggsci)
library(scales)

get_col <- function(one){
  one_list <- names(table(one))
  all_col <- pal_nejm()(10)[4:8][1:length(one_list)]
  names(all_col) <- one_list
  res <- c()
  for(i in one){
    res <- c(res, all_col[i])
  }
  return(res)
}

G <- list()
set.seed(3)

#Pseudo-p-valued matrix
exp1 <- read.delim(paste0(work_dir,'t1.txt'), row.names = 1, sep = '\t', check.names = FALSE)
pvals1 <- read.delim(paste0(work_dir,'boot1/pvals.two_sided.txt'), row.names = 1, sep = '\t', check.names = FALSE)

corm1 <- read.delim(paste0(work_dir,'cor_sp_t1.txt'), row.names = 1, sep = '\t', check.names = FALSE)
cnet1 <- graph_from_adjacency_matrix(as.matrix(corm1),weight=T)
e_cnet1 <- as_data_frame(cnet1,'edge')
w_list1 <- e_cnet1$weight
names(w_list1) <- paste0(e_cnet1[,1],'-',e_cnet1[,2])

corm1[corm1 <= 0] <- -1
corm1[corm1 > 0] <- 1
corm1[corm1 == -1] <- 0
corm1 <- 1

pvals1[pvals1 >= 0.05] <- -1
pvals1[pvals1 < 0.05 & pvals1 >= 0] <- 1
pvals1[pvals1 == -1] <- 0
adj1 <- pvals1 * corm1
diag(adj1) <- 0

exp2 <- read.delim(paste0(work_dir,'t2.txt'), row.names = 1, sep = '\t', check.names = FALSE)
pvals2 <- read.delim(paste0(work_dir,'boot2/pvals.two_sided.txt'), row.names = 1, sep = '\t', check.names = FALSE)

corm2 <- read.delim(paste0(work_dir,'cor_sp_t2.txt'), row.names = 1, sep = '\t', check.names = FALSE)
cnet2 <- graph_from_adjacency_matrix(as.matrix(corm2),weight=T)
e_cnet2 <- as_data_frame(cnet2,'edge')
w_list2 <- e_cnet2$weight
names(w_list2) <- paste0(e_cnet2[,1],'-',e_cnet2[,2])

corm2[corm2 <= 0] <- -1
corm2[corm2 > 0] <- 1
corm2[corm2 == -1] <- 0
corm2 <- 1

pvals2[pvals2 >= 0.05] <- -1
pvals2[pvals2 < 0.05 & pvals2 >= 0] <- 1
pvals2[pvals2 == -1] <- 0
adj2 <- pvals2 * corm2
diag(adj2) <- 0

adj <- adj1 + adj2

net <- graph_from_adjacency_matrix(as.matrix(adj), mode = 'undirected')
l <- layout_with_graphopt(net,charge=0.1,niter = 5000)

net <- graph_from_adjacency_matrix(as.matrix(adj1), mode = 'undirected')
clu <- components(net)
conn <- groups(clu)

big <- -Inf
pick <- 1
for(i in 1:length(conn)){
  if(length(conn[[i]]) > big){
    big <- length(conn[[i]])
    pick <- i
  }
}

V(net)$vertex.color <- '#E18727FF'

edge_table <- as_data_frame(net, what = c("edges"))
edge_table$cor <- unlist(lapply(as.list(w_list1[paste0(edge_table[,1],'-',edge_table[,2])]),function(x) {if(x>0){'+'}else{'-'}} ))
write.table(edge_table,file=paste0(work_dir,'edge_table.1.txt'),quote=F,sep='\t')

E(net)$edge.color <- unlist(lapply(as.list(w_list1[paste0(edge_table[,1],'-',edge_table[,2])]),function(x) {if(x>0){'#B33B2A'}else{'#106CA9'}} ))

abundance_score_1 <- as.numeric(rowMeans(exp1))
abundance_score_mat_1 <- data.frame(
  value = abundance_score_1,
  name = 'abundance_score_1'
)

hub_score_1 <- hub_score(net,scale=F)$vector
score_mat_1 <- data.frame(
  value = hub_score_1,
  name = 'hub_score_1'
)

# scale_hub_score_1 <- (hub_score_1 - min(hub_score_1))/(max(hub_score_1)-min(hub_score_1))
# ran <- quantile(as.numeric(hub_score_1+1), probs = c(0.05,0.95))
# V(net)$vertex.color <- circlize::colorRamp2(c(ran[1],ran[2]),c("white","#B33B2A"))(hub_score_1+1)

G[[1]] <- net

net <- graph_from_adjacency_matrix(as.matrix(adj2), mode = 'undirected')
clu <- components(net)
conn <- groups(clu)

big <- -Inf
pick <- 1
for(i in 1:length(conn)){
  if(length(conn[[i]]) > big){
    big <- length(conn[[i]])
    pick <- i
  }
}

V(net)$vertex.color <- '#E18727FF'

edge_table <- as_data_frame(net, what = c("edges"))
edge_table$cor <- unlist(lapply(as.list(w_list2[paste0(edge_table[,1],'-',edge_table[,2])]),function(x) {if(x>0){'+'}else{'-'}} ))
write.table(edge_table,file=paste0(work_dir,'edge_table.2.txt'),quote=F,sep='\t')

E(net)$edge.color <- unlist(lapply(as.list(w_list2[paste0(edge_table[,1],'-',edge_table[,2])]),function(x) {if(x>0){'#B33B2A'}else{'#106CA9'}} ))

abundance_score_2 <- as.numeric(rowMeans(exp2))
abundance_score_mat_2 <- data.frame(
  value = abundance_score_2,
  name = 'abundance_score_2'
)

hub_score_2 <- hub_score(net,scale=F)$vector
score_mat_2 <- data.frame(
  value = hub_score_2,
  name = 'hub_score_2'
)

G[[2]] <- net

abundance_all_score_mat <- rbind(abundance_score_mat_1,abundance_score_mat_2)
abundance_all_score_mat[,1] <- (abundance_all_score_mat[,1]-min(abundance_all_score_mat[,1]))/(max(abundance_all_score_mat[,1])-min(abundance_all_score_mat[,1]))
abundance_all_score_mat[,1][which(abundance_all_score_mat[,1] < 0.05)] <- 0.05
abundance_score_1 <- abundance_all_score_mat[abundance_all_score_mat[,2] == 'abundance_score_1',][,1]
abundance_score_2 <- abundance_all_score_mat[abundance_all_score_mat[,2] == 'abundance_score_2',][,1]

all_score_mat <- rbind(score_mat_1,score_mat_2)
all_score_mat[,1] <- (all_score_mat[,1]-min(all_score_mat[,1]))/(max(all_score_mat[,1])-min(all_score_mat[,1]))
write.table(all_score_mat,file=paste0(work_dir,'all_score_mat.txt'),quote=F,sep='\t')

ran <- quantile(as.numeric(abundance_score_1), probs = c(0.01,0.99))
all_score_mat[,1][which(all_score_mat[,1] < 0.2)] <- 0.2
hub_score_1 <- all_score_mat[all_score_mat[,2] == 'hub_score_1',][,1]
hub_score_2 <- all_score_mat[all_score_mat[,2] == 'hub_score_2',][,1]

V(G[[1]])$vertex.color <- circlize::colorRamp2(c(ran[1],ran[2]),c("pink","#B33B2A"))(abundance_score_1*2)
V(G[[2]])$vertex.color <- circlize::colorRamp2(c(ran[1],ran[2]),c("pink","#B33B2A"))(abundance_score_2*2)

pdf(paste0(work_dir,'network.all.pdf'),width = 20,height = 10)
par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(G[[1]],layout=l, #layout_with_fr
     edge.color = E(G[[1]])$edge.color,
     edge.width = 3,
     vertex.color = V(G[[1]])$vertex.color,
     vertex.frame.color = V(G[[1]])$vertex.color,
     vertex.label.color = 'black', # 3B3B3B
     vertex.label.cex = 1.2,
     vertex.size = hub_score_1*1.5*15)
text(x = -1, y = 1.1, 'A', cex=2.5)
legend("bottomleft", legend = c("      0.9","      0.6","      0.3","      0"),title.adj=0.5,cex=1,bty='o',col='#E18727FF',y.intersp = 3.2,pch = 19, pt.cex = (1*11)*seq(1,0.2,-0.2),pt.lwd = 1,title = "Hub centrality score", trace=F)
plot(G[[2]],layout=l, #layout_with_fr
     edge.color = E(G[[2]])$edge.color,
     edge.width = 3,
     vertex.color = V(G[[2]])$vertex.color,
     vertex.frame.color = V(G[[2]])$vertex.color,
     vertex.label.color = 'black',
     vertex.label.cex = 1.2,
     vertex.size = hub_score_2*1.5*15)
text(x = -1, y = 1.1, 'B', cex=2.5)
legend("bottomright", legend = c("Positive correlation","Negative correlation"),lty = 1,lwd=3.5,pt.cex=3.5,y.intersp = 2,col=c('#B33B2A','#106CA9'))
dev.off()

pdf(paste0(work_dir,'colorlegend.pdf'),width = 5,height = 10)
color3 <- colorRampPalette(c("pink","#B33B2A"))(200)
plot(y = c(1, (length(color3) + 1)), x = c(1, 2),xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n',type = 'n', ann = F)
rect(xleft = rep(0, (length(color3) + 1)),ybottom = 1:length(color3),xright = rep(2, (length(color3) + 1)),ytop = 2:(length(color3) + 1),col = color3,border = color3)
dev.off()
