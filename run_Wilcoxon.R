##############update 2023/5/24
####Wilcoxon rank-sum test-Difference analysis
setwd("/path/test/")
dat <- read.table('pro_matrix.txt',sep = '\t',header = T,row.names = 1)
group <- read.table('group.txt',header = T)
group <- group[order(group$id),]
dat <- t(dat)
dat <- dat[order(rownames(dat)),]
dat <- as.data.frame(dat)
dat$group <- group$group
wl <- c()
pl <- c()
LFC <- c()
T_mean<- c()
N_mean<- c()
for (i in colnames(dat)[1:(ncol(dat)-1)]) {
 x <- dat[dat$group=='T',i]
 y <- dat[dat$group=='NAT',i]
 xid <- which(!is.na(x))
 yid <- which(!is.na(y))
 se <- intersect(xid,yid)
 if (sum(se) != 0 & length(se) >=3) {
   x <- x[se]
   y <- y[se]
   LFC <- c(LFC,mean(x)-mean(y))
   wt <-  wilcox.test(x,y,paired = T)
   w <- wt$statistic
   p <- wt$p.value
   T_mean<- c(T_mean,mean(x))
   N_mean<- c(N_mean,mean(y))
   wl <- c(wl,w)
   pl <- c(pl,p)
 }else{wl <- c(wl,NA) 
       pl <- c(pl,NA)
       LFC <- c(LFC,NA)
       T_mean<- c(T_mean,NA)
       N_mean<- c(N_mean,NA)}
}

p_just <- p.adjust(pl, "BH")
output <- data.frame(gene=colnames(dat)[-ncol(dat)],V=wl,pvalue=pl,p_just=p_just,LogFC=LFC,T_mean=T_mean,N_mean=N_mean)#Difference analysis results 


