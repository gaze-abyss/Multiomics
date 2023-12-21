###############date 2023/5/24
####correlation analysis for RNA-PRO
setwd("/path/")

dat <- read.table('RNA-PRO-merge.txt',header=T,stringsAsFactors=F)
rn <- dat[,1]
dat <- dat[,-1]
rownames(dat) <- rn
dat <- dat[order(rownames(dat)),]
dat2 <- t(dat)
genename <- strsplit(x = rn,split = '|',fixed = T)
genename <- sapply(genename,'[[',1)
genename <- unique(genename)

p <- c()
cor_value <- c()
genelist <- c()
for (gene in genename) {
  exp <- paste0(gene,'|PRO')
  methy <- paste0(gene,'|RNA')
  if (exp %in% rn & methy %in% rn) {
    genelist <- c(genelist,gene)
    dat_cor <- cbind(dat2[,exp],dat2[,methy])
    colnames(dat_cor) <- c(exp,methy)
    res <- rcorr(as.matrix(dat_cor),type = 'spearman')
    cor_value <- c(cor_value,res$r[1,2])
    p <- c(p,res$P[1,2])
  }
}
out <- data.frame(gene=genelist,cor=cor_value,p_value=p)  #correlation results
