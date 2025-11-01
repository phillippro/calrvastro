library(fields)
target <-'HD210193'
load(paste0('results/',target,'/',target,'_PFS_white_MA_AR_GP011_Esd0.2_ofac2_xi10_N1e+06.Robj'))
pdf('paper_MP.pdf',8,8)
size <- 2
par(cex.axis=size,cex.lab=size,cex=size)
source('MPplot.R')
dev.off()
