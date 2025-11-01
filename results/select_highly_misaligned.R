args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    mi <- as.numeric(args[1])
    m0 <- as.numeric(args[2])
}else{
    mi <- 10
    m0 <- 0.6
}
psi0 <- 90
tab <- read.table('misaligned_starpar.txt',header=TRUE)
cn <- colnames(tab)
cn <- gsub('plus','upper',cn)
cn <- gsub('minus','lower',cn)
colnames(tab) <- cn
ind0 <- grep('psi.map',cn)
ind1 <- grep('psi.lower',cn)
ind2 <- grep('psi.upper',cn)
ind.per <- grep('Per.map',cn)
candidates2 <- candidates1 <- c()
index <- c()
dpsis <- psis <- c()
dpsis1 <- psis1 <- c()
dpsis2 <- psis2 <- c()
pp2 <- pp1 <- pp <- c()
for(j in 1:nrow(tab)){
    ii <- which(!is.na(tab[j,ind0]))
    psi.opt <- unlist(lapply(ii, function(i) as.numeric(unlist(strsplit(as.character(tab[j,ind0[i]]),'d')))))
    psi.lower <- unlist(lapply(ii, function(i) as.numeric(unlist(strsplit(as.character(tab[j,ind1[i]]),'d')))))
    psi.upper <- unlist(lapply(ii, function(i) as.numeric(unlist(strsplit(as.character(tab[j,ind2[i]]),'d')))))
    per.opt <- unlist(lapply(ii, function(i) rep(tab[j,ind.per[i]],length(unlist(strsplit(as.character(tab[j,ind0[i]]),'d'))))))
    dpsi <- 0.5*abs(psi.upper-psi.lower)
    psis <- c(psis,psi.opt)
    dpsis <- c(dpsis,dpsi)
    pp <- c(pp,per.opt)
    star <- tab[j,1]
    mstar <- tab[j,'Mstar']
    emstar <- tab[j,'eMstar']
    emstar <- tab[j,'eMstar']
    f <- tab[j,'file']
    if(mstar<m0){
        if(any(dpsi<mi & psi.opt>psi0)){
            kk <- which(dpsi<mi & psi.opt>psi0)
            for(k in kk){
                if(!any(psis1==psi.opt[k])){
                    psis1 <- c(psis1,psi.opt[k])
                    dpsis1 <- c(dpsis1,dpsi[k])
                    pp1 <- c(pp1,per.opt[k])
                }
                candidates1 <- c(candidates1,star)
                index <- c(index,j)
            }
        }
        if(any(dpsi<mi)){
            jj <- which(dpsi<mi)
            for(j1 in jj){
                if(!any(psis2==psi.opt[j1])){
                    psis2 <- c(psis2,psi.opt[j1])
                    dpsis2 <- c(dpsis2,dpsi[j1])
                    pp2 <- c(pp2,per.opt[j1])
                }
                candidates2 <- c(candidates2,star)
            }
            cat(star,';Mstar=',mstar,'Msun; psi=',psi.opt[jj],'deg; dpsi=',dpsi[jj],'deg\n',f,'\n\n')
        }
    }
}
fout <- 'misaglined_system.txt'
write.table(candidates2,file=fout,quote=FALSE,row.names=FALSE)
#tab[index,c(1,ind0,ind1,ind2)]
#ii <- which(psis1>90)
fpdf <- paste0('misalignment_',mi,'deg.pdf')
cat(fpdf,'\n')
pdf(fpdf,8,12)
par(mfrow=c(3,2))
ii <- which(duplicated(pp))
pp <- pp[-ii]
psis <- psis[-ii]
plot(exp(pp),psis,xlab='Period [d]',ylab=expression(psi*' [deg]'),log='x')
points(exp(pp2),psis2,col='blue')
points(exp(pp1),psis1,col='red')

plot(exp(pp),cos(psis/180*pi),xlab='Period [d]',ylab=expression(cos*psi),log='x')
points(exp(pp2),cos(psis2/180*pi),col='blue')
points(exp(pp1),cos(psis1/180*pi),col='red')

p <- hist(psis,xlab=expression(psi*' [deg]'),ylab='Freq.',main='')
hist(psis2,col='blue',add=TRUE,breaks=p$breaks)
hist(psis1,col='red',add=TRUE,breaks=p$breaks)

p <- hist(cos(psis/180*pi),xlab=expression(cos*psi),ylab='Freq.',main='')
hist(cos(psis2/180*pi),col='blue',add=TRUE,breaks=p$breaks)
hist(cos(psis1/180*pi),col='red',add=TRUE,breaks=p$breaks)

p <- hist(psis2,xlab=expression(psi*' [deg]'),ylab='Freq.',main='Conservative sample')
hist(psis1,col='red',add=TRUE,breaks=p$breaks)

p <- hist(cos(psis2/180*pi),xlab=expression(cos*psi),ylab='Freq.',main='Conservative sample')
hist(cos(psis1/180*pi),col='red',add=TRUE,breaks=p$breaks)

title(main=paste('Error upper limit:',mi,'deg;Ncons=',length(unique(candidates2)),'(',length(unique(candidates1)),')'),xpd=NA,outer=TRUE,line=-2)
dev.off()
