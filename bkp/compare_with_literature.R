source('mcmc_func.R')
library(magicaxis)
#dr2 <- read.table('dr2_par0.txt',header=TRUE)
#hg3 <- read.table('hg3_par0.txt',header=TRUE)
#hgca <- read.table('hgca_par0.txt',header=TRUE)
stars <- c('HD22049','HD209100')
#par0 <- tab[index,c('Companion','ImFlag','Host','mc','e_mc','E_mc','Per-d','e_Per-d','E_Per-d','a','e_a','E_a','e','e_e','E_e','I','e_I','E_I','reference')]
par0 <- data.frame(rbind(c('HD22049','b','HD22049',0.78,0.12,0.38,2692,26,26,3.48,0.02,0.02,0.07,0.05,0.06,89,42,42,'M19'),
              c('HD22049','b','HD22049',0.66,0.09,0.12,2671,23,17,3.53,0.04,0.03,0.07,0.05,0.07,78.81,22.41,29.34,'L21'),
              c('HD22049','b','HD22049',0.63,0.04,0.12,2775,5,5,3.56,0.06,0.12,0.16,0.01,0.01,45,8,8,'B22'),
              c('HD209100','b','HD209100',3.25,0.65,0.39,16509.0,1742.0,2096.0,11.55,-0.71,4.18,0.26,0.03,0.07,64.25,6.09,13.8,'F19'),
              c('HD209100','b','HD209100',3.0,0.1,0.1,10932,228,266,8.8,-0.1,0.2,0.48,0.01,0.01,91,5,4,'P23')
              ))
colnames(par0) <- c('Companion','ImFlag','Host','m','e_m','E_m','Per-d','e_Per-d','E_Per-d','a','e_a','E_a','e','e_e','E_e','I','e_I','E_I','reference')
par1 <- read.table('companion_mpa.txt',header=TRUE)
cc <- rainbow(nrow(par1))
###repeat and order
cols <- kk <- jj <- c()
for(j in 1:nrow(par1)){
    ii <- which(par0[,1]==par1[j,1] & par0[,2]==par1[j,2])
    jj <- c(jj,ii)
    kk <- c(kk,rep(j,length(ii)))
    cols <- c(cols,rep(cc[j],length(ii)))
}
par0 <- par0[jj,]
par1 <- par1[kk,]

for(j in 1:(ncol(par0)-1)){
    if(j>3) par0[,j] <- as.numeric(par0[,j])
}
for(j in 1:ncol(par1)){
    if(j>2) par1[,j] <- as.numeric(par1[,j])
}

fpdf <- 'compare_solution_mpa.pdf'
cat(fpdf,'\n')
#cols <- rainbow(nrow(par0))

pdf(fpdf,15,10)
size <- 1
par(mfrow=c(2,3),cex.lab=size,cex.axis=size,cex=size,mar=c(5,5,1,1))

ltys <- c(2, 1)
pchs <- c(1, 20)
for(j in 1:4){
if(j==1){
xlab0 <- c('Per-d','e_Per-d','E_Per-d')
xlab1 <- c('p','e_p','E_p')
xlab <- 'P [day]'
}else if(j==2){
xlab0 <- c('a','e_a','E_a')
xlab1 <- c('a','e_a','E_a')
xlab <- 'a [au]'
}else if(j==3){
xlab1 <- xlab0 <- c('e','e_e','E_e')
xlab <- 'eccentricity'
}else if(j==4){
xlab1 <- xlab0 <- c('I','e_I','E_I')
xlab <- 'Inclination [deg]'
}
ylab1 <- ylab0 <- c('m','e_m','E_m')
ylab <- 'mass [Mjup]'

xlim <- range(par0[,xlab0[1]]-par0[,xlab0[2]],par0[,xlab0[1]]+par0[,xlab0[3]],par1[,xlab1[1]]-par1[,xlab1[2]],par1[,xlab1[1]]+par1[,xlab1[3]])
ylim <- range(par0[,ylab0[1]]-par0[,ylab0[2]],par0[,ylab0[1]]+par0[,ylab0[3]],par1[,ylab1[1]]-par1[,ylab1[2]],par1[,ylab1[1]]+par1[,ylab1[3]])

plot(par0[,xlab0[1]],par0[,ylab0[1]],xlab=xlab,ylab=ylab,col=cols,cex=2,log='xy',xlim=xlim,ylim=ylim,pch=pchs[1],axes=FALSE)
magaxis(1:2,frame.plot=TRUE,grid=TRUE,minorn=10)
magaxis(3:4,labels=FALSE)
#texts <- paste0(par0[,'Companion'],';',par0[,'reference'])
texts1 <- par0[,'Companion']
texts2 <- par0[,'reference']
text(par0[,xlab0[1]],par0[,ylab0[1]],labels=texts2,col=cols,pos=4,xpd=NA)

points(par1[,xlab1[1]],par1[,ylab1[1]],col=cols,cex=2,pch=pchs[2])
text(par1[,xlab1[1]],par1[,ylab0[1]],labels=texts1,col=cols,pos=4,xpd=NA)

arrows(par0[,xlab0[1]],par0[,ylab0[1]]-par0[,ylab0[2]],par0[,xlab0[1]],par0[,ylab0[1]]+par0[,ylab0[3]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[1])
arrows(par0[,xlab0[1]]-par0[,xlab0[2]],par0[,ylab0[1]],par0[,xlab0[1]]+par0[,xlab0[3]],par0[,ylab0[1]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[1])

arrows(par1[,xlab1[1]],par1[,ylab1[1]]-par1[,ylab1[2]],par1[,xlab1[1]],par1[,ylab1[1]]+par1[,ylab1[3]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[2])
arrows(par1[,xlab1[1]]-par1[,xlab1[2]],par1[,ylab1[1]],par1[,xlab1[1]]+par1[,xlab1[3]],par1[,ylab1[1]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[2])

legend('topleft',legend=c('Old solution','New solution'),pch=pchs,lty=ltys,bty='n')

plot(par0[,xlab0[1]],par1[,xlab1[1]],xlab=xlab,ylab=xlab,log='xy',axes=FALSE)
magaxis(1:2)
magaxis(3:4,labels=FALSE)
abline(a=0,b=1,col='grey')

arrows(par0[,xlab0[1]]-par0[,xlab0[2]],par1[,xlab1[1]],par0[,xlab0[1]]+par0[,xlab0[3]],par1[,xlab1[1]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[1])
arrows(par0[,xlab0[1]],par1[,xlab1[1]]-par1[,xlab1[2]],par0[,xlab0[1]],par1[,xlab1[1]]+par1[,xlab1[3]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[2])
text(par0[,xlab0[1]],par1[,xlab1[1]],labels=texts2,col=cols,pos=1,xpd=NA)
#}

plot(par0[,ylab0[1]],par1[,ylab1[1]],xlab=ylab,ylab=ylab,log='xy',axes=FALSE)
magaxis(1:2)
magaxis(3:4,labels=FALSE)
abline(a=0,b=1,col='grey')

arrows(par0[,ylab0[1]]-par0[,ylab0[2]],par1[,ylab1[1]],par0[,ylab0[1]]+par0[,ylab0[3]],par1[,ylab1[1]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[1])
arrows(par0[,ylab0[1]],par1[,ylab1[1]]-par1[,ylab1[2]],par0[,ylab0[1]],par1[,ylab1[1]]+par1[,ylab1[3]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[2])
text(par0[,ylab0[1]],par1[,ylab1[1]],labels=texts2,col=cols,pos=1,xpd=NA)
}
dev.off()

