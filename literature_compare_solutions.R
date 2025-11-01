source('mcmc_func.R')
library(magicaxis)
#dr2 <- read.table('dr2_par0.txt',header=TRUE)
#hg3 <- read.table('hg3_par0.txt',header=TRUE)
#hgca <- read.table('hgca_par0.txt',header=TRUE)
#stars <- c('GL676A','HD145675','HD149806','HD62364','HD115404A','HD42581')
stars <- c('HD209100','HD22049')
tab <- read.csv('vizier_167companion_pars.csv',sep='|',check.names=FALSE)
star0 <- gsub(' ','',tab[,3])
stars[stars=='GL676A'] <- 'CD-5110924'
stars[stars=='HD145675'] <- '14Her'
stars <- gsub('GL','GJ',stars)
stars <- gsub('A$','',stars)
index <- c()
for(star in stars){
    ii <- which(star0==star & tab[,'Per-d']>100)
#    ii <- ii[sort(tab[ii,'Per-d'],index.return=TRUE)$ix]
    index <- c(index,ii)
}
par0 <- tab[index,c('Companion','ImFlag','Host','mc','e_mc','E_mc','Per-d','e_Per-d','E_Per-d','a','e_a','E_a')]
par1 <- read.table('companion_mpa.txt',header=TRUE)

fpdf <- 'compare_solution_mpa.pdf'
cat(fpdf,'\n')
cols <- rainbow(nrow(par0))
pdf(fpdf,16,16)
size <- 1.2
par(mfrow=c(2,2),cex.lab=size,cex.axis=size,cex=size,mar=c(5,5,2,2))

ltys <- c(2, 1)
pchs <- c(1, 20)

for(j in 1:2){
if(j==1){
xlab0 <- c('Per-d','e_Per-d','E_Per-d')
xlab1 <- c('pc','e_pc','E_pc')
xlab <- 'P [day]'
}else if(j==2){
xlab0 <- c('a','e_a','E_a')
xlab1 <- c('ac','e_ac','E_ac')
xlab <- 'a [au]'
}
ylab1 <- ylab0 <- c('mc','e_mc','E_mc')
ylab <- 'mass [Mjup]'

xlim <- range(par0[,xlab0[1]]-par0[,xlab0[2]],par0[,xlab0[1]]+par0[,xlab0[3]],par1[,xlab1[1]]-par1[,xlab1[2]],par1[,xlab1[1]]+par1[,xlab1[3]])
ylim <- range(par0[,ylab0[1]]-par0[,ylab0[2]],par0[,ylab0[1]]+par0[,ylab0[3]],par1[,ylab1[1]]-par1[,ylab1[2]],par1[,ylab1[1]]+par1[,ylab1[3]])

plot(par0[,xlab0[1]],par0[,ylab0[1]],xlab=xlab,ylab=ylab,col=cols,cex=2,log='xy',xlim=xlim,ylim=ylim,pch=pchs[1],axes=FALSE)
magaxis(1:2,frame.plot=TRUE,grid=TRUE,minorn=10)
magaxis(3:4)
text(par0[,xlab0[1]],par0[,ylab0[1]],labels=par0[,'Companion'],col=cols,pos=4,xpd=NA)

points(par1[,xlab1[1]],par1[,ylab1[1]],col=cols,cex=2,pch=pchs[2])

arrows(par0[,xlab0[1]],par0[,ylab0[1]]-par0[,ylab0[2]],par0[,xlab0[1]],par0[,ylab0[1]]+par0[,ylab0[3]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[1])
arrows(par0[,xlab0[1]]-par0[,xlab0[2]],par0[,ylab0[1]],par0[,xlab0[1]]+par0[,xlab0[3]],par0[,ylab0[1]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[1])

arrows(par1[,xlab1[1]],par1[,ylab1[1]]-par1[,ylab1[2]],par1[,xlab1[1]],par1[,ylab1[1]]+par1[,ylab1[3]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[2])
arrows(par1[,xlab1[1]]-par1[,xlab1[2]],par1[,ylab1[1]],par1[,xlab1[1]]+par1[,xlab1[3]],par1[,ylab1[1]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[2])

legend('topleft',legend=c('Old solution','New solution'),pch=pchs,lty=ltys,bty='n')

plot(par0[,xlab0[1]],par1[,xlab1[1]],xlab=xlab,ylab=xlab,log='xy',axes=FALSE)
magaxis(1:4)
abline(a=0,b=1,col='grey')

arrows(par0[,xlab0[1]]-par0[,xlab0[2]],par1[,xlab1[1]],par0[,xlab0[1]]+par0[,xlab0[3]],par1[,xlab1[1]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[1])
arrows(par0[,xlab0[1]],par1[,xlab1[1]]-par1[,xlab1[2]],par0[,xlab0[1]],par1[,xlab1[1]]+par1[,xlab1[3]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[2])
text(par0[,xlab0[1]],par1[,xlab1[1]],labels=par0[,'Companion'],col=cols,pos=1,xpd=NA)
}

plot(par0[,ylab0[1]],par1[,ylab1[1]],xlab=ylab,ylab=ylab,log='xy',axes=FALSE)
magaxis(1:4)
abline(a=0,b=1,col='grey')

arrows(par0[,ylab0[1]]-par0[,ylab0[2]],par1[,ylab1[1]],par0[,ylab0[1]]+par0[,ylab0[3]],par1[,ylab1[1]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[1])
arrows(par0[,ylab0[1]],par1[,ylab1[1]]-par1[,ylab1[2]],par0[,ylab0[1]],par1[,ylab1[1]]+par1[,ylab1[3]],length=0.001,angle=0,code=1,col=tcol(cols,20),lty=ltys[2])
text(par0[,ylab0[1]],par1[,ylab1[1]],labels=par0[,'Companion'],col=cols,pos=1,xpd=NA)

dev.off()

