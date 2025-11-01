source('mcmc_func.R')
tab <- read.table('companion_files.txt')[,1]
stars <- c('HD39091')
ss <- gsub('\\/.+','',tab)
ind <- match(stars,ss)
pdf.file <- tab[ind]
robj.file <- paste0('results/',gsub('.pdf','.Robj',pdf.file))
if(!exists(out)) load(robj.file)
mlow <- 13
mup <- 75
###select BDs
for(j in 1:Nsig){

}

fout <- 'plot_misalignment.pdf'
cat(fout,'\n')
pdf(fout,6,6)
#plot(x2,y2,log='x',ylim=c(0,150))
#plot(x2[index],y2[index],ylim=c(0,150),xlim=c(0.5,2))
#plot(x2[index],y2[index],ylim=c(0,180))
plot(x2[index],y2[index],xlab='log(a2/a1)',ylab='Misalignment [deg]')
arrows(x2[index],y2[index]-ey2[index],x2[index],y2[index]+ey2[index],length=0.05,angle=90,code=3,col='grey')

if(length(x3)>1){
    points(x3,y3,col='red')
    arrows(x3,y3-ey3,x3,y3+ey3,length=0.05,angle=90,code=3,col='red')
}
#points(exp(tmp[,1]),tmp[,2],pch=20)
if(length(index)>20){
    points(tmp[,1],tmp[,2],pch=20)
    arrows(tmp[,1],tmp[,2]-tmp[,3],tmp[,1],tmp[,2]+tmp[,3],length=0.05,angle=90,code=3,col='black')
}

dev.off()
