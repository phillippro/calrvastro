fpdf <- paste0('results/',target,'/',target,'_res.pdf')
cat(fpdf,'\n')
pdf(fpdf,6,6)
tmp  <- RV.kepler(par.opt)
oT <- cbind(tmp$tauT,out$res.comb.all[[paste0('sig',Nsig)]],out$all$eRV.all)
if(grepl('HD128620|HD128621',target)){
                                        #   ind <- which(abs(oT[,2])<20)
    ind <- which(abs(oT[,2]-mean(oT[,2]))<1*sd(oT[,2]))
    oT <- oT[ind,]

###binary RV residual
    resC <- out$data.binary[,'RV']-tmp$rvC
    tauC <- tmp$tauC
    ervC <- out$data.binary[,'eRV']
    oC <- cbind(tauC,resC,ervC)
    ind <- which(abs(oC[,2]-mean(oC[,2]))<1*sd(oC[,2]))
    oC <- oC[ind,]
    colnames(oC) <-  c('BJD','RV','eRV')
    binary <- gsub('HD128620','HD128621',target)
    fC <- paste0('results/',target,'/',binary,'_res.dat')
    cat(fC,'\n')
    write.table(oC,file=fC,quote=FALSE,row.names=FALSE)
    plot(oC[,1],oC[,2],xlab='Emission Time [JD]',ylab='O-C [m/s]',main=binary)
    arrows(oC[,1],oC[,2]-oC[,3],oC[,1],oC[,2]+oC[,3],length=0.01,angle=90,code=3,col='grey')
}
colnames(oT) <- c('BJD','RV','eRV')
fT <- paste0('results/',target,'/',target,'_res.dat')
cat(fT,'\n')
write.table(oT,file=fT,quote=FALSE,row.names=FALSE)
plot(oT[,1],oT[,2],xlab='Emission Time [JD]',ylab='O-C [m/s]',main=target)
arrows(oT[,1],oT[,2]-oT[,3],oT[,1],oT[,2]+oT[,3],length=0.01,angle=90,code=3,col='grey')
dev.off()
