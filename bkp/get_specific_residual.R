j3 <- 1
psig <- Popt[j3]
Popt <- extract.par(par.opt,Np=Nsig,bases=bases)$P
tmp <- cal.residual(par.opt,bases=bases)
rvms.list <- tmp$res.sig# rv-noise
##calculate residual after subtraction of single signal
par1 <- rep(0,length(par.opt))
if(prior.type!='e0'){
    par1[(j3-1)*5+1:5] <-par.opt[(j3-1)*5+1:5]
}else{
    par1[(j3-1)*3+1:3] <-par.opt[(j3-1)*3+1:3]
}
names(par1) <- names(par.opt)
rv.kep <- cal.residual(par1,bases=bases)$rv.sig# residual
for(i in ins){
    t <- out[[i]]$RV[,1]
    ey <- out[[i]]$RV[,3]
    y <- rvms.list[[i]]+rv.kep[[i]]
    dat <- cbind(t,y,ey)
    names(dat) <- c('BJD','RV','eRV')
    fout <- paste0('../data/combined/',target,'/',target,'_',i,'.dat')
    cat(fout,'\n')
    write.table(dat,file=fout,quote=FALSE,row.names=FALSE)
}
