args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    file <- as.character(args[1])
    ms <- as.numeric(args[2])
    dms.low <- as.numeric(args[3])
    dms.up <- as.numeric(args[4])
}else{
#    file <- 'keppure_priormt_poly11_Ndata144_quantifyTRUE_1per1_Nw1_LHS1140_ind0_2planet_ARMA00_Nsamp20000000_tem1_acc0.35_pretem1P92.3d24.6d_negLmax399'
#    file <- 'keppure_priormt_poly11_Ndata144_quantifyTRUE_1per1_Nw1_LHS1140_ind0_2planet_GPqp_Nsamp4000000_tem1_acc0.67_pretem1P93.1d24.5d_negLmax400'
#    file <- 'keppure_priormt_poly11_Ndata144_quantifyTRUE_1per1_Nw1_LHS1140_ind0_2planet_ARMA01_Nsamp10000000_tem1_acc0.43_pretem1P92.1d24.6d_negLmax398'
#file <- 'keppure_priormt_poly11_Ndata138_quantifyTRUE_1per1_Nw1_LHS1140_HARPS_ind7_3planet_ARMA00_Nsamp6000000_tem1_acc0.19_pretem0.02P3.8d24.7d91.6d_negLmax383'
#file <- 'keppure_priormt_poly11_Ndata138_quantifyTRUE_1per1_Nw1_LHS1140_HARPS_ind7_3planet_ARMA01_Nsamp6000000_tem1_acc0.22_pretem0.02P3.8d24.6d91.4d_negLmax384'
#file <- 'keppure_priormt_poly11_Ndata138_quantifyTRUE_1per1_Nw1_LHS1140_HARPS_ind7_3planet_ARMA02_Nsamp6000000_tem1_acc0.26_pretem0.02P3.8d24.6d92d_negLmax382'
file <- 'keppure_priormt_poly11_Ndata138_quantifyTRUE_1per1_Nw1_LHS1140_HARPS_ind7_3planet_GPqp_Nsamp6000000_tem1_acc0.35_pretem0.02P3.8d24.7d91.5d_negLmax383'
    ms <- 0.146 
    dms.low <- 0.019
    dms.up <- 0.019
}
uncertainty <- 'q1'#; sd or q10 or q1
###note: dms is one sigma
###convert dms.low and dms.up to 1% and 99% confidence level uncertainty which are corresponding to -2.326sigma and 2.326sigma. 
if(uncertainty=='sd'){
    dms.low <- dms.low
    dms.up <- dms.up
}else if(uncertainty=='q1'){
    dms.low <- dms.low*2.3270
    dms.up <- dms.up*2.3270
}else if(uncertainty=='q10'){
    dms.low <- dms.low*1.2814
    dms.up <- dms.up*1.2814
}
if(!grepl('Robj',file)){
    file <- paste0(file,'.Robj')
}
f1 <- gsub('_ind.+','',file)
f1 <- gsub('_HARPS','',f1)
f1 <- gsub('.+Nw\\d_','',f1)
target <- gsub('.+_','',f1)
cat('target=',target,'\n')
folder <- paste0('../output/',target,'/')
file <- paste0(folder,file)
load(file)
load(gsub('.Robj','_optpar.Robj',file))
source('mcmc_func.R')
Ps <- exp(mcmc.out[,grepl('per',names(par.opt.post))])
Ks <- mcmc.out[,grepl('K',names(par.opt.post))]
es <- mcmc.out[,grepl('e\\d',names(par.opt.post))]
Ps.opt <- exp(par.opt.post[grepl('per',names(par.opt.post))])
Ks.opt <- par.opt.post[grepl('K',names(par.opt.post))]
es.opt <- par.opt.post[grepl('e\\d',names(par.opt.post))]
omega.opt <- par.opt.post[grepl('omega\\d',names(par.opt.post))]
Mo.opt <- par.opt.post[grepl('Mo\\d',names(par.opt.post))]
if(uncertainty=='q1'){
    qq1 <- 'x1per'
    qq2 <- 'x99per'
}else if(uncertainty=='q10'){
    qq1 <- 'x10per'
    qq2 <- 'x90per'
}else{
    qq1 <- 'xminus.1sig'
    qq2 <- 'xplus.1sig'
}
Ps.low <- par.stat[qq1,grep('per',colnames(par.stat))]
Ps.up <- par.stat[qq2,grep('per',colnames(par.stat))]
Ks.low <- par.stat[qq1,grep('K',colnames(par.stat))]
Ks.up <- par.stat[qq2,grep('K',colnames(par.stat))]
es.low <- par.stat[qq1,grep('e\\d',colnames(par.stat))]
es.up <- par.stat[qq2,grep('e\\d',colnames(par.stat))]
omega.low <- par.stat[qq1,grep('omega\\d',colnames(par.stat))]
omega.up <- par.stat[qq2,grep('omega\\d',colnames(par.stat))]
Mo.low <- par.stat[qq1,grep('Mo\\d',colnames(par.stat))]
Mo.up <- par.stat[qq2,grep('Mo\\d',colnames(par.stat))]
###jupiter to Earth mass ratio
Mj2e <- 317.83
mma.opt <- K2msini(Ks.opt,Ps.opt,es.opt,ms)
mj.opt <- mma.opt$mj
me.opt <- mma.opt$me
ap.opt <- mma.opt$a
mls <- c()
mus <- c()
als <- c()
aus <- c()
as <- ap.opt
mes <- me.opt
dms.low.rel <- dms.low/ms##used to calculate the uncertainty intervals; total uncertainty=quarterly sum 2/3*msini*dms.rel and dmsini caused by other parameters; 
dms.up.rel <- dms.up/ms
if(!is.matrix(Ps)){
    Ps <- matrix(Ps,ncol=1)
    Ks <- matrix(Ks,ncol=1)
    es <- matrix(es,ncol=1)
}
for(j in 1:ncol(Ps)){
#    mma <- K2msini(Ks[,j],Ps[,j],es[,j],rnorm(length(es[,j]),ms,dms))
    mma <- K2msini(Ks[,j],Ps[,j],es[,j],ms)
    mj <- mma$mj
    me <- mma$me
    ap <- mma$a
    pj <- data.distr(mj,xlab='x',ylab='y',plotf=FALSE)
    if(uncertainty=='sd'){
        dmj.low <- sqrt((mj.opt[j]-pj['xminus.1sig'])^2+(mj.opt[j]*2/3*dms.low.rel)^2)
        dmj.up <- sqrt((mj.opt[j]-pj['xplus.1sig'])^2+(mj.opt[j]*2/3*dms.up.rel)^2)
    }else if(uncertainty=='q1'){
        dmj.low <- sqrt((mj.opt[j]-pj['x1per'])^2+(mj.opt[j]*2/3*dms.low.rel)^2)
        dmj.up <- sqrt((mj.opt[j]-pj['x99per'])^2+(mj.opt[j]*2/3*dms.up.rel)^2)
    }else if(uncertainty=='q10'){
        dmj.low <- sqrt((mj.opt[j]-pj['x10per'])^2+(mj.opt[j]*2/3*dms.low.rel)^2)
        dmj.up <- sqrt((mj.opt[j]-pj['x90per'])^2+(mj.opt[j]*2/3*dms.up.rel)^2)
    }
    mj.low <- mj.opt[j]-dmj.low
    mj.up <- mj.opt[j]+dmj.up
    me.low <- mj.low*Mj2e
    me.up <- mj.up*Mj2e
    mls <- c(mls,me.low)
    mus <- c(mus,me.up)
    q <- data.distr(ap,xlab='x',ylab='y',plotf=FALSE)
    da.low <- sqrt((ap.opt[j]-q['xminus.1sig'])^2+(ap.opt[j]*1/3*dms.low.rel)^2)
    da.up <- sqrt((ap.opt[j]-q['xplus.1sig'])^2+(ap.opt[j]*1/3*dms.low.rel)^2)
    ap.low <- ap.opt[j]-da.low
    ap.up <- ap.opt[j]+da.up
    als <- c(als,ap.low)
    aus <- c(aus,ap.up)
    cat('\nsignal',j,' P=',Ps.opt[j],'d:\n')
    cat(paste0('msini(Mj)=\n',formatC(mj.opt[j],format='f',digits=2),' [',formatC(mj.low,format='f',digits=2),', ',formatC(mj.up,format='f',digits=2),']'),'\n')
    cat(paste0('msini(Me)=\n',formatC(me.opt[j],format='f',digits=2),' [',formatC(me.low,format='f',digits=2),', ',formatC(me.up,format='f',digits=2),']'),'\n')
    cat(paste0('ap(au)=\n',formatC(ap.opt[j],format='f',digits=3),' [',formatC(ap.low,format='f',digits=3),', ',formatC(ap.up,format='f',digits=3),']'),'\n\n')
}

pars <- rbind(Ps.opt, Ps.low,Ps.up,Ks.opt,Ks.low,Ks.up,es.opt,es.low,es.up,omega.opt,omega.low,omega.up,Mo.opt,Mo.low,Mo.up,mes,mls,mus,as,als,aus)
ind <- sort(Popt,index.return=TRUE)$ix
pars <- pars[,ind]
#ind2 <- c(2,4,5,6,1,3)
#pars <- pars[,ind2]
for(j in 1:(nrow(pars)/3)){
    if(j==7){
        dig <- 3
    }else{
        dig <- 2
    }
#    lapply(1:ncol(pars),function(k) cat(paste0(formatC(pars[3*(j-1)+1,k],format='f',digits=dig),' [',formatC(pars[3*(j-1)+2,k],format='f',digits=dig),', ',formatC(pars[3*(j-1)+3,k],format='f',digits=dig),']&')))
    lapply(1:ncol(pars),function(k) cat(paste0(formatC(pars[3*(j-1)+1,k],format='f',digits=dig),' [',formatC(pars[3*(j-1)+2,k],format='f',digits=dig),', ',formatC(pars[3*(j-1)+3,k],format='f',digits=dig),']&')))
cat('\n')
}

for(j in ((Np*5)+1):ncol(par.stat)){
    dig <- 2
    cat('par:',colnames(par.stat)[j],'\n')
#    cat(paste0(formatC(par.stat[1,j],format='f',digits=dig),' [',formatC(par.stat[4,j],format='f',digits=dig),', ',formatC(par.stat[5,j],format='f',digits=dig),']&'))
    cat(paste0(formatC(par.stat[1,j],format='f',digits=dig),' [',formatC(par.stat[4,j],format='f',digits=dig),', ',formatC(par.stat[5,j],format='f',digits=dig),']&'))
cat('\n')
}
