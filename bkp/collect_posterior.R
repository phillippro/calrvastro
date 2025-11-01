library(RColorBrewer)
source('mcmc_func.R')
options(scipen=9)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    version <- as.integer(args)
}else{
    version <- 2
}
mstars <- emstars <- c()
yr2d <- 365.25
Ns <- 1e5
value.type <- c('ms','mq')#mean+sd and map+quantile
#value.type <- c('ms')
type <- 'format'
xup <- 'x99per'
xlow <- 'x1per'
#xup <- 'xplus.1sig'
#xlow <- 'xminus.1sig'
show.digit <- function(x,Ndig,type='format'){
    if(type=='format'){
        format(round(x,Ndig),nsmall=Ndig)
    }else if(type=='signif'){
        signif(x,Ndig)
    }
}

###UVES list
#name.alt <- c('HIP 73408 B','GJ 105 C','$\\mu$ Herculis Ab','GJ 777 b','GJ 779 B','$\\chi^1$ Orionis B','GJ 36 B','GJ 22 C')
#fs <- read.table('companion_files.txt')[,1]
feng <- read.table('companion_files_sc.txt')[,1]
xiao <- read.csv('feng_cross_SWEET_new3.csv')
stars <- gsub('\\/.+','',fs)
fs <- read.table('companion_files_sc.txt')[,1]
#fs <- gsub('\\+','\\\\+',fs)
name.alt <- gsub('\\/.+','',fs)#ineed to be changed
if(version==2){
    ii <- which(is.na(match(xiao[,1],name.alt)))
    missing <- xiao[ii,1]
    fs <- read.table('companion_files.txt')[,1]
    name.alt <- gsub('\\/.+','',fs)
    ind <- match(missing,name.alt)
    fs <- fs[ind]
    name.alt <- name.alt[ind]
}
#fs <- gsub('pdf','Robj',gsub('.+\\/','',fs))
fs <- paste0('results/',gsub('pdf','Robj',fs))
#psigs <- gsub('.+\\_P|_acc.+','',fs)
#fs <- fs[grep('d',psigs)]
#bds <- read.table('companion_multi.txt',header=TRUE)[,1]
#stars <- gsub('_.+','',fs)
#ind <- match(stars,bds)
#fs <- fs[which(!is.na(ind))]
cat('length(fs)=',length(fs),'\n')
###checking whether files exist
n <- c('per','K','e','omega','Mo','Inc','Omega')
ns <- as.vector(outer(n,c(1:4),'paste0'))
Nsamp <- 1e4
NN <- Nsamp+1
Nrow <- NN*length(fs)
#fs <- fs[1:5]
out <- data.frame(array(NA,dim=c(Nrow,length(ns)+5+8)))
mcJs <- paste0('mcJ',1:4)
acs <- paste0('ac',1:4)
colnames(out) <- c('star','mass','emass1','emass2','tmin',mcJs,acs,ns)
N1 <- 1
N2 <- NN
set.seed(9999)
for(j3 in 1:length(fs)){
    f <- fs[j3]
    star <- gsub('_.+','',f)
    if(!file.exists(f)){
        f <- gsub('results','../rvastro/results',f)
    }
    cat('read:',f,'\n')
    load(f,env=e0 <- new.env())
    mc <- e0$out$mcmc.opt[[paste0('sig',e0$Nsig)]]
    ll <- mc[,'loglike']
    lp <- mc[,'logpost']
    mc <- mc[,1:(ncol(mc)-2)]
    par.opt <- mc[which.max(ll),]
#    inds <- sort(sample(1:nrow(mc),Nsamp))
#    mc <- mc[inds,]
    mc <- mc[sample(1:nrow(mc),Nsamp),]
    pp <- rbind(par.opt,mc)
    nn <- intersect(ns,colnames(mc))
    ind1 <- match(nn,colnames(out))
    ind2 <- match(nn,colnames(mc))
    out[N1:N2,ind1] <- pp[,ind2]
    Mstar <- e0$out$Mstar
    eMstar <- e0$out$eMstar
    if(eMstar<0){
        Mstar <- e0$out$astrometry[nrow(e0$out$astrometry),'mass']
        eMstar <- 0.5*(e0$out$astrometry[nrow(e0$out$astrometry),'mass.lower']+e0$out$astrometry[nrow(e0$out$astrometry),'mass.upper'])
    }
    Ms <- rnorm(2*Nsamp,Mstar,eMstar)
    Ms <- Ms[Ms>0][1:Nsamp]
    Ms <- c(Mstar,Ms)
    out[N1:N2,1] <- e0$target
    out[N1:N2,2] <- e0$out$astrometry[1,'mass']
    out[N1:N2,3] <- e0$out$astrometry[1,'mass.lower']
    out[N1:N2,4] <- e0$out$astrometry[1,'mass.upper']
    out[N1:N2,5] <- e0$tmin
    tmp <- k2m(pp[,'K1'],exp(pp[,'per1']),pp[,'e1'],Ms,Inc=pp[,'Inc1'],more=TRUE)
    out[N1:N2,6] <- tmp$mj
    out[N1:N2,10] <- tmp$a
    if(e0$Nsig>1){
        tmp <- k2m(pp[,'K2'],exp(pp[,'per2']),pp[,'e2'],Ms,Inc=pp[,'Inc2'],more=TRUE)
        out[N1:N2,7] <- tmp$mj
        out[N1:N2,11] <- tmp$a
    }
    if(e0$Nsig>2){
        tmp <- k2m(pp[,'K3'],exp(pp[,'per3']),pp[,'e3'],Ms,Inc=pp[,'Inc3'],more=TRUE)
        out[N1:N2,8] <- tmp$mj
        out[N1:N2,12] <- tmp$a
    }
    if(e0$Nsig>3){
        tmp <- k2m(pp[,'K4'],exp(pp[,'per4']),pp[,'e4'],Ms,Inc=pp[,'Inc4'],more=TRUE)
        out[N1:N2,9] <- tmp$mj
        out[N1:N2,13] <- tmp$a
    }
#t(replicate(c(e0$target,mstar,emstar1,emstar2,e0$tmin),NN))
    N1 <- N1+NN
    N2 <- N2+NN
}
fout <- paste0('companion_sc_post_v',version,'.txt')
cat(fout,'\n')
write.table(out,file=fout,quote=FALSE,row.names=FALSE)

