library(RColorBrewer)
source('mcmc_func.R')
options(scipen=9)
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

planets <- read.csv('../data/code/SigClassify/ranking_complex1.csv')
n0 <- gsub(' ','',as.character(planets[,'Name']))
n1 <- gsub(' ','',as.character(planets[,'ID']))
n2 <- gsub(' ','',as.character(planets[,'StarKnown']))

###UVES list
#name.alt <- c('HIP 73408 B','GJ 105 C','$\\mu$ Herculis Ab','GJ 777 b','GJ 779 B','$\\chi^1$ Orionis B','GJ 36 B','GJ 22 C')
fs <- read.table('all_companions.txt',header=TRUE)[,1]
stars <- gsub('_.+','',fs)
cat('length(fs)=',length(fs),'\n')
###checking whether files exist
n <- c('per','K','e','omega','Mo','Inc','Omega')
ns <- as.vector(outer(n,c(1:4),'paste0'))
Nsamp <- 1e4
NN <- Nsamp+1
Nrow <- NN*length(fs)
out <- data.frame(array(NA,dim=c(Nrow,length(ns)+5)))
colnames(out) <- c('star','mass','emass1','emass2','tmin',ns)
N1 <- 1
N2 <- NN
for(j3 in 1:length(fs)){
    f <- fs[j3]
    cat('read:',f,'\n')
    star <- gsub('_.+','',f)
    load(f,env=e0 <- new.env())
    mc <- e0$out$mcmc.opt[[paste0('sig',e0$Nsig)]]
    ll <- mc[,'loglike']
    lp <- mc[,'logpost']
    mc <- mc[,1:(ncol(mc)-2)]
    par.opt <- mc[which.max(ll),]
    inds <- sort(sample(1:nrow(mc),Nsamp))
    mc <- mc[inds,]
    pp <- rbind(par.opt,mc)
    nn <- intersect(ns,colnames(mc))
    ind1 <- match(nn,colnames(out))
    ind2 <- match(nn,colnames(mc))
    out[N1:N2,ind1] <- pp[,ind2]
    mstar <- e0$out$astrometry[1,'mass']
    emstar1 <- e0$out$astrometry[1,'mass.lower']
    emstar2 <- e0$out$astrometry[1,'mass.upper']
    out[N1:N2,1] <- e0$target
    out[N1:N2,2] <- mstar
    out[N1:N2,3] <- emstar1
    out[N1:N2,4] <- emstar2
    out[N1:N2,5] <- e0$tmin
#t(replicate(c(e0$target,mstar,emstar1,emstar2,e0$tmin),NN))
    N1 <- N1+NN
    N2 <- N2+NN
}
fout <- 'all_companion_post.txt'
cat(fout,'\n')
write.table(out,file=fout,quote=FALSE,row.names=FALSE)

