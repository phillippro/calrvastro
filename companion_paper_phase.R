library(RColorBrewer)
source('mcmc_func.R')
options(scipen=9)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    target.type <- args[1]
}else{
#    target.type <- 'bd'
    target.type <- 'companion'
}
version <- 0
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

planets <- read.csv('../data/code/SigClassify/ranking_complex1.csv')
n0 <- gsub(' ','',as.character(planets[,'Name']))
n1 <- gsub(' ','',as.character(planets[,'ID']))
n2 <- gsub(' ','',as.character(planets[,'StarKnown']))

###UVES list
#name.alt <- c('HIP 73408 B','GJ 105 C','$\\mu$ Herculis Ab','GJ 777 b','GJ 779 B','$\\chi^1$ Orionis B','GJ 36 B','GJ 22 C')
fs <- read.table(paste0(target.type,'_files.txt'))[,1]
#fs <- fs[1:10]
#fs <- fs[grep('HD2151_',fs):length(fs)]
name.alt <- gsub('\\/.+','',fs)#ineed to be changed
fs <- gsub('\\.pdf|.+\\/','',fs)
###checking whether files exist
f1 <- c()
ind <- 1
#ind <- grep('HD127506',fs)
#ind <- grep('HD153557',fs)
for(j3 in ind:length(fs)){
#for(j3 in (ind-1):(ind+1)){
#for(j3 in ind){
    f <- fs[j3]
    star <- gsub('_.+','',f)
    if(!grepl('Robj',f)){
        ff <- list.files(paste0('results/',star),pattern=paste0(f,'.+Robj'),full.name=TRUE)
    }else{
        ff <- list.files(paste0('results/',star),pattern=f,full.name=TRUE)
    }
    cat('input file:',ff,'\n')
    f5 <- list.files(paste0('results/',star),pattern=paste0(f,'.+Robj'),full.name=TRUE)
    f3 <- list.files(paste0('../data/combined/',star),pattern='hg3$')
    if(length(f5)>0 & length(f3)>0){
        f1 <- c(f1,f5[1])
    }
}
fs <- f1
targets0 <- targets <- stars <- gsub('.+\\/|_.+','',fs)

####get catalog parameters for the targets
inds <- c()
for(s in targets){
    s1 <- gsub('GJ','GL',s)
    ind <- which(n0==s | n1==s | n2==s | n0==s1 | n1==s1 | n2==s1)
    if(length(ind)==0){
        if(s=='GJ9349') ind <- which(n0=='HIP54532')
        if(s=='GJ251') ind <- which(n0=='HD265866')
        if(s=='GJ880') ind <- which(n0=='HD216899')
        if(s=='GJ4070') ind <- which(n0=='HIP91699')
        if(s=='GJ2056') ind <- which(n0=='HIP34785')
        if(s=='GJ9066') ind <- which(n0=='GL83.1')
        if(s=='GJ480') ind <- which(n0=='HIP61706')
        if(s=='GJ3072') ind <- which(n0=='HIP4845')
    }
    inds <- c(inds,ind)
}
tab <- planets[inds,]

###modify star name
stars <- gsub('GL','GL ',stars)
stars <- gsub('GJ','GJ ',stars)
stars <- gsub('HIP','HIP ',stars)
stars <- gsub('HD','HD ',stars)
stars0 <- stars
####loading data
fobj <- paste0('hg3_phase',version,'.Robj')
#if(file.exists(fobj)){
if(FALSE){
    if(!exists('hg3_out.phase')){
        load(fobj)
    }
    reload <- FALSE
}else{
    reload <- TRUE
    out.phase <- list()
    ins.all <- c()
}
inss <- c()
lnBF3s <- c()
lnBF5s <- c()
#exception <- 'HIP22762'
#exception <- 'HIP19165'
exception <- ''
#exception <- c('HD42581','HIP80268')
#exception <- c('HIP117886')
prior.type <- 'mt'
time.unit <- 365.25
#colors <- c('HARPSpre'='black','HARPSpost'='darkgrey','UVES'='blue','PFSpre'='green','PFSpost'='darkgreen','KECK'='purple','APF'='cyan','AAT'='brown','HARPSN'='yellow','LC'='orange','VLC'='pink','CRIRES'='steelblue')
#colors <- c('HARPSpre'='black','HARPSpost'='brown','UVES'='blue','PFSpre'='green','PFSpost'='orange','KECK'='purple')
#sets <- ins.all <- sort(c('UVES','PFS','KECK','HARPSpre','HARPSpost','HARPN','AAT','APF','SOPHIE','CAR','HET','LICK','ELODIE'))
#sets <- ins.all <- sort(c('UVES','PFS','KECK','HARPSpre','HARPSpost','HARPN','AAT','APF','SOPHIE','CAR','HET','LICK'))
sets <- ins.all <- sort(c('UVES','PFS','KECK','HARPSpre','HARPSpost','HARPN','APF','SOPHIE','CAR','HET','LICK'))
nc <- length(sets)
colors <- c(brewer.pal(9, "Set1"),brewer.pal(n = 8, name = "Set2"))
#colors <- rainbow(nc)
#colors[c(5,6)] <- c('cyan','black')
names(colors) <- sets
colors['UVES'] <- 'blue'
colors['LICK'] <- 'green'
colors['HET'] <- 'tan'
colors['PFS'] <- 'red'
colors['HARPSpre'] <- 'black'
colors['HARPSpost'] <- 'gold'
colors['SOPHIE'] <- 'cyan'

###
par.out <- c()
table.out <- c()
star.report <- c()
Ndig0 <- Ndig <- 2
#pp <- c('b','c','d','e','f','g','h','i')
#pp <- c('C','Ab','B','b','B','B','B','C')
pp <- c('B','C','Ab','b','B','B','B','C')

ind.out <- c()

fout <- paste0(target.type,'_phase_coverage',version,'.pdf')
cat('output pdf:\n',fout,'\n')
pch.type <- 20
pch.size <- 1
Ncol <- 4
Nrow <- 4
pdf(fout,height=4*Nrow,width=4*Ncol)
#pic.lay <- rbind(rep(1,4),matrix((1:(Nrow*Ncol*2))+1,ncol=Ncol,byrow=FALSE))
pic.lay <- rbind(rep(1,Ncol),matrix((1:(Nrow*Ncol*2))+1,ncol=Ncol,byrow=FALSE))
layout(pic.lay,heights=as.numeric(c(0.2,replicate(6,c(2,1)))))
par(oma=c(4,4,4,4),mar=rep(0,4))
np <- 0
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#legend('top',inset=-1,legend=ins.all,col=colors[ins.all],pch=pch.type,xpd=NA,horiz=TRUE,bty='n',cex=2)
ind1 <- 1:7
ind2 <- 8:length(ins.all)
legend('top',inset=-1.5,legend=ins.all[ind1],col=colors[ins.all[ind1]],pch=pch.type,xpd=NA,horiz=TRUE,bty='n',cex=2)
legend('top',inset=-0.5,legend=ins.all[ind2],col=colors[ins.all[ind2]],pch=pch.type,xpd=NA,horiz=TRUE,bty='n',cex=2)
#dev.off()
#stop()
for(k5 in 1:length(fs)){
    ff <- fs[k5]
    cat('\nff:',ff,'\n')
    cat(k5,'\n')
    target0 <- targets0[k5]
    target.name <- stars[k5]
#    if(any(target0==exception)){
    if(reload){
        load(ff,env=e0 <- new.env())
        Nastro <- 2
        astro <- astrometry <- 8
        if(!exists('basis')) basis <- 'natural'
        if(!exists('bases')) bases <- rep('natural',10)
        cat('star=',target.name,'\n')
        cat('ins=',e0$ins,'\n')
        cat('ns=',unlist(lapply(e0$ins,function(i) nrow(e0$out[[i]][['RV']]))),'\n')
        cat('Tp=',unlist(lapply(e0$ins,function(i) max(e0$out[[i]][['RV']][,1])-min(e0$out[[i]][['RV']][,1]))),'\n')
        inss <- c(inss,e0$ins)
        offset <- TRUE
        Nsig <- e0$Nsig
        if(target0=='GJP9066' | target0=='GJ300' | target0=='GJ880' | target0=='GJ251') Nsig <- 1
        trv.all <- e0$trv.all
        basis <- e0$basis
        ins <- e0$ins
        out <- e0$out
        tmin <- e0$tmin
#        par.opt <- e0$out$par.stat[[paste0('sig',Nsig)]]['xopt',]
        par.opt <- e0$par.opt
        mcmc <- e0$out$mcmc.opt
        lnBF5 <- c()
        lnBF3 <- c()
        for(kk in 1:Nsig){
            base.model <- paste0('sig',kk-1)
            target.model <- paste0('sig',kk)
            nn <- names(e0$out$mcmc.out)
            if(!any(nn==base.model) & !any(nn==target.model)){
                lnBF5 <- c(lnBF5,NA)
                lnBF3 <- c(lnBF3,NA)
            }else{
                lnBF5 <- c(lnBF5,max(e0$out$mcmc.opt[[paste0('sig',kk)]][,'loglike'])-max(e0$out$mcmc.opt[[paste0('sig',kk-1)]][,'loglike'])-2.5*log(nrow(e0$out$all)))
                lnBF3 <- c(lnBF3,max(e0$out$mcmc.opt[[paste0('sig',kk)]][,'loglike'])-max(e0$out$mcmc.opt[[paste0('sig',kk-1)]][,'loglike'])-1.5*log(nrow(e0$out$all)))
            }
        }
        e0$out$mcmc.all <- NULL
        e0$out$mcmc.opt <- NULL
        out.phase[[ff]] <- list(par.opt=par.opt,ins=ins,Nsig=Nsig,trv.all=trv.all,basis=basis,par.stat=e0$out$par.stat)
        rv.list <- list()
        for(i in ins){
            rv.list[[i]] <- out[[i]]$RV
        }
        out.phase[[ff]]$rv.list <- rv.list
    }else{
        par.opt <- out.phase[[ff]]$par.opt
        ins <- as.character(out.phase[[ff]]$ins)
        rv.list <- out.phase[[ff]]$rv.list
        Nsig <- out.phase[[ff]]$Nsig
        trv.all <- out.phase[[ff]]$trv.all
        basis <- out.phase[[ff]]$basis
    }
    ins.name <- ins
    if(any(ins=='HARPS')){
        th <- rv.list$HARPS[,1]
        if(min(th)%%2400000>57174.5) ins.name[ins.name=='HARPS'] <- 'HARPSpost'
        if(max(th)%%2400000<=57174.5) ins.name[ins.name=='HARPS'] <- 'HARPSpre'
    }
    if(any(ins=='PFS') & FALSE){
        th <- rv.list$PFS[,1]
        cat('range(th)=',range(th),'\n')
        ins.name[ins=='PFS'] <- 'PFSpre'
        if(min(th)%%2400000>58157) ins.name[ins=='PFS'] <- 'PFSpost'
        if(max(th)%%2400000<=58157) ins.name[ins=='PFS'] <- 'PFSpre'
        cat('ins.name=',ins.name,'\n')
    }
    if(any(target0==exception)) ins.all <- c(ins.all,ins.name)
#    if(Nsig>0 & target0!='HIP86214') source('plot_phase.R')
#    source('plot_phase.R')
    source('generate_table.R')
    lnBF3s[[ff]] <- lnBF3
    lnBF5s[[ff]] <- lnBF5
}
inss <- unique(inss)
dev.off()
if(length(fs)>1) source('latex_table.R')
#colnames(par.out) <- c('Star','P','Plow','Pup','K','Klow','Kup','e','elow','eup','omega','omega.low','omega.up','Omega','Omega.low','Omega.up','mp','mp.lower','mp.upper')
fout <- paste0(target.type,'_par',version,'.txt')
par.out <- cbind(stars,mstars,emstars,par.out)
par.out[,1] <- gsub(' ','',par.out[,1])
cat(fout,'\n')
write.table(par.out,file=fout,quote=FALSE,row.names=FALSE)

