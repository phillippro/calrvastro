library(RColorBrewer)
source('mcmc_func.R')
options(scipen=9)
args <- commandArgs(trailingOnly=TRUE)
if(length(args>0)){
    f <- args
}else{
    f <- 'HD182488_fix1_relativityFALSE_Niter4000000_Ncores8_ofac2_Nset8_Esd1_transit0_P15344d36902_acc0.35_lnlmax'
}
star <- gsub('_.+','',f)
fin <- list.files(paste0('results/',star),pattern=paste0(f,'.+Robj'),full.name=TRUE)[1]
mstars <- emstars <- c()
yr2d <- 365.25
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
name.alt <- '15 Sge'
###checking whether files exist
fs <- fin
target0 <- targets0 <- targets <- stars <- star

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
out.phase <- list()
ins.all <- c()
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
#pp <- c('B','C','Ab','b','B','B','B','C')

ind.out <- c()
target <- gsub('.+\\/|_.+','',fin)
pch.type <- 20
pch.size <- 1
Ncol <- 4
Nrow <- 4
ff <- fin
cat('ff:\n',ff,'\n')
target.name <- stars
    load(ff,env=e0 <- new.env())
    Nastro <- 2
    astro <- astrometry <- 8
    if(!exists('basis')) basis <- 'natural'
    if(!exists('bases')) bases <- rep('natural',10)
    cat('\nstar=',target.name,'\n')
    cat('ins=',e0$ins,'\n')
    cat('ns=',unlist(lapply(e0$ins,function(i) nrow(e0$out[[i]][['RV']]))),'\n')
    cat('Tp=',unlist(lapply(e0$ins,function(i) max(e0$out[[i]][['RV']][,1])-min(e0$out[[i]][['RV']][,1]))),'\n')
    inss <- c(inss,e0$ins)
    offset <- TRUE
    Nsig <- e0$Nsig
pp <- c('b','c','d','e','f','g','h')[1:Nsig]
    if(target0=='GJP9066' | target0=='GJ300' | target0=='GJ880' | target0=='GJ251') Nsig <- 1
    trv.all <- e0$trv.all
    basis <- e0$basis
    ins <- e0$ins
    out <- e0$out
    tmin <- e0$tmin
    par.opt <- e0$out$par.stat[[paste0('sig',Nsig)]]['xopt',]
    mcmc <- e0$out$mcmc.opt
    lnBF5 <- c()
    lnBF3 <- c()
    for(kk in 1:Nsig){
        lnBF5 <- c(lnBF5,max(e0$out$mcmc.opt[[paste0('sig',kk)]][,'loglike'])-max(e0$out$mcmc.opt[[paste0('sig',kk-1)]][,'loglike'])-2.5*log(nrow(e0$out$all)))
        lnBF3 <- c(lnBF3,max(e0$out$mcmc.opt[[paste0('sig',kk)]][,'loglike'])-max(e0$out$mcmc.opt[[paste0('sig',kk-1)]][,'loglike'])-1.5*log(nrow(e0$out$all)))
    }
    e0$out$mcmc.all <- NULL
    e0$out$mcmc.opt <- NULL
    out.phase[[ff]] <- list(par.opt=par.opt,ins=ins,Nsig=Nsig,trv.all=trv.all,basis=basis,par.stat=e0$out$par.stat)
    rv.list <- list()
    for(i in ins){
        rv.list[[i]] <- out[[i]]$RV
    }
    out.phase[[ff]]$rv.list <- rv.list
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
source('generate_table.R')
lnBF3s[[ff]] <- lnBF3
lnBF5s[[ff]] <- lnBF5
inss <- unique(inss)
fout <- paste0(target,'_par.txt')
par.out[,1] <- gsub(' ','',par.out[,1])
cat(fout,'\n')
write.table(par.out,file=fout,quote=FALSE,row.names=FALSE)
