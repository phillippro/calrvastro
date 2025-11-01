library(RColorBrewer)
source('mcmc_func.R')
options(scipen=9)
version <- 0
yr2d <- 365.25
Ns <- 1e4
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
#   star <- args[1]
   comptype <- args[1]
}else{
#   star <- 'GJ9163'
   comptype <- 'companion'
}
cat('comptype=',comptype,'\n')

value.type <- c('ms','mq')#mean+sd and map+quantile
#value.type <- c('ms')
type <- 'format'
xup <- 'x99per'
xlow <- 'x1per'
opt <- 'med'
#opt <- 'map'
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
fs <- paste0('results/',gsub('pdf','Robj',read.table(paste0(comptype,'_files.txt'))[,1]))
#fs <- fs[grep('HD75393',fs)]
fs <- fs[670:680]
#fs <- fs[500:]
Nf <- length(fs)
name.alt <- gsub('\\/.+','',fs)#ineed to be changed
ind.rm <- c()

#fs <- fs[grep(paste0(star,'_'),fs)][1]
###checking whether files exist
targets0 <- targets <- stars <- gsub('.+\\/|_.+','',fs)

Ncores <- min(48,length(fs))
#Ncores <- 4
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}

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
fobj <- paste0(comptype,'_phase',version,'.Robj')
#if(file.exists(fobj)){
reload <- TRUE
out.phase <- list()
ins.all <- c()
inss <- c()
#exception <- 'HIP22762'
#exception <- 'HIP19165'
exception <- ''
#exception <- c('HD42581','HIP80268')
#exception <- c('HIP117886')
prior.type <- 'mt'
time.unit <- 365.25

###
table.out <- c()
star.report <- c()
Ndig0 <- Ndig <- 2
pp <- c('B','C','Ab','b','B','B','B','C')

ind.out <- c()

par.all <- foreach(k5 = 1:length(fs), .combine='list', .errorhandling='pass') %dopar% {
result <- tryCatch({
    object_that_doesnt_exist[i]}, 
     warning = function(war) {
         return('a warning')}, 
     error = function(err) {
         return('an error')}, 
     finally = {
         return('other things')
     }) # END tryCatch

  return(result) 
mstars <- emstars <- par.out <- c()
#for(k5 in 1:length(fs)){
    ff <- fs[k5]
#    cat('ff:\n',ff,'\n')
    target0 <- targets0[k5]
    target.name <- stars[k5]
     if(!file.exists(ff)){
            cmd <- paste0('scp tdlffb@data.hpc.sjtu.edu.cn:Documents/projects/dwarfs/rvastro/',ff,' results/',target0,'/')
	    cat(cmd,'\n')
	    system(cmd)
     }
#    if(any(target0==exception)){
     if(file.exists(ff)){
    load(ff,env=e0 <- new.env())
#    pp <- any(names(e0$out)=='astrometry')	
    nastro <- e0$out$Nastro
    if(is.null(nastro)) nastro <- 2
    cat('\nwhether astrometry is used in ',ff,'?',nastro>0,'!\n')
    if(nastro>0){
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
        par.opt <- e0$out$par.stat[[paste0('sig',Nsig)]]['xopt',]
        mcmc <- e0$out$mcmc.opt
        e0$out$mcmc.all <- NULL
#        e0$out$mcmc.opt <- NULL
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
    y <- try(source("generate_table.R",local=TRUE),TRUE)
#    y <- try(source("generate_table.R",local=TRUE),FALSE)
    if(class(y)=='try-error'){
       cat('y=',y,';\ntarget with error:',target0,'\n')
       stop()
    }
#    source('generate_table.R',local=TRUE)
#    source('generate_table.R')
inss <- unique(inss)
#colnames(par.out) <- c('Star','P','Plow','Pup','K','Klow','Kup','e','elow','eup','omega','omega.low','omega.up','Omega','Omega.low','Omega.up','mp','mp.lower','mp.upper')
par.out <- cbind(target0,mstars,emstars,par.out)
    par.out[,1] <- gsub(' ','',par.out[,1])
    par.out
}
}
}

par.out <- par.all[[1]]
if(length(par.all)>1){
for(j in 2:length(par.all)){
par <- par.all[[j]][1,]
 par.name <- names(par)	
   out.name <- colnames(par.out)
   ind.par <- match(out.name,par.name)
   ind.out <- match(par.name,out.name)
   if(!any(is.na(ind.out))){
     par.out <- rbind(par.out,par[ind.par])
   }else{
     par.name <- names(par)
     out.name <- colnames(par.out)
     ind.out <- match(par.name,out.name)
     if(any(is.na(ind.out))){
         array.add <- array(NA,dim=c(nrow(par.out),length(which(is.na(ind.out)))))
         colnames(array.add) <- par.name[is.na(ind.out)]
         par.out <- cbind(par.out,array.add)
         ind.par <- match(colnames(par.out),names(par))
         par.out <- rbind(par.out,par[ind.par])
     }
  }
}
}

fout <- paste0(comptype,'_par',version,'_',opt,'.txt')
cat(fout,'\n')
write.table(par.out,file=fout,quote=FALSE,row.names=FALSE)
