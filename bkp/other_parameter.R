library(magicaxis)
library(graphics)
library(RColorBrewer)
source('mcmc_func.R')
options(scipen = 0)
library("plotrix")
source('OrbitFunction.R')
source('timing_function.R')
source('general_function.R')
source('sofa_function.R')
set.seed(9999)
Ndig <- 2
modifypar <- function(par){
    Nsig <- length(grep('^per',names(par)))
    ps <- exp(par[paste0('per',1:Nsig)])
    ind <- which(ps<1000)
    if(length(ind)>0){
        par[paste0('K',ind)] <- 0
    }
    par
}
m22 <- read.table('../data/combined/mamajek22.txt',header=TRUE)
ms <- as.numeric(m22[,'Msun'])
mg <- as.numeric(m22[,'M_G'])
ind <- which(!is.na(ms) & !is.na(mg))
mrl.m22 <- approxfun(ms[ind],mg[ind])
mlow.m22 <- min(ms[ind])
mup.m22 <- max(ms[ind])
eta0 <- rep(NA,10)

#Nmc <- 0
Nmc <- 1000
#Nmc <- 100
comptype <- 'companion'
version <- 1
set.seed(100)
ff <- read.table('1companion_files.txt')[,1]
#ff <- read.table('2companion_files.txt')[,1]

#nn <- c('GL676A','HD145675','HD149806','HD62364','HD115404A','HD42581')
#nn <- 'G3425'
nn <- 'UCAC4569-026385'
#nn <- c('HD209100','HD22049')
#nn <- c('HD222237')
#nn <- c('HR8799')
stars <- gsub('\\/.+','',ff)
ii <- match(nn,stars)
#ii <- as.integer(sapply(nn,function(i) ))
ff <- ff[ii]
stars <- stars[ii]
host <- read.table('host_par.txt',header=TRUE)
rds <- mpacs <- c()

#stars <- gsub('\\/.+','',ff)
#stars <- c('HD42581','HD39060','HD182488','HD215257')
#stars <- c('HD39060','HD182488','HD215257','HD39091')
#stars <- c('HD39060','HD182488')
#stars <- c('HD7449')
Nstar <- length(stars)
etaH <- etaG <- 1
###modifiy input file names
#Nsig <- length(fs)
planets <- read.csv('../data/code/SigClassify/ranking_complex2.csv')
n0 <- gsub(' ','',as.character(planets[,'Name']))
n1 <- gsub(' ','',as.character(planets[,'ID']))
n2 <- gsub(' ','',as.character(planets[,'StarKnown']))
nfig <- 0
npdf <- 0

####pdf figure
for(j in 1:Nstar){
    ss <- c()
    star <- stars[j]
    star.name <- gsub('HD','HD ',star)
    star.name <- gsub('HIP','HIP ',star.name)
    star.name <- gsub('GJ|GL','GJ ',star.name)
###check whether to show
    ind <- which(host[,1]==gsub('GL','GJ',star))
    plx <- host[ind,'plxG']
#    out.orbit <- rbind(out.orbit,c(pp[pp[,1]==star,],host[ind,]))
    ii <- grep(paste0('^',star,'\\/'),ff)
    fobj <- gsub('pdf','Robj',ff[ii])
    fs <- paste0('results/',fobj)
    cat('load ',fs,'\n')
    ind <- which(stars==star)
#    if(!exists('out')) load(fs)
    load(fs)
    mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
    par.name <- names(par.opt)
#    if(star=='HD209100') pars <- read.table('par_hd209100.txt',header=TRUE)
#    if(star=='HD22049') pars <- read.table('par_hd22049.txt',header=TRUE)
#    if(star=='HD222237') pars <- read.table('par_hd222237.txt',header=TRUE)
#    par.min <- sprintf("%.0f", pars[2,])
#    par.max <- sprintf("%.0f", pars[3,])
#   names(par.min) <- names(par.max) <- colnames(pars)

    Np <- length(grep('^per',names(par.opt)))
    source('initial_condition.R')
     for(pn in par.name){
         ps <- quantile(mc[,pn],c(0.16,0.5,0.84))
        if(!grepl('^per|^K|^e\\d|^omega|^Omega|^Inc|^Mo|dra|ddec|dplx|dpmra|dpmdec',pn)){
            nn <- pn
            if(grepl('^b_',pn)){
                ins <- gsub('b_','',pn)
                ins[ins=='APF'] <- 'APFp'
                ins[ins=='APFj'] <- 'APFh'
#                nn <- paste0("$\\gamma^{\\rm ",ins,'}$')
                unit <- 'm\\,s$^{-1}$'
                if(star=="UCAC4569-026385"){
                    ps <- ps/1e3#m/s -> km/s
                    unit <- 'km\\,s$^{-1}$'
#                    nn <- paste0("$b^{\\rm ",ins,'}$')
                    nn <- paste0("$b^{\\rm ",toupper(ins),'}$')
                    par.min[pn] <- par.min[pn]/1e3
                    par.max[pn] <- par.max[pn]/1e3
                }
                meaning <- paste('RV offset for',toupper(ins))
            }
            if(grepl('^J_',pn)){
                ins <- gsub('J_','',pn)
                ins[ins=='APF'] <- 'APFp'
                ins[ins=='APFj'] <- 'APFh'
#                nn <- paste0('$J^{\\rm ',ins,'}$')
                nn <- paste0('$J^{\\rm ',toupper(ins),'}$')
                unit <- 'm\\,s$^{-1}$'
                if(star=="UCAC4569-026385"){
                    ps <- ps/1e3#m/s -> km/s
                    unit <- 'km\\,s$^{-1}$'
                    par.min[pn] <- par.min[pn]/1e3
                    par.max[pn] <- par.max[pn]/1e3
                }
                meaning <- paste('RV jitter for',toupper(ins))
            }
            if(grepl('^logJ_',pn)){
                ins <- gsub('logJ_','',pn)
                if(grepl('hip2',ins)) ins <- 'hip'
#                nn <- paste0('ln$J_{\\rm ',ins,'}$')
                nn <- paste0('ln$J_{\\rm ',toupper(ins),'}$')
                unit <- '---'
#                meaning <- paste('Logarithmic jitter for',ins)
                meaning <- paste('Logarithmic jitter for',toupper(ins))
            }
            if(grepl('^w\\d_',pn)){
                order <- gsub('^w|_.+','',pn)
                ins <- gsub('.+_','',pn)
                ins[ins=='APF'] <- 'APFp'
                ins[ins=='APFj'] <- 'APFh'
#                nn <- paste0('$w_{',order,'}^{\\rm ',ins,'}$')
                nn <- paste0('$w_{',order,'}^{\\rm ',toupper(ins),'}$')
                unit <- '---'
                q <- max(as.numeric(gsub('w|_.+','',par.name[grepl('^w\\d_',par.name) & grepl(paste0('_',ins),par.name)])))
#                meaning <- paste0('Amplitude of component ',order,' of MA(',q,') for ',ins)
                meaning <- paste0('Amplitude of component ',order,' of MA(',q,') for ',toupper(ins))
            }
            if(grepl('^beta_',pn)){
                ins <- gsub('.+_','',pn)
                ins[ins=='APF'] <- 'APFp'
                ins[ins=='APFj'] <- 'APFh'
#                nn <- paste0('ln$\\tau^{\\rm ',ins,'}$')
                nn <- paste0('ln$\\tau^{\\rm ',toupper(ins),'}$')
                unit <- '---'
                q <- max(as.numeric(gsub('w|_.+','',par.name[grepl('^w\\d_',par.name) & grepl(paste0('_',ins),par.name)])))
#                meaning <- paste0('Logarithmic time scale of MA(',q,') for ',ins)
                meaning <- paste0('Logarithmic time scale of MA(',q,') for ',toupper(ins))
            }
            pp <- sprintf("%.1f",c(ps[2],ps[2]-ps[1],ps[3]-ps[2]))
            pm <- sprintf("%.1f",par.min[pn])
            pm[pm=='-1000000'] <- '$-10^6$'
            pa <- sprintf("%.1f",par.max[pn])
            pa[pa=='1000000'] <- '$10^6$'
            tmp <- paste0(nn,'&',unit,'&',meaning,'&$',pp[1],'_{-',pp[2],'}^{+',pp[3],'}$','&Uniform&',pm,'&',pa,'\\\\')
            ss <- c(ss,tmp)
        }
    }
#    fout <- paste0(star,'_other_par.tex')
    fout <- paste0('/Users/ffeng/Documents/projects/dwarfs/paper/nearest_jupiters/v3/',star,'_other_par.tex')
    cat(fout,'\n')
    write.table(ss,file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
}
