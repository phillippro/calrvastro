library(RColorBrewer)
source('mcmc_func.R')
options(scipen=9)
show.digit <- function(x,Ndig,type='format'){
    if(type=='format'){
        format(round(x,Ndig),nsmall=Ndig)
    }else if(type=='signif'){
        signif(x,Ndig)
    }
}
Ndig0 <- Ndig <- 2
inss <- c()
lnBF3s <- c()
lnBF5s <- c()
par.out <- c()
table.out <- c()
star.report <- c()
yr2d <- 365.25
type <- 'format'
#xup <- 'xup'
xup <- 'xplus.1sig'
#xlow <- 'xlow'
xlow <- 'xminus.1sig'
Ns <- 1e5
fs <- read.table('all_companions.txt',header=TRUE)[,1]
fs <- fs[!grepl('HD128620|HD128621',fs)]
name.alt <- gsub('\\/.+','',fs)#ineed to be changed
targets <- stars <- gsub('_.+','',gsub('.+\\/','',fs))
stars <- gsub('GL','GL ',stars)
stars <- gsub('GJ','GJ ',stars)
stars <- gsub('HIP','HIP ',stars)
stars <- gsub('HD','HD ',stars)
opt <- 'med'
pp <- c('B','C','Ab','b','B','B','B','C')
value.type <- c('ms','mq')#mean+sd and map+quantile

mstars <- emstars <- c()
for(k5 in 1:length(fs)){
#for(k5 in 1076:length(fs)){
#for(k5 in 1:10){
    ff <- fs[k5]
    cat('\nff:',ff,'\n')
    cat(k5,'\n')
    target0 <- target <- targets[k5]
    target.name <- stars[k5]
    load(ff,env=e0 <- new.env())
    tmin <- e0$tmin
    out <- e0$out
    Nsig <- e0$Nsig
    mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
    source('generate_table.R')
}
inss <- unique(inss)
#if(length(fs)>1) source('latex_table.R')
#colnames(par.out) <- c('Star','P','Plow','Pup','K','Klow','Kup','e','elow','eup','omega','omega.low','omega.up','Omega','Omega.low','Omega.up','mp','mp.lower','mp.upper')
fout <- paste0('all_par.txt')
par.out <- cbind(stars,mstars,emstars,par.out)
par.out[,1] <- gsub(' ','',par.out[,1])
cat(fout,'\n')
write.table(par.out,file=fout,quote=FALSE,row.names=FALSE)

