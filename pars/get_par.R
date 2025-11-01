library(foreach)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  fs <- args
}else{
  fs <- '2MASSJ12451043+1217401_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1e+05_jitteruniform_240715_Nset1_hg123+_Nsig1_P740_Esd1_astro5TRUE_acc7.8_rvcbarycentric_priorf1_lnpm'
}
for(f in fs){
if(grepl('txt',f)){
  ff <- read.table(f)[,1]
}else{
  ff <- f
}
#ff <- gsub('\\+','\\\\+',ff)
for(f in ff){
    star <- gsub('_.+','',f)
    if(grepl('fix|astro|lpspm',f)){
        dir <- paste0('/malibu/ffeng/astro_output/',star)
        if(!file.exists(dir)){
            dir <- paste0('../results/',star)
        }
        astroPar <- TRUE
    }else{
        dir <- paste0('../../red3/results/',star)
        if(!file.exists(dir)){
            dir <- paste0('/malibu/ffeng/dwarfs_output/',star)
        }
        astroPar <- FALSE
    }
    f <- gsub('\\+','\\\\+',f)
    ff <- list.files(dir,pattern=paste0(f,'.+Robj'),full.name=TRUE)
    if(length(ff)>0){
    cat('\nload:',ff,'\n')
    cat(f,'\n')
    load(ff)
    fout <- paste0(star,'.par')
    cat(fout,'\n')
    write.table(par.opt,file=fout,quote=FALSE,col.names=FALSE)
    }
}
}
