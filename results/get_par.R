args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  f <- args
}else{
  f <- 'HD42581_natural_astro8_Niter10000000_Ncores4_ofac2_Nset5_Esd1_transit0_P39195_acc0.039_lnlmax'
}
if(grepl('txt',f)){
  ff <- read.table(paste0('results/',f))[,1]
}else{
  ff <- f
}
ff <- gsub('\\+','\\\\+',ff)
for(f in ff){
    star <- gsub('_.+','',f)
    dir1 <- paste0('~/Documents/projects/dwarfs/calrvastro/results/',star)
    dir2 <- paste0('~/Documents/projects/dwarfs/rvastro/results/',star)
    dir3 <- paste0('~/Documents/projects/dwarfs/red3/results/',star)
    ff <- list.files(dir1,pattern=paste0(f,'.+Robj'),full.name=TRUE)
    if(length(ff)==0){
    ff <- list.files(dir2,pattern=paste0(f,'.+Robj'),full.name=TRUE)
    }
    if(length(ff)==0){
    ff <- list.files(dir3,pattern=paste0(f,'.+Robj'),full.name=TRUE)
    }
    cat('\nload:',ff,'\n')
    cat(f,'\n')
    load(ff)
    fout <- paste0('pars/',star,'.par')
    cat(fout,'\n')
    write.table(par.opt,file=fout,quote=FALSE,col.names=FALSE)
}
