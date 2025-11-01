args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  f <- args
}else{
  f <- 'GJ551_lpspm_photovaryFALSE_relativityFALSE_Niter2200000_Ncores16_ofac2_Nset11_hg123+_transit0_P11d1901_acc0.06_sinI_lnlmax-3693'
}
f <- gsub('.+\\/','',f)
f <- gsub('.pdf','',f)
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
    dd <- '~/Documents/projects/dwarfs/calrvastro/pars/'
    fout <- paste0('~/Documents/projects/dwarfs/calrvastro/pars/',target,'.par')
    if(star=='HR8799'){
       if(!exists('coplanar')){
          coplanar <- FALSE
       } 
       if(Esd==1){
          lowe <- FALSE
       }else{
          lowe <- TRUE
       }
       if(!exists('resonance')){
          resonance <- FALSE
       } 
       fout <- paste0(dd,'HR8799_c',coplanar,'_l',lowe,'_r',resonance,'_n',Nsig,'.par')
    }
    cat(fout,'\n')
    write.table(t(t(par.opt)),file=fout,quote=FALSE,col.names=FALSE)
#    write.table(par.opt,file=fout,quote=FALSE,col.names=FALSE)
}
