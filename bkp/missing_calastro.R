tab <- read.table('all_targets.txt')#[,1]
stars <- unique(tab[,1])
missing <- c()
missing.wRobj <- c()
missing.woRobj <- c()
missing.woHG <- c()
fpars <- list.files('optpar',pattern='par')
for(star in stars){
    fpar <- paste0('pars/',star,'.par')
    ind <- which(tab[,1]==star)
# &
    if(!any(tab[ind,2]=='calrvastro')){
        ind1 <- ind[which(tab[ind,2]=='red3')]
        ind2 <- ind[which(tab[ind,2]=='rvastro')]
        if(length(ind2)>0){
            f <- tab[ind2,3]
            folder <- tab[ind2,2]
        }else{
            f <- tab[ind1,3]
            folder <- tab[ind1,2]
        }
        p <- as.numeric(unlist(strsplit(gsub('.+_P|_acc.+','',f),split='d')))

        if(any(p>100,na.rm=TRUE)){
            missing <- c(missing,star)
            f123 <- list.files(paste0('../data/combined/',star),pattern='hg123')
            if(length(f123)==0){
                missing.woHG <- c(missing.woHG,star)
            }
            if(!file.exists(fpar)){
                fobj <- paste0('../',folder,'/results/',star,gsub('pdf','Robj',f))
                if(file.exists(fobj)){
                    ff <- paste0('../',folder,'/results/',star,gsub('\\.pdf','',f))
                    cmd <- paste('Rscript get_par.R',ff)
                    cat(cmd,'\n')
                    missing.wRobj <- c(missing.wRobj,star)
                }else if(length(ind.par)>0){
                    fin <- paste0('optpar/',fpars[ind.par])
                    par <- t(read.csv(fin)[1,-1])
                    write.table(par,file=paste0('pars/',star,'.par'),quote=FALSE,col.names=FALSE)
                    missing.wRobj <- c(missing.wRobj,star)
                }else{
                    missing.woRobj <- c(missing.woRobj,star)
                }
            }
        }
    }
#
#            if(!file.exists(fobj) & length(ind1)>0){
#                fobj <- paste0('../red3/results/',star,gsub('pdf','Robj',tab[ind1,3]))
#                if(!file.exists(fobj))
#            }
#        }


#        fs <- list.files(paste0('results/',star),pattern='_hg123_.+pdf')
#        if(length(fs)==0){
#            f1 <- list.files(paste0('../rvastro/results/',star),pattern='_astro_|_hg3\\+_|_hg3_|_fix1_|astro')
#            if(length(f1)==0){
#            }
#        }
#    }
#    cat('star:',star,'\n')
#    cat(length(fs),'\n')
#    if(length(fs)==0){
#          missing <- c(missing,star)
#    }
}
f1 <- 'missing.txt'
f2 <- 'missing_wRobj.txt'
f3 <- 'missing_woRobj.txt'
f4 <- 'missing_woHG.txt'
cat(f1,'\n')
cat(f2,'\n')
cat(f3,'\n')
cat(f4,'\n')
write.table(missing,file=f1,quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(missing.wRobj,file=f2,quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(missing.woRobj,file=f3,quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(missing.woHG,file=f4,quote=FALSE,row.names=FALSE,col.names=FALSE)
