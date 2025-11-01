fs <- readLines('candidates.txt')
fs <- fs[grep('pdf$',fs)]
stars <- gsub('\\/.+','',fs)
targets <- unique(gsub('\\/.+','',fs))
###select the best solution
for(target in targets){
    ii <- which(target==stars)
    lls <- as.numeric(gsub('.+lnlmax|\\.pdf','',fs[ii]))
    jj <- which.max(lls)
    fpdf <- paste0('results/',fs[ii[jj]])
    if(file.exists(fpdf)){
        cmd <- paste('rsync -azP',fpdf,'/home/share/ffeng-group/kepler/results/')
        cat(cmd,'\n')
        system(cmd)
        f1 <- gsub(paste0(target,'\\/'),'',fpdf)

        cmd <- paste('Rscript get_par.R',f1)
        cat(cmd,'\n')
        system(cmd)

        cmd <- paste0('mv pars/',target,'.par',' /home/share/ffeng-group/kepler/results/pars/')
        cat(cmd,'\n')
        system(cmd)
        stop()
    }
}
