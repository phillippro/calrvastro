fs <- readLines('candidates.txt')
fs <- fs[grepl('pdf$',fs) & grepl('lnlmax',fs)]
stars <- gsub('\\/.+','',fs)
targets <- unique(stars)
#dir.out <- '/home/share/ffeng-group/kepler/'
dir.out <- '~/eb/'
###select the best solution
for(target in targets){
    ii <- which(target==stars)
    lls <- as.numeric(gsub('.+lnlmax|\\.pdf','',fs[ii]))
    jj <- which.max(lls)
    fpdf <- paste0('results/',fs[ii[jj]])
    if(file.exists(fpdf)){
        cmd1 <- paste('cp ',fpdf,paste0(dir.out,'results/'))
        cat(cmd1,'\n')
        system(cmd1)

        f1 <- gsub(paste0(target,'\\/'),'',fpdf)
        cmd2 <- paste('Rscript get_par.R',f1)
        cat(cmd2,'\n')
        system(cmd2)

        cmd3 <- paste0('cp pars/',target,'.par ',dir.out,'pars/')
        cat(cmd3,'\n')
        system(cmd3)
    }
}
