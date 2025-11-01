args  <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    type <- args[1]
}else{
#    type <- 'jupiter'
    type <- 'companion'
}
if(type=='BD'){
    f <- 'new_bd.txt'
    pat <- 'BD|bd'
    on <- 'BD'
}
if(type=='jupiter'){
    f <- 'new_jupiter.txt'
    pat <- 'jupiter|Jupiter'
    on <- 'jupiter'
}
if(type=='red'){
    f <- 'new_red.txt'
    pat <- 'binary|stellar|Binary'
    on <- 'red'
}

if(type=='companion'){
    f <- 'new_companion.txt'
    pat <- 'binary|stellar|jupiter|Jupiter|neptune|Neptune|saturn|Saturn|BD|bd|Binary'
    on <- 'companion'
}

tab <- readLines(f)
stars <- stars0 <- c()

####get stars
inds <- c()
for(j in 1:length(tab)){
    s  <- tab[j]
    if(grepl(pat,s)){
	ss <- unlist(strsplit(s,')| '))
	ss <- ss[ss!='']
	ss <- ss[grepl('HD|GJ|GL|HIP',ss)]
        ss <- ss[grepl('[0-9]',ss)]
        ss <- unique(ss)[1]
        ss <- gsub('.+\\/|_.+','',ss)
        if(!is.na(ss)){
            stars0 <- c(stars0,ss)
###        cat(ss,'\n')
            inds <- c(inds,j)
        }
    }
}
fout <- paste0(on,'_hosts.txt')
cat(fout,'\n')
write.table(stars0,file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)

inds <- c(inds,length(tab)+1)
####get files
r <- FALSE
fs <- c()
types <- c()
ll0 <- ll1 <- -1e6
for(j in 1:(length(inds)-1)){
###the lines where the file names are given
    k0  <- inds[j]+1
    k1  <- inds[j+1]-1

####initial setting
    f3 <- tab[k0:k1]
#    ff0 <- f3[grepl('fix\\d_rel|_astro',f3)]
    frv <- f3[grepl('_natural_',f3)]
    fdr2 <- f3[grepl('_astro|_relativity',f3) & !grepl('hg3|hg123',f3)]
    fhg3 <- f3[grepl('_hg3_',f3)]
    fhg3plus <- f3[grepl('_hg3\\+_',f3)]
    fhg123 <- f3[grepl('_hg123_',f3)]
    fhg123plus <- f3[grepl('_hg123\\+_',f3)]
    addf <- TRUE
    if(length(fhg123plus)>0){
#        cat('hg3!\n')
        types <- c(types,'hg123plus')
        ff <- fhg123plus
    }else if(length(fhg123)>0){
        types <- c(types,'hg123')
        ff <- fhg123
    }else if(length(fhg3plus)>0){
        types <- c(types,'hg3plus')
        ff <- fhg3plus
    }else if(length(fhg3)>0){
        types <- c(types,'hg3')
        ff <- fhg3
    }else if(length(fdr2)>0){
#        cat('dr2!\n')
        types <- c(types,'dr2')
        ff <- fdr2
    }else if(length(frv)>0){
#        cat(stars0[j],':rv!\n')
        types <- c(types,'rv')
        ff <- frv
    }else{
#        cat('none!\n')
        cat('No mcmc files for ',stars0[j],' is found\n')
        addf <- FALSE
    }
    if(addf){
        ind <- grep('lpspm',ff)
	if(length(ind)>0){
	     ff <- ff[ind]
	}
        lls <- as.numeric(gsub('.+lnlmax|\\.pdf.+|\\.pdf','',ff))
        f <- ff[which.max(lls)][1]
        stars <- c(stars,gsub('.+\\/|_.+','',f))
        fs <- c(fs,f)
    }
}
out <- cbind(stars,types,fs)
colnames(out) <- c('star','type','file')
##more info for targets
fout <- paste0(on,'_info.txt')
cat(fout,'\n')
write.table(out,file=fout,quote=FALSE,row.names=FALSE)

##only files
fout <- paste0(on,'_files.txt')
if(file.exists(fout)){
    fs0 <- read.table(fout)[,1]
    df <- setdiff(fs,fs0)
    if(length(df)>0){
        cat('The following files differ from previous ones:',df,'\n')
	ids <- gsub('\\/.+','',df)
	fout1 <- 'update_hosts.txt'
	cat(fout1,'\n')
	write.table(ids,file=fout1,row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
}
cat(fout,'\n')
write.table(fs,file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
