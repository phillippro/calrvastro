tab <- readLines('new_bd.txt')
stars <- c()

####get stars
inds <- c()
for(j in 1:length(tab)){
    s  <- tab[j]
    if(grepl('BD',s)){
	ss <- unlist(strsplit(s,')| '))
	ss <- ss[ss!='']
	ss <- ss[grepl('HD|GJ|GL|HIP',ss)]
        ss <- ss[grepl('[0-9]',ss)]
        ss <- unique(ss)[1]
#        cat(ss,'\n')
        stars <- c(stars,ss)
        inds <- c(inds,j)
    }
}
fout <- 'BD_hosts.txt'
cat(fout,'\n')
write.table(stars,file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)

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
    ff <- tab[k0:k1]
    ff <- ff[!grepl('natural',ff)]
    ind <- grep('hg3',ff)
    if(length(ind)>0){
        ff <- ff[ind]
        types <- c(types,'hg3')
    }else{
        types <- c(types,'dr2')
    }
    ind <- grep('lpspm',ff)
    if(length(ind)>0) ff <- ff[ind]
    lls <- as.numeric(gsub('.+lnlmax|\\.pdf','',ff))
    fs <- c(fs,ff[which.max(lls)])
}
out <- cbind(stars,types,fs)
colnames(out) <- c('star','type','file')
##more info for targets
fout <- 'BD_info.txt'
cat(fout,'\n')
write.table(out,file=fout,quote=FALSE,row.names=FALSE)

##only files
fout <- 'BD_files.txt'
cat(fout,'\n')
write.table(fs,file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
