#tab <- read.table('stars_with_hg123_simbad.csv',sep='|',header=TRUE)
tab <- read.csv('stars_with_hg123_simbad.csv',sep='|')
ind <- grep('^M',tab[,'spec..type'])
f1 <- 'stars_Mtype.txt'
cat(f1,'\n')
write.table(gsub(' ','',tab[ind,2]),file=f1,quote=FALSE,row.names=FALSE,col.names=FALSE)

f2 <- 'stars_nearby.txt'
cat(f2,'\n')
ind <- which(as.numeric(tab[,'plx'])>100)#10pc sample
write.table(gsub(' ','',tab[ind,2]),file=f2,quote=FALSE,row.names=FALSE,col.names=FALSE)
