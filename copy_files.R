fs <- read.table('all_companions.txt',header=FALSE)[,1]
f1 <- gsub('.+\\/','',fs)
fo <- list.files('/home/share/companion/robj_file')
fp <- list.files('/home/share/companion/pdf_file')
ii <- match(f1,fo)
jj <- (1:length(f1))[which(is.na(ii))]
for(i in jj){
   cmd1 <- paste('cp',fs[i],'/home/share/companion/robj_file/')
   cat(cmd1)
   system(cmd1)

   f <- gsub('Robj','pdf',fs[i])
   cmd2 <- paste('cp',f,'/home/share/companion/pdf_file/')
   system(cmd2)
}
