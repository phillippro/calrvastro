fs <- list.files('WD1202−232',full.name=TRUE)
fs1 <- gsub('WD1202−232','WD1202-232',fs)
for(j in 1:length(fs)){
  f <- fs[j]
  cmd <- paste('mv',f,fs1[j])
  cat(cmd,'\n')
  system(cmd)
}
