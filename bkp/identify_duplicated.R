tab <- read.table('companion_files.txt',header=TRUE)[,1]
stars <- gsub('\\/.+','',tab)
ii <- which(duplicated(stars))
print(stars[ii])