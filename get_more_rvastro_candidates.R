ss <- readLines('more_systems.txt')
out2 <- out1 <- c()
for(s in ss){
    if(grepl('\\/',s) & grepl('^HD|^GJ|^HIP|^TYC|^WASP|^XO|^GL|^BD',s)){
        pp <- gsub('.+_P|_acc.+','',s)
        if(pp!='NA'){
            ps <- as.numeric(unlist(strsplit(pp,split='d')))
            target <- gsub('\\/.+','',s)
            if(grepl('fix|astro',s)){
                out1 <- rbind(out1,c(target,pp))
            }
            if(any(ps>100)){
                out2 <- rbind(out2,c(target,pp))
            }
        }
    }
}
star1 <- unique(out1[,1])
star2 <- unique(out2[,1])
stars <- setdiff(star2,star1)

fs <- list.files('results/',recursive=TRUE)
targets <- unique(gsub('\\/.+','',fs))
systems <- setdiff(stars,targets)
missing  <- read.table('candidates1.txt',header=TRUE)[,1]
candidate <- setdiff(systems,missing)
write.table(candidate,file='candidates_more.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
