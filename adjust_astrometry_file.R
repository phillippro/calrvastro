if(!exists('astro.type')){
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    star <- args[1]
    type <- args[2]
}else{
    star <- 'HD42581'
    type <- 'abs'
}
}else{
    star <- target
    type <- astro.type
}
dir.in <- paste0('../data/combined/',star,'/')
dir.out <- '../data/combined/bkp/'
fs <- list.files(dir.in,full.name=TRUE)
fs <- fs[grepl('rel|astro',fs)]
if(type=='abs'){
    if(!any(grepl('hipgaia',fs))){
        cmd <- paste0('mv ',dir.out,star,'_hipgaia.astro ',dir.in)
        cat(cmd,'\n')
        system(cmd)
    }
    if(any(grepl('\\.rel',fs))){
        ff <- fs[grepl('\\.rel',fs)]
        for(f in fs){
            cmd <- paste('mv',f,dir.out)
            cat(cmd,'\n')
            system(cmd)
        }
    }
}

if(grepl('rel',type)){
    if(!any(grepl(type,fs))){
        cmd <- paste0('mv ',dir.out,star,'_astrometry.',type,' ',dir.in)
        cat(cmd,'\n')
        system(cmd)
    }
    ind <- grep('\\.rel',fs)
    ind1 <- which(!grepl(type,fs[ind]))
    if(length(ind1)>0){
        for(f in fs[ind[ind1]]){
            cmd <- paste('mv',f,dir.out)
            cat(cmd,'\n')
            system(cmd)
        }
    }

    if(any(grepl('\\.astro',fs))){
        ff <- fs[grep('\\.astro',fs)]
        for(f in ff){
            cmd <- paste('mv',f,dir.out)
            cat(cmd,'\n')
            system(cmd)
        }
    }
}

if(grepl('both',type)){
    ip <- as.integer(gsub('both','',type))
    n <- paste0('rel',ip)
    if(!any(grepl('hipgaia',fs))){
        cmd <- paste0('mv ',dir.out,star,'_hipgaia.astro ',dir.in)
        cat(cmd,'\n')
        system(cmd)
    }

    index <- grep(n,fs)
    if(!any(grepl(n,fs))){
        cmd <- paste0('mv ',dir.out,star,'_astrometry.',n,' ',dir.in)
        cat(cmd,'\n')
        system(cmd)
    }

    index <- which(grepl('rel',fs) & !grepl(n,fs))
    if(length(index)>0){
        for(f in fs[index]){
            cmd <- paste('mv',f,dir.out)
            cat(cmd,'\n')
            system(cmd)
        }
    }
}
fs <- list.files(dir.in,full.name=TRUE)
cat('remain files:',paste(fs,collapse='\n'),'\n')
