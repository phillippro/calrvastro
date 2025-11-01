##this file is called by calmag.R
bands <- c('J','H','K','L')
pars <- read.csv('par_companion.csv',header=TRUE)
#nn <- as.character(outer(paste(bands,'band;'),paste(ages,'Gyr'),'paste0'))
nn <- c('$\\rho$ [mas]',as.character(outer(bands,ages,'paste0')))
#ii <- grep('L band',nn)
#nn[ii] <- gsub('L band','Luminosity',nn[ii])
#tex <- c('\\startlongtable\\begin{longtable}','\\centering','\\caption{Flux ratio to host stars for the companions identified in this work. }','\\label{tab:contrast}',
tex <- c('\\startlongtable',paste0('\\begin{deluxetable*}{',paste(rep('l',length(bands)*length(ages)+1),collapse=''),'r}'),'\\tablecaption{Mean separation and logarithmic flux ratio of companions to host stars. \\label{tab:contrast}}',
         paste('\\tablehead{\\colhead{Companion}&',paste(paste0('\\colhead{',nn,'}'),collapse='&'),'}'),'\\startdata')
cn <- as.character(outer(paste0('Contrast',bands),ages,'paste0'))
for(j in 1:nrow(out)){
    ind <- which(pars[,1]==out[j,1] | pars[,1]==gsub('GL','GJ',out[j,1]))
    jj <- ind[which.min(abs(pars[ind,'per']-out[j,'per.opt']))]
    comp <- pars[jj,'companion']
    cat(comp,'\n')
    if(comp=='') stop()
    tex <- c(tex,paste(comp,'&',round(out[j,'rho.opt']),'&',paste(format(log10(out[j,cn]),digits=2,nsmall=1),collapse='&'),'\\\\'))
}
tex <- c(tex,'\\enddata','\\tablecomments{The contrasts in different bands are denoted by "BandAge" where "Band" is H, J, or K and "Age" is 0.1, 1, or 10 Gyrs. The contrast or flux ratio is in a base-10 log scale. This table only lists the mean values of the separation $\\rho$ and the contrast in different bands for each companion. The separation of the companion to the host vary over time and the exact values for different epochs are shown in Fig. \\ref{fig:orbit}. The luminosity is either given by the Gaia DR2 or by the mass-luminosity relation provided by \\cite{eker15}. The J, H, and K magntiudes of each star are obtained from the Simbad database \\citep{wenger00}. The J, H, and K magnitudes as well as the luminosity for each companion  are derived by the ATMO 2020 cooling models \\citep{phillips20}. }','\\end{deluxetable*}')
fout <- '../paper/bd/v2/table_contrast.tex'
cat(fout,'\n')
write.table(tex,file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
