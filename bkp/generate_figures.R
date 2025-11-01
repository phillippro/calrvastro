###settings
colors <- c('black','blue','orange','cyan','brown','green','steelblue','purple','pink',rainbow(10))
fpdf <- paste0(fname,'.pdf')
cat(fpdf,'\n')
pdf(fpdf,16,16)
mar0 <- c(5,5,5,1)
size0 <- 1.5
par(mfrow=c(4,4),mar=mar0,cex.axis=size0,cex.lab=size0,oma=c(1,1,4,1))

##plot raw RV data
cat('\nPlot raw data and activity indices!\n')
if(out$Nrv>0){
    source('plot_data.R')
}

##plot phase curve
cat('\nPlot phase curve!\n')
if(Nsig>0){
    source('phase_plot.R')
}

###astrometry plot
cat('\nPlot orbital motion!\n')
if(length(out$astrometry)>0){
    source('ellipse_plot.R')
}

if(!is.null(out$timing)){
cat('\nPlot Timing model and data!\n')
    source('timing_plot.R')
}

##plot BFPs
cat('\nPlot BFP!\n')
if(out$Nrv>0){
    source('plot_BFP.R')
}

##trace plot
cat('\nTrace plot!\n')
source('trace_plot.R')

##pair plot
cat('\nPair plot!\n')
try(source('pairplot.R'),TRUE)
#try(source('pairplot_wopdf.R'),TRUE)

##plot MPs
cat('\nPlot MP!\n')
if(Nsig>0 & out$Nrv>0){
    source('MPplot.R')
}

##plot WPs
cat('\nPlot WP!\n')
if(Nsig>0 & out$Nrv>0){
   try(source('WPplot.R'),TRUE)
}

##report statistics
if(Nsig>0){
try(source('report.R'),TRUE)
}

##end
dev.off()
