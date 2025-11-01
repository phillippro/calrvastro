#library(system)
fs <- read.table('results/misaligned_ma_dpsi50.txt')[,1]
for(f in fs){
    cmd  <-  paste0('cp ',f,' results/misaligned/')
    system(cmd)
}
