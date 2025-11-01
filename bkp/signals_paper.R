library(magicaxis)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    file <- as.character(args[1])
}else{
#     file <- 'keppure_priormt_poly11_Ndata144_quantifyTRUE_1per1_Nw1_LHS1140_ind0_2planet_GPqp_Nsamp4000000_tem1_acc0.67_pretem1P93.1d24.5d_negLmax400'
#     file <- 'keppure_priormt_poly11_Ndata144_quantifyTRUE_1per1_Nw1_LHS1140_ind0_2planet_ARMA01_Nsamp10000000_tem1_acc0.43_pretem1P92.1d24.6d_negLmax398'
#     file <- 'keppure_priormt_poly11_Ndata144_quantifyTRUE_1per1_Nw1_LHS1140_ind0_2planet_ARMA00_Nsamp20000000_tem1_acc0.35_pretem1P92.3d24.6d_negLmax399'
#    file <- 'keppure_priormt_poly11_Ndata144_quantifyTRUE_1per1_Nw1_LHS1140_ind0_2planet_ARMA02_Nsamp10000000_tem1_acc0.48_pretem1P92.3d24.6d_negLmax399'
    file <- 'keppure_priormt_poly11_Ndata138_quantifyTRUE_1per1_Nw1_LHS1140_HARPS_ind7_3planet_ARMA01_Nsamp6000000_tem1_acc0.22_pretem0.02P3.8d24.6d91.4d_negLmax384'
}
if(!grepl('Robj',file)){
    file <- paste0(file,'.Robj')
}
f1 <- gsub('_ind.+','',file)
f1 <- gsub('_HARPS','',f1)
target <- gsub('.+Nw\\d_','',f1)
cat('target=',target,'\n')
folder <- paste0('/car-data/ffeng/dwarfs/',target,'output/')
if(!file.exists(folder)){
    folder <- paste0('../output/',target,'/')
}
solutions <- file1 <- paste0(folder,file)
load(file1)
source('paper_phase_fold.R')
