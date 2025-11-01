###This routine choose the activity data to determine the GP parameters
if(conservative){
    Prot <- Protation
    activity.name <- rot.name
}else{
    Prot <- mean(Protations)
    Prot.min <- max(Prot-3*sd(Protations),0)
    Prot.max <- Prot+3*sd(Protations)
    n1 <- sapply(1:nrow(Pact),function(k) any(abs(Pact[k,]-Prot)/Prot<0.1))
    activity.name <- rownames(Pact)[n1]
}
if(any(activity.name=='Halpha') & Halpha){
    ind.act <- ihh
    activity <- hh
    Dactivity <- dhh
}else if(any(activity.name=='Sindex') & Sindex){
    ind.act <- iss
    activity <- ss
    Dactivity <- dss
}else if(Halpha){
    ind.act <- ihh
    activity <- hh
    Dactivity <- dhh
}else if(Sindex){
    ind.act <- iss
    activity <- ss
    Dactivity <- dss
}else{
    ind.act <- ipp
    activity <- pp
    Dactivity <- dpp
}
