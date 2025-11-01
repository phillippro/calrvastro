p <- system('pbsnodes -a',intern = TRUE)
ind3 <- grep('properties',p)
ind2 <- ind3-1
ind1 <- ind3-2
ind0 <- ind3-3
node <- p[ind0]
state <- gsub('state|=| ','',p[ind1])
Np <- gsub('np|=| ','',p[ind2])
out <- cbind(node,state,Np)
ind <- which(state=='free')
write.table(out[ind,1][-1],'valide_nodes.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(out,'all_node_info.txt',quote=FALSE,row.names=FALSE,col.names=c('node','state','Nprocessors'))