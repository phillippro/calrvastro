ahg3 <- read.table('hg3_par0.txt',header=TRUE)
bhg3 <- read.table('hg3_test0.txt',header=TRUE)
ahg2 <- read.table('hg2_par0.txt',header=TRUE)
bhg2 <- read.table('hg2_test0.txt',header=TRUE)
ahgca3 <- read.table('hgca3_par0.txt',header=TRUE)
bhgca3 <- read.table('hgca3_test0.txt',header=TRUE)
ahgca2 <- read.table('hgca2_par0.txt',header=TRUE)
bhgca2 <- read.table('hgca2_test0.txt',header=TRUE)

eta.hg2 <- bhg2[,'mpJ1.sd']/ahg2[,'mpJ1.sd']
eta.hg3 <- bhg3[,'mpJ1.sd']/ahg3[,'mpJ1.sd']
eta.hgca2 <- bhgca2[,'mpJ1.sd']/ahgca2[,'mpJ1.sd']
eta.hgca3 <- bhgca3[,'mpJ1.sd']/ahgca3[,'mpJ1.sd']
cat('hg2 stellar mass error contribution:',(1-eta.hg2)*100,'\n')
cat('hg3 stellar mass error contribution:',(1-eta.hg3)*100,'\n')
cat('hgca2 stellar mass error contribution:',(1-eta.hgca2)*100,'\n')
cat('hgca3 stellar mass error contribution:',(1-eta.hgca3)*100,'\n')

eta.hg2 <- bhg2[,'per1.sd']/ahg2[,'per1.sd']
eta.hg3 <- bhg3[,'per1.sd']/ahg3[,'per1.sd']
eta.hgca2 <- bhgca2[,'per1.sd']/ahgca2[,'per1.sd']
eta.hgca3 <- bhgca3[,'per1.sd']/ahgca3[,'per1.sd']
cat('hg2 stellar mass error contribution:',(1-eta.hg2)*100,'\n')
cat('hg3 stellar mass error contribution:',(1-eta.hg3)*100,'\n')
cat('hgca2 stellar mass error contribution:',(1-eta.hgca2)*100,'\n')
cat('hgca3 stellar mass error contribution:',(1-eta.hgca3)*100,'\n')

