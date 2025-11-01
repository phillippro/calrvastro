ii <- 1:16
strs <- c(
    '$P$&day&Orbital period&',
    '$K$&\\ms&RV semi-amplitude&',
    '$e$&---&Eccentricity&',
    '$\\omega$$^{b}$&deg&Argument of periapsis&',
    '$I$&deg&Inclination&',
    '$\\Omega$&deg&Longitude of ascending node&',
    '$M_0$&deg&Mean anomaly at the reference epoch&',
                                        #  '\\hline'
    '$\\Delta \\alpha$&mas&$\\alpha$ offset&',
    '$\\Delta \\delta$&mas&$\\delta$ offset&',
    '$\\Delta \\varpi$&mas&$\\varpi$ offset&',
    '$\\Delta \\mu^*_\\alpha$&mas\\,yr$^{-1}$&$\\mu^*_\\alpha$ offset &',
    '$\\Delta \\mu_\\delta$&mas\\,yr$^{-1}$&$\\mu_\\delta$ offset&',
                                        #  \\hline
    '$P$&yr&Orbital period &',
    '$a$&au&Semi-major axis &',
    '$m_p$&$M_{\\rm Jup}$&Planet mass &',
    '$T_p-2400000$&JD&Periapsis epoch&')
strs2 <- c(
    'Log-Uniform&-1&16',
    'Uniform&$10^{-6}$&$10^{6}$',
    'Uniform&0&1',
    'Uniform&0&2$\\pi$',
    'CosI-Uniform&0&1',
    'Uniform&0&2$\\pi$',
    'Uniform&0&2$\\pi$',

    'Uniform&$-10^6$&$10^6$',
    'Uniform&$-10^6$&$10^6$',
    'Uniform&$-10^6$&$10^6$',
    'Uniform&$-10^6$&$10^6$',
    'Uniform&$-10^6$&$10^6$',

    '---&---&---',
    '---&---&---',
    '---&---&---',
    '---&---&---')

out <- c()
for(i in ii){
    i1 <- 3*i
    i2 <- 3*i+1
    i3 <- 3*i+2
    tmp <- strs[i]
    for(j in 1:nrow(dat)){
        dd <- dat[j,c(i1,i2,i3)]
        dd <- sprintf("%.2f", round(as.numeric(dd),2))
        tmp <- paste0(tmp,'$',dd[1],'_{-',dd[2],'}^{+',dd[3],'}$')
        if(j< nrow(dat)){
            tmp <- paste0(tmp,'&')
        }
    }
    tmp <- paste(tmp,'&',strs2[i])
    if(i==7 | i==12){
        tmp <- paste0(tmp,'\\\\\\hline')
    }else{
        tmp <- paste0(tmp,'\\\\')
    }
#    cat(tmp,'\n')
    out <- c(out,tmp)
}
cat(out,'\n\n')
if(star=='HD22049' | star=='HD209100') fout <- '/Users/ffeng/Documents/projects/dwarfs/paper/nearest_jupiters/v3/keppar.tex'
if(star=='UCAC4569-026385') fout <- '/Users/ffeng/Documents/projects/dwarfs/paper/BH/keppar.tex'
cat(fout,'\n')
write.table(out,file=fout,quote=FALSE,col.names=FALSE,row.names=FALSE)
