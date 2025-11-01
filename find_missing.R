stars <- read.table('companion_hosts.txt')[,1]
targets <- gsub('_.+','',list.files('companion_plot',pattern='pdf'))
miss <- setdiff(stars,targets)
miss <- miss[grepl('HD\\d|HIP\\d|')]
print(miss)

