

set.seed(1)
data <- matrix(sample(c(0,1),100*100,replace=TRUE,prob=c(0.9,0.1)),nrow=100,ncol=100)
data[1:10,1:10] <- 1 # BC1
data[11:20,11:20] <- 1 # BC2
data[21:30,21:30] <- 1 # BC3
data <- data[sample(1:nrow(data),nrow(data)),sample(1:ncol(data),ncol(data))]



result <- bibit2(data,minr=5,minc=5,noise=0.2,extend_columns = "naive",extend_mincol=1,extend_limitcol=1)
result

result <- bibit2(data,minr=5,minc=5,noise=0.2,extend_columns = "recursive",extend_mincol=1,extend_limitcol=1)
result