#!/usr/local/apps/R/3.6/3.6.1/lib64/R/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Two arguments must be supplied: input_file output_file", call.=FALSE)
}

options(scipen = 999)

# loading data
resdata<-read.table(args[1],header=TRUE)

# Cheking lineraity between square roots of 20th percentile and of valid unique pairs at 20kb window
r20kb<-cor((resdata[which(resdata$Window==20000),2])^(1/3),(resdata[which(resdata$Window==20000),4])^(1/3))
test_20kb<-ifelse(r20kb<0.98,"error","ok")

# Cheking lineraity between square roots of 20th percentile and of window size at 20M valid unique pairs
r20M<-cor((resdata[c(2,5,8),3])^(1/3),(resdata[c(2,5,8),4])^(1/3))
test_20M<-ifelse(r20M<0.98,"error","ok")

# Solving 20th percentile = f(valid unique pairs, window size) 
equamatrix<-matrix(c((as.numeric(resdata[1,3]))^(1/3)*(as.numeric(resdata[1,2]))^(1/3),(resdata[1,3])^(1/3),(resdata[1,2])^(1/3),1,(as.numeric(resdata[2,3]))^(1/3)*(as.numeric(resdata[2,2]))^(1/3),(resdata[2,3])^(1/3),(resdata[2,2])^(1/3),1,(as.numeric(resdata[4,3]))^(1/3)*(as.numeric(resdata[4,2]))^(1/3),(resdata[4,3])^(1/3),(resdata[4,2])^(1/3),1,(as.numeric(resdata[5,3]))^(1/3)*(as.numeric(resdata[5,2]))^(1/3),(resdata[5,3])^(1/3),(resdata[5,2])^(1/3),1),byrow=TRUE,nrow=4)
equasol<-(resdata[c(1,2,4,5),4])^(1/3)
coefequa<-solve(equamatrix,equasol)

# outputing result
output<-file(args[2])
if(test_20kb=="error" | test_20M=="error"){
		      writeLines("Error: non linear data.", output)
} else {
       writeLines(c("Resolution = (( (1000)^(1/3) - d - c * (Unique valid pairs)^(1/3) / ( a * (Unique valid pairs)^(1/3) + b)) ^ 3",paste0("a = ",coefequa[1]),paste0("b = ",coefequa[2]), paste0("c = ",coefequa[3]), paste0("d = ",coefequa[4])), output)
}
close(output)
