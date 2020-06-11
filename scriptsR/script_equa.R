args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Two arguments must be supplied: input_file output_file", call.=FALSE)
}

options(scipen = 999)

# loading data
resdata<-read.table(args[1],header=TRUE)

# Cheking lineraity between 20th percentile and valid unique pairs at 20kb window
r20kb<-cor(resdata[which(resdata$Window==20000),2],resdata[which(resdata$Window==20000),4])
test_20kb<-ifelse(r20kb<0.95,"error","ok")

# Cheking lineraity between 20th percentile and	window size at 20M valid unique pairs
r20M<-cor(resdata[which(resdata$SizeU<=20000000),3],resdata[which(resdata$SizeU<=20000000),4])
test_20M<-ifelse(r20M<0.95,"error","ok")

# Solving 20th percentile = f(valid unique pairs, window size) 
equamatrix<-matrix(c(as.numeric(resdata[1,3])*as.numeric(resdata[1,2]),resdata[1,3],resdata[1,2],1,as.numeric(resdata[2,3])*as.numeric(resdata[2,2]),resdata[2,3],resdata[2,2],1,as.numeric(resdata[4,3])*as.numeric(resdata[4,2]),resdata[4,3],resdata[4,2],1,as.numeric(resdata[5,3])*as.numeric(resdata[5,2]),resdata[5,3],resdata[5,2],1),byrow=TRUE,nrow=4)
equasol<-resdata[c(1,2,4,5),4]
coefequa<-solve(equamatrix,equasol)

# outputing result
output<-file(args[2])
if(test_20kb=="error" | test_20M=="error"){
		      writeLines("Error: data distribution is not linear.", output)
} else {
	writeLines(c("Resolution = ( 1000 - d - c * (Unique valid pairs)) / ( a * (Unique valid pairs) - b )",paste0("a = ",coefequa[1]),paste0("b = ",coefequa[2]), paste0("c = ",coefequa[3]), paste0("d = ",coefequa[4])), output)
}
close(output)