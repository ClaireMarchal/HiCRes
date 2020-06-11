args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
  stop("One argument must be supplied: input_file", call.=FALSE)
}

bgdata<-read.table(args[1],header=FALSE)
val<-quantile(bgdata[,4],0.2)
cat(val,"\n")

