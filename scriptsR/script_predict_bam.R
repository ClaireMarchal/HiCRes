args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("5 arguments must be supplied: a b c d output_file", call.=FALSE)
}

options(scipen = 999)

# Importing values
coefequa<-c(as.numeric(args[1]),as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]))

# Generating a table
preseqout<-data.frame(Valid_Reads=matrix(rep(NA,999),ncol=1))
n=0; for(i in seq(1000000, 1000000000, 1000000)){preseqout[n,1]<-i; n=n+1}

# Calculating predicted resolution
preseqout$Predicted_Resolution<-(((1000)^(1/3) - coefequa[4] - coefequa[3] * (preseqout$Valid_Reads)^(1/3)) / (coefequa[1] * (preseqout$Valid_Reads)^(1/3) + coefequa[2]))^3

# Removing bad lines:
preseqout[preseqout<=0]<-NA
preseqout<-preseqout[complete.cases(preseqout),]

# Exporting results in a text format
write.table(preseqout,args[5],sep="\t",quote=FALSE,row.names=FALSE)

# Exporting results in graph format 
png(gsub(".txt",".png",args[5]),width=2000,height=2000,pointsize=20,res=200)
plot(data.frame(preseqout$Valid_Reads/1000000,preseqout$Predicted_Resolution),xlim=c(0,1000),xlab="Valid read pairs (M)",ylab="Predicted resolution (bp)",main="Resolution prediction from valid interactions",type="l",log="y"); dev.off()
