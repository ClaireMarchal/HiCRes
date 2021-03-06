#!/usr/local/apps/R/3.6/3.6.1/lib64/R/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=9) {
  stop("9 arguments must be supplied: a b c d %map %valid %cis %cis_far output_file", call.=FALSE)
}

options(scipen = 999)

# Importing data
coefequa<-c(as.numeric(args[1]),as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]))
percmap<-as.numeric(args[5])
percval<-as.numeric(args[6])
perccis<-as.numeric(args[7])
perccisfar<-as.numeric(args[8])

# Generating a table
preseqout<-data.frame(Uniq_Sequenced_Pairs=matrix(rep(NA,999),ncol=1))
n=0; for(i in seq(1000000, 1000000000, 1000000)){preseqout[n,1]<-i; n=n+1}


# Calculating intermediate read pairs number
preseqout$Uniq_Mapped_Pairs<-preseqout$Uniq_Sequenced_Pairs*percmap/100
preseqout$Predicted_Unique_Valid_Pairs<-preseqout$Uniq_Mapped_Pairs*percval/100
preseqout$Predicted_Unique_Valid_Cis_Pairs<-preseqout$Predicted_Unique_Valid_Pairs*perccis/100
preseqout$Predicted_Unique_Valid_Cis_Far_Pairs<-preseqout$Predicted_Unique_Valid_Pairs*perccisfar/100

# Calculating predicted resolution
preseqout$Predicted_resolution_all_interactions<-(((1000)^(1/3) - coefequa[4] - coefequa[3] * (preseqout$Predicted_Unique_Valid_Pairs)^(1/3)) / (coefequa[1] * (preseqout$Predicted_Unique_Valid_Pairs)^(1/3) + coefequa[2]))^3
preseqout$Predicted_resolution_cis_interactions_only<-(((1000)^(1/3) - coefequa[4] - coefequa[3] * (preseqout$Predicted_Unique_Valid_Cis_Pairs)^(1/3)) / (coefequa[1] * (preseqout$Predicted_Unique_Valid_Cis_Pairs)^(1/3) + coefequa[2]))^3
preseqout$Predicted_resolution_cis_far_interactions_only<-(((1000)^(1/3) - coefequa[4] - coefequa[3] * (preseqout$Predicted_Unique_Valid_Cis_Far_Pairs)^(1/3)) / (coefequa[1] * (preseqout$Predicted_Unique_Valid_Cis_Far_Pairs)^(1/3) + coefequa[2]))^3


# Removing bad lines:
preseqout[preseqout<=0]<-NA
preseqout<-preseqout[complete.cases(preseqout),]

# Exporting results in a text format
write.table(preseqout,args[9],sep="\t",quote=FALSE,row.names=FALSE)

# Exporting results in graph format 
png(gsub(".txt",".png",args[9]),width=2000,height=2000,pointsize=20,res=200)
plot(data.frame(preseqout$Uniq_Sequenced_Pairs/1000000,preseqout$Predicted_resolution_all_interactions),ylim=c(min(preseqout$Predicted_resolution_all_interactions),max(preseqout$Predicted_resolution_cis_far_interactions_only)),xlim=c(0,1000),xlab="Uniquely sequenced read pairs (M)",ylab="Predicted resolution (bp)",main="Resolution prediction from sequencing depth",type="l",col="RED",log="y"); par(new=TRUE); plot(data.frame(preseqout$Uniq_Sequenced_Pairs/1000000,preseqout$Predicted_resolution_cis_interactions_only),ylim=c(min(preseqout$Predicted_resolution_all_interactions),max(preseqout$Predicted_resolution_cis_far_interactions_only)),xlim=c(0,1000),xlab="",ylab="",type="l",col="BLUE",log="y"); par(new=TRUE); plot(data.frame(preseqout$Uniq_Sequenced_Pairs/1000000,preseqout$Predicted_resolution_cis_far_interactions_only),ylim=c(min(preseqout$Predicted_resolution_all_interactions),max(preseqout$Predicted_resolution_cis_far_interactions_only)),xlim=c(0,1000),xlab="",ylab="",type="l",col="BLUE",lty=2,log="y"); legend("topright",legend=c("All interactions","Cis interactions only","Cis far (>10kb) interactions only"),col=c("RED","BLUE","BLUE"),lty=c(1,1,2))
dev.off()
