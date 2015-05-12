#Script GetTheModel.R
#Author: Chi Yang
#Date: 2014/11/13
#License: MIT licence
library(glmnet);
args <- commandArgs(trailingOnly = TRUE);
#R CMD BATCH --vanilla --slave '--args ../example/TrainSeq/PurR_training.matrix ../example/TrainSeq/PurR.elrm 0.5 0.00015' ./GetTheModel.R /dev/null";
trainingFile = args[1];
outputFile = args[2];
alpha = args[3];
s = args[4];

a = as.matrix(read.table(trainingFile, header=TRUE, sep="\t"));
colnames(a) = gsub("\\.", ":",colnames(a), perl=TRUE );
x = as.matrix(a[,2:dim(a)[2]]);
y = a[,1];
if(is.na(alpha)){
	alpha = 0.5;
}

if(is.na(s)){
	cv.fit= cv.glmnet(x, y, family="binomial",  nfold = 10, alpha = alpha);
	Coefficients = coef(cv.fit, s="lambda.min");
}else{
	s = as.numeric(s);
	cv.fit= cv.glmnet(x, y, family="binomial",  nfold = 10, alpha = alpha);
	Coefficients = coef(cv.fit, s=s);
}
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index,]
coef= as.matrix(Coefficients);
if(is.na(s)){
	coef = rbind(coef, cv.fit$lambda.min);
}else{
	coef = rbind(coef, s);
}
rownames(coef)[dim(coef)[1]] = "s";
write.table(coef, outputFile, quote=FALSE, sep="\t", col.names=FALSE);
