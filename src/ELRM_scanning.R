#Script ELRM_scanning.R
#Author: Chi Yang
#Date: 2014/11/13
#License: MIT licence
library(glmnet);
args <- commandArgs(trailingOnly = TRUE);
#R CMD BATCH --vanilla --slave '--args ../example/TrainSeq/PurR_training.matrix ../example/TestSeq/PurR_testing.matrix ../example/ScanningResults/PurR_results.txt 0.5 0.8' ./ELRM_scanning.R /dev/null";
trainingFile = args[1];
testingFile = args[2];
outputFile = args[3];
alpha = args[4];
threshold = args[5];
s = args[6];


a = as.matrix(read.table(trainingFile, header=TRUE, sep="\t"));
colnames(a) = gsub("\\.", ":",colnames(a), perl=TRUE );
x = as.matrix(a[,2:dim(a)[2]]);
y = a[,1];
if(is.na(alpha)){
	alpha = 0.5;
}

testData = read.table(testingFile, header=TRUE, sep="\t");
colnames(testData) = gsub("\\.", ":",colnames(testData), perl=TRUE );
testX = as.matrix(testData[,6:dim(testData)[2]]);

if(is.na(s)){
	cv.fit= cv.glmnet(x, y, family="binomial",  nfold = 10, alpha = alpha);
	predict1 = predict(cv.fit, testX, s="lambda.min", type="response")
}else{
	s = as.numeric(s);
	fit = glmnet(x, y, family="binomial", alpha = alpha);
	predict1 = predict(fit, testX, s=s, type="response")
}
testResults = cbind(testData[,1:5], predict1);
colnames(testResults) = c("Title","Start","End", "Strand","Motif" ,"Probability");

temp = sort(testResults[,6], index.return=T, method="sh",decreasing =T);
testResults = testResults[temp$ix,]


if(!is.na(threshold)){
	threshold = as.numeric(threshold);
	testResults = testResults[which(testResults[,"Probability"] >= threshold),];
	#write.table(testResults, file=outputFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE);
}
write.table(testResults, file=outputFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE);
