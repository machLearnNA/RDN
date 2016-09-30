#' @name getRDN
#' @title Compute the RDN applicability domain for increasing neighbourhoods.
#' @description Function that automates the characterization of the RDN applicability domain, by computing the neighbourhood thresholds
#' for an increasing number of nearest neighbours, and for each set of thresholds collecting the respective in-domain accuracy
#' and outlying test instances.
#'
#' @param performance Data frame of dimensions (N,1) with 1 for correctly predicted instances in testSet, and 0 otherwise.
#' @param trainingSet Data frame with the training descriptors (raw) used to calculate \code{NNthreshold}.
#' @param testSet Data frame with descriptors (raw) of the new instances to be tested against \code{trainingSet}.
#' @param agreementInput A data frame (N,1) with the agreement measure from M ensemble models, for each instance in \code{trainingSet}.
#' @param STDInput A data frame (N,1) with the ensemble standard deviation for each instance in  \code{trainingSet}
#' @param stepLimit Number of domain expansion iterations to be computed (optional)
#' @param initialcompression Integer setting the iteration limit up to (but not including) which the threshold distances are compressed to a third of their original values; from initialcompression + 1, threshold values get decompressed to half of their original distance values.
#' @param decompression  Integer setting the starting iteration at which threshold distances get fully decompressed (Euclidean distances used as is)
#'
#' @return a matrix called \code{resultSummary} which stores a column of outlying instance count and another with the in-domain accuracy, for each calculation step.
#'
#' @details This function calls \code{getEDmatrix}, \code{getThreshold} and \code{TestInTrain}  to compute the applicability domain accuracy and in-domain test instances at iteratively larger threshold distances. It will do so by calling \code{getEDmatrix} once and storing it under \code{Dij.sort} (implicit variable),
#' and then call \code{getThreshold} and \code{TestInTrain} for an iteratively larger number of \code{k} nearest neighbours; at each such step the in-domain accuracy and the outlying instance count will be stored in \code{resultSummary}.
#'
#' @examples
#' library(randomForest)
#' library(mlbench)
#' data("BreastCancer")
#' # remove first Col (IDs)
#' bcdata <- BreastCancer[,c(2:ncol(BreastCancer))]
#' bcdata <- bcdata[complete.cases(bcdata),]
#'
#' #sample 70% for training
#' trainID <- sample(1:nrow(bcdata),round(nrow(bcdata)*0.7))
#'
#' #QSAR model; gather external predictions
#' # ntree was set to was purposefully to challenge the RDN calculation with more error.
#' m<-randomForest(Class~.,data=bcdata[c(trainID),], kepp.forest=TRUE, ntree=1, norm.votes=TRUE)
#' p<-predict(m, bcdata[-trainID,], predict.all=TRUE, type="class")
#' test_pred <- p$aggregate #class pred external
#' performance <- as.double(test_pred==bcdata[-trainID,"Class"])
#' performance <- data.frame(performance)
#'
#' ensembleClass.t <- data.frame(bcdata$Class[trainID])
#' ensembleProb.t <- data.frame(bcdata$Class[trainID])
#'
#' # train ensemble for AD calculation
#' for (i in 2:11){
#'   #sampling
#'   samp <- sample(trainID, round(length(trainID)*0.8))
#'   # train ensemble
#'   m<-randomForest(Class~.,data=bcdata[samp,], kepp.forest=TRUE, ntree=100, norm.votes=TRUE)
#'   pred<-predict(m, bcdata[trainID,], predict.all=TRUE, type="Prob")
#'   pred_class<-predict(m, bcdata[trainID,], predict.all=TRUE, type="class")
#'   P_AD <- pred$aggregate[,1] #prob train
#'  class_AD <-pred_class$aggregate #class train
#'  ensembleClass.t[,i]=class_AD
#'  ensembleProb.t[,i]=P_AD
#'
#'}
#'
#' # compute agreement for TRAIN
#' agree <- data.frame(trainID)
#' for (i in 1:length(trainID)){
#'   agree[i,2] <- sum(ensembleClass.t[i,-1]==toString(ensembleClass.t[i,1]))/10
#' }
#' agree <- data.frame(agree[,-1])
#'
#' #compute std for TRAIN
#' std <- apply(ensembleProb.t[,-1],1,sd)
#' std <- data.frame(std)
#'
#' # Prepare descriptors to be passed into getRDN
#' train <- data.matrix(bcdata[trainID,1:ncol(bcdata)-1])
#' test <- data.matrix(bcdata[-trainID,1:ncol(bcdata)-1])
#'
#' # Compute RDN; this will take care of scaling trainingSet and testSet internally before any use.
#' resultSummary <- getRDN(performance=performance, trainingSet=train,
#'                         testSet=test, agreementInput=agree, STDInput=std)
#' # The results saved in resultSummary show a decreasing overall quality of predictions as
#' # AD gets expanded (i.e. instances out of AD decrease).
#' #> resultSummary
#' #NNout ACC in AD
#' #[1,]    108 0.9793814
#' #[2,]     97 0.9722222
#' #[3,]     92 0.9734513
#' #[4,]     89 0.9741379
#' #[5,]     88 0.9743590
#' #[6,]     88 0.9743590
#' #[7,]     88 0.9743590
#' #[8,]     87 0.9745763
#' #[9,]     84 0.9752066
#' #[10,]     83 0.9754098
#' #[11,]     80 0.9760000
#' #[12,]     78 0.9763780
#' #[13,]     77 0.9765625
#' #[14,]     75 0.9769231
#' #[15,]     73 0.9772727
#' #[16,]     71 0.9776119
#' #[17,]     67 0.9710145
#' #[18,]     67 0.9710145
#' #[19,]     67 0.9710145
#' #[20,]     64 0.9645390
#' #[21,]     64 0.9645390
#' #[22,]     64 0.9645390
#' #[23,]     64 0.9645390
#' #[24,]     61 0.9652778
#' #[25,]     60 0.9655172
#' #[26,]     59 0.9657534
#' #[27,]     58 0.9659864
#' #[28,]     58 0.9659864
#' #[29,]     58 0.9659864
#' #[30,]     58 0.9659864
#' #[31,]     14 0.9476440
#' #[32,]     14 0.9476440
#' #[33,]     14 0.9476440
#' #[34,]     14 0.9476440
#' #[35,]     12 0.9430052
#' #[36,]     12 0.9430052
#' #[37,]     12 0.9430052
#' #[38,]     12 0.9430052
#' #[39,]     12 0.9430052
#' #[40,]     12 0.9430052
#' #[41,]      0 0.9365854
#' #[42,]      0 0.9365854
#' #...          ...
#' #[65,]      0 0.9365854
#'
#'
#' @export
getRDN <- function(performance, trainingSet, testSet, agreementInput, STDInput, stepLimit = 65, initialcompression = 31, decompression = 41){
  # Compute de distance matrix ONCE
  Dij.sort <<- getEDmatrix(trainingSet, trainingSet)

  listOut <- c()
  AccIn <- c()
  trainingSet <- scale(trainingSet, center = mins, scale = maxs - mins)
  testSet <- scale(testSet, center = mins, scale = maxs - mins)  # mins and maxs were implicitly created when getEDmatrix is called above
  for (i in 1:stepLimit){
    NNthreshold <- getThreshold(agreementInput,STDInput,trainingSet,i)
    if (i < initialcompression){
      NNthreshold   # D/3
    }
    if (i >= initialcompression && i < decompression){
      NNthreshold <- NNthreshold*3/2 #D/2
    }
    if (i >= decompression){
      NNthreshold <- NNthreshold*3   #D
    }
    resultOutput <- TestInTrain(performance,testSet,trainingSet,NNthreshold)
    TE.NNcount <- resultOutput[[1]]
    outADcount <- resultOutput[[2]]
    Acc <- resultOutput[[3]]
    listOut <- append(listOut, outADcount)
    AccIn <- append(AccIn, Acc)
  }
  resultSummary <- cbind(as.matrix(listOut),as.matrix(AccIn))
  colnames(resultSummary) <- list("#NNout","ACC in AD")
  return(resultSummary)
}
