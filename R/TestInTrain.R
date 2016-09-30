#' @name TestInTrain
#' @title Determine if new queries fall inside training neighbourhood
#' @description Function that searches for training instances whose neighbourhood radius covers new (test) instances (i.e., is equal to, or larger than the distance to new instances). Test instances falling in at least a training neighbourhood are taken to calculate in-domain accuracy.
#' @param performance Data frame of dimensions (N,1) with 1 for correctly predicted instances in testSet, and 0 otherwise.
#' @param testSet Data frame with the scaled descritptors of the new instances to be tested against trainingSet. scaling should be done with scale(testSet, center = mins, scale = maxs - mins), where mins and maxs are inherited from previously calling \code{getEDmatrix}
#' @param trainingSet Data frame with the scaled descritptors of the training instances used to calculate NNthreshold. scaling should be done with scale(trainingSet, center = mins, scale = maxs - mins), where mins and maxs are inherited from previously calling \code{getEDmatrix}
#' @param NNthreshold Data frame with threshold distances output by getThreshold.
#'
#' @return \code{TestInTrain} outputs a list with three elements:
#' \itemize{
#'           \item \code{TE.NNcount} A vector with the number of training nearest neighbours found within the threshold distance of each test instances
#'           \item \code{outAD} Number of test instances outside all threshold distances (i.e., outside of the applicability domain)
#'           \item \code{Acc} Within-applicability domain accuracy; sum of in-domain performances divided over count of in-domain test instances
#'           }
#' @details This function computes the Euclidean distance between each test instance and each of all training instances, and determines how many training instances are within threshold distance of each test instance. Test instances that have 1 or more training neighbours are considered
#' in-domain and only those will be used for computing in-applicability domain accuracy (Acc). This function calls the output of getThreshold.
#' @examples
#' library(randomForest)
#' library(mlbench)
#' data("BreastCancer")
#' #remove ID col
#' bcdata <- BreastCancer[,c(2:ncol(BreastCancer))]
#' bcdata <- bcdata[complete.cases(bcdata),]
#'
#' #sample 70% for training
#' trainID <- sample(1:nrow(bcdata),round(nrow(bcdata)*0.7))
#'
#' #QSAR model; gather external predictions
#' # ntree was set to 1 to purposefully to challenge the RDN calculation with more error.
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
#' # Prepare descriptors to be passed into getRDN; dataframe should be numerical not in levels.
#' train <- data.matrix(bcdata[trainID,1:ncol(bcdata)-1])
#' test <- data.matrix(bcdata[-trainID,1:ncol(bcdata)-1])
#'
#' ## calculate sorted distance matrix to training neighbours
#' Dij.sort <- getEDmatrix(train, train)
#'
#' # data needs to be scaled as all functions rely on the Euclidean distance;
#' # maxs and mins have been implicitly
#' # produced from the previous line.
#' train <- scale(train, center = mins, scale = maxs - mins)
#' test <- scale(test, center = mins, scale = maxs - mins)
#'
#' ## Compute the coverage threshold corresponding to k=3 nearest neighbours
#' NNthreshold <- getThreshold(agree, std, train, k=3)
#'
#' # Place test into train neighbouhoods
#' resultOutput <- TestInTrain(performance=performance, testSet=test, trainingSet=train, NNthreshold)
#'
#' # Colect results
#' test.NNcount <- resultOutput[[1]]
#' outADcount <- resultOutput[[2]]
#' Acc <- resultOutput[[3]]
#'
#' @export
#testset sould be scaled previous to being passed to this functions using mins and maxs are implicitly created at getEDmatrix
TestInTrain <- function(performance,testSet,trainingSet,NNthreshold){
  cat("Computing neighborhood for the set k value", "\n")
  TE.NNcount <- c()
  Acc = 0
  for (i in 1:nrow(testSet)){
    NNcount = 0
    for (j in 1:nrow(trainingSet)){
      ED = sqrt(sum((testSet[i,]-trainingSet[j,])^2))         # testset need to be scaled as per trainset
      if (ED <= NNthreshold[j]){
        NNcount = NNcount + 1
      }
    }
    if (NNcount > 0){
      Acc = Acc + performance[i,1]
    }
    TE.NNcount = append(TE.NNcount, NNcount)
  }
  outAD = sum(TE.NNcount==0)
  Acc = Acc /(nrow(testSet)-outAD)
  resultOutput <- list(TE.NNcount, outAD, Acc)
  return(resultOutput)
}
