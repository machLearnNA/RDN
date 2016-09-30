#' @name getThreshold
#' @title Calculate the training neighbourhood threshold distances
#' @description Function that calculates the neighbourhood threshold distance associated to each training instance. The individual thresholds are computed from the average distances of the nearest neighbours within the overall average distance of all training instances to their respective kth nearest neighbour, where k is selected by the user.
#' @param agreementInput A data frame (N,1) with the agreement measure from M ensemble models, for each instance in \code{trainingSet}.
#' @param STDInput A data frame (N,1) with the ensemble standard deviation for each instance in  \code{trainingSet}.
#' @param trainingSet Data frame with the scaled descritptors of the training instances used to calculate \code{NNthreshold}. scaling should be done with
#' \code{scale(trainingSet, center = mins, scale = maxs - mins)}, where \code{mins} and \code{maxs} are inherited from previously calling \code{getEDmatrix}
#' @param k Number of nearest neighbours to account for when computing the neighbourhood distance.
#' @return \code{getThreshold} returns a data frame with dimensions (N,1) with the threshold neighbourhood distances for each input instance.
#' @details Agreement is calculated from the amount of matching observed and predicted responses in an ensemble of models, divided by the total number of models in the ensemble \eqn{Agreement = \frac{|Obs\cap Pred|}{M}}. The ensemble standard deviation is calculated according to Tetko et al [2]. Both \code{agreementInput} and \code{STDInput} data frames should be provided with dimensions (N,1).
#' This function implicitly uses the object \code{Dij.sort} output by \code{getEDmatrix}, which consists of a data frame with the distances between each \code{trainingSet} instance and its training neighbours (sorted in ascending order of distance).
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
#' # Prepare descriptors to be used
#' train <- data.matrix(bcdata[trainID,1:ncol(bcdata)-1])
#' test <- data.matrix(bcdata[-trainID,1:ncol(bcdata)-1])
#'
#'
#' ## compute distance to neighbours
#' Dij.sort <- getEDmatrix(train, train)
#'
#' #data needs to be scaled as all functions rely on the Euclidean distance
#' train <- scale(train, center = mins, scale = maxs - mins)
#'
#' ## Compute the threshold corresponding to k=3 nearest neighbours, using scaled data
#' NNthreshold <- getThreshold(agree, std, train, k=3)
#'
#' @references [2] IV Tetko, I Sushko, et al. Critical Assessment of QSAR models of
#' environmental toxicity against Tetrahymena pyriformis: focusing on applicability
#' domain and overfitting by variable selection. 2008. 48(9):1733-46. doi:10.1021/ci800151m
#'
#' @export
# train needs previous scaling
getThreshold <- function(agreementInput,STDInput,trainingSet,k){ # k has to start @ 1
  cat("Computing threshold for", toString(k), "nearest neighbors","\n")
  Dk <- data.frame(Dij.sort[1:nrow(Dij.sort), 2:(k+1)]) #assumes Dij.sort is TRvsTR

  Dk <- rowSums(Dk)/k

  Q3 <- unname(quantile(Dk,0.75))
  Q1 <- unname(quantile(Dk,0.25))

  RefVal <- Q3 + 1.5*(Q3-Q1)

  NNthreshold <- c()
  for (i in 1:nrow(Dij.sort)) {
    NNdist <- Dij.sort[i,-1][Dij.sort[i,-1]<=RefVal]
    if (length(NNdist)>0){
      threshold.DistNN <- sum(NNdist)/length(NNdist)
      NNthreshold <- append(NNthreshold, threshold.DistNN)
      #if (threshold.DistNN<0){print("NEGATIVE")}
      #print(threshold.DistNN)
    } else {
      NNthreshold <- append(NNthreshold, "x") # ATTENTION: this will coerse the rest of the numbers into str
    }
  }
  # After NNthreshold is complete for current k, replace X for the smallest threshold in the list
  minT = min(NNthreshold)  # min is able to correctly calculate with "x" in the list
  NNthreshold[NNthreshold=="x"]=minT
  NNthreshold <- as.numeric(NNthreshold)

  # After the distance thresholds are found for every compound, this will be corrected with the P(class) taken as input
  agree = unname(as.matrix(agreementInput))
  STD = unname(as.matrix(STDInput))
  correct = (1-STD)*agree
  minCorrect = min(correct[correct>0])
  correct[correct==0]=minCorrect
  NNthreshold <- NNthreshold * correct
  NNthreshold = NNthreshold/3
  return(NNthreshold)
}
