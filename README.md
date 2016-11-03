# RDN (README UNDER CONSTRUCTION)
## Description
Reliability Density Neighbourhood is an Applicability Domain (AD) technique meant to estimate the reliability of the predictions from a QSAR (classification) model. This method scans the chemical space, starting from the locations of training instances, taking into account local density, and local bias and precision. After the chemical space has been mapped, the established RDN AD can be used to sort new (external) predictions according to their reliability. Even though the RDN mapping is calculated using `getRDN`, the different tasks that this entails are separately available through the remaining functions in the package, which are listed below. However, despite being available for use, functions should ideally not be called isolated, and the user should use `getRDN` directly instead.

`getRDN` is the automated implementation of the following workflow:

**STEP 1**: Calculation of an Euclidean Distance matrix of the training set through `getEDmatrix`. This matrix will contain the distance between each training instance and each of its training neighbours, sorted in ascending order of distance.

**STEP 2**: Calculation of individual average distance to the k-th nearest neighbours through `getThreshold`. This distance will be used as coverage threshold around each training instance.

**STEP 3**: Place new queries onto the established coverage map using `TestInTrain`. If an instance is located within the radius of coverage around any training instances, it will be deemed as covered by the AD.

This workflow is fully automated in getRDN which runs these steps iteratively for a range of k values, which allows scanning chemical space from the near vicinity around training instances outwards. The full details on the theoretical background of this algorithm are available in the literature. (Aniceto et al. J Cheminf 2016. Submitted.)

## Calculating the RDN AD profile for a given QSAR Classification model
The following example demonstrates how to use `getRDN` to obtain an AD profile for an example dataset (BreastCancer)

### getRDN
Before calculating the RDN mapping, the data needs to be prepared. The data to be passed to the function should only contain descriptors and the dependent variable (response), and an external set should be set aside.

```
library(randomForest)
library(mlbench)
data("BreastCancer")

# remove first Col (IDs)
bcdata <- BreastCancer[,c(2:ncol(BreastCancer))]
bcdata <- bcdata[complete.cases(bcdata),]

#sample 70% for training
trainID <- sample(1:nrow(bcdata),round(nrow(bcdata)*0.7))
```

Now the data will be used to build a QSAR model and gather the external predictions for using later when testing the RDN map.

```
#QSAR model; gather external predictions
# ntree was set to was purposefully to challenge the RDN calculation with more error.
m <- randomForest(Class~.,data=bcdata[c(trainID),], kepp.forest=TRUE, ntree=1, norm.votes=TRUE)
p <- predict(m, bcdata[-trainID,], predict.all=TRUE, type="class")
test_pred <- p$aggregate #class pred external
performance <- as.double(test_pred==bcdata[-trainID,"Class"])
performance <- data.frame(performance)

ensembleClass.t <- data.frame(bcdata$Class[trainID])
ensembleProb.t <- data.frame(bcdata$Class[trainID])
```

Then, an ensemble model should be created (ideally independently from the QSAR model). This will be used to obtain the bias an precision associated with each training compound. Note that **_only_** the training set is used for this, as the external test set will be later on placed onto the reliability-density training matrix across chemical space, built by `getRDN`.

```
# train ensemble for AD calculation
for (i in 2:11){
  #sampling
  samp <- sample(trainID, round(length(trainID)*0.8))
  # train ensemble
  m <- randomForest(Class~.,data=bcdata[samp,], kepp.forest=TRUE, ntree=100, norm.votes=TRUE)
  pred <- predict(m, bcdata[trainID,], predict.all=TRUE, type="Prob")
  pred_class <- predict(m, bcdata[trainID,], predict.all=TRUE, type="class")
  P_AD <- pred$aggregate[,1] #prob train
class_AD <-pred_class$aggregate #class train
ensembleClass.t[,i]=class_AD
ensembleProb.t[,i]=P_AD
}

# compute agreement for TRAIN
agree <- data.frame(trainID)
for (i in 1:length(trainID)){
  agree[i,2] <- sum(ensembleClass.t[i,-1]==toString(ensembleClass.t[i,1]))/10
}

agree <- data.frame(agree[,-1])
#compute std for TRAIN
std <- apply(ensembleProb.t[,-1],1,sd)
std <- data.frame(std)
# Prepare descriptors to be passed into getRDN
train <- data.matrix(bcdata[trainID,1:ncol(bcdata)-1])
test <- data.matrix(bcdata[-trainID,1:ncol(bcdata)-1])
# Compute RDN; this will take care of scaling trainingSet and testSet internally before any use.
resultSummary <- getRDN(performance=performance, trainingSet=train,
testSet=test, agreementInput=agree, STDInput=std)
# The results saved in resultSummary show a decreasing overall quality of predictions as
# AD gets expanded (i.e. instances out of AD decrease).
```
