#' RDN: Reliability Density Neighborhood for Applicability Domain characterization.
#'
#' The RDN package provides a straightforward way of computing a QSAR model's applicability domain (AD),
#' being currently only applicable for classification models.
#' This method scans the chemical space, starting from the locations of training instances,
#' taking into account local density, and local bias and precision. After the chemical space
#' has been mapped, the established RDN AD can be used to sort new (external) predictions
#' according to their reliability.
#' Even though the RDN mapping is calculated using \code{getRDN}, the different tasks that this entails
#' are separately available through the remaining functions in the package, which are listed below. However,
#' Despite being available for use functions should ideally not be called isolated, and the user should use
#' \code{getRDN} directly instead.
#'
#' The AD will be established according to the following workflow:
#' \itemize{
#'    \item STEP #1: Calculation of an Euclidean Distance matrix of the training set through \code{getEDmatrix}.
#'    This matrix will contain the distance between each training instance and each of its training neighbours, sorted
#'    in ascending order of distance.
#'    \item STEP #2: Calculation of individual average distance to the k-th nearest neighbours through \code{getThreshold}.
#'    This distance will be used as coverage threshold around each training instance.
#'    \item STEP #3: Place new queries onto the established coverage map using \code{TestInTrain}. If an instance is
#'    located within the radius of coverage around any training instances, it will be deemed as covered by the AD.
#' }
#'
#' This workflow is fully automated in \code{getRDN} which runs these steps iteratively for a range of k values, which allows
#' scanning chemical space from the near vicinity around training instances outwards.
#' The full details on the theoretical background of this algorithm are available in the literature.[1]
#'
#' @references [1] N Aniceto, AA Freitas, et al. A Novel Applicability Domain Technique for Mapping Predictive Reliability Accross the Chemical
#' Space of a QSAR: Reliability-Density Neighbourhood. J Cheminf. 2016. Submitted.
#'
#' @docType package
#' @name RDN-package
#' @aliases RDN
#' @importFrom randomForest randomForest
NULL
