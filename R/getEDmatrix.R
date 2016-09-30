#' @name getEDmatrix
#' @title Calculate the Euclidean distances between two datasets
#'
#' @description Function that calculates a matrix of Euclidean distances between each pair of instances from two datasets.
#'
#' @param set1 a data frame containing only the molecular features meant for the calculation of the Euclidean Distance
#' @param set2 a data frame containing only the molecular features meant for the calculation of the Euclidean Distance
#'
#' @return a \code{getEDmatrix} object which consists of \code{set1} vs \code{set2} data frame of Euclidean distances, where the values in each row are sorted in
#' ascending order. As a consequence columns have no meaning on their own. \code{getEDmatrix} also implicitly creates two
#' variables, \code{maxs} and \code{mins}, which are automatically saved under such names and do not need explicitly variable
#' assignment. They are created for later data scaling.
#'
#' @details No \code{NA} values are accepted, so either the respective instance is previously removed or empty values should be replaced
#' (eg., with the respective column median or average). For the purpose of using this package, \code{set1} and \code{set2} should be the same dataset.
#' All columns present at the data frames are used in the calculation of the Euclidean distances, i.e. Euclidean distance between set1[rowi] and set2[rowj].
#' Prior to calculating the Euclidean distances, the datasets will be scaled using \code{scale}, which applies \eqn{\frac{x_{ij}-min_j}{max_j-min_j}}
#' to each instance \eqn{x_i} under column (feature) j.
#'
#' @examples
#' train <- matrix(1:9,nrow=3, ncol=3)
#' a <- getEDmatrix(train, train)
#'
#' @export
getEDmatrix <- function(set1,set2) {  #store under var; will scale values before calc ED; Descriptors only!
  # in case set2 != set1, set2 will be scaled to set1(max,min)
  cat("computing the distance matrix...","\n")
  # check NA values
  if(sum(is.na(set1))) stop("Cannot handle NA values; remove instance or replace value from set1", call.=FALSE)
  if(sum(is.na(set2))) stop("Cannot handle NA values; remove instance or replace value from set2", call.=FALSE)
  EDmatrix <- matrix(nrow=nrow(set1), ncol=nrow(set2))
  maxs <<- apply(set1, 2, max)
  mins <<- apply(set1, 2, min)
  set1 <- scale(set1, center = mins, scale = maxs - mins)

  set2 <- scale(set2, center = mins, scale = maxs - mins)

  for (i in 1:nrow(set1)){
    for (j in 1:nrow(set2)){
      EDmatrix[i,j]=sqrt(sum((set1[i,]-set2[j,])**2))
    }
    EDmatrix[i,]=sort(EDmatrix[i,])
  }
  Dij.sort <- EDmatrix
  return(Dij.sort)
}
