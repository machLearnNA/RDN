# RDN (README UNDER CONSTRUCTION)
## Description
Reliability Density Neighbourhood is an Applicability Domain (AD) technique meant to estimate the reliability of the predictions from a QSAR model. This method scans the chemical space, starting from the locations of training instances, taking into account local density, and local bias and precision. After the chemical space has been mapped, the established RDN AD can be used to sort new (external) predictions according to their reliability. Even though the RDN mapping is calculated using `getRDN`, the different tasks that this entails are separately available through the remaining functions in the package, which are listed below. However, Despite being available for use functions should ideally not be called isolated, andthe user should use `getRDN` directly instead.

This function is the automated implementation of the following workflow:

**STEP 1**: Calculation of an Euclidean Distance matrix of the training set through `getEDmatrix`. This matrix will contain the distance to each of the training neighbours of each training instances, sorted in ascending oderder of distance.

**STEP 2**: Calculation of individual average distance to the k-th nearest neighbours through `getThreshold`. This distance will be used as coverage threshold around each training instance.

**STEP 3**: Place new queries onto the established coverage map using `TestInTrain`. If an instance is located within the radius.

