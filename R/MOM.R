#' Calculate the moment matrix of an order-of-addition design
#'
#' \code{PWO} returns the moment matrix of an order-of-addition design
#'
#' @param X A matrix object. \code{X} must be an order-of-addition design (OofA) matrix.
#'
#' @return If the input is logical, then the output will be a \code{1+{k \\choose 2}} by \code{1+{k \\choose 2}} matrix, where \code{k} is the factor size of the input OofA.
#'
#' @references Lin, D. K., and Peng, J. (2019) Order-of-addition experiments: A review and some new thoughts. \emph{Quality Engineering}, \strong{31}, 49-59.
#'
#' @examples
#' #create a full OofA with 3 factors.
#' toy=rOofA(k=3);toy
#'
#' #Calculate the moment matrix of toy 
#' MOM(X=toy)
#' @export

MOM=function(X){
  
  #X must be an OofA design matrix
  
  Z=PWO(X)
  
  X=cbind(1,Z)
  
  M=t(X)%*%X/nrow(X)
  
  M
}
