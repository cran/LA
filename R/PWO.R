#' Calculate the pairwise-order (PWO) of an order-of-addition design
#'
#' \code{PWO} returns the pairwise-order matrix of an order-of-addition design
#'
#' @param X A matrix object. \code{X} must be an order-of-addition design (OofA) matrix.
#'
#' @return If the input is logical, then the output will be a \code{n} by \code{{k \\choose 2}} matrix, where \code{n} and \code{k} are the run size and factor size of the input OofA. Each column of PWO represents one distinct pair of elements from the input OofA. For example, an OofA with 4 factors has the following pairs: \{(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)\}.
#'
#' @references Van Nostrand, R. C. (1995) Design of experiments where the order of addition is important. \emph{ASA proceedings of the Section on Physical and Engineering Sciences}, 155-160.
#'
#' @examples
#' #create a full OofA with 3 factors.
#' toy=rOofA(k=3);toy
#'
#' #Calculate the pairwise-order of toy 
#' PWO(X=toy)
#' @export

PWO=function(X){
  #Pair-wise Order
  
  #X must be an OofA design matrix
  
  m=ncol(X)
  
  Y=matrix(0,ncol=choose(m,2),nrow=nrow(X))
  
  for (k in 1:nrow(Y)) {
    index=1
    
    for (i in 1:(m-1)) {
      for (j in (i+1):m) {
        
        e1=which(X[k,]==i)
        
        e2=which(X[k,]==j)
        
        if(e1<e2){
          Y[k,index]=1
        }
        
        if(e1>e2){
          Y[k,index]=-1
        }
        
        
        index=index+1
      }
      
    }
    
  }
  
  Y
}
