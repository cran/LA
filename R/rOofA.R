#' Generate a random order-of-addition design (OofA)
#'
#' \code{rOofA} returns a random \code{n} by \code{k} order-of-addition design matrix
#'
#' @param n A positive integer, which stands for the number of rows (or run size). The default setting of \code{n} is \code{k} factorial, which yields a full order-of-addition design matrix. Note that the maximum of \code{n} cannot be greater than \code{k} factorial.
#' @param k A positive integer, which stands for the number of columns (or factor size).
#'
#' @return If all inputs are positive integer, then the output will be a \code{n} by \code{k} design matrix.
#'
#' @examples
#' #generate a full OofA with 4 factors.
#' toy=rOofA(k=4);toy
#'
#' #generate a 12-run random OofA with 4 factors.
#' toy=rOofA(n=12,k=4);toy
#' @export


#Generate a random OofA
rOofA=function(n=factorial(k),k){

  if(n>factorial(k)){
    stop("Run size, n, must be no greater than k factorial")
  }

  if(k==1){
    return(k)
  }

  else if(2<=k&&k<=8){

    X=NULL

    for(i in 1:k){

      X=rbind(X,cbind(i,(rOofA(k=(k-1))+i)))

    }

    for (i in 1:nrow(X)) {
      for (j in 1:ncol(X)) {

        if(X[i,j]>k){
          X[i,j]=X[i,j]%%k
        }

      }
    }

    colnames(X)=NULL

    if(n<factorial(k)){

      rrow=sample(1:factorial(k),n,replace=FALSE)

      X=X[rrow,]

    }

    return(X)
  }

  else {

    if(n==factorial(k)){

      X=NULL

      for(i in 1:k){

        X=rbind(X,cbind(i,(rOofA(k=(k-1))+i)))

      }

      for (i in 1:nrow(X)) {
        for (j in 1:ncol(X)) {

          if(X[i,j]>k){
            X[i,j]=X[i,j]%%k
          }

        }
      }

      colnames(X)=NULL

      return(X)
    }

    if(n<factorial(k)){

      X=NULL

      for (i in 1:n) {
        X=rbind(X,sample(x=seq(1,k),k,replace=F))
      }

      return(X)
    }

  }

}
