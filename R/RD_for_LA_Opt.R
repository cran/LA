#' Lioness Algorithm for experimental designs with continuous input
#'
#' \code{LA_OptC} returns optimal designs with continuous input
#'
#' @param n A positive integer, which stands for the number of rows (or run size) for a design.
#' @param lb A vector contains the lower bounds of all the input variables. For example, if there are 3 input variables whose lower bounds are 0, 5, and 15, \code{lb} should be \code{lb=c(0,5,15)}.
#' @param ub A vector contains the upper bounds of all the input variables. For example, if there are 3 input variables whose upper bounds are 10, 15, and 25, \code{ub} should be \code{lb=c(10,15,25)}.
#' @param m A positive integer, which stands for the number of starting lionesses agents. The default is set to be 100.
#' @param N A positive integer, which stands for the number of iterations. The default is set to be 5000. A large value of \code{N} will result a high CPU time.
#' @param OC An optimality criterion. The default setting is "D-optimality". It could be one of the following: "D", "A", and "G", which stands for "D-optimality", "A-optimality", and "G-optimality", respectively.
#' @param alpha A tuning parameter in algorithm for controlling how big the change would be when updating elements in the step of avoiding local optimum. The default is set to be 0.1, which is the recommended value.
#'
#' @return If all inputs are logical, then the output will be either a \code{n} by \code{length(lb)} optimal design. Here, the \code{length(lb)} is assumed to be at least 2.
#'
#' @examples
#' #Assume in a simple linear regression model, we want to find a D-optimal
#' #20-run design, where the input variable takes values between 0 and 24.
#' #In theory, we know the optimal design is the following:
#' #matrix(c(rep(1,20),rep(0,10),rep(24,10)),ncol=2,nrow=20,byrow=FALSE)
#'
#' #Use LA with default setting to find the optimal design for above problem.
#' try=LA_OptC(n=20,lb=c(1,0),ub=c(1,24))
#' round(try,8)
#' @export

LA_OptC <- function(n, lb, ub, m = 100L, N = 5000L, OC = "D", alpha = 0.1) {
  .Call(`_LA_LA_OptC`, n, lb, ub, m, N, OC, alpha)
}

#' @rdname LA_OptC
#' 
#' @param X A matrix object. In general, \code{X} stands for the design matrix.
#'
D <- function(X) {
  .Call(`_LA_D`, X)
}

#' @rdname LA_OptC
#' 
#' @param X A matrix object. In general, \code{X} stands for the design matrix.
#'
A <- function(X) {
  .Call(`_LA_A`, X)
}

#' @rdname LA_OptC
#' 
#' @param X A matrix object. In general, \code{X} stands for the design matrix.
#' @param x is a vector.
#'
GscoreC <- function(X, x) {
  .Call(`_LA_GscoreC`, X, x)
}

#' @rdname LA_OptC
#' 
#' @param m A positive integer.
#'
rSign <- function(m = 2L) {
  .Call(`_LA_rSign`, m)
}

#' @rdname LA_OptC
#' 
#' @param Y A matrix object. In general, \code{Y} stands for the design matrix.
#'
G <- function(Y) {
  .Call(`_LA_G`, Y)
}

