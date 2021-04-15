#' Lioness Algorithm for experimental designs with continuous input and function optimization
#'
#' \code{LA_opt} returns optimal designs with continuous input or optimal solutions for function optimization
#'
#' @param of An objective function to be evaluated. If user is seeking for optimal design rather than function optimization, \code{of} should be left with \code{NULL}, which is the default setting.
#' @param n A positive integer, which stands for the number of rows (or run size) for a design. If user is seeking for optimal design, \code{n} should be the number of rows of the design matrix. If user is seeking for function optimization, \code{n} should be left with \code{NULL}, which is the default setting.
#' @param lb A vector contains the lower bounds of all the input variables. For example, if there are 3 input variables whose lower bounds are 0, 5, and 15, \code{lb} should be \code{lb=c(0,5,15)}.
#' @param ub A vector contains the upper bounds of all the input variables. For example, if there are 3 input variables whose upper bounds are 10, 15, and 25, \code{ub} should be \code{lb=c(10,15,25)}.
#' @param N A positive integer, which stands for the number of iterations. The default is set to be 10. A large value of \code{N} will result a high CPU time, and it is recommended to be no greater than 500.
#' @param OC An optimality criterion. If user is seeking for optimal design, \code{OC} should be an optimality criterion for how to evaluate the design matrix. If user is seeking for function optimization, \code{OC} should be left with \code{NULL}, which is the default setting.
#' @param type A logic input argument, which indicates the type of optimization. If \code{type} is \code{mini} (the default setting), minimization will be implemented in the algorithm. If \code{type} is \code{maxi},  maximization will be implemented in the algorithm.
#' @param maxtime A positive number, which indicates the expected maximum CPU time given by user, and it is measured by minutes. For example, maxtime=3.5 indicates the CPU time will be no greater than three and half minutes. The default is set to be 5.
#'
#' @return If all inputs are logical, then the output will be either a \code{n} by \code{length(lb)} optimal design or a 1 by \code{length(lb)} vector of optimal solutions for function optimization.
#'
#' @examples
#' #We start with function optimization
#' #Now define an objective function: Sum of Different Powers
#' SDP=function(x){i=1:length(x);y=sum(abs(x)^(i=1));return(y)}
#'
#' #Use LA to find the optimal solution under 20-dimensional setting
#' #for SDP function with 10 iterations.
#' try=LA_opt(of=SDP,lb=rep(-1,20),ub=rep(1,20),N=10,type="mini")
#' SDP(try)  #Note that the true global optimum is 0, but we only have 10 iterations
#'
#' #Another example
#' #Define an objective function: Cross-in-Tray
#' CiT=function(x){x1=x[1];x2=x[2];y=-0.0001*(abs(sin(x1)*sin(x2)*
#' exp(abs(100-sqrt(x1^2+x2^2)/pi)))+1)^0.1;return(y)}
#'
#' #Use LA to find the optimal solution under 2-dimensional setting
#' #for CiT function with 10 iterations.
#' try2=LA_opt(of=CiT,lb=rep(-10,2),ub=rep(10,2),N=10,type="mini")
#' CiT(try2) #Note that the true global optimum is -2.06261, but we only have 10 iterations
#'
#' #Next we introduce the optimal design part
#' #Assume in a simple linear regression model, we want to find a D-optimal
#' #20-run design, where the input variable takes values between 0 and 24.
#' #In theory, we know the optimal design is the following:
#' #matrix(c(rep(1,20),rep(0,10),rep(24,10)),ncol=2,nrow=20,byrow=FALSE)
#' #Let's see if LA is able to identify that.
#'
#' #Define the D-optimality criterion in simple linear regression model:
#' D=function(x){IM=t(x)%*%x;return(det(IM))}
#'
#' #Use LA to find the optimal solution for above problem.
#' #We want to maximize the determinant of the information matrix.
#' try3=LA_opt(n=20,lb=c(1,0),ub=c(1,24),N=10,OC=D,type="maxi")
#' try3 #with more iterations, LA would return even better result.
#' @export

LA_opt=function(of=NULL,n=NULL,lb,ub,N=10,OC=NULL,type="mini",maxtime=5){

  if(length(lb)!=length(ub)){
    stop("The number of lower bounds is not equal to the number of upper bounds, please check.")
  }

  maxtime=maxtime*60  #convert minutes to seconds
  timeALL=NULL        #record all cpu time

  if(is.null(of)==FALSE){
    n=1

    C=1  #Initialize counter index

    k=length(lb)    #count how many input variables

    m=6*k+3         #calculate how many lionesses

    #step 1 starts
    S=array(0,dim=c(n,k,m))        #solution matrix

    for (i in 1:n) {
      for (j in 1:k) {
        for (l in 1:m) {

          S[i,j,l]=stats::runif(1,min=lb[j],max=ub[j])

        }
      }
    }
    #step 1 ends

    if(type=="mini"){

    #step 2 starts
    result=rep(0,m)

    for (i in 1:m) {
      result[i]=of(S[,,i])
    }

    #step 2 ends

    progressbar = utils::txtProgressBar(min = 0, max = N, style = 3)

    while (C<=N) {   #step 3

      time0=Sys.time()

      temp=cbind(result,1:m)

      temp=temp[order(temp[,1]),]

      #step 4: determine the top 3 agents
      centre=S[,,temp[1,2]]

      LW=S[,,temp[2,2]]

      RW=S[,,temp[3,2]]

      S[,,1]=centre
      S[,,2]=LW
      S[,,3]=RW

      #centre troop
      index=4

      for (j in 1:k) {

        S[,,index]=centre
        S[,j,index]=LW[j]
        index=index+1

        S[,,index]=centre
        S[,j,index]=RW[j]
        index=index+1

      }

      #LW troop

      for (j in 1:k) {

        S[,,index]=LW
        S[,j,index]=centre[j]
        index=index+1

        S[,,index]=LW
        S[,j,index]=RW[j]
        index=index+1

      }

      #RW troop

      for (j in 1:k) {

        S[,,index]=RW
        S[,j,index]=centre[j]
        index=index+1

        S[,,index]=RW
        S[,j,index]=LW[j]
        index=index+1

      }

      #Until here, S has been fully updated

      for (l in 2:m) {
        for (i in 1:n) {
          for (j in 1:k) {

            z=stats::runif(1,0,1)

            if (z<=((N-C+1)/N)){

              S[i,j,l]=S[i,j,l]+sample(c(-1,1),1,prob=c(0.5,0.5))*0.1*S[i,j,l]

            }

            if(S[i,j,l]>ub[j]){S[i,j,l]=ub[j]}

            if(S[i,j,l]<lb[j]){S[i,j,l]=lb[j]}

          }
        }
      }

      #update result for all
      result=rep(0,m)

      for (i in 1:m) {
        result[i]=of(S[,,i])
      }

      time1=Sys.time()
      timediff=time1-time0
      timeALL=c(timeALL,timediff)

      ##########progress bar codes
      utils::setTxtProgressBar(progressbar, C)
      ##########

      if(as.numeric(sum(timeALL)+timediff)<=maxtime){C=C+1}
      if(as.numeric(sum(timeALL)+timediff)>maxtime){C=N+1}
    }

    temp=cbind(result,1:m)

    temp=temp[order(temp[,1]),]

    centre=S[,,temp[1,2]]

    }

    if(type=="maxi"){

    #step 2 starts
    result=rep(0,m)

    for (i in 1:m) {
      result[i]=of(S[,,i])
    }

    #step 2 ends

    progressbar = utils::txtProgressBar(min = 0, max = N, style = 3)

    while (C<=N) {   #step 3

      time0=Sys.time()

      temp=cbind(result,1:m)

      temp=temp[order(temp[,1],decreasing=TRUE),]

      #step 4: determine the top 3 agents
      centre=S[,,temp[1,2]]

      LW=S[,,temp[2,2]]

      RW=S[,,temp[3,2]]

      S[,,1]=centre
      S[,,2]=LW
      S[,,3]=RW

      #centre troop
      index=4

      for (j in 1:k) {

        S[,,index]=centre
        S[,j,index]=LW[j]
        index=index+1

        S[,,index]=centre
        S[,j,index]=RW[j]
        index=index+1

      }

      #LW troop

      for (j in 1:k) {

        S[,,index]=LW
        S[,j,index]=centre[j]
        index=index+1

        S[,,index]=LW
        S[,j,index]=RW[j]
        index=index+1

      }

      #RW troop

      for (j in 1:k) {

        S[,,index]=RW
        S[,j,index]=centre[j]
        index=index+1

        S[,,index]=RW
        S[,j,index]=LW[j]
        index=index+1

      }

      #Until here, S has been fully updated

      for (l in 2:m) {
        for (i in 1:n) {
          for (j in 1:k) {

            z=stats::runif(1,0,1)

            if (z<=((N-C+1)/N)){

              S[i,j,l]=S[i,j,l]+sample(c(-1,1),1,prob=c(0.5,0.5))*0.1*S[i,j,l]

            }

            if(S[i,j,l]>ub[j]){S[i,j,l]=ub[j]}

            if(S[i,j,l]<lb[j]){S[i,j,l]=lb[j]}

          }
        }
      }

      #update result for all
      result=rep(0,m)

      for (i in 1:m) {
        result[i]=of(S[,,i])
      }

      time1=Sys.time()
      timediff=time1-time0
      timeALL=c(timeALL,timediff)

      ##########progress bar codes
      utils::setTxtProgressBar(progressbar, C)
      ##########

      if(as.numeric(sum(timeALL)+timediff)<=maxtime){C=C+1}
      if(as.numeric(sum(timeALL)+timediff)>maxtime){C=N+1}
    }

    temp=cbind(result,1:m)

    temp=temp[order(temp[,1],decreasing=TRUE),]

    centre=S[,,temp[1,2]]

    }

  }

  if(is.null(of)==TRUE){

    if(is.null(n)==TRUE){
      stop("Design run size, n, needs to be specified.")
    }

    if(is.null(OC)==TRUE){
      stop("An optimality criterion, OC, needs to be specified.")
    }

    C=1  #Initialize counter index

    k=length(lb)    #count how many input variables

    m=6*k+3         #calculate how many lionesses

    #step 1 starts
    S=array(0,dim=c(n,k,m))        #solution matrix

    for (i in 1:n) {
      for (j in 1:k) {
        for (l in 1:m) {

          S[i,j,l]=stats::runif(1,min=lb[j],max=ub[j])

        }
      }
    }
    #step 1 ends

    if(type=="mini"){

      #step 2 starts
      result=rep(0,m)

      for (i in 1:m) {
        result[i]=OC(S[,,i])
      }

      #step 2 ends

      progressbar = utils::txtProgressBar(min = 0, max = N, style = 3)

      while (C<=N) {   #step 3

        time0=Sys.time()

        temp=cbind(result,1:m)

        temp=temp[order(temp[,1]),]

        #step 4: determine the top 3 agents
        centre=S[,,temp[1,2]]

        LW=S[,,temp[2,2]]

        RW=S[,,temp[3,2]]

        S[,,1]=centre
        S[,,2]=LW
        S[,,3]=RW

        #centre troop
        index=4

        for (j in 1:k) {

          S[,,index]=centre
          S[,j,index]=LW[j]
          index=index+1

          S[,,index]=centre
          S[,j,index]=RW[j]
          index=index+1

        }

        #LW troop

        for (j in 1:k) {

          S[,,index]=LW
          S[,j,index]=centre[j]
          index=index+1

          S[,,index]=LW
          S[,j,index]=RW[j]
          index=index+1

        }

        #RW troop

        for (j in 1:k) {

          S[,,index]=RW
          S[,j,index]=centre[j]
          index=index+1

          S[,,index]=RW
          S[,j,index]=LW[j]
          index=index+1

        }

        #Until here, S has been fully updated

        for (l in 2:m) {
          for (i in 1:n) {
            for (j in 1:k) {

              z=stats::runif(1,0,1)

              if (z<=((N-C+1)/N)){

                S[i,j,l]=S[i,j,l]+sample(c(-1,1),1,prob=c(0.5,0.5))*0.1*S[i,j,l]

              }

              if(S[i,j,l]>ub[j]){S[i,j,l]=ub[j]}

              if(S[i,j,l]<lb[j]){S[i,j,l]=lb[j]}

            }
          }
        }

        #update result for all
        result=rep(0,m)

        for (i in 1:m) {
          result[i]=OC(S[,,i])
        }

        time1=Sys.time()
        timediff=time1-time0
        timeALL=c(timeALL,timediff)

        ##########progress bar codes
        utils::setTxtProgressBar(progressbar, C)
        ##########

        if(as.numeric(sum(timeALL)+timediff)<=maxtime){C=C+1}
        if(as.numeric(sum(timeALL)+timediff)>maxtime){C=N+1}
      }

      temp=cbind(result,1:m)

      temp=temp[order(temp[,1]),]

      centre=S[,,temp[1,2]]

    }

    if(type=="maxi"){

      #step 2 starts
      result=rep(0,m)

      for (i in 1:m) {
        result[i]=OC(S[,,i])
      }

      #step 2 ends

      progressbar = utils::txtProgressBar(min = 0, max = N, style = 3)

      while (C<=N) {   #step 3

        time0=Sys.time()

        temp=cbind(result,1:m)

        temp=temp[order(temp[,1],decreasing=TRUE),]

        #step 4: determine the top 3 agents
        centre=S[,,temp[1,2]]

        LW=S[,,temp[2,2]]

        RW=S[,,temp[3,2]]

        S[,,1]=centre
        S[,,2]=LW
        S[,,3]=RW

        #centre troop
        index=4

        for (j in 1:k) {

          S[,,index]=centre
          S[,j,index]=LW[j]
          index=index+1

          S[,,index]=centre
          S[,j,index]=RW[j]
          index=index+1

        }

        #LW troop

        for (j in 1:k) {

          S[,,index]=LW
          S[,j,index]=centre[j]
          index=index+1

          S[,,index]=LW
          S[,j,index]=RW[j]
          index=index+1

        }

        #RW troop

        for (j in 1:k) {

          S[,,index]=RW
          S[,j,index]=centre[j]
          index=index+1

          S[,,index]=RW
          S[,j,index]=LW[j]
          index=index+1

        }

        #Until here, S has been fully updated

        for (l in 2:m) {
          for (i in 1:n) {
            for (j in 1:k) {

              z=stats::runif(1,0,1)

              if (z<=((N-C+1)/N)){

                S[i,j,l]=S[i,j,l]+sample(c(-1,1),1,prob=c(0.5,0.5))*0.1*S[i,j,l]

              }

              if(S[i,j,l]>ub[j]){S[i,j,l]=ub[j]}

              if(S[i,j,l]<lb[j]){S[i,j,l]=lb[j]}

            }
          }
        }

        #update result for all
        result=rep(0,m)

        for (i in 1:m) {
          result[i]=OC(S[,,i])
        }

        time1=Sys.time()
        timediff=time1-time0
        timeALL=c(timeALL,timediff)

        ##########progress bar codes
        utils::setTxtProgressBar(progressbar, C)
        ##########

        if(as.numeric(sum(timeALL)+timediff)<=maxtime){C=C+1}
        if(as.numeric(sum(timeALL)+timediff)>maxtime){C=N+1}
      }

      temp=cbind(result,1:m)

      temp=temp[order(temp[,1],decreasing=TRUE),]

      centre=S[,,temp[1,2]]

    }

  }


  avgtime=round(mean(timeALL),2)
  iterations=length(timeALL)

  close(progressbar)
  print(paste0("average CPU time per iteration is: ", avgtime, " seconds"))
  print(paste0("the number of iterations completed is: ", iterations))


  centre
}
