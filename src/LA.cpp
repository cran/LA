#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec permuC(arma::vec x){
  
  int n=x.size();
  
  for (int i=0;i<=(n-2);i++){
    
    int j=i+floor(unif_rand()*(n-i));
    
    std::swap(x[i],x[j]);
    
  }
  return x;
}

// [[Rcpp::export]]
arma::vec seqC(unsigned a, int b){
  
  int n=b-a+1;
  arma::vec seq(n);
  
  seq[0]=a;
  
  for(int i=0;i<n;i++){
    
    seq[i]=seq[0]+i;
    
  }
  
  return seq;
  
}

// [[Rcpp::export]]
arma::mat rLHDC(int n, int k){
  arma::mat X(n,k);
  
  for(int j=0;j<k;j++){
    
    X.col(j)=permuC(seqC(1,n));
    
  }
  
  return X;
}

// [[Rcpp::export]]
double dijC(arma::mat X, int i, int j, int q=1){
  
  return pow(sum(pow(abs(X.row(i-1)-X.row(j-1)),q)),1.0/q);
  
}

// [[Rcpp::export]]
double phi_pC(arma::mat X,int p=15,int q=1){
  
  int n=X.n_rows;
  
  double result=0;
  
  for(int i=0;i<n-1;i++){
    
    for (int j=i+1;j<n;j++){
      
      result=result+pow(pow(sum(pow(abs(X.row(i)-X.row(j)),q)),1.0/q),-p);
      
    }
    
  }
  
  return pow(result,1.0/p);
}

// [[Rcpp::export]]
double MaxProCriterionC(arma::mat X){
  
  int n=X.n_rows;
  int k=X.n_cols;
  
  double result=0;
  
  double denom=1;
  
  for(int i=0;i<n-1;i++){
    
    for (int j=i+1;j<n;j++){
      
      for (int l=0;l<k;l++){
        
        denom=denom*pow(X(i,l)-X(j,l),2.0);
        
      }
      
      result=result+1.0/denom;
      
      denom=1;
      
    }
    
  }
  
  return pow(2.0/(n*(n-1))*result,1.0/k);
}

// [[Rcpp::export]]
double corC(arma::vec x, arma::vec y){
  
  int n=x.size();
  int m=y.size();
  double meanx=sum(x)/n;
  double meany=sum(y)/m;
  double term1=0;
  double term2=0;
  double term3=0;
  
  for(int i=0;i<n;i++){
    
    term1=term1+(x[i]-meanx)*(y[i]-meany);
    
    term2=term2+pow((x[i]-meanx),2.0);
    
    term3=term3+pow((y[i]-meany),2.0);
  }
  
  return term1/(pow(term2,0.5)*pow(term3,0.5));
  
}

// [[Rcpp::export]]
double MaxAbsCorC(arma::mat X){
  
  int k=X.n_cols;
  
  double maxcor=0;
  
  for (int i=0;i<k-1;i++){
    for (int j=i+1;j<k;j++){
      
      if(std::abs(corC(X.col(i),X.col(j)))>maxcor){maxcor=std::abs(corC(X.col(i),X.col(j)));}
      
    }
  }
  
  return maxcor;
  
}

// [[Rcpp::export]]
double AvgAbsCorC(arma::mat X){
  
  int k=X.n_cols;
  
  double totalcor=0;
  
  int C=0;
  
  for (int i=0;i<k-1;i++){
    for (int j=i+1;j<k;j++){
      
      totalcor=totalcor+std::abs(corC(X.col(i),X.col(j)));
      
      C=C+1;
      
    }
  }
  
  return totalcor/C;
  
}

// [[Rcpp::export]]
arma::mat exchangeC(arma::mat X, int j, String type="col"){
  
  if(type=="col"){
    
    int n=X.n_rows;
    
    arma::uvec location=arma::randperm(n,2);
    
    std::swap(X(location[0],j-1),X(location[1],j-1));
    
  }
  
  if(type=="row"){
    
    int k=X.n_cols;
    
    arma::uvec location=arma::randperm(k,2);
    
    std::swap(X(j-1,location[0]),X(j-1,location[1]));
    
  }
  
  return X;
}

// [[Rcpp::export]]
arma::mat LA_LHDC(int n, int k, int m=100, int N=5000,  
                  String OC="phi_p",int p=15, int q=1){
  
  double prun=1.0/(k-1);
  
  int C=1;
  
  arma::cube X(n,k,m);
  
  for (int i=0;i<m;i++){
    
    X.slice(i)=rLHDC(n,k);
    
  }
  
  arma::vec result(m);
  
  arma::uvec indices;
  
  arma::mat centre;
  
  arma::mat LW;
  
  arma::mat RW;
  
  int mnew=6*k+3;
  
  arma::cube Xnew(n,k,mnew);
  
  int index;
  
  double z;
  
  if(OC=="phi_p"){
    
    for (int i=0;i<m;i++){
      
      result(i,0)=phi_pC(X.slice(i),p,q);
      
    }
    
  }
  
  if(OC=="MaxProCriterion"){
    
    for (int i=0;i<m;i++){
      
      result(i,0)=MaxProCriterionC(X.slice(i));
      
    }
    
  }
  
  if(OC=="MaxAbsCor"){
    
    for (int i=0;i<m;i++){
      
      result(i,0)=MaxAbsCorC(X.slice(i));
      
    }
    
  }
  
  if(OC=="AvgAbsCor"){
    
    for (int i=0;i<m;i++){
      
      result(i,0)=AvgAbsCorC(X.slice(i));
      
    }
    
  }
  
  while(C<=N){
    
    indices=sort_index(result.col(0),"ascend");
    
    centre=X.slice(indices(0,0));
    
    LW=X.slice(indices(1,0));
    
    RW=X.slice(indices(2,0));
    
    Xnew.slice(0)=centre;
    
    Xnew.slice(1)=LW;
    
    Xnew.slice(2)=RW;
    
    index=3;
    
    for (int j=0;j<k;j++){
      
      Xnew.slice(index)=centre;
      Xnew.slice(index).col(j)=LW.col(j);
      index=index+1;  
      
      Xnew.slice(index)=centre;
      Xnew.slice(index).col(j)=RW.col(j);
      index=index+1;
      
    }
    
    for (int j=0;j<k;j++){
      
      Xnew.slice(index)=LW;
      Xnew.slice(index).col(j)=centre.col(j);
      index=index+1;  
      
      Xnew.slice(index)=LW;
      Xnew.slice(index).col(j)=RW.col(j);
      index=index+1;
      
    }
    
    for (int j=0;j<k;j++){
      
      Xnew.slice(index)=RW;
      Xnew.slice(index).col(j)=centre.col(j);
      index=index+1;  
      
      Xnew.slice(index)=RW;
      Xnew.slice(index).col(j)=LW.col(j);
      index=index+1;
      
    }
    //Until here, Xnew has been fully updated
    
    X=Xnew;
    
    for (int i=1;i<mnew;i++){
      
      for (int j=0;j<k;j++){
        
        z=unif_rand();
        
        if(z<=prun){
          
          X.slice(i)=exchangeC(X.slice(i),j+1);
          
        }
        
      }
      
    }
    
    result.set_size(mnew);
    
    if(OC=="phi_p"){
      
      for (int i=0;i<mnew;i++){
        
        result(i,0)=phi_pC(X.slice(i),p,q);
        
      }
      
    }
    
    if(OC=="MaxProCriterion"){
      
      for (int i=0;i<mnew;i++){
        
        result(i,0)=MaxProCriterionC(X.slice(i));
        
      }
      
    }
    
    if(OC=="MaxAbsCor"){
      
      for (int i=0;i<mnew;i++){
        
        result(i,0)=MaxAbsCorC(X.slice(i));
        
      }
      
    }
    
    if(OC=="AvgAbsCor"){
      
      for (int i=0;i<mnew;i++){
        
        result(i,0)=AvgAbsCorC(X.slice(i));
        
      }
      
    }
    
    C=C+1;
    
  }
  
  indices=sort_index(result.col(0),"ascend");
  
  centre=X.slice(indices(0,0));
  
  return centre;
}


// [[Rcpp::export]]
double factorialC(int x){
  
  return(R::gammafn(x+1.0));
  
}

// [[Rcpp::export]]
int modC(int a,int b){
  
  return a-floor(a/b)*b;
  
}

// [[Rcpp::export]]
arma::mat rOofAC(int n, int k){
  
  if(n>factorialC(k)){
    stop("Run size, n, must be no greater than k factorial");
  }
  
  int C=2;
  
  arma::mat X;
  
  arma::mat Xnew;
  
  if(k==1){
    
    X.set_size(1,1);
    
    X(0,0)=1;
    
  }
  
  else{
    
    X.set_size(1,1);
    
    X(0,0)=1;
    
    while(C<=k){
      
      Xnew.set_size(factorialC(C),C);
      
      for(int i=0;i<C;i++){
        
        Xnew.submat(factorialC(C-1)*i,0,factorialC(C-1)*(i+1)-1,0)=arma::ones(factorialC(C-1))+i;
        
        Xnew.submat(factorialC(C-1)*i,1,factorialC(C-1)*(i+1)-1,C-1)=X+arma::ones(factorialC(C-1),C-1)+i;
        
      }
      
      for(int i=0;i<factorialC(C);i++){
        
        for(int j=0;j<C;j++){
          
          if(Xnew(i,j)>C){Xnew(i,j)=modC(Xnew(i,j),C);}
          
        }
        
      }
      
      X=Xnew;
      C=C+1;
    }
    
  }
  
  if(n<factorialC(k)){
    
    X=X.rows(arma::randperm(factorialC(k),n));
    
  }
  
  return X;
}

// [[Rcpp::export]]
arma::mat PWOC(arma::mat X){
  
  int n=X.n_rows;
  int k=X.n_cols;
  
  int m=k*(k-1)/2.0;
  arma::mat result(n,m);
  
  int C;
  
  arma::uvec indices;
  
  for(int h=0;h<n;h++){
    C=0;
    
    for(int i=1;i<k;i++){
      
      for(int j=i+1;j<k+1;j++){
        
        indices=sort_index(X.row(h),"ascend");
        
        if(indices(i-1,0)<indices(j-1,0)){
          
          result(h,C)=1;
          
        } else{
          
          result(h,C)=-1;
          
        }
        
        C=C+1;
      }
      
    }
    
    
  }
  
  
  return result;
}

// [[Rcpp::export]]
arma::mat TC(arma::mat X){
  
  int n=X.n_rows;
  int k=X.n_cols;
  
  arma::mat Y(k,n);
  
  for(int i=0;i<k;i++){
    
    for(int j=0;j<n;j++){
      
      Y(i,j)=X(j,i);
      
    }
    
  }
  
  return Y;
  
}

// [[Rcpp::export]]
arma::mat MOMC(arma::mat X){
  
  arma::mat Z=PWOC(X);
  
  int n=Z.n_rows;
  int k=Z.n_cols;
  
  arma::mat Y(n,1+k);
  
  Y.col(0)=arma::ones(n);
  Y.cols(1,k)=Z;
  
  return TC(Y)*Y/n;
}

// [[Rcpp::export]]
arma::mat LA_OofAC(int n, int k, int m=100, int N=5000){
  
  double prun=1.0/(n-1);
  
  int C=1;
  
  arma::cube X(n,k,m);
  
  for (int i=0;i<m;i++){
    
    X.slice(i)=rOofAC(n,k);
    
  }
  
  arma::vec result(m);
  
  arma::uvec indices;
  
  arma::mat centre;
  
  arma::mat LW;
  
  arma::mat RW;
  
  int mnew=6*n+3;
  
  arma::cube Xnew(n,k,mnew);
  
  int index;
  
  double z;
  
  //if(OC=="D"){
  
  for (int i=0;i<m;i++){
    
    result(i,0)=arma::det(MOMC(X.slice(i)));
    
  }
  
  //}
  
  //if(OC=="A"){
  
  //for (int i=0;i<m;i++){
  
  //result(i,0)=arma::trace(arma::inv(MOMC(X.slice(i))));
  
  //}
  
  //}
  
  while(C<=N){
    
    indices=sort_index(result.col(0),"descend");
    
    centre=X.slice(indices(0,0));
    
    LW=X.slice(indices(1,0));
    
    RW=X.slice(indices(2,0));
    
    Xnew.slice(0)=centre;
    
    Xnew.slice(1)=LW;
    
    Xnew.slice(2)=RW;
    
    index=3;
    
    for (int i=0;i<n;i++){
      
      Xnew.slice(index)=centre;
      Xnew.slice(index).row(i)=LW.row(i);
      index=index+1;
      
      Xnew.slice(index)=centre;
      Xnew.slice(index).row(i)=RW.row(i);
      index=index+1;
      
    }
    
    for (int i=0;i<n;i++){
      
      Xnew.slice(index)=LW;
      Xnew.slice(index).row(i)=centre.row(i);
      index=index+1;
      
      Xnew.slice(index)=LW;
      Xnew.slice(index).row(i)=RW.row(i);
      index=index+1;
      
    }
    
    for (int i=0;i<n;i++){
      
      Xnew.slice(index)=RW;
      Xnew.slice(index).row(i)=centre.row(i);
      index=index+1;
      
      Xnew.slice(index)=RW;
      Xnew.slice(index).row(i)=LW.row(i);
      index=index+1;
      
    }
    //Until here, Xnew has been fully updated
    
    X=Xnew;
    
    for (int i=1;i<mnew;i++){
      
      for (int j=0;j<n;j++){
        
        z=unif_rand();
        
        if(z<=prun){
          
          X.slice(i)=exchangeC(X.slice(i),j+1,"row");
          
        }
        
      }
      
    }
    
    result.set_size(mnew);
    
    //if(OC=="D"){
    
    for (int i=0;i<mnew;i++){
      
      result(i,0)=arma::det(MOMC(X.slice(i)));
      
    }
    
    //}
    
    //if(OC=="A"){
    
    //for (int i=0;i<mnew;i++){
    
    //result(i,0)=arma::trace(arma::inv(MOMC(X.slice(i))));
    
    //}
    
    //}
    
    C=C+1;
    
  }
  
  indices=sort_index(result.col(0),"descend");
  
  centre=X.slice(indices(0,0));
  
  return centre;
}

// [[Rcpp::export]]
double D(arma::mat X){
  
  return arma::det(arma::inv(X.t()*X));
  
}

// [[Rcpp::export]]
double A(arma::mat X){
  
  return arma::trace(arma::inv(X.t()*X));
  
}

// [[Rcpp::export]]
double GscoreC(arma::mat X, arma::vec x){
  
  int nX=X.n_rows;
  int kX=X.n_cols;
  int pX=(kX+1)*(kX+2)/2.0;
  
  arma::mat Fx(nX,pX);
  
  Fx.col(0)=arma::ones(nX);
  
  for(int i=0;i<kX;i++){
    
    Fx.col(i+1)=X.col(i);
    
  }
  
  int z=kX+1;
  
  if(kX>=2){
    
    for(int i=0;i<kX-1;i++){
      
      for(int j=i+1;j<kX;j++){
        
        for(int l=0;l<nX;l++){
          
          Fx(l,z)=X(l,i)*X(l,j);
          
        }
        z=z+1;
        
      }
      
    }
    
  }
  
  for(int i=0;i<kX;i++){
    
    Fx.col(z)=pow(X.col(i),2.0);
    z=z+1;
    
  }
  
  arma::mat fx(pX,1);
  
  fx(0,0)=1;
  
  for(int i=0;i<kX;i++){
    
    fx(i+1,0)=x(i,0);
    
  }
  
  z=kX+1;
  
  if(kX>=2){
    
    for(int i=0;i<kX-1;i++){
      
      for(int j=i+1;j<kX;j++){
        
        fx(z,0)=x(i,0)*x(j,0);
        z=z+1;
        
      }
      
    }
    
  }
  
  for(int i=0;i<kX;i++){
    
    fx(z,0)=pow(x(i,0),2.0);
    z=z+1;
    
  }
  
  arma::mat res=nX*(fx.t()*arma::inv((Fx.t()*Fx))*fx);
  
  return res(0,0);
  
}

// [[Rcpp::export]]
double rSign(int s=2){
  
  arma::vec x(s);
  
  x(0,0)=-1.0;
  x(1,0)=1.0;
  
  arma::uvec r=arma::randperm(2,1);
  
  return x(r(0,0),0);
  
}

// [[Rcpp::export]]
double G(arma::mat Y){
  
  int k=Y.n_cols;
  
  arma::vec lb=-1*arma::ones(k);
  arma::vec ub=arma::ones(k);
  
  int n=1;
  int m=100;
  int N=500;
  double alpha=0.1;
  
  int C=1;
  
  arma::cube X(n,k,m);
  
  for (int i=0;i<m;i++){
    
    for (int j=0;j<k;j++){
      
      for (int l=0;l<n;l++){
        
        X(l,j,i)=unif_rand()*(ub(j,0)-lb(j,0))+lb(j,0);
        
      }
      
    }
    
  }
  
  arma::vec result(m);
  
  arma::uvec indices;
  
  arma::mat centre;
  
  arma::mat LW;
  
  arma::mat RW;
  
  int mnew=6*k+3;
  
  arma::cube Xnew(n,k,mnew);
  
  int index;
  
  double z;
  
  for (int i=0;i<m;i++){
    
    result(i,0)=GscoreC(Y,X.slice(i).t());
    
  }
  
  while(C<=N){
    
    indices=sort_index(result.col(0),"descend");
    
    centre=X.slice(indices(0,0));
    
    LW=X.slice(indices(1,0));
    
    RW=X.slice(indices(2,0));
    
    Xnew.slice(0)=centre;
    
    Xnew.slice(1)=LW;
    
    Xnew.slice(2)=RW;
    
    index=3;
    
    for (int j=0;j<k;j++){
      
      Xnew.slice(index)=centre;
      Xnew.slice(index).col(j)=LW.col(j);
      index=index+1;  
      
      Xnew.slice(index)=centre;
      Xnew.slice(index).col(j)=RW.col(j);
      index=index+1;
      
    }
    
    for (int j=0;j<k;j++){
      
      Xnew.slice(index)=LW;
      Xnew.slice(index).col(j)=centre.col(j);
      index=index+1;  
      
      Xnew.slice(index)=LW;
      Xnew.slice(index).col(j)=RW.col(j);
      index=index+1;
      
    }
    
    for (int j=0;j<k;j++){
      
      Xnew.slice(index)=RW;
      Xnew.slice(index).col(j)=centre.col(j);
      index=index+1;  
      
      Xnew.slice(index)=RW;
      Xnew.slice(index).col(j)=LW.col(j);
      index=index+1;
      
    }
    //Until here, Xnew has been fully updated
    
    X=Xnew;
    
    for (int i=1;i<mnew;i++){
      
      for (int j=0;j<k;j++){
        
        for (int l=0;l<n;l++){
          
          z=unif_rand();
          
          if(z<=((N-C+1.0)/N)){
            
            X(l,j,i)=X(l,j,i)*(1.0+rSign(2)*alpha);
            
          }
          
          if(X(l,j,i)>ub(j,0)){X(l,j,i)=ub(j,0);}
          if(X(l,j,i)<lb(j,0)){X(l,j,i)=lb(j,0);}
          
        }
        
      }
      
    }
    
    result.set_size(mnew);
    
    for (int i=0;i<mnew;i++){
      
      result(i,0)=GscoreC(Y,X.slice(i).t());
      
    }
    
    C=C+1;
    
  }
  
  indices=sort_index(result.col(0),"descend");
  
  return result(indices(0,0),0);
  
}

// [[Rcpp::export]]
arma::mat LA_OptC(int n, arma::vec lb, arma::vec ub, int m=100, int N=5000,
                  String OC="D", double alpha=0.1){
  
  int klb=lb.n_rows;
  int kub=ub.n_rows;
  
  if(klb!=kub){
    
    stop("The number of lower bounds is not equal to the number of upper bounds, please check.");
    
  }
  
  if(klb<=1){
    
    stop("The number of input variables needs to be at least two.");
    
  }
  
  int C=1;
  
  arma::cube X(n,klb,m);
  
  for (int i=0;i<m;i++){
    
    for (int j=0;j<klb;j++){
      
      for (int l=0;l<n;l++){
        
        X(l,j,i)=unif_rand()*(ub(j,0)-lb(j,0))+lb(j,0);
        
      }
      
    }
    
  }
  
  arma::vec result(m);
  
  arma::uvec indices;
  
  arma::mat centre;
  
  arma::mat LW;
  
  arma::mat RW;
  
  int mnew=6*klb+3;
  
  arma::cube Xnew(n,klb,mnew);
  
  int index;
  
  double z;
  
  if(OC=="A"){
    
    for (int i=0;i<m;i++){
      
      result(i,0)=A(X.slice(i));
      
    }
    
  }
  
  if(OC=="D"){
    
    for (int i=0;i<m;i++){
      
      result(i,0)=D(X.slice(i));
      
    }
    
  }
  
  if(OC=="G"){
    
    for (int i=0;i<m;i++){
      
      result(i,0)=G(X.slice(i));
      
    }
    
  }
  
  
  while(C<=N){
    
    indices=sort_index(result.col(0),"ascend");
    
    centre=X.slice(indices(0,0));
    
    LW=X.slice(indices(1,0));
    
    RW=X.slice(indices(2,0));
    
    Xnew.slice(0)=centre;
    
    Xnew.slice(1)=LW;
    
    Xnew.slice(2)=RW;
    
    index=3;
    
    for (int j=0;j<klb;j++){
      
      Xnew.slice(index)=centre;
      Xnew.slice(index).col(j)=LW.col(j);
      index=index+1;  
      
      Xnew.slice(index)=centre;
      Xnew.slice(index).col(j)=RW.col(j);
      index=index+1;
      
    }
    
    for (int j=0;j<klb;j++){
      
      Xnew.slice(index)=LW;
      Xnew.slice(index).col(j)=centre.col(j);
      index=index+1;  
      
      Xnew.slice(index)=LW;
      Xnew.slice(index).col(j)=RW.col(j);
      index=index+1;
      
    }
    
    for (int j=0;j<klb;j++){
      
      Xnew.slice(index)=RW;
      Xnew.slice(index).col(j)=centre.col(j);
      index=index+1;  
      
      Xnew.slice(index)=RW;
      Xnew.slice(index).col(j)=LW.col(j);
      index=index+1;
      
    }
    //Until here, Xnew has been fully updated
    
    X=Xnew;
    
    for (int i=1;i<mnew;i++){
      
      for (int j=0;j<klb;j++){
        
        for (int l=0;l<n;l++){
          
          z=unif_rand();
          
          if(z<=((N-C+1.0)/N)){
            
            X(l,j,i)=X(l,j,i)*(1.0+rSign(2)*alpha);
            
          }
          
          if(X(l,j,i)>ub(j,0)){X(l,j,i)=ub(j,0);}
          if(X(l,j,i)<lb(j,0)){X(l,j,i)=lb(j,0);}
          
        }
        
      }
      
    }
    
    result.set_size(mnew);
    
    if(OC=="A"){
      
      for (int i=0;i<mnew;i++){
        
        result(i,0)=A(X.slice(i));
        
      }
      
    }
    
    if(OC=="D"){
      
      for (int i=0;i<mnew;i++){
        
        result(i,0)=D(X.slice(i));
        
      }
      
    }
    
    if(OC=="G"){
      
      for (int i=0;i<mnew;i++){
        
        result(i,0)=G(X.slice(i));
        
      }
      
    }
    
    C=C+1;
    
  }
  
  indices=sort_index(result.col(0),"ascend");
  
  centre=X.slice(indices(0,0));
  
  return centre;
}
