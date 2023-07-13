#' Find test statistic value from a given data and corresponding critical value
#'
#' More detailed description
#'
#' @param data a real matrix
#' @param alpha real number between 0 and 1 called significance level
#' @param k positive integer
#' @return numeric vector
#'
#' @examples
#' k=3;N=c(20,10,10);S=c(1,2,6)
#' g=NULL
#' for(i in 1:k)
#' {
#'  g[[i]]=rnorm(N[i],0,sqrt(S[i]))
#' }
#' data=g
#' onewayAMinT(data,k,0.05)
#' @export
onewayAMinT<-function(data,k,alpha)
{
  fun1<-function(data,k)
  {
    N=unlist(rbind(lapply(data,length)))
    M=unlist(rbind(lapply(data,mean)))
    S=unlist(rbind(lapply(data,var)))
    #k=NROW(data)
    V=NULL
    T=NULL
    for(i in 1:k-1)
    {
      V[i]=sqrt((S[i]/N[i])+(S[i+1]/N[i+1]))
      T[i]=(M[i+1]-M[i])/V[i]
    }
    value=min(T)
    return(value)
  }
  fun2<-function(S,N,nu,d,k)
  {
    x=NULL
    for ( i in 1:k-2)
    {
      temp=-sqrt(nu[i+2]*nu[i])*S[i+1]
      x[i]=temp
    }
    m <- diag(d^2)
    m[row(m)-col(m)== 1]<-x
    m[row(m)-col(m)==-1]<-x
    D=diag(d)
    R=solve(D)%*%m%*%solve(D)
    mu=rep(0,k-1)
    g=MASS::mvrnorm(1,mu,R)
    T1=min(g)
    return(T1)
  }
  fun3<-function(S,N,nu,d,k,alpha)
  {
    x<-replicate(5000,fun2(S,N,nu,d,k))
    y<-sort(x,decreasing=FALSE)
    m=(1-alpha)*5000
    c<-y[m]
    return(c)
  }
  data1<-lapply(data, function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(data1,length)))
  M=unlist(rbind(lapply(data1,mean)))
  S=unlist(rbind(lapply(data1,var)))
  #k=NROW(data1)
  Ntot=sum(N)
  nu=N/Ntot
  d=NULL
  for(i in 1:k-1)
  {
    temp=sqrt((nu[i]*S[i+1])+(nu[i+1]*S[i]))
    d[i]=temp
  }
  set.seed(490)
  statistic_value<-fun1(data1,k)
  crit_value<-fun3(S,N,nu,d,k,alpha)
  result<-c(statistic_value, crit_value)
  print("test statistic value and critical value")
  print(result)
  r1=result[1];r2=result[2]
  if(r1>r2)
  {
    print("Null hypothesis is rejected")
  }
  else
  {
    print("Null hypothesis is not rejected")
  }
}
