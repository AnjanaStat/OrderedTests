#' Find test statistic value and corresponding critical value from a given data set and level of significance
#'
#' More detailed description
#'
#' @param data a real matrix
#' @param alpha real number between 0 and 1 called significance level
#' @param k positive integer
#' @return numeric vector
#'
#' @examples
#' k=3;N=c(40,20,10);S=c(1,3,5)
#' g=NULL
#' for(i in 1:k)
#' {
#'  g[[i]]=rnorm(N[i],0,sqrt(S[i]))
#' }
#' data=g
#' onewayMaxT(data,k,0.05)
#' @export
onewayMaxT<-function(data,k,alpha)
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
    value=max(T)
    return(value)
  }
  fun2<-function(S,N,k)
  {
    g=NULL
    for(i in 1:k)
    {
      g[[i]]=rnorm(N[i],0,sqrt(S[i]))
    }
    data=g
    T1=fun1(data,k)
    return(T1)
  }
  fun3<-function(S,N,k,alpha)
  {
    x<-replicate(5000,fun2(S,N,k))
    y<-sort(x,decreasing=FALSE)
    m=(1-alpha)*5000
    c<-y[m]
    return(c)
  }
  data1<-lapply(data, function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(data1,length)))
  S=unlist(rbind(lapply(data1,var)))
  #k=NROW(data1)
  set.seed(439)
  statistic_value<-fun1(data1,k)
  crit_value<-fun3(S,N,k,alpha)
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

