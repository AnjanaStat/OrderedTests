#' Find test statistic value from a given data and corresponding critical value
#'
#' More detailed description
#'
#' @param data a real data set considering only numerical values
#' @param a number of levels of the row factor
#' @param b number of levels of the column factor
#' @param alpha real number between 0 and 1 called significance level
#'
#' @return numeric vector
#'
#' @examples
#' a=3;b=2;N=c(20,30,20,30,40,60);S=c(1,1,2,1,4,2)
#' g=NULL
#' for(i in 1:(a*b))
#' {
#'  g[[i]]=rnorm(N[i],0,sqrt(S[i]))
#' }
#' data=g
#' twowayMaxT(data,3,2,0.05)
#' @export
twowayMaxT<-function(data,a,b,alpha)
{
  fun1<-function(M,S,a,b)
  {
    T=NULL
    for(i in 1:a-1)
    {
      result1=sqrt(S[i]+S[i+1])
      result2=(M[i+1]-M[i])/result1
      T[i]=result2
    }
    value=max(T)
    return(value)
  }
  fun2<-function(S,N,a,b)
  {
    M=NULL;v=NULL
    r=S/N
    k=a*b
    for (i in 1:k)
    {
      data=(rnorm(N[i],0,sqrt(S[i])))
      M[i]=mean(data)
      v[i]=var(data)
    }
    m=NULL
    for(i in 1:a)
    {
      result=sum(M[(b*(i-1)+1):(b*(i-1)+b)])/b
      m[i]=result
    }
    s=NULL
    for(i in 1:a)
    {
      result=sum(r[(b*(i-1)+1):(b*(i-1)+b)])/(b^2)
      s[i]= result
    }
    value=fun1(m,s,a,b)
    return(value)
  }
  fun3<-function(S,N,a,b,alpha)
  {
    x<-replicate(5000,fun2(S,N,a,b))
    y<-sort(x,decreasing=FALSE)
    m=(1-alpha)*5000
    c<-y[m]
    return(c)
  }
  data1<-lapply(data, function(col)col[!is.na(col)])
  #Var_data1<-lapply(Var_data,function(col)col[!is.na(col)] )
  N=unlist(rbind(lapply(data1,length)))
  M=unlist(rbind(lapply(data1,mean)))
  S=unlist(rbind(lapply(data1,var)))
  r=S/N
  v=NULL
  for(i in 1:a)
  {
    result=sum(r[(b*(i-1)+1):(b*(i-1)+b)])/(b^2)
    v[i]= result
  }
  m=NULL
  for(i in 1:a)
  {
    result=sum(M[(b*(i-1)+1):(b*(i-1)+b)])/b
    m[i]=result
  }
  set.seed(40)
  statistic_value<-fun1(m, v, a, b)
  crit_value<-fun3(S,N,a,b,alpha)
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
