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
#' k=4;N=c(20,30,20,30);S=c(1,1,2,1)
#' g=NULL
#' for(i in 1:k)
#' {
#'  g[[i]]=rnorm(N[i],0,sqrt(S[i]))
#' }
#' data=g
#' onewayLRT(data,k,0.05)
#'
#' @export
onewayLRT<-function(data,k,alpha)
{
  fun1<-function(data,k)
  {
    #data1<-lapply(data, function(col)col[!is.na(col)])
    g=data
    N=unlist(rbind(lapply(data,length)))
    M=unlist(rbind(lapply(data,mean)))
    S=unlist(rbind(lapply(data,var)))
    mu0=mean(M)
    s0=S;x1=mu0;s1=S;x2=M
    repeat
    {
      u=NULL
      for(i in 1:k)
      {
        tmp1=N[i]/s0[i]
        u[i]=tmp1
      }
      neu=NULL
      for(i in 1:k)
      {
        tmp2=u[i]*M[i]
        neu[i]=tmp2
      }
      mun=sum(neu)/sum(u)
      S0=NULL
      for(i in 1:k)
      {
        tmp4=sum((g[[i]]-mun)^2)/N[i]
        S0[i]=tmp4
      }
      diff=abs(mun-x1)
      if(diff<=0.00001)
      {
        break
      }
      x1=mun
      s0=S0
    }
    repeat
    {
      w=NULL
      for(i in 1:k)
      {
        tmp5=N[i]/s1[i]
        w[i]=tmp5
      }
      tmp6=Iso::pava(M,w)
      mu=tmp6
      S1=NULL
      for(i in 1:k)
      {
        tmp7=sum((g[[i]]-mu[i])^2)/N[i]
        S1[i]=tmp7
      }
      diff2=max(abs(mu-x2))
      if(diff2<=0.00001)
      {
        break
      }
      x2=mu;s1=S1
    }
    rto=S1/S0;ratio=NULL
    for(i in 1:k)
    {
      rslt=(rto[i])^(N[i]/2)
      ratio[i]=rslt
    }
    value=prod(ratio)
    return(value)
  }
  fun2<-function(N,S,k)
  {
    g=NULL
    for(i in 1:k)
    {
      g[[i]]=rnorm(N[i],0,sqrt(S[i]))
    }
    value=fun1(g,k)
    return(value)
  }
  fun3<-function(N,S,k,alpha)
  {
    x<-replicate(5000,fun2(N,S,k))
    y<-sort(x,decreasing=FALSE)
    m=(alpha)*5000
    c<-y[m]
    return(c)
  }
  data1<-lapply(data, function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(data1,length)))
  S=unlist(rbind(lapply(data1,var)))
  set.seed(911)
  test_statistic<-fun1(data1,k)
  crit_value<-fun3(N,S,k,alpha)
  result<-c(test_statistic,crit_value)
  print("test statistic value and critical value")
  print(result)
  r1=result[1];r2=result[2]
  if(r1<r2)
  {
    print("Null hypothesis is rejected")
  }
  else
  {
    print("Null hypothesis is not rejected")
  }
}
