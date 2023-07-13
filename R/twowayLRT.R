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
#' twowayLRT(data,3,2, 0.05)
#' @export
twowayLRT<-function(data,a,b,alpha)
{
  fun1<-function(data,a,b)
  {
    g=data
    N=unlist(rbind(lapply(data,length)))
    M=unlist(rbind(lapply(data,mean)))
    S=unlist(rbind(lapply(data,var)))
    m=mean(M)
    M1=NULL
    for(i in 1:a)
    {
      rst1=sum(M[(b*(i-1)+1):(b*(i-1)+b)])/b
      M1[i]=rst1
    }
    M2=NULL
    rl=NULL
    for(j in 1:b)
    {
      for(i in 1:a)
      {
        rst1=M[(i-1)*b+j]
        rl[i]=rst1
      }
      M2[j]=sum(rl)/a
    }
    x14=M1;x15=M2-m;S1=S
    repeat
    {
      U=N/S1
      um=U*M
      ub=NULL;ubb=NULL;tmp4=NULL;tmp5=NULL
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          rst2=U[(i-1)*b+j]*x15[j]
          rst3=U[(i-1)*b+b]*x15[j]
          tmp4[j]=rst2
          tmp5[j]=rst3
        }
        ub[i]=sum(tmp4[1:b-1])
        ubb[i]=sum(tmp5[1:b-1])
      }
      m1=NULL
      for(i in 1:a)
      {
        rst1=sum(um[(b*(i-1)+1):(b*(i-1)+b)])
        rst5=sum(U[(b*(i-1)+1):(b*(i-1)+b)])
        m1[i]=(rst1-ub[i]+ubb[i])/rst5
      }
      w=NULL
      for(i in 1:a)
      {
        rst1=sum(U[(b*(i-1)+1):(b*(i-1)+b)])
        w[i]=rst1
      }
      alpha1=Iso::pava(m1,w)
      tmp1=NULL;ua=NULL
      tmp2=NULL;Um=NULL
      for(j in 1:b-1)
      {
        for( i in 1:a)
        {
          rst1=U[(i-1)*b+j]*alpha1[i]
          rst2=um[(i-1)*b+j]
          tmp1[i]=rst1
          tmp2[i]=rst2
        }
        ua[j]=sum(tmp1)
        Um[j]=sum(tmp2)
      }
      q=NULL
      for(i in 1:a)
      {
        rst1=(M[(i-1)*b+b]-alpha1[i])*(U[(i-1)*b+b])
        q[i]=rst1
      }
      v=NULL
      for(j in 1:b-1)
      {
        rst1=Um[j]-ua[j]-sum(q)
        v[j]=rst1
      }
      u.=NULL
      tmp3=NULL
      for(j in 1:b)
      {
        for(i in 1:a)
        {
          rst1=U[(i-1)*b+j]
          tmp3[i]=rst1
        }
        u.[j]=sum(tmp3)
      }
      d=NULL
      for(j in 1:(b-1))
      {
        rst1=u.[j]
        d[j]=rst1
      }
      Unit=u.[b]*matrix(1,nrow=b-1,ncol=b-1)
      Q=(diag(d,nrow=b-1,ncol=b-1)+Unit)
      B=solve(Q)%*%v
      Bb=-sum(B[1:b-1])
      beta=c(B,Bb)
      s=NULL
      sum((g[[1]]-alpha1[1]-beta[1])^2)/N[1]
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          result=sum((g[[(i-1)*b+j]]-alpha1[i]-beta[j])^2)/N[(i-1)*b+j]
          s[(i-1)*b+j]=result
        }
      }
      diff1=max(abs(x14-alpha1));diff2=max(abs(x15-beta))
      if(diff1<=0.0000001 & diff2<=0.0000001)
      {
        break
      }
      S1=s
      x14=alpha1;x15=beta
    }


    x16=M2;S0=S
    repeat
    {
      U0=N/S0
      um0=U0*M
      tmp1=NULL;tmp2=NULL;Um0=NULL;u0=NULL;beta0=NULL
      for(j in 1:b)
      {
        for( i in 1:a)
        {
          rst1=U0[(i-1)*b+j]
          rst2=um0[(i-1)*b+j]
          tmp1[i]=rst1
          tmp2[i]=rst2
        }
        Um0[j]=sum(tmp2)
        u0[j]=sum(tmp1)
        beta0[j]=(Um0[j])/u0[j]
      }
      s0=NULL
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          result=sum((g[[(i-1)*b+j]]-beta0[j])^2)/N[(i-1)*b+j]
          s0[(i-1)*b+j]=result
        }
      }
      diff=max(abs(x16-beta0))
      if(diff<=0.0000001)
      {
        break
      }
      sum((g[[6]]-beta0[2])^2)/N[6]
      S0=s0;x16=beta0
    }
    #((g[[1]]-beta0[1])^2)
    #sum((g[[1]]-beta0[1])^2)/N[1]
    rto=S1/S0;ratio=NULL
    for(i in 1:(a*b))
    {
      rslt=(rto[i])^(N[i]/2)
      ratio[i]=rslt
    }
    value=prod(ratio)
    return(value)
  }
  fun2<-function(N,S,a,b)
  {
    g=NULL
    for(i in 1:(a*b))
    {
      g[[i]]=rnorm(N[i],0,sqrt(S[i]))
    }
    value=fun1(g,a,b)
    return(value)
  }
  fun3<-function(N,S,a,b,alpha)
  {
    x<-replicate(5000,fun2(N,S,a,b))
    y<-sort(x,decreasing=FALSE)
    m=(alpha)*5000
    c<-y[m]
    return(c)
  }
  data1<-lapply(data, function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(data1,length)))
  S=unlist(rbind(lapply(data1,var)))
  set.seed(419)
  test_statistic<-fun1(data1,a,b)
  crit_value<-fun3(N,S,a,b,alpha)
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
