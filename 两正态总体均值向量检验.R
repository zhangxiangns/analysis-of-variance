source('/Users/zhangxiang/GitHub/协方差检验/协方差检验.R')

### Population: N(mu1,sigma),N(mu2,sigma)
### Hypothesis test: H0:mu1=mu2 H1:mu1!=mu2
### Attention: the sigmas are considered as the same.So, the first step is to test whether
###            the sigmas are same or not.If the sigmas are not the same, the function'll
###            warn.
### input:
### 1)my_data--my_data is a data frame such as:
####### X1   X2   X3   X4   group
####### x11  x12  x13  x14  1
####### x21  x22  x23  x24  1
####### ...  ...  ...  ...  ...
####### xn1  xn2  xn3  xn4  2
### 2)multiple--is the population tested a multiple normal population or a single one.
### 3)alpha--the confidence level to decide whether the cov.test pass to show that the
###   sigmas are considered as the same.
bi.mean.test=function(my_data,multiple=TRUE,alpha=0.05){
  groups = unique(my_data$group)
  data1 = subset(my_data,group==groups[1],-group)
  data2 = subset(my_data,group==groups[2],-group)
  data1 = as.matrix(data1)
  data2 = as.matrix(data2)
  n = nrow(data1)
  m = nrow(data2)
  p = ncol(data1)
  warn = 0
  
  if (multiple==FALSE){
    x_bar = mean(data1)
    y_bar = mean(data2)
    sx = var(data1)
    sy = var(data2)
    t = ((x_bar-y_bar)/sqrt(1/n+1/m))/sqrt(((n-1)*sx+(m-1)*sy)/(n+m-2))
    p_cov = var.test(data1,data2,alternative = 'two.sided')$p.value
    if (p_cov<alpha){
      cat('warning:the sigmas are not the same!The p value of cov test is',p_cov,'\n')
    }
    if (p_cov>=alpha){
      cat('Cov test pass.The sigmas could be considered as the same.The p value of cov test is',p_cov,'\n')
    }
    p_value = 1 - pt(t,n+m-2)
    statistics = t
  }
  if (multiple==TRUE){
    X_bar = apply(data1,2,mean)
    Y_bar = apply(data2,2,mean)
    A1 = (n-1)*var(data1)
    A2 = (m-1)*var(data2)
    T_square = (n+m-2)*n*m/(n+m)*t(X_bar-Y_bar)%*%solve(A1+A2)%*%(X_bar-Y_bar)
    F = (n+m-2-p+1)/(n+m-2)/p*T_square
    p_cov = cov.test.bi(data1,data2)[['p_value']]
    if (p_cov<alpha){
      cat('warning:the sigmas are not the same!The p value of cov test is',p_cov,'\n')
    }
    if (p_cov>=alpha){
      cat('Cov test pass.The sigmas could be considered as the same.The p value of cov test is',p_cov,'\n')
    }
    p_value = 1 - pf(F,p,n+m-p-1)
    statistics = F
  }
  show = list()
  show['statistic'] = statistics
  show['p_value'] = p_value
  return(show)
}
