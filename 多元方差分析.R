source('/Users/zhangxiang/GitHub/协方差检验/协方差检验.R')


### Population:N(mu1,sigma),N(mu2,sigma),...,N(muk,sigma).
### Hypothesis test:
### H0:mu1=mu2=...=muk   H1:mu1,mu2,...,muk are all not the same.
### note:the sigma of populations is considered as the same.
### input:
### 1)my_data--is a data frame such as
#####################################
####### X1   X2   X3   X4   group ###
####### x11  x12  x13  x14  1     ###
####### x21  x22  x23  x24  1     ###
####### ...  ...  ...  ...  ...   ###
####### xn1  xn2  xn3  xn4  k     ###
#####################################
### Attention: the classification column needs to be 'group'!
### 2)multiple--is the population tested a multiple normal population or a single one.
### 3)alpha--the confidence level to decide whether the cov.test pass to show that the
###   sigmas are considered as the same.

analysis_variance  = function(my_data,multiple=TRUE,alpha=0.05){
  groups = unique(my_data$group)
  k = length(groups)
  data = list()
  n = list()
  N = nrow(my_data)
  p = ncol(my_data)-1
  for (i in 1:k){
    label = groups[i]
    sub_data = subset(my_data,group==label,-group)
    sub_data = as.matrix(sub_data)
    data[[i]] = sub_data
    n[[i]] = nrow(sub_data)
  }
  if (multiple==FALSE){
   x_sum = 0
   ssa = 0
   sse = 0
   for (i in 1:k){
     x_sum = x_sum + sum(data[[i]])
     sse = sse + var(data[[i]])*(n[[i]]-1)
     }
   x_bar = x_sum/N
   for (i in 1:k){
     ssa = ssa + n[[i]]*(mean(data[[i]])-x_bar)^2
   }
   F = (ssa/(k-1))/(sse/(n-k))
   p_value = pf(F,k-1,n-k)
   statistics = F
  }
  if (multiple==TRUE){
    A = 0
    X_sum = 0
    for (i in 1:k){
      A = A + (n[[i]]-1)*var(data[[i]])
      X_sum = X_sum + apply(data[[i]],2,sum)
    }
    X_bar = X_sum/N
    B = 0
    for (i in 1:k){
      x_bar = apply(data[[i]],2,mean)
      B = B + n[[i]]*(x_bar-X_bar)%*%t(x_bar-X_bar)
    }
    V = det(A)/det(A+B)
    
    r = (N-k)-1/2*(p-k+2)
    Chi_square = -r*log(V)
    p_value =1- pchisq(Chi_square,p*(k-1))
    statistics = Chi_square
    p_cov = cov.test.multi(my_data)[['p_value']]
    if (p_cov<alpha){
      cat('warning:the sigmas are not the same!The p value of cov test is',p_cov,'\n')
    }
    if (p_cov>=alpha){
      cat('Cov test pass.The sigmas could be considered as the same.The p value of cov test is'
          ,p_cov,'\n')
    }
    
  }
  show = list()
  show['statistic'] = statistics
  show['p_value'] = p_value
  return(show)
  
}