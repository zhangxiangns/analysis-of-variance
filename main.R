setwd('/Users/zhangxiang/GitHub/两正态总体均值向量检验&&方差分析')
source('两正态总体均值向量检验.R')
data = read.csv('Japan.csv')
bi.mean.test(data)

source('多元方差分析.R')
data = read.csv('example3_3_2.csv')
analysis_variance(data)
