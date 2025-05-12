---
layout: single
title: "Statistical Matching using Optimal Transport"
author: "Raphaël Jauslin"
date: "2021-10-19"
output: html_document
permalink: /matching/
excerpt: Example of statistical matching using optimal transport.
# classes: wide
---

## Introduction  

In this vignette, we explore how key functions from the package can be used to estimate a contingency table. Our analysis is based on the `eusilc` dataset from the `laeken` package. Each function discussed here is thoroughly explained in the manuscript by Raphaël Jauslin and Yves Tillé (2021), available on [doi:10.1016/j.jspi.2022.12.003](https://doi.org/10.1016/j.jspi.2022.12.003).  

## Contingency Table  

To construct the contingency table, we examine the factor variable `pl030`, which represents economic status, in combination with a discretized version of the equivalized household income, `eqIncome`. The discretization process involves calculating specific percentiles (0.15, 0.30, 0.45, 0.60, 0.75, 0.90) of `eqIncome` and defining categorical intervals based on these values.  

 

{% highlight r %}
library(laeken)
library(sampling)
library(StratifiedSampling)
{% endhighlight %}



<!-- {% highlight text %}
#> Warning: le package 'StratifiedSampling' a été compilé avec la version R 4.2.0
{% endhighlight %} -->



<!-- {% highlight text %}
#> Le chargement a nécessité le package : Matrix
{% endhighlight %}
 -->


{% highlight r %}

data("eusilc")
eusilc <- na.omit(eusilc)
N <- nrow(eusilc)


# Xm are the matching variables and id are identity of the units
Xm <- eusilc[,c("hsize","db040","age","rb090","pb220a")]
Xmcat <- do.call(cbind,apply(Xm[,c(2,4,5)],MARGIN = 2,FUN = disjunctive))
Xm <- cbind(Xmcat,Xm[,-c(2,4,5)])
id <- eusilc$rb030


# categorial income splitted by the percentile
c_income  <- eusilc$eqIncome
q <- quantile(eusilc$eqIncome, probs = seq(0, 1, 0.15))
c_income[which(eusilc$eqIncome <= q[2])] <- "(0,15]"
c_income[which(q[2] < eusilc$eqIncome & eusilc$eqIncome <= q[3])] <- "(15,30]"
c_income[which(q[3] < eusilc$eqIncome & eusilc$eqIncome <= q[4])] <- "(30,45]"
c_income[which(q[4] < eusilc$eqIncome & eusilc$eqIncome <= q[5])] <- "(45,60]"
c_income[which(q[5] < eusilc$eqIncome & eusilc$eqIncome <= q[6])] <- "(60,75]"
c_income[which(q[6] < eusilc$eqIncome & eusilc$eqIncome <= q[7])] <- "(75,90]"
c_income[which(  eusilc$eqIncome > q[7] )] <- "(90,100]"

# variable of interests
Y <- data.frame(ecostat = eusilc$pl030)
Z <- data.frame(c_income = c_income)

# put same rownames
rownames(Xm) <- rownames(Y) <- rownames(Z)<- id

YZ <- table(cbind(Y,Z))
addmargins(YZ)
{% endhighlight %}



{% highlight text %}
#>        c_income
#> ecostat (0,15] (15,30] (30,45] (45,60] (60,75] (75,90] (90,100]   Sum
#>     1      409     616     722     807     935    1025      648  5162
#>     2      189     181     205     184     165     154       82  1160
#>     3      137      90      72      75      59      52       33   518
#>     4      210     159     103      95      74      49       46   736
#>     5      470     462     492     477     459     435      351  3146
#>     6       57      25      28      30      17      11       10   178
#>     7      344     283     194     149     106      91       40  1207
#>     Sum   1816    1816    1816    1817    1815    1817     1210 12107
{% endhighlight %}


## Sampling schemes

Here we set up the sampling designs and define all the quantities we will need for the rest of the vignette. The sample is selected with simple random sampling without replacement and the weights are equal to the inverse of the inclusion probabilities.


{% highlight r %}

# size of sample
n1 <- 1000
n2 <- 500

# samples
s1 <- srswor(n1,N)
s2 <- srswor(n2,N)
  
# extract matching units
X1 <- Xm[s1 == 1,]
X2 <- Xm[s2 == 1,]
  
# extract variable of interest
Y1 <- data.frame(Y[s1 == 1,])
colnames(Y1) <- colnames(Y)
Z2 <- as.data.frame(Z[s2 == 1,])
colnames(Z2) <- colnames(Z)
  
# extract correct identities
id1 <- id[s1 == 1]
id2 <- id[s2 == 1]
  
# put correct rownames
rownames(Y1) <- id1
rownames(Z2) <- id2
  
# here weights are inverse of inclusion probabilities
d1 <- rep(N/n1,n1)
d2 <- rep(N/n2,n2)
  
# disjunctive form
Y_dis <- sampling::disjunctive(as.matrix(Y))
Z_dis <- sampling::disjunctive(as.matrix(Z))
  
Y1_dis <- Y_dis[s1 ==1,]
Z2_dis <- Z_dis[s2 ==1,]

{% endhighlight %}


## Harmonization

Then the harmonization step must be performed. The `harmonize` function returns the harmonized weights. If the true population totals are known, it is possible to use these instead of the estimate made within the function.


{% highlight r %}


re <- harmonize(X1,d1,id1,X2,d2,id2)  

# if we want to use the population totals to harmonize we can use 
re <- harmonize(X1,d1,id1,X2,d2,id2,totals = c(N,colSums(Xm)))

w1 <- re$w1
w2 <- re$w2

colSums(Xm)
{% endhighlight %}



{% highlight text %}
#>      1      2      3      4      5      6      7      8      9     10     11 
#>    476    887   2340    763   1880   1021   2244   1938    558   6263   5844 
#>     12     13     14  hsize    age 
#>  11073    283    751  36380 559915
{% endhighlight %}



{% highlight r %}
colSums(w1*X1)
{% endhighlight %}



{% highlight text %}
#>      1      2      3      4      5      6      7      8      9     10     11 
#>    476    887   2340    763   1880   1021   2244   1938    558   6263   5844 
#>     12     13     14  hsize    age 
#>  11073    283    751  36380 559915
{% endhighlight %}



{% highlight r %}
colSums(w2*X2)
{% endhighlight %}



{% highlight text %}
#>      1      2      3      4      5      6      7      8      9     10     11 
#>    476    887   2340    763   1880   1021   2244   1938    558   6263   5844 
#>     12     13     14  hsize    age 
#>  11073    283    751  36380 559915
{% endhighlight %}



## Optimal transport matching

The statistical matching is done by using the `otmatch` function. The estimation of the contingency table is calculated by extracting the `id1` units (respectively `id2` units) and by using the function `tapply` with the correct weights.


{% highlight r %}

# Optimal transport matching
object <- otmatch(X1,id1,X2,id2,w1,w2)
head(object[,1:3])
{% endhighlight %}



{% highlight text %}
#>         id1    id2    weight
#> 702     702 251702 11.509002
#> 1      1401   1401 13.550397
#> 2506   2506 194004  8.315938
#> 2506.1 2506 324205  2.013395
#> 3001   3001 494702 10.976034
#> 3602   3602 503002 12.651816
{% endhighlight %}



{% highlight r %}

Y1_ot <- cbind(X1[as.character(object$id1),],y = Y1[as.character(object$id1),])
Z2_ot <- cbind(X2[as.character(object$id2),],z = Z2[as.character(object$id2),])
YZ_ot <- tapply(object$weight,list(Y1_ot$y,Z2_ot$z),sum)

# transform NA into 0
YZ_ot[is.na(YZ_ot)] <- 0

# result
round(addmargins(YZ_ot),3)
{% endhighlight %}



{% highlight text %}
#>       (0,15]  (15,30]  (30,45]  (45,60]  (60,75]  (75,90] (90,100]       Sum
#> 1    908.206  732.717  739.153  768.961  886.744  804.436  505.966  5346.183
#> 2    229.633  157.397  125.376  231.835  232.178  166.663  105.011  1248.094
#> 3    111.164  105.834   47.015   68.783   81.041   51.124   51.106   516.067
#> 4     60.987   66.797  104.289  210.126   76.667   92.290   12.085   623.241
#> 5    549.988  566.912  482.875  446.948  400.627  362.297  356.626  3166.273
#> 6      8.577   37.881   14.943   51.063    0.000   10.138    0.000   122.602
#> 7    176.177  164.798  193.780  190.152  119.938  166.566   73.129  1084.540
#> Sum 2044.732 1832.336 1707.432 1967.869 1797.195 1653.514 1103.922 12107.000
{% endhighlight %}



## Balanced sampling

As you can see from the previous section, the optimal transport results generally do not have a one-to-one match, meaning that for every unit in sample 1, we have more than one unit with weights not equal to 0 in sample 2.  The `bsmatch` function creates a one-to-one match by selecting a balanced stratified sampling to obtain a data.frame where each unit in sample 1 has only one imputed unit from sample 2. 



{% highlight r %}

# Balanced Sampling 
BS <- bsmatch(object)
head(BS$object[,1:3])
{% endhighlight %}



{% highlight text %}
#>         id1    id2    weight
#> 702     702 251702 11.509002
#> 1      1401   1401 13.550397
#> 2506   2506 194004  8.315938
#> 3001   3001 494702 10.976034
#> 3602   3602 503002 12.651816
#> 3901.1 3901 294601 10.970414
{% endhighlight %}



{% highlight r %}


Y1_bs <- cbind(X1[as.character(BS$object$id1),],y = Y1[as.character(BS$object$id1),])
Z2_bs <- cbind(X2[as.character(BS$object$id2),],z = Z2[as.character(BS$object$id2),])
YZ_bs <- tapply(BS$object$weight/BS$q,list(Y1_bs$y,Z2_bs$z),sum)
YZ_bs[is.na(YZ_bs)] <- 0
round(addmargins(YZ_bs),3)
{% endhighlight %}



{% highlight text %}
#>       (0,15]  (15,30]  (30,45]  (45,60]  (60,75]  (75,90] (90,100]       Sum
#> 1    950.180  747.323  753.780  706.298  903.800  786.384  498.417  5346.183
#> 2    202.833  138.314  153.513  220.030  237.620  175.001  120.782  1248.094
#> 3     93.911   92.113   52.996   78.145   85.264   42.861   70.776   516.067
#> 4     69.117   58.966  102.611  227.049   68.017   85.771   11.710   623.241
#> 5    516.973  554.095  495.790  464.549  395.341  331.096  408.429  3166.273
#> 6      8.367   37.881   14.943   51.274    0.000   10.138    0.000   122.602
#> 7    171.059  177.551  181.066  213.189  102.669  172.524   66.482  1084.540
#> Sum 2012.440 1806.244 1754.699 1960.535 1792.710 1603.775 1176.597 12107.000
{% endhighlight %}



{% highlight r %}

# With Z2 as auxiliary information for stratified balanced sampling.
BS <- bsmatch(object,Z2)

Y1_bs <- cbind(X1[as.character(BS$object$id1),],y = Y1[as.character(BS$object$id1),])
Z2_bs <- cbind(X2[as.character(BS$object$id2),],z = Z2[as.character(BS$object$id2),])
YZ_bs <- tapply(BS$object$weight/BS$q,list(Y1_bs$y,Z2_bs$z),sum)
YZ_bs[is.na(YZ_bs)] <- 0
round(addmargins(YZ_bs),3)
{% endhighlight %}



{% highlight text %}
#>       (0,15]  (15,30]  (30,45]  (45,60]  (60,75]  (75,90] (90,100]       Sum
#> 1    916.607  733.295  727.348  783.917  893.807  804.205  487.003  5346.183
#> 2    215.298  139.840  115.908  246.175  238.891  195.369   96.613  1248.094
#> 3    105.798  103.158   75.489   55.427   72.337   42.861   60.997   516.067
#> 4     46.193   70.368  114.037  190.878   91.443   98.613   11.710   623.241
#> 5    571.876  569.674  459.916  460.539  378.955  356.515  368.799  3166.273
#> 6      8.367   37.881   14.943   51.274    0.000   10.138    0.000   122.602
#> 7    186.803  174.473  191.122  180.442  130.024  144.693   76.984  1084.540
#> Sum 2050.942 1828.688 1698.763 1968.652 1805.457 1652.393 1102.105 12107.000
{% endhighlight %}


## Prediction  



{% highlight r %}

# split the weight by id1
q_l <- split(object$weight,f = object$id1)
# normalize in each id1
q_l <- lapply(q_l, function(x){x/sum(x)})
q <- as.numeric(do.call(c,q_l))
  
Z_pred <- t(q*disjunctive(object$id1))%*%disjunctive(Z2[as.character(object$id2),])
colnames(Z_pred) <- levels(factor(Z2$c_income))
head(Z_pred)
{% endhighlight %}



{% highlight text %}
#>         (0,15]   (15,30]    (30,45] (45,60] (60,75] (75,90] (90,100]
#> [1,] 0.0000000 0.0000000 1.00000000       0       0       0        0
#> [2,] 0.0000000 0.0000000 0.00000000       1       0       0        0
#> [3,] 0.1949201 0.8050799 0.00000000       0       0       0        0
#> [4,] 0.0000000 0.0000000 0.00000000       0       1       0        0
#> [5,] 0.0000000 0.0000000 1.00000000       0       0       0        0
#> [6,] 0.7749145 0.1455486 0.07953691       0       0       0        0
{% endhighlight %}

