---
layout: single
title: "Sequential balanced sampling"
author: "Raphaël Jauslin"
date: "2022-06-06"
output: html_document
bibliography: ref.bib
permalink: /balseq/
excerpt: Example of sequential spatially balanced sampling.
---


## Introduction

Balanced sampling plays a crucial role in applied statistics. In this vignette, we explain how to use the `balseq` function to select a balanced and spatially distributed sample. For a detailed explanation of the method, refer to [doi:10.1002/env.2776](https://doi.org/10.1002/env.2776).

## Loading Data

We will use the `belgianmunicipalities` dataset from the `sampling` package, which does not contain spatial coordinates. Fortunately, a `GEOjson` file is available on [ArcGIS Hub](https://hub.arcgis.com/datasets/esribeluxdata::belgium-municipalities-1/about). We transform it into an `sf` object and compute the municipalities' centroids to distribute the sample across space.


{% highlight r %}
# to load data
library(sampling)
library(geojsonio)
library(ggplot2) 
library(viridis)
library(rgeos) 
library(sf)
library(rmapshaper)

data("belgianmunicipalities")

# load geojson directly from the url
belg  <- geojson_read("https://opendata.arcgis.com/datasets/9589d9e5e5904f1ea8d245b54f51b4fd_0.geojson",what = "sp")

# simplify the variable and transform it into a sf object
belg <- rmapshaper::ms_simplify(input = belg,keep = 0.01) %>%
  st_as_sf()

coord <- gCentroid(as(belg, "Spatial"), byid = TRUE)

# concatenated file
Belgium <- cbind(belg,belgianmunicipalities,coord)
head(Belgium)
{% endhighlight %}



{% highlight text %}
#> Simple feature collection with 6 features and 26 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 4.216306 ymin: 51.08066 xmax: 4.564865 ymax: 51.37764
#> Geodetic CRS:  WGS 84
#>   OBJECTID   ADMUNAFR   ADMUNADU   ADMUNAGE   Communes CODE_INS arrond    Commune   INS Province Arrondiss
#> 1        1 AARTSELAAR AARTSELAAR AARTSELAAR Aartselaar    11001     11 Aartselaar 11001        1        11
#> 2        2     ANVERS  ANTWERPEN  ANTWERPEN  Antwerpen    11002     11     Anvers 11002        1        11
#> 3        3   BOECHOUT   BOECHOUT   BOECHOUT   Boechout    11004     11   Boechout 11004        1        11
#> 4        4       BOOM       BOOM       BOOM       Boom    11005     11       Boom 11005        1        11
#> 5        5   BORSBEEK   BORSBEEK   BORSBEEK   Borsbeek    11007     11   Borsbeek 11007        1        11
#> 6        6 BRASSCHAAT BRASSCHAAT BRASSCHAAT Brasschaat    11008     11 Brasschaat 11008        1        11
#>    Men04 Women04  Tot04  Men03 Women03  Tot03 Diffmen Diffwom DiffTOT TaxableIncome Totaltaxation averageincome
#> 1   6971    7169  14140   7010    7243  14253     -39     -74    -113     242104077      74976114         33809
#> 2 223677  233642 457319 221767  232405 454172    1910    1237    3147    5416418842    1423715652         22072
#> 3   6027    5927  11954   6005    5942  11947      22     -15       7     167616996      50739035         29453
#> 4   7640    8066  15706   7535    7952  15487     105     114     219     186075961      46636930         21907
#> 5   4948    5328  10276   4951    5322  10273      -3       6       3     143225590      40564374         26632
#> 6  18142   18916  37058  18217   18903  37120     -75      13     -62     533368826     153629397         30574
#>   medianincome        x        y                       geometry
#> 1        23901 4.382005 51.13223 MULTIPOLYGON (((4.400451 51...
#> 2        17226 4.369578 51.26067 MULTIPOLYGON (((4.368136 51...
#> 3        21613 4.516463 51.16576 MULTIPOLYGON (((4.530071 51...
#> 4        17537 4.371434 51.09387 MULTIPOLYGON (((4.385267 51...
#> 5        20739 4.488162 51.19138 MULTIPOLYGON (((4.509002 51...
#> 6        21523 4.500281 51.30940 MULTIPOLYGON (((4.540815 51...
{% endhighlight %}

## Data visualization

We visualize Belgian municipalities with their average income using ggplot2.

{% highlight r %}
p <- ggplot()+
  geom_sf(data = Belgium,aes(fill = averageincome),size = 0.1)+
  scale_fill_viridis_c(option = "G")
p
{% endhighlight %}

<img src="/figs/sequential_balanced/unnamed-chunk-3-1.png" title="center" alt="center" width="100%" />

## Inclusion probabilites

A good sample should maintain the population's characteristics. By defining proportional inclusion probabilities, we ensure better representativity. We set up here the inclusion probabilities equal with sum equal to 50. i.e. the sample will contain 50 units.



{% highlight r %}

N <- nrow(Belgium) # population total
n <- 50 # sample size

# variable of interest
y <- belgianmunicipalities$averageincome

# auxiliary variables
Xaux <- cbind(belgianmunicipalities$Tot04,
              belgianmunicipalities$Women04,
              belgianmunicipalities$TaxableIncome,
              belgianmunicipalities$Diffmen,
              belgianmunicipalities$Diffwom)


# inclusion probabilities
pik <- rep(n/N,N)

Xaux <- cbind(pik,Xaux) # add pik to fixed sample size
Xspread <- coord
{% endhighlight %}


## Balanced sampling

We compare balanced sampling `balseq` with two other methods:  `samplecube` and simple random sampling `srswor`. The percentage deviation from auxiliary totals helps evaluate the balance quality.


{% highlight r %}

library(StratifiedSampling)

s <- balseq(pik,Xaux)
s_cube <- samplecube(Xaux,pik,comment = FALSE)
s_srs <- srswor(n,N)


TOT <- colSums(Xaux)
EST1 <- colSums(Xaux[s,]/pik[s])
EST2 <- colSums(Xaux[s_cube == 1,]/pik[s_cube == 1])
EST3 <- colSums(Xaux[s_srs == 1,]/pik[s_srs == 1])


100* (EST1 - TOT)/TOT
{% endhighlight %}



{% highlight text %}
#>       pik                                                   
#>  0.000000  4.908396  4.890552  6.873691 14.225213  6.440032
{% endhighlight %}



{% highlight r %}
100* (EST2 - TOT)/TOT
{% endhighlight %}



{% highlight text %}
#>        pik                                                        
#>   0.000000  -4.068363  -4.153367  -8.105537 -19.212958 -14.508554
{% endhighlight %}



{% highlight r %}
100* (EST3 - TOT)/TOT
{% endhighlight %}



{% highlight text %}
#>       pik                                                   
#>   0.00000  12.99101  12.52005  12.72183 -13.28123 -21.77427
{% endhighlight %}


## Spread sampling

To incorporate spatial distribution, we use geographic coordinates as a matrix to the argument `Xspread` of the function. Here `coord` is an output of the function `gCentroid` which is by construction an `S4` object. To get the `data.frame` that are encapsulated inside, we simply use the `@coords` operator.


{% highlight r %}

s <-  balseq(pik,
             Xaux,
             Xspread = as.matrix(Xspread@coords))
p <- p + 
  geom_point(data = Belgium,aes(x = x,y = y),shape = 1,alpha = 0.9)+
  geom_point(data = Belgium[s,],aes(x = x,y = y),colour = "red")+
  scale_fill_viridis_c(option = "G")
{% endhighlight %}


{% highlight r %}
p
{% endhighlight %}

<img src="/figs/sequential_balanced/unnamed-chunk-6-1.png" title="center" alt="center" width="100%" />

