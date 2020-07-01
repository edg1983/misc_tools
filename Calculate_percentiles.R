#Suppose a data frame mydf with numeric value in column V4
perc.rank <- function(x) trunc(rank(x))/length(x)
mydf <- within(mydf, perc_dist <- perc.rank(V4))
