zscore <- function(x, na.rm=TRUE){
  (x-mean(x, na.rm=na.rm))/sd(x, na.rm=na.rm)
}


