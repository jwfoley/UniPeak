# This script defines some functions that are useful for dealing with data in UniPeak regions format
# It doesn't actually do anything by itself, but use  source("functions.R")  to load it for whatever you want to do


# read regions file
read.unipeak <- function(file = "", ...) read.table(file, sep = "\t", header = TRUE, row.names = 1, ...)

# write regions file
write.unipeak <- function(x, file = "", ...) write.table(format(x, scientific = FALSE), file, quote = FALSE, sep = "\t", col.names = NA, ...)

# find which column has kurtosis (last column before expression values)
kurtosis.col <- function(x) {
  result <- which(colnames(x) == "kurtosis")
  if (length(result) > 1) {
    stop("kurtosis defined twice")
  } else if (length(result) == 0) {
    stop("no kurtosis column")
  }
  result
}

# isolate the counts from a table
values.unipeak <- function(x) as.matrix(x[,-(1:kurtosis.col(x))])

# isolate everything else
annotations.unipeak <- function(x) as.matrix(x[,1:kurtosis.col(x)])

# apply a matrix function without breaking annotation/kurtosis columns
apply.unipeak <- function(x, f) { # x is table, f is function
  this.kurtosis.col <- kurtosis.col(x)
  cbind(x[,1:this.kurtosis.col], f(x[,(this.kurtosis.col + 1):ncol(x)]))
}

