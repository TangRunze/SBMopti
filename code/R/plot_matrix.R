setwd("/Users/Runze/Dropbox/GitHub/SBMopti_neg")

require(Matrix)

fname <- "matrix_sim.txt"
data <- as.matrix(read.table(fname))

image(Matrix(data))

