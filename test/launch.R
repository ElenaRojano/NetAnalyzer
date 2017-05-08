#! /usr/bin/env Rscript

data <- read.table('adjacency_matrix.txt')
pcc <- cos(data)
write.table(pcc, file = 'pcc_adjacency_results.txt', sep = "\t")
