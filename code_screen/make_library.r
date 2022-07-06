#create fasta library
# Orgininal Author Nicholas J. Neill

#input file should have 3 columns (shRNA, Gene, Sequence)
#$1 is path to library (no extension)

args = commandArgs(trailingOnly = TRUE)
library = as.matrix(read.table(paste(args[1], '.txt', sep = ''), sep = '\t', stringsAsFactors = F))
library[,1] = paste(library[,1], library[,2], sep = '-')
library = library[,c(1,3)]
library[,1] = sapply(library[,1], function(x) paste('>', x, sep = ''))
library = unlist(split(library, row(library)))
writeLines(library, con = paste(args[1], '.fasta', sep = ''))
