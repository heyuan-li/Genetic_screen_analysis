setwd("/Users/Heyuan/Documents/Local/PublicData/DRIVE/RSA")
library(plyr)
library(dplyr)
library(fastmatch)
library(edgeR)

# Note: this analysis is a replicate analysis of the project DRIVE paper (Robert et al, Cell, 2017). It starts with shRNA level count data "DriveCountData.RDS" from the original paper.

## generating inputs: dataframes for each of the cell lines for RSA ##

# reading and mapping the shRNA to genes #
data=readRDS('DriveCountData.RDS')
colnames(data) = c('shrna', 'sequence', 'plasmid_count', 'read_count', 'experiment_id', 'pool', 'cell_line')
shrna_annotations = read.table('drive_shrna_table_annotated.txt',header=T,sep='\t')[,c(1,3)]
merge_data = left_join(data, shrna_annotations, by = 'shrna')
rows = which(!is.na(merge_data$plasmid_count))
merge_data = merge_data[rows,]

# split data to each cellline and each pool #
data_list = split(merge_data, merge_data$cell_line)
data_list = lapply(seq(data_list), function(x){
  split(data_list[[x]], data_list[[x]]$pool)
})

# TMM normalization of read-plasmid count pair and fitting negative binomial for log fold change#
data_fit = list()
for (i in seq(length(data_list))) {
  res = lapply(seq(data_list[[i]]), function(x){
    dat = subset(data_list[[i]][[x]], select = c('plasmid_count','read_count'))
    d = DGEList(
      counts = dat,
      norm.factors = calcNormFactors(dat, method = "TMM"),
      genes = subset(data_list[[i]][[x]], select = c("shrna","gene")),
      group = c("plasmid_count","read_count")
    )
    fit = exactTest(d, pair = c("plasmid_count","read_count"), dispersion = 0.2, prior.count = 12)
    cbind(fit$genes, fit$table)
  })
  data_fit[[i]] = rbind(res[[1]], res[[2]], res[[3]])
  print(i)
}

test = aggregate(data_fit[[x]][3:5], list(data_fit[[x]]$shrna, data_fit[[x]]$gene), mean)
a = aggregate(orderd[3:5], list(orderd$shrna, orderd$gene), mean)

# Taking averages of the repilcated lines #
data_ave = lapply(seq(data_fit), function(x){
  if (nrow(data_fit[[x]]) > nrow(shrna_annotations)){
    ave = aggregate(data_fit[[x]][3:5], list(data_fit[[x]]$shrna, data_fit[[x]]$gene), mean)
    colnames(ave) = c('shrna','gene',colnames(ave)[3:5])
    ave
    }
  else {data_fit[[x]]}
})

# quantile normalization across cell lines #

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

# extract logFC and do quantile normalization
data_logFC = lapply(seq(data_ave), function(x){
  subset(data_ave[[x]], select = c('shrna','gene','logFC'))})
data_logFC = Reduce(function(x, y) {left_join(x,y,by = c('shrna','gene'))}, data_logFC)
data_norm = quantile_normalisation(data_logFC[3:ncol(data_logFC)])


# filter the genes that has 5 or more hairpins and unique mapping
genelist = as.character(subset(as.data.frame(table(data_norm$gene)), Freq >= 5)$Var1)
#genelist = genelist[-grep(",",genelist)]
data_filtered = subset(data_norm, gene %in% genelist)

# formating for running RSA pakage #
colnames(sub_data) = c('Gene_ID', 'Well_ID', 'Score', 'cell_line')
data_list = split(sub_data, sub_data$cell_line)

# writing .csv files for RSA.r to read #
for (i in 1:length(data_list))
{
  write.csv(data_list[[i]], paste("input_", i, ".csv", sep = ""),
            row.names = FALSE, quote = FALSE)
}

## formatting outputs: data 

output = lapply(1:length(data_list), function(x){
  dat = read.csv(paste("output_", x, ".csv", sep = ""), header = T, stringsAsFactors = F)[,c(1,4,5)]
  dat = dat[!duplicated(dat),]
  #dat = dat[match(dat$Gene_ID, genelist),]
  colnames(dat) = c("gene", "cell_line", paste(dat[2,2]))
  dat
})

rsa = join_all(output, by = "gene")
rsa = rsa[-grep("cell_line", colnames(rsa))]

write.table(rsa, "RSA_output_compiled.txt", sep = "\t", quote = F, row.names = F)

#plotting
plot(density(as.numeric(rsa[126,-1])))


