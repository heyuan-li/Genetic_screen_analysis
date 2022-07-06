## Testing the Enrichment score for project DRIVE ##
setwd("/Users/Heyuan/Documents/Local/PublicData/DRIVE/RSA/Post_RSA")
library(plyr)
library(ggplot2)
library(ggrepel)

# read data
genelist = read.csv("genelist.csv", header = F, stringsAsFactors = F)[,1]
data_RSA = read.table("Compiled_DRIVE_candidates.txt", header = T, stringsAsFactors = F)
anno = read.csv('CCLE_sample_info_file_2012-10-18.csv', header = T)

rownames(data_RSA) = gsub(pattern = "_.*", replacement = "", x = data_RSA$name)
rownames(rsa) = rsa$gene
colnames(rsa) = toupper(colnames(rsa))
data_rsa = as.data.frame(t(rsa[,-1]))
data_rsa$name = rownames(data_rsa)

namelist = intersect(rownames(data_RSA), rownames(data_rsa))
index = apply(data_rsa,2,function(x){sum(is.na(x)) == 0})
data_rsa = subset(data_rsa, name %in% namelist)
data_RSA$name = rownames(data_RSA)
data_RSA = subset(data_RSA, name %in% namelist)
dat = data.frame(join(data_RSA, data_rsa, by = "name"))

# waterful plot for sensitivity score #
gene = "U2AF1"
p_score = ggplot(data_rsa, aes(x=reorder(name, get(gene)), y = get(gene))) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = c('grey', 'red')) +
  #ylim(-30,0) +
  ylab(paste(gene,'RSA logP', sep = ' ')) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
p_score

# correlation test w/ RSA score for each gene #
cor_score = lapply(genelist,function(x){
  p = ggplot(dat, aes(x = get(x), y = get(paste(x,".1",sep = "")))) + 
    geom_point() +
    ylab(paste(x,"custom RSA logP",sep = " ")) +
    xlab(paste(x, "RSA score", sep = " "))
  cor = cor.test(dat[[x]], dat[[paste(x,".1",sep = "")]], method = "pearson")
  list(p,cor)
})
rvr = sapply(seq(cor_score), function(x){cor_score[[x]][[2]]$estimate})
plot = hist(rvr, breaks = 30)


# correlation of MYC correlation #
cor_MYC = t(sapply(genelist,function(x){
  cor_1 = cor.test(dat$MYC, dat[[paste(x)]], method = "pearson")
  cor_2 = cor.test(dat$MYC.1, dat[[paste(x,".1",sep = "")]], method = "pearson")
  c(cor_1$estimate, cor_1$p.value, cor_2$estimate, cor_2$p.value)
}))
colnames(cor_MYC) = c("dr_rho","dr_pval","cr_rho","cr_pval")

p = ggplot(as.data.frame(cor_MYC), aes(x = dr_rho, y = cr_rho)) + 
  geom_point() +
  geom_label_repel(label = rownames(cor_MYC)) +
  xlim(-0.25,0.3) +
  ylim(-0.25,0.3) +
  ylab("downloaded MYC_cor") +
  xlab("custom MYC_cor") 
cor = cor.test(cor_MYC[,1], cor_MYC[,3], method = "pearson")

# gene correlation test #

res = t(sapply(colnames(data_rsa)[1:ncol(data_rsa)-1], function(x){
    t = cor.test(data_rsa$MYC, data_rsa[[x]], na.action = na.omit, method = "pearson")
    c(t$estimate, t$p.value)
}))
colnames(res) = c("rho","p_val")
res_sorted = res[order(as.data.frame(res)$rho, decreasing = T),]
plot = hist(res_sorted[,1], breaks = 50)

write.table(res_sorted, "DRIVE_MYC_correlation.txt", quote = F, sep = "\t", row.names = T)
