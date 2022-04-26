
suppressMessages(library("DESeq2"))

# take all samples htseq-count
raw_count <- "htseq_raw.txt"
cnt_raw <- read.table(raw_count, sep = "\t", header = T, row.names = 1)

# remove genes in lowly expression
cnt_expressed <- cnt_raw[rowMeans(cnt_raw > 10) > 0.2, ]
expressed_genes <- row.names(cnt_expressed)

# library-size normalization using DESeq2
htseq_mat <- as.matrix(cnt_raw)
col <- data.frame(Samples = colnames(cnt_raw))
dds <- DESeqDataSetFromMatrix(htseq_mat, colData = col, design = ~1) # ~1 means no design
dds <- estimateSizeFactors(dds)
norm_cts <- counts(dds, normalized = T)
# write normalized read-count to a new file without the lowly expressed genes
norm_expression <- norm_cts[expressed_genes, ]
write.table(norm_expression, file = "DESeq_norm_exp.txt", sep = "\t", quote = F)

# quantile normalization (mean: 0, standard deviation: 1, individual gene across all samples)
exp_mat <- as.matrix(norm_clean)
exp_mat <- t(apply(exp_mat, 1, rank, ties.method = "average"))
htseq_quantile <- qnorm(exp_mat/(ncol(exp_mat)+1))
# write quantile normalized read-count to a new file 
write.table(htseq_quantile, file = "quantile_norm_exp.txt", sep = "\t", quote = F)

