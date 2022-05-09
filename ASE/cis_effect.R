
suppressMessages(library("argparse"))

parse_arg <- function() {
    # accept arguments from terminal
    parser <- ArgumentParser(description = "Usage: Rscript cis_effect.R -ASE <ASE.txt> -O <prefix> ")
    parser$add_argument("-ASE", "--ASE", help = "Allele-specific expression on gene-basis for all sample. ")
    parser$add_argument("-remove_homo", "--remove_homo", action = "store_true", help = "Remove genes within homozygous genotype region. ")
    parser$add_argument("-O", "--output", help = "the prefix for the output name.")
    parser$parse_args()
}

args <- parse_arg()
ASE <- args$ASE
output <- args$output

###
gase <- read.table(ASE, sep = "\t", header = T)
### 
gase <- gase[(gase$alt1_read != 0) & (gase$alt2_read != 0), ]
if (args$remove_homo == TRUE) {
    gase$alt1_perc <- gase$alt1_read/(gase$alt1_read + gase$alt2_read)*100
    gase <- gase[(gase$alt1_perc > 5) & (gase$alt1_perc < 95), ]
}

gase$log2alt1_vs_alt2 <- log2(gase$alt1_read/gase$alt2_read)

gid <- unique(gase$gene_id)

gname <- c() # gene id
parcont <- c() # 
pcont <- c()
CI1cont <- c()
CI2cont <- c()
mmcont <- c()
errcont <- c()
i <- 0
for (g in gid) {
    #print(g)
    sub_gase <- gase[gase$gene_id == g, ]
    if (nrow(sub_gase) > 80) { # at least 80 sample with heterozygous genotype on that gene region
        i <- i+1
        gname <- c(gname, g)
        res <- t.test(sub_gase$log2alt1_vs_alt2, mu = 0, alternative = "two.sided")
        par <- as.numeric(res$parameter)
        parcont <- c(parcont, par+1)
        p <- res$p.value
        pcont <- c(pcont, p)
        ci1 <- res$conf.int[1]
        ci2 <- res$conf.int[2]
        CI1cont <- c(CI1cont, ci1)
        CI2cont <- c(CI2cont, ci2)
        mm <- as.numeric(res$estimate)
        mmcont <- c(mmcont, mm)
        stderr <- res$stderr
        errcont <- c(errcont, stderr)
    }
}

df <- data.frame(gene_id = character(i), sample_size = integer(i), p_value = numeric(i), CIL = numeric(i), CIR = numeric(i), alt1_vs_alt2_mean = numeric(i), stderr = numeric(i))

df$gene_id <- gname
df$sample_size <- parcont
df$p_value <- pcont
df$CIL <- CI1cont
df$CIR <- CI2cont
df$alt1_vs_alt2_mean <- mmcont
df$stderr <- errcont

df$p_adj <- p.adjust(df$p_value, method = "BH")
df <- df[with(df, order(p_value, p_adj, decreasing = F)), ]
write.table(df, file = paste0(output, ".txt"), sep = "\t", quote = F, row.names = F)
