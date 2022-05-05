
suppressMessages(library("argparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("plyr"))

parse_arg <- function() {
    parser <- ArgumentParser(description = "Usage: Rscript block2bin.R -genodir <block_directory> -chrLen <chrlen.txt> -SNP <SNP_pos.txt> ")
    parser$add_argument("-genodir", "--genodir", help = 'directory of genotype block files for each sample. ')
    parser$add_argument("-chrLen", "--chrLen", help = "chromosome length file. ")
    parser$add_argument("-SNP", "--SNP_loc", help = "SNP location file to get the representative SNP. ")
    parser$parse_args()
}

args <- parse_arg()
genodir <- args$genodir
chrlen <- args$chrLen
snp <- args$SNP_loc

fis <- sort(list.files(genodir, pattern = ".txt$"))
df_list <- list()

for (f in fis) {
    ###
    pre <- strsplit(f, ".txt")[[1]][1]
    f.df <- read.table(file.path(genodir, f), sep = '\t', header = T, stringsAsFactors = F)
    for (chr_ret in c('chromosome_1', 'chromosome_2', 'chromosome_3')) {
        f.sub <- f.df[f.df$chromosome == chr_ret, ]
        rownames(f.sub) <- seq(1, nrow(f.sub))
        ###
        geno_1 <- ''
        chr_vec <- c()
        start_vec <- c()
        end_vec <- c()
        ###
        for (i in seq(1, nrow(f.sub))) {
            if (geno_1 == "") {
                geno_1 <- f.sub[i, ]$genotype
                start_snp <- as.numeric(strsplit(f.sub[i, ]$transition_site, "-")[[1]][1])
                end_snp <- as.numeric(strsplit(f.sub[i, ]$transition_site, "-")[[1]][2])
            } else {
                geno_new <- f.sub[i, ]$genotype
                start_new <- as.numeric(strsplit(f.sub[i, ]$transition_site, "-")[[1]][1])
                end_new <- as.numeric(strsplit(f.sub[i, ]$transition_site, "-")[[1]][2])
                if (geno_1 != geno_new) {
                    chr_vec <- c(chr_vec, chr_ret)
                    start_vec <- c(start_vec, start_snp)
                    end_vec <- c(end_vec, end_snp)
                    geno_1 <- geno_new
                    start_snp <- start_new
                    end_snp <- end_new
                } else {
                    start_snp <- as.numeric(strsplit(f.sub[i, ]$transition_site, "-")[[1]][1])
                    end_snp <- as.numeric(strsplit(f.sub[i, ]$transition_site, "-")[[1]][2])
                }
            }
        }
        t <- length(chr_vec)
        geno_df <- data.frame(sample = character(t), chromosome = character(t), start_info = integer(t), end_info = integer(t))
        ###
        geno_df$sample <- c(rep(pre, t))
        geno_df$chromosome <- chr_vec
        geno_df$start_info <- start_vec
        geno_df$end_info <- end_vec
        df_list[[paste0(pre, chr_ret)]] <- geno_df
    }
}

###
df_chr <- rbind.fill(df_list)
df_chr['midpoint'] <- (df_chr['start_info'] + df_chr['end_info'])%/%2
print("Write out breakpoint for all samples. ")
write.table(df_chr, file = "chr_breakpoint.txt", sep = "\t", quote = F, row.names = F)

### chromosome length file
len_df <- read.table(chrlen, sep = "\t", header = T)

### SNP location
snp_df <- read.table(snp, sep = '\t', header = T, stringsAsFactors = F)

### assign bin and get the representative SNP
bin_start_vec <- c()
bin_end_vec <- c()
chr_vec <- c()

for (i in c("chromosome_1", "chromosome_2", "chromosome_3")) {
    print(i)
    bin_start <- 1
    chr_len <- len_df[len_df$chromosome == i, 'length']
    sub_df <- df_chr[df_chr$chromosome == i, ]
    sub_snp <- snp_df[snp_df$chromosome == i, ]
    sub_df_order <- sub_df[order(sub_df$midpoint), ]
    sub_snp_order <- sub_snp[order(sub_snp$position), ]
    min_mid <- min(sub_df_order$midpoint)
    max_mid <- max(sub_df_order$midpoint)
    rownames(sub_df_order) <- seq(1, nrow(sub_df_order))
    write.table(sub_df_order, file = paste0(i,"_break.txt"), sep = '\t', quote = F, row.names = F)
    ###
    for (n in seq(1, nrow(sub_df_order))) {
        bin_end <- sub_df_order[n, 'midpoint']
        if (bin_end == min_mid) {
            bin_snp <- sub_snp_order[sub_snp_order$position >= bin_start & sub_snp_order$position < bin_end, ]
            if (length(row.names(bin_snp)) > 0) {
                rownames(bin_snp) <- seq(1, nrow(bin_snp))
                bin_start_vec <- c(bin_start_vec, bin_start)
                bin_end_vec <- c(bin_end_vec, bin_end)
                chr_vec <- c(chr_vec, i)
                bin_start <- bin_end
            }
        } else if (bin_end == max_mid) {
            bin_snp <- sub_snp_order[sub_snp_order$positioin >= bin_start & sub_snp_order$position < bin_end, ]
            bin_snp_end <- sub_snp_order[sub_snp_order$position >= bin_end & sub_snp_order$position < chr_len, ]
            ##
            if (length(row.names(bin_snp)) > 0) {
                rownames(bin_snp) <- seq(1, nrow(bin_snp))
                bin_start_vec <- c(bin_start_vec, bin_start)
                bin_end_vec <- c(bin_end_vec, bin_end)
                chr_vec <- c(chr_vec, i)
                bin_start <- bin_end
            }
            if (length(row.names(bin_snp_end)) > 0) {
                rownames(bin_snp_end) <- seq(1, nrow(bin_snp_end))
                bin_start_vec <- c(bin_start_vec, bin_start)
                bin_end_vec <- c(bin_end_vec, chr_len)
                chr_vec <- c(chr_vec, i)
            }
        } else {
            bin_snp <- sub_snp_order[sub_snp_order$position >= bin_start & sub_snp_order$position < bin_end, ]
            ##
            if (length(row.names(bin_snp)) > 0) {
                rownames(bin_snp) <- seq(1, nrow(bin_snp))
                bin_start_vec <- c(bin_start_vec, bin_start)
                bin_end_vec <- c(bin_end_vec, bin_end)
                chr_vec <- c(chr_vec, i)
                bin_start <- bin_end
            }
        }
    }
}

t <- length(chr_vec)
bin_df <- data.frame(chromosome = character(t), bin_start = integer(t), bin_end = integer(t))
###
bin_df$chromosome <- chr_vec
bin_df$bin_start <- bin_start_vec
bin_df$bin_end <- bin_end_vec

print("Write out bins for chromosomes. ")
for (i in c("chromosome_1", "chromosome_2", "chromosome_3")) {
    chr_bin <- bin_df[bin_df$chromosome == i, ]
    print(paste(i, nrow(chr_bin)))
}

bin_df['mid_bin'] <- (bin_df$bin_start + bin_df$bin_end)%/%2

bin_df['sudo_SNP'] <- paste0(bin_df$chromosome,":", bin_df$mid_bin)
#bin_df <- bin_df %>% select(-mid_bin)

len_df <- read.table(chrlen, sep = "\t", header = T)
scaffold_df <- len_df[!len_df$chromosome %in% c("chromosome_1", "chromosome_2", "chromosome_3"), ]
scaffold_df['bin_start'] <- 1
scaffold_df['bin_end'] <- scaffold_df['length']
scaffold_df['mid_bin'] <- (scaffold_df['bin_start'] + scaffold_df['bin_end'])%/%2
scaffold_df['sudo_SNP'] <- paste(scaffold_df$chromosome, scaffold_df$mid_bin, sep = ":")
scaffold_df <- scaffold_df[, c("chromosome", "bin_start", "bin_end", "mid_bin", "sudo_SNP")]

bin_df <- rbind(bin_df, scaffold_df) 
write.table(bin_df, file = "chr_bin.txt", sep = "\t", quote = F, row.names = F)