suppressMessages(library("argparse"))

parse_arg <- function() {
    # accept arguments from terminal
    parser <- ArgumentParser(description = "Usage: Rscript clean_count.R -raw <raw.txt> -bad <bad.txt> -O <prefix> ")
    parser$add_argument("-raw", "--raw_count", help = "raw allele-specific read count table ")
    parser$add_argument("-bad", "--bad", help = "bad SNPs which should be filtered out. ")
    parser$add_argument("-O", "--output", default = "default", help = "the prefix for the output name. If not provided, then the file will be named after the input prefix with .new extension. ")
    parser$parse_args()
}

args <- parse_arg()
raw <- args$raw_count
bad <- args$bad
out <- args$output
###
raw.df <- read.table(raw, sep = "\t", header = T)
bad.df <- read.table(bad, sep = "\t", header = T, row.names = 1)

rownames(raw.df) <- with(raw.df, paste(chromosome, position, sep = ":"))
remove <- rownames(bad.df)
gen.new <- raw.df[!rownames(raw.df) %in% remove, ]

if (out == "default") {
    raw <- basename(raw)
    pre <- strsplit(raw, ".txt")[[1]]
} else {
    pre <- out
}

write.table(gen.new, file = paste0(pre, ".txt"), sep = "\t", row.names = F, quote = F)