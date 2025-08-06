#### This script takes the multicoverage file containing the reads in every peak detected by MACS on both replicates, select only
#### those that are up-regulated compared to the control. Then of these select those that can be part of a deDSB and those that
#### represent seDSB. It also produces the BED files for each category of peak and for each sample and a .csv file for the peak
#### number of peaks.

# Set working directory
setwd("C:/Users/user/WORK/ENDseq")

## Import multicoverage data keeping duplicates
CPT10_fwd <- read.delim("multicov/CPT10_fwd.multicov")
CPT10_rev <- read.delim("multicov/CPT10_rev.multicov")
CPT20_fwd <- read.delim("multicov/CPT20_fwd.multicov")
CPT20_rev <- read.delim("multicov/CPT20_rev.multicov")

multicov_list_keep_dup<-list(CPT10_fwd, CPT10_rev, CPT20_fwd, CPT20_rev)

# Set table names
names(multicov_list)<-c("CPT10_fwd", "CPT10_rev", "CPT20_fwd", "CPT20_rev") 

# Select only peaks localizations, peak ids and read counts columns
library(tidyverse)
multicov_list<-lapply(multicov_list, function(x){
  x <- x |> select(1:4, 11:22)
  return(x)
})

# Remove intermediate tables
rm(CPT10_fwd, CPT10_rev, CPT20_fwd, CPT20_rev)
gc()

multicov_list <- lapply(multicov_list, function(x){
  colnames(x) <- c("chr", "start", "end", "name", 
                   "CPT10_a_fwd",	"CPT10_a_rev",	"CPT10_b_fwd",	"CPT10_b_rev",	"CPT20_a_fwd",	
                   "CPT20_a_rev",	"CPT20_b_fwd",	"CPT20_b_rev",	"NT_a_fwd",	"NT_a_rev",	
                   "NT_b_fwd",	"NT_b_rev")
  return(x)
})

# Check table headings
head(multicov_list[["CPT10_fwd"]])

# Set column names
multicov_list[["CPT10_fwd"]] <- multicov_list[["CPT10_fwd"]] |> select("chr", "start", "end", "name", "CPT10_a_fwd","CPT10_b_fwd", 
                                                                                         "NT_a_fwd","NT_b_fwd")

multicov_list[["CPT10_rev"]] <- multicov_list[["CPT10_rev"]] |> select("chr", "start", "end", "name", "CPT10_a_rev","CPT10_b_rev", 
                                                                                         "NT_a_rev","NT_b_rev")

multicov_list[["CPT20_fwd"]] <- multicov_list[["CPT20_fwd"]] |> select("chr", "start", "end", "name", "CPT20_a_fwd","CPT20_b_fwd", 
                                                                                         "NT_a_fwd","NT_b_fwd")

multicov_list[["CPT20_rev"]] <- multicov_list[["CPT20_rev"]] |> select("chr", "start", "end", "name", "CPT20_a_rev","CPT20_b_rev", 
                                                                                         "NT_a_rev","NT_b_rev")

# Pivot of each table in the list
multicov_list_long <- lapply(multicov_list, function(x){
  return(pivot_longer(
    x, 
    !c("chr", "start", "end", "name"),
    cols_vary = "fastest",
    names_to = "sample",
    values_to = "reads"
  ))
})

# Check table structure
head(multicov_list_long[["CPT10_fwd"]])

# Remove intermediate tables
rm(multicov_list)
gc()

# Import scaling matrix.
scaling_matrix <- read.csv2("normalization_matrices/final_scaling_matrix_lib_size_spikein.csv")
scaling_matrix <- scaling_matrix |> mutate(sample = paste0(condition, "_", replicate, "_", strand))

# Add the scalign factors to the counts table
multicov_list_long_sf <- lapply(multicov_list_long, function(x){
  return(
    merge(x, scaling_matrix, by = "sample")
  )})

# Remove duplicated peaks
multicov_list_long_sf <- lapply(multicov_list_long_sf, function(x){
  x <- x |> distinct(sample, chr, start, end, name, .keep_all = T)
  return(x)
})

# Normalize with the scaling factor
multicov_list_long_sf <- lapply(multicov_list_long_sf, function(x){
  x$reads_norm <- x$reads * x$scaling_factor
  return(x)
})

# Create a peak ID column merging Chr, Start, End
multicov_list_long_sf <- lapply(multicov_list_long_sf, function(x){
  x <- x |> mutate(peakID = paste0(x$name, "_", x$chr, "_", x$start, "_", x$end))
  return(x)
})

# Pivot wider to perform differential analysis
multicov_list_wide <- lapply(multicov_list_long_sf, function(x){
  x <- x %>% 
    select(c(peakID, sample, reads_norm)) %>%
    pivot_wider(
      id_cols = peakID,
      names_from = sample,
      values_from = reads_norm
    )
  return(x)
})

# Check tables structure
head(multicov_list_keep_dup_wide[["CPT10_rev"]])

# Remove intermediate tables
rm(multicov_list_keep_dup_long, multicov_list_keep_dup_long_sf)
gc()

# Set the peak_ID column as row names
multicov_list_wide <- lapply(multicov_list_wide, function(x){
  x <- column_to_rownames(x, var = "peakID")
  return(x)
}) 

# Remove rows containing NAs
multicov_list_wide <- lapply(multicov_list_wide, function(x){
  x <- na.omit(x)
  return(x)
})

# Keep only those peaks that have more that one read in the treated and at least one read in the control samples
multicov_list_wide_filt <- lapply(multicov_list_wide, function(x){
  x <- x |> filter(x[,1] > 1 & x[,2] > 1, x[,3] > 0 & x[,4] > 0)
  return(x)
})

# Remove intermediate tables
rm(multicov_list_keep_dup_wide)


# Differential expression analysis ----------------------------------------------------------------------------------------------
# 
library(limma)
library(edgeR)

DGE_list <- lapply(multicov_list_wide_filt, function(x){
  d0 <- DGEList(x, lib.size = rep(1,4), norm.factors = rep(1,4))
  return(d0)
})

# Create a group which contains experiment information from the sample names. in this case, we should have two factor levels:
# CPT and CTRL and a vector with four elements, since we have four replicates. 
group <- interaction(c("CPT", "CPT", "CTRL", "CTRL"))
# Multidimensional scaling (MDS) plot
plotMDS(DGE_list[["CPT10_fwd"]], col = as.numeric(group))
plotMDS(DGE_list[["CPT10_rev"]], col = as.numeric(group))
plotMDS(DGE_list[["CPT20_fwd"]], col = as.numeric(group))
plotMDS(DGE_list[["CPT20_rev"]], col = as.numeric(group))


# Specify a model to be fitted. This one specifies a model where each coefficient corresponds to a group mean. 
mm <- model.matrix(~0 + group)
# Voom transformation and calculation of variance weights
voom_list <- lapply(DGE_list, function(x){
  y <- voom(x, mm)
  return(y)
})

# lmFit fits a linear model using weighted least squares for each gene
fit_list <- lapply(voom_list, function(x){
  fit <- lmFit(x, mm)
  return(fit)
})

# Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models.
# Specify which groups to compare:
contr_list <- lapply(fit_list, function(x){
  contr <- makeContrasts(groupCPT - groupCTRL, levels = colnames(coef(x)))
  return(contr)
})

# Estimate contrast for each gene
contrasts_fit <- function(x, y){
  tmp <- list(
    contrasts.fit(x, y)
  )
  return(tmp)
}

contrast_fit_list <- mapply(contrasts_fit, fit_list, contr_list)

# Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or smaller than those 
# from other genes towards the average standard error)
contrast_fit_list <- lapply(contrast_fit_list, function(x){
  return(eBayes(x))
})

# What peaks are most differentially expressed?
top_table_list <- lapply(contrast_fit_list, function(x){
  top_table <- topTable(x, sort.by = "P", n = Inf)
  return(top_table)
})

# Remove intermediate tables
rm(contrast_fit_list, DGE_list, fit_list, voom_list)
gc()

# Calculating log10 of the p value and adj p value
top_table_volcano_list <- lapply(top_table_list, function(x){
  x <- x %>% mutate(
    `-log_p_value`= -log10(P.Value), 
    `-logdj_p_value`= -log10(adj.P.Val)
  )
  return(x)
})

# Add a tag for the up-regulation and down-regulation
top_table_volcano_list <- lapply(top_table_volcano_list, function(x){
  x$diffexpressed_sig[x$logFC > 1 & x$P.Value < 0.05] <- "UP_sig"
  x$diffexpressed_sig[x$logFC < -1 & x$P.Value < 0.05] <- "DOWN_sig"
  return(x)
})

# Print the number of up-regulated and down-regulated peaks for every condition and strand
nrow(top_table_volcano_list[["CPT10_fwd"]] |> filter(diffexpressed_sig == "UP_sig"))
nrow(top_table_volcano_list[["CPT10_rev"]] |> filter(diffexpressed_sig == "UP_sig"))
nrow(top_table_volcano_list[["CPT20_fwd"]] |> filter(diffexpressed_sig == "UP_sig"))
nrow(top_table_volcano_list[["CPT20_rev"]] |> filter(diffexpressed_sig == "UP_sig"))

nrow(top_table_volcano_list[["CPT10_fwd"]] |> filter(diffexpressed_sig == "DOWN_sig"))
nrow(top_table_volcano_list[["CPT10_rev"]] |> filter(diffexpressed_sig == "DOWN_sig"))
nrow(top_table_volcano_list[["CPT20_fwd"]] |> filter(diffexpressed_sig == "DOWN_sig"))
nrow(top_table_volcano_list[["CPT20_rev"]] |> filter(diffexpressed_sig == "DOWN_sig"))

## Volcano plots ----------------------------------------------------------------------
# Create a directory for volcano plots
mkdir("multicov/volcanoes")

library(ggrepel)
library(EnhancedVolcano)
EnhancedVolcano(top_table_volcano_list[["CPT10_fwd"]],
                title = "10' CPT fwd vs NT",
                lab = NA,
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-4,4),
                ylim = c(0, 4),
                pointSize = 0.5,
                subtitle = "",
                axisLabSize = 19,
                legendPosition = "bottom",
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", expression(p - value ~ and
                                                                                      ~ log[2] ~ FC)),)
ggsave("CPT10_fwd_volcano.pdf", 
       plot = last_plot(),
       path = "multicov/volcanoes",
       width = 6.5,
       height = 6)

EnhancedVolcano(top_table_volcano_list[["CPT10_rev"]],
                title = "10' CPT rev vs NT",
                lab = NA,
                x = 'logFC',
                y = 'P.Value',
                xlim = c(-4,4),
                ylim = c(0, 4),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 0.5,
                subtitle = "",
                axisLabSize = 19,
                legendPosition = "bottom",
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", expression(p - value ~ and
                                                                                      ~ log[2] ~ FC))
)
ggsave("CPT10_rev_volcano.pdf", 
       plot = last_plot(),
       path = "multicov/volcanoes",
       width = 6.5,
       height = 6)

EnhancedVolcano(top_table_volcano_list[["CPT20_fwd"]],
                title = "20' CPT fwd vs NT",
                lab = NA,
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-4,4),
                ylim = c(0, 4),
                pointSize = 0.5,
                subtitle = "",
                axisLabSize = 19,
                legendPosition = "bottom",
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", expression(p - value ~ and
                                                                                      ~ log[2] ~ FC))
)
ggsave("CPT20_fwd_volcano.pdf", 
       plot = last_plot(),
       path = "multicov/volcanoes",
       width = 6.5,
       height = 6)

EnhancedVolcano(top_table_volcano_list[["CPT20_rev"]],
                title = "20' CPT rev vs NT",
                lab = NA,
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-4,4),
                ylim = c(0, 4),
                pointSize = 0.5,
                subtitle = "",
                axisLabSize = 19,
                legendPosition = "bottom",
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", expression(p - value ~ and
                                                                                      ~ log[2] ~ FC))
)
ggsave("CPT20_rev_volcano.pdf", 
       plot = last_plot(),
       path = "multicov/volcanoes",
       width = 6.5,
       height = 6)


### Produce the BED files for up-regulated and down-regulated peaks ------------------------------------------

# Prepare the tables
up_sig_list <- lapply(top_table_volcano_list, function(x){
  x <- x |> filter(diffexpressed_sig == "UP_sig")
  y <- tibble(up_sig_peaks = row.names(x))
  y$chr <- gsub(".*_(chr.|chr..)_.*", "\\1", y$up_sig_peaks)
  y$start <- as.numeric(str_extract(y$up_sig_peaks, "(\\d+)(?=_\\d+$)"))
  y$end <- as.numeric((str_extract(y$up_sig_peaks, "(\\d+$)")))
  return(y)
})

down_sig_list <- lapply(top_table_volcano_list, function(x){
  x <- x |> filter(diffexpressed_sig == "DOWN_sig")
  y <- tibble(up_sig_peaks = row.names(x))
  y$chr <- gsub(".*_(chr.|chr..)_.*", "\\1", y$up_sig_peaks)
  y$start <- as.numeric(str_extract(y$up_sig_peaks, "(\\d+)(?=_\\d+$)"))
  y$end <- as.numeric((str_extract(y$up_sig_peaks, "(\\d+$)")))
  return(y)
})

# Select only chromosome start and end
BED_upsig_list <- lapply(up_sig_list, function(x){
  x <- x |> select(chr, start, end)
  return(x)
})

BED_downsig_list <- lapply(down_sig_list, function(x){
  x <- x |> select(chr, start, end)
  return(x)
})

names <- c("CPT10_fwd", "CPT10_rev", "CPT20_fwd",  "CPT20_rev")

# Generate .bed files
lapply(1:length(BED_list), function(i) 
  write.table(BED_list[[i]],
              file = paste0("bed/", names[i], "_upsig.bed"),
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)
)

lapply(1:length(BED_downsig_list), function(i) 
  write.table(BED_downsig_list[[i]],
              file = paste0("bed/", names[i], "_downsig.bed"),
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)
)