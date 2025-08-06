# Set working directory
setwd("C:/Users/user/WORK/ENDseq")

# Import the data obtained by the bedtools closest tool
CPT10_H1 <- read.delim("closest/CPT10.closest", header=FALSE)
CPT20_H1 <- read.delim("closest/CPT20.closest", header=FALSE)

closest_list <- list(CPT10_H1, CPT20_H1)
names(closest_list) <- c("CPT10_closest_peaks", "CPT20_closest_peaks")

# Remove peaks in which the distance is not valid
closest_list <- lapply(closest_list, function(x){
  x <- x |> filter(V4 != ".")
  return(x)
})


# deDSB detection ----------------------------------------------------------

# Start of the reverse must be before the start of the forward and the distance must be lower that 150bp
deDSB_list <- lapply(closest_list, function(x){
  x <- x |> filter(V7 > -150 & V7 <= 0) |> 
    filter(V5 < V2)
  return(x)
})

# Create the table with the peak reconstruction relative to the peak in which both strands are upregulated
entire_deDSB_both_strands_upsig_list <- lapply(deDSB_list, function(x){
  y <- tibble(
    chr = x$V1,
    start = x$V5,
    end = x$V3,
    peak_id = paste0(x$V1, "_", x$V5, "_", x$V3),
    zeros = rep(0, nrow(x)),
    plus = rep(".", nrow(x))
  )
  return(y)
}) 

# Produce the BED files containing the entire peak relative to the DSB in which both strands are upregulated
mkdir("bed/deDSBs")

names <- c("CPT10", "CPT20")

lapply(1:length(entire_deDSB_both_strands_upsig_list), function(i) 
  write.table(entire_deDSB_both_strands_upsig_list[[i]],
              file = paste0("bed/deDSBs/", names[i], "_deDSB_both_upsig.bed"),
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)
)


# seDSB detection ----------------------------------------------------------

# Make a directory for seDSBs files
mkdir("bed/seDSBs")

# Detect fwd seDSBs
closest_list_seDSB <- lapply(closest_list, function(x){
  x <- x |> filter(V7 > 150 | V7 < -150)
  return(x)
}) 

seDSB_CPT10_fwd <- tibble(chr = closest_list_seDSB[["CPT10_fwd_closest_peaks"]]$V1,
                          start = closest_list_seDSB[["CPT10_fwd_closest_peaks"]]$V2,
                          end = closest_list_seDSB_fwd[["CPT10_fwd_closest_peaks"]]$V3,
                          peak_id = paste0(closest_list_seDSB[["CPT10_fwd_closest_peaks"]]$V1,
                                           "_", 
                                           closest_list_seDSB[["CPT10_fwd_closest_peaks"]]$V2,
                                           "_",
                                           closest_list_seDSB[["CPT10_fwd_closest_peaks"]]$V3),
                          zeros = rep(0, nrow(closest_list_seDSB[["CPT10_fwd_closest_peaks"]])),
                          plus = rep("+", nrow(closest_list_seDSB[["CPT10_fwd_closest_peaks"]]))
)|>
  distinct(chr, start, end, .keep_all = T)

seDSB_CPT20_fwd <- tibble(chr = closest_list_seDSB[["CPT20_fwd_closest_peaks"]]$V1,
                          start = closest_list_seDSB[["CPT20_fwd_closest_peaks"]]$V2,
                          end = closest_list_seDSB[["CPT20_fwd_closest_peaks"]]$V3,
                          peak_id = paste0(closest_list_seDSB[["CPT20_fwd_closest_peaks"]]$V1,
                                           "_", 
                                           closest_list_seDSB[["CPT20_fwd_closest_peaks"]]$V2,
                                           "_",
                                           closest_list_seDSB[["CPT20_fwd_closest_peaks"]]$V3),
                          zeros = rep(0, nrow(closest_list_seDSB[["CPT20_fwd_closest_peaks"]])),
                          plus = rep("+", nrow(closest_list_seDSB[["CPT20_fwd_closest_peaks"]]))
) |>
  distinct(chr, start, end, .keep_all = T)

# Detect rev seDSBs
seDSB_CPT10_rev <- tibble(chr = closest_list_seDSB[["CPT10_rev_closest_peaks"]]$V4,
                          start = closest_list_seDSB[["CPT10_rev_closest_peaks"]]$V5,
                          end = closest_list_seDSB[["CPT10_rev_closest_peaks"]]$V6,
                          peak_id = paste0(closest_list_seDSB[["CPT10_rev_closest_peaks"]]$V4,
                                           "_", 
                                           closest_list_seDSB[["CPT10_rev_closest_peaks"]]$V5,
                                           "_",
                                           closest_list_seDSB[["CPT10_rev_closest_peaks"]]$V6),
                          zeros = rep(0, nrow(closest_list_seDSB[["CPT10_rev_closest_peaks"]])),
                          plus = rep("-", nrow(closest_list_seDSBv[["CPT10_rev_closest_peaks"]]))
) |>
  distinct(chr, start, end, .keep_all = T)

seDSB_CPT20_rev <- tibble(chr = closest_list_seDSB[["CPT20_rev_closest_peaks"]]$V4,
                          start = closest_list_seDSB[["CPT20_rev_closest_peaks"]]$V5,
                          end = closest_list_seDSB[["CPT20_rev_closest_peaks"]]$V6,
                          peak_id = paste0(closest_list_seDSB[["CPT20_rev_closest_peaks"]]$V4,
                                           "_", 
                                           closest_list_seDSB[["CPT20_rev_closest_peaks"]]$V5,
                                           "_",
                                           closest_list_seDSB[["CPT20_rev_closest_peaks"]]$V6),
                          zeros = rep(0, nrow(closest_list_seDSB[["CPT20_rev_closest_peaks"]])),
                          plus = rep("-", nrow(closest_list_seDSB[["CPT20_rev_closest_peaks"]]))
) |>
  distinct(chr, start, end, .keep_all = T)

# Create a list with all the seDSBs
seDSB_list <- list(seDSB_CPT10_fwd, seDSB_CPT10_rev, seDSB_CPT20_fwd, seDSB_CPT20_rev)
names(seDSB_list) <- c("seDSB_CPT10_fwd", "seDSB_CPT10_rev", "seDSB_CPT20_fwd", "seDSB_CPT20_rev")

# Produce .bed files for the seDSBS
names <- c("CPT10_fwd", "CPT10_rev", "CPT20_fwd", "CPT20_rev")
lapply(1:length(seDSB_list), function(i) 
  write.table(seDSB_list[[i]],
              file = paste0("bed/seDSBs/", names[i], "_seDSB_upsig.bed"),
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)
)

