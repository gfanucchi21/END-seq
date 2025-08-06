library(tidyverse)

# Set working directory
setwd("C:/Users/user/WORK/ENDseq")

# Create the spike-in peak count matrix-------------------------------

# Import data relative to the reads in spike in peak
spike_in <- read.delim("multicov/spike_in_peak.multicov")

# Remove Chr, Start, End columns
spike_in <- select(spike_in, -c("chr", "start", "end"))

# Create the spike-in counts matrix 
spike_in <- spike_in |> pivot_longer(
  cols = everything(), 
  cols_vary = "slowest", 
  names_to = "sample", 
  values_to = "reads_spikein"
)

# Add columns relative to condition, replicate and strand
spike_in$Condition <- gsub("(CPT10|CPT20|NT)_.*", "\\1", spike_in$sample)
spike_in$Replicate <- gsub(".*_(rep1|rep2)_.*", "\\1", spike_in$sample)
spike_in$Strand <- gsub(".*_(fwd|rev)", "\\1", spike_in$sample)
colnames(spike_in) <- c("Sample", "Reads", "Condition", "Replicate", "Strand")

write.csv(spike_in,
          file = "multicov/spike_in_counts_temp.csv",
          row.names = FALSE
)

# Create the mouse library size----------------------------------------------------------
multiqc_samtools_flagstat <- read.delim("bam/bam_stats/multiqc_samtools_flagstat.txt")
multiqc_samtools_flagstat <- multiqc_samtools_flagstat |> select(Sample, mapped_passed)
multiqc_samtools_flagstat$genome <- gsub(".*_(hg|mm)_.*", "\\1", multiqc_samtools_flagstat$Sample)
library_size_mm <- multiqc_samtools_flagstat |> filter (genome == "mm") |> select(Sample, mapped_passed)
library_size_mm$Condition <- gsub("(CPT10|CPT20|NT)_.*", "\\1", library_size_mm10_H1$Sample)
library_size_mm$Replicate <- gsub(".*_(rep1|rep2)_.*", "\\1", library_size_mm10_H1$Sample)
library_size_mm$Strand <- gsub(".*_(fwd|rev)", "\\1", library_size_mm10_H1$Sample)
library_size_mm1$Sample <- gsub("mm_", "", library_size_mm10_H1$Sample)

mkdir("normalization_matrices")

write.csv(library_size_mm,
          file = "normalization_matrices/library_size_mm.csv",
          row.names = FALSE
)

# Remove genome column from the mouse library size matrix
library_size_mm <- library_size_mm |> select(!Genome)

# Normalization of the reads count in the spike-in peak over the total number of reads mapped on the mouse genome---------------------
spike_in <- merge(spike_in, library_size_mm) |> select("Sample", "Condition", "Strand", "Replicate", "reads_spikein", "mapped_passed")
spike_in <- spike_in |> mutate(reads_norm = reads_spikein / mapped_passed)

# Normalize every value using the smallest for each strand.
spike_in_norm_fwd <- spike_in |> 
  filter(Strand == "fwd") |>
  mutate(
    reads_norm_norm = reads_norm / min(reads_norm)
  )

spike_in_norm_rev <- spike_in |> 
  filter(Strand == "rev") |>
  mutate(
    reads_norm_norm = reads_norm / min(reads_norm)
  )

spike_in_norm_fwd_rev <- rbind(spike_in_norm_fwd,
                               spike_in_norm_rev
)

write.csv(spike_in_norm_fwd_rev,
          file = "normalization_matrices/spike_in_matrix_norm_fwd_rev.csv",
          row.names = FALSE
)

# Generating human library scaling matrices ---------------------------------------

# Import data relative to the library size
multiqc_samtools_stats_human <- read.delim("bam/bam_stats/multiqc_samtools_stats.txt")

# Select only the data relative to the number of reads
multiqc_samtools_stats_human <- multiqc_samtools_stats_human |> select(Sample, reads_mapped)
multiqc_samtools_stats_human$Condition <- gsub("(CPT10|CPT20|NT)_.*", "\\1", multiqc_samtools_stats_human$Sample)
multiqc_samtools_stats_human$Replicate <- gsub(".*_(rep1|rep2)_.*", "\\1", multiqc_samtools_stats_human$Sample)
multiqc_samtools_stats_human$Strand <- gsub(".*_(fwd|rev)", "\\1", multiqc_samtools_stats_human$Sample)
multiqc_samtools_stats_human$Genome <- gsub(".*_(hg|mm)_.*", "\\1", multiqc_samtools_stats_human$Sample)
multiqc_samtools_stats_human <- multiqc_samtools_stats_human |> filter(Genome == "hg")

write.csv(multiqc_samtools_stats_human,
          file = "library_size_hg.csv",
          row.names = FALSE
)

library_size_hg <- read.csv("normalization_matrices/library_size_hg.csv")

# Normalize every value using the smallest between each strand.
lib_size_matrix_fwd <- library_size_hg |>
  filter(Strand == "fwd") |>
  mutate(
    lib_size_scaling_factor =  reads_mapped / min(reads_mapped)
  )

lib_size_matrix_rev <- library_size_hg |> 
  filter(Strand == "rev") |>
  mutate(
    lib_size_scaling_factor =  reads_mapped / min(reads_mapped)
  )

# Merge the two tables
lib_size_matrix_fwd_rev <- full_join(lib_size_matrix_fwd, lib_size_matrix_rev)

# Create a scaling matrix that take into account only the number of reads mapped on the human genome
lib_size_scaling_matrix <- data.frame(
  condition = lib_size_matrix_fwd_rev$Condition,
  strand = lib_size_matrix_fwd_rev$Strand,
  replicate = lib_size_matrix_fwd_rev$Replicate,
  scaling_factor = lib_size_matrix_fwd_rev$lib_size_scaling_factor
)

lib_size_scaling_matrix <- lib_size_scaling_matrix |> 
  mutate(
    scaling_factor = 1 / scaling_factor
  )

write.csv(lib_size_scaling_matrix,
          file = "normalization_matrices/lib_size_scaling_matrix.csv",
          row.names = FALSE
)


# Generating final scaling matrix -----------------------------------------

final_scaling_matrix <- data.frame(
  condition = spike_in_norm_fwd_rev$Condition,
  strand = spike_in_norm_fwd_rev$Strand,
  replicate = spike_in_norm_fwd_rev$Replicate,
  scaling_factor = spike_in_norm_fwd_rev$reads_norm_norm * lib_size_matrix_fwd_rev$lib_size_scaling_factor
)

# The end seq signal have to be multiplied by the scaling factor
final_scaling_matrix <- final_scaling_matrix |> 
  mutate(
    scaling_factor = 1 / scaling_factor)


# Create final scaling matrix
write.csv(final_scaling_matrix,
          file = "normalization_matrices/final_scaling_matrix_lib_size_spikein.csv",
          row.names = FALSE,
          col.names = FALSE
)

