# END-seq Pipeline
END-seq is a next-generation sequencing technique that allows the quantitative mapping of endogenous DNA double-strand breaks (DSBs) at nucleotide resolution.
Process END-seq data with this pipeline.

## Project description
This project has been created to analyse the END-seq data.
The experimental design consisted of the treatment of H1 stem cells with Camptothecin, followed by the END-seq analysis at two time points: 10 minutes and 20 minutes. 
The mapping of DSBs allowed us to obtain insights into the role of RNA polymerase and DNA polymerase in the formation of DSBs.

## Methods
A spike-in internal control has been used in the experiment. Each sample sequenced contained a constant number of murine cells harbouring an induced DSB in a known genomic location. Given that the efficiency of the formation of the spike-in DSB is nearly 100%, the intensity of the END-seq signal over the spike-in peak can be used to normalize the signal of the experimental samples. 
Paired-end sequencing was performed.
As the sequenced samples contained both murine and human cells, the reference genome used to align the reads was a merge of the hg38 and mm10 genomes.

## Analysis steps
- Reads are trimmed and aligned on the hybrid hg38-mm10 reference genome
- Reads are filtered based on the alignment quality
- Statistical analysis is performed both before and after filtering
- Reads aligning on the forward and reverse strands are split
- The normalization factors are computed based on the END-seq signal over the spike-in peak and the total library sizes
- BigWig coverage files are generated
- Strand-specific peak-calling was performed using MACS2
- Only peaks overlapping in both replicates were retained
- Differential analysis is performed using the limma package and edgeR
- Double-end DSBs (deDSBs) and single-end DSBs (seDSBs) are detected based on the relative distances between forward and reverse END-seq peaks
- Transient, late, and persistent DSBs are detected based on their formation time
- DSB validation has been performed by studying the END-seq signal intensity over the detected DSBs using deeptools.

## Running instructions
- END_seq_analysis_main_script.sh contains all the bash commands to be executed through the Linux terminal or WSL. It also specifies when to execute the R scripts (compute_scaling_factors.r, differential_analysis.r, closest_analysis.R) on R Studio.
- compute_scaling_factors.r takes the number of reads aligning on the spike-in peak, the murine and human library sizes to compute a scaling factor for each sample and for each DNA strand.
- differential_analysis.r perform a differential analysis to detect which END-seq peaks detected by MACS2 result up-regulated compared to the controls.
- closest_analysis.r takes the output from bedtools closest to produce the bed files for the deDSBs and seDSBs. 
