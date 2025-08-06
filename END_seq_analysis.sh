#!/bin/bash
conda activate bioinfo
# Folder containing .fastq.gz files
cd /mnt/c/Users/user/WORK/ENDseq
# Set the number of threads to use
THR=5

# Generate ids file
# It is important that the control filename start with NT and replicates are listed as rep1, rep2 ecc
parallel -j 1 echo {1}{2} ::: CPT10_ CPT20_ NT_ ::: rep1 rep2 > ids

# Read all sample names into an array
samples=($(cat ids))

# Extract all unique condition names
conditions=($(for s in "${samples[@]}"; do echo "$s" | cut -d'_' -f1; done | sort -u))

# Extract all unique replicate names
replicates=($(for s in "${samples[@]}"; do echo "$s" | cut -d'_' -f2; done | sort -u))

# Define the strands
strands=("fwd" "rev")

### Check files
md5sum --check md5.txt

# Perform fastqc + multiqc
parallel fastqc {} ::: *.fastq.gz
multiqc . -n multiqc_fastqc_ENDseq

## Trimming
mkdir trim_galore
cat ids | parallel "trim_galore {}_R1.fastq.gz --illumina --fastqc -o trim_galore/"
cat ids | parallel "trim_galore {}_R2.fastq.gz --illumina --fastqc -o trim_galore/"

## Alignment
# Define the bwa genome index.
BWA_index="/mnt/c/Users/user/WORK/ENDseq/hybrid_hg19_mm10.fa"
mkdir bam
# Run Alignment
parallel -j2 "bwa mem -t 4 $BWA_index {1} {2} | samtools sort -@ 4 -o bam/{1/_R}.bam" ::: trim_galore/*_R1.fastq.gz :::+ trim_galore/*_R2.fastq.gz

# Sort alignments depending on the genome on which they align
hg_chrs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrR"
mm_chrs="mm10_chr10 mm10_chr11 mm10_chr12 mm10_chr13 mm10_chr14 mm10_chr15 mm10_chr16 mm10_chr17 mm10_chr18 mm10_chr19 mm10_chr1 mm10_chr2 mm10_chr3 mm10_chr4 mm10_chr5 mm10_chr6 mm10_chr7 mm10_chr8 mm10_chr9 mm10_chrX mm10_chrY"
for i in $(cat ids)
    samtools view -b -@ $THR bam/"$i".bam $hg38_chrs -o bam/"$i"_hg.bam
    samtools view -b -@ $THR bam/"$i".bam $mm10_chrs -o bam/"$i"_mm.bam
done

# Generate a new ids file accounting for the genome
parallel -j 1 echo {1}{2}{3} ::: CPT10_ CPT20_ NT_ ::: rep1_ rep2_ ::: hg mm > ids_genome

# Pre-filtering statistics
mkdir bam/bam_stats
for i in $(cat ids_genome)
    picard MarkDuplicates -I bam/"$i".bam -O bam/"$i"_md.bam -M bam/bam_stats/"$i"_md.txt
    samtools stats -@ $THR bam/$i.bam > bam/bam_stats/$i.stats
    samtools flagstat -@ $THR bam/$i.bam > bam/bam_stats/$i.flagstat
    samtools idxstats -@ $THR bam/$i.bam > bam/bam_stats/$i.idxstat
done
multiqc bam/bam_stats

# Filter by genome reference and alignment quality
for i in $(cat ids_genome)
    samtools view -b -q 30 -f PROPER_PAIR -@ $THR bam/"$i"_md.bam -o bam/"$i"_md_filt.bam
done

# Post-filtering quality control
for i in $(cat ids_genome)
    samtools stats -@ $THR bam/"$i"_md_filt.bam > bam/bam_stats/"$i"_md_filt.stats
    samtools flagstat -@ $THR bam/"$i"_md_filt.bam > bam/bam_stats/"$i"_md_filt.flagstat
    samtools idxstats -@ $THR bam/"$i"_md_filt.bam > bam/bam_stats/"$i"_md_filt.idxstat
done
multiqc bam/bam_stats/*filt*

## Split bam in fwd e rev
for i in $(cat ids_genome)
do
    FLAG_FW1=$(echo fwd99) 
    FLAG_FW2=$(echo fwd147) 
    FLAG_REV1=$(echo rev83) 
    FLAG_REV2=$(echo rev163) 
    filename="$sample"_hg19
    samtools view -f 99 -@ $THR bam/"$i"_md_filt.bam -o bam/"$i"_"$FLAG_FW1".bam 
    samtools view -f 147 -@ $THR bam/"$i"_md_filt.bam -o bam/"$i"_"$FLAG_FW2".bam 
    samtools view -f 83 -@ $THR bam/"$i"_md_filt.bam -o bam/"$i"_"$FLAG_REV1".bam 
    samtools view -f 163 -@ $THR bam/"$i"_md_filt.bam -o bam/"$i"_"$FLAG_REV2".bam 
    samtools merge -f -@ $THR bam/"$i"_fwd.bam bam/"$i"_"$FLAG_FW1".bam bam/"$i"_"$FLAG_FW2".bam 
    samtools merge -f -@ $THR bam/"$i"_rev.bam bam/"$i"_"$FLAG_REV1".bam bam/"$i"_"$FLAG_REV2".bam 
    samtools index -b -M -@ $THR bam/"$i"_fwd.bam bam/"$i"_rev.bam
done

rm bam/*99.bam bam/*147.bam bam/*83.bam bam/*163.bam

## Generate a new ids file accounting for the strand
parallel -j 1 echo {1}{2}{3} ::: CPT10_ CPT20_ NT_ ::: rep1_ rep2_ ::: hg mm ::: fwd rev > ids_genome_strand


## To normalize the number of reads in the peaks for the immunoprecipitation efficiency we 
## evaluated the number of reads in the spike-in peak (chr6:41554859-41555099) within each 
## condition and replicate and divided this number for the total number of reads aligned on the murine 
## m10 genome within each condition and replicate. The resulting values are then divided by the 
## minimum resulting value within each replicate. We refer to these values as the spike-in scaling factors.  

## To normalize the library size, we evaluated the total number of quality-filtered reads for each 
## replicate and then divided by the minimum resulting value within each replicate. We refer to these 
## values as the library size scaling factors. 

## The final scaling factor is calculated by multiplying the spike-in scaling factor for the library size 
## scaling factor.

## To calculate the scaling factors you need to run bedtools multicov on the spike-in peak.
## A bed file specifying the coordinates of the spike-in peak is needed
mkdir multicov
spikein_bed="bed/spikein_peak.bed"
output="multicov/spike_in_peak.multicov"
bedtools multicov -D -bams $(cat ids_genome_strand | grep mm | awk '{print "bam/"$1".bam"}' | tr '\n' ' ') -bed $spikein_bed > $output

### Run compute_scaling_factors.r in Rstudio

# Import scaling vector. The order of the scaling factors must correspond to the order of samples in input in ids_strand
scaling_vector=($(awk -F',' '{print $4}' normalization_matrices/final_scaling_matrix_lib_size_spikein.csv))
echo "${scaling_vector[@]}"

## Generate a new ids file torefer only to the human genome bam files
parallel -j 1 echo {1}{2}{3} ::: CPT10_ CPT20_ NT_ ::: rep1_ rep2_ ::: hg ::: fwd rev > ids_hg_strand

## Generate bigwig coverage file normalized on the spike-in and library size
conda activate deeptools
mkdir bigwig
i=0
for name in $(cat ids_hg_strand)
do
    SCALE_FACTOR=${scaling_vector[i]}
    bamCoverage -b bam/"$name".bam --scaleFactor $SCALE_FACTOR -o BigWig/"$name"_spikein.bw --binSize 10 -p $THR
    i=$((i+1))
done

## Average BigWigs
# Loop through each condition and strand
for condition in "${conditions[@]}"; do
  for strand in "${strands[@]}"; do
    rep1="bigwig/${condition}_rep1_${strand}_spikein.bw"
    rep2="bigwig/${condition}_rep2_${strand}_spikein.bw"
    mean="bigwig/${condition}_mean_${strand}_spikein.bw"
    # Check if both files exist in the sample list
    if [[ " ${samples[@]} " =~ " $rep1 " && " ${samples[@]} " =~ " $rep2 " ]]; then
      echo "Running Bigwig compare for condition=$condition, strand=$strand with $rep1 and $rep2"
      # Run bigwig compare
      bigwigCompare --operation "mean" -p $THR -b1 $rep1 -b2 $rep2 -o $mean
    else
      echo "Skipping condition=$condition, strand=$strand — one or both replicates missing"
    fi
  done
done

# Create the directory for bed files
mkdir bed

## Perform strand-specific peak-calling
for condition in "${conditions[@]}"; do
    for replicate in "${replicates[@]}"; do
        for strand in "${strands[@]}"; do
            if [[" ${condition[@]} " != "NT"]]
            macs2 callpeak -t bam/${condition}_${replicate}_hg_${strand}.bam -c bam/NT_${replicate}_hg_${strand}.bam -g hs --nomodel --nolambda -q 0.05 --keep-dup all -n ${condition}_${replicate}_${strand} -f BAMPE --outdir bed
            fi
        done 
    done
done

## Remove all peaks that overlap blacklisted regions
for condition in "${conditions[@]}"; do
for replicate in "${replicates[@]}"; do
for strand in "${strands[@]}"; do
# Input filename
bed="bed/${condition}_${replicate}_${strand}.narrowPeak"+
# Specify path to blacklist
blacklist=""
# Output filename
output="bed/${condition}_${replicate}_${strand}_blacklist.bed"
# -v produces all peaks in A that do not overlap with B
bedtools intersect -v -a $bed -b $blacklist > $output

## Keep only those peaks that overlap in both replicates and mantain them as separete entities
rep1="bed/${condition}_rep1_${strand}_blacklist.bed"
rep2="bed/${condition}_rep2_${strand}_blacklist.bed"
bedtools intersect -wa -u -a $rep1 -b $rep2 > bed/${condition}_keep_${strand}_blacklist.bed
bedtools intersect -wa -u -a $rep2 -b $rep1 >> bed/${condition}_keep_${strand}_blacklist.bed


### Differential analysis between the number of reads in each peak treated vs ctrl -----------------------------

## Evaluation of the number of reads in each peak
for condition in "${conditions[@]}"; do
    for strand in "${strands[@]}"; do
        if [[" ${condition[@]} " != "NT"]]
        treated_rep1="bam/${condition}_rep1_${strand}.bam"
        treated_rep2="bam/${condition}_rep2_${strand}.bam"
        NT_rep1="bam/NT_rep1_${strand}.bam"
        NT_rep2="bam/NT_rep2_${strand}.bam"
        BED="bed/${condition}_keep_${strand}_blacklist.bed"
        output="multicov/${condition}_${strand}.multicov"
        # -D keep duplicates
        bedtools multicov -D -bams $treated_rep1 $treated_rep2 $NT_rep1 $NT_rep2 -bed $BED > $output
        fi
    done
done

### Run differential_analysis.r in Rstudio to detect upregulated peaks compared to the controls

# Run bedtools closest to find the closest peaks between fwd and rev strands
mkdir closest
for condition in ${conditions[@]};do 
    fwd="bed/"$condition"_fwd_upsig.bed"
    rev="bed/"$condition"_rev_upsig.bed"
    output="closest/"$condition".closest"
    bedtools closest -D ref -a $fwd -b $rev > $output
done

### Run closest_analysis.r in Rstudio to detect the deDSBs and seDBSs

# Deduplicate deDSBs
for condition in ${conditions[@]};do 
    input="bed/deDSBs/"$condition"_deDSB.bed"
    output="bed/deDSBs/"$condition"_deDSB_merged.bed"
    bedtools sort -i $input |  bedtools merge -i - > $output
    n_deDSBs=$(wc -l $output)
    echo "Number of deDSBs detected=$n_deDSBs"
done

# Deduplicate seDSBs
for condition in ${conditions[@]};do 
for strand in "${strands[@]}"; do
    input="bed/seDSBs/"$condition"_"$strand"_seDSB.bed"
    output="bed/seDSBs/"$condition"_"$strand"_seDSB_merged.bed"
    bedtools sort -i $input | 
    bedtools merge -i - |
    bedtools intersect -wa -v -a - -b $input > $output
    n_seDSBs=$(wc -l $output)
    echo "Number of seDSBs detected=$n_seDSBs"
done
done

# Find transient late persistent deDSBs and seDSBs -----------------------------------------------------------

### The experimental setting we used leveraged two time points, 10 minutes and 20 minutes after CPT treatment, 
### to distinguish transient, late and persistent deDSBs and seDSBs. Transient are revealed only at 10 minutes,
### late are revealed only at 20 minutes, while persistent are revealed at both time points.

## Find transient, late and persistent deDSBs
CPT10="bed/deDSBs/CPT10_deDSB_merged.bed"
CPT20="bed/deDSBs/CPT20_deDSB_merged.bed"
# -v Only report those entries in A that have no overlap in B. 
bedtools intersect -v -a $CPT10 -b $CPT20 > bed/deDSBs/transient_deDSB.bed

bedtools intersect -v -a $CPT20 -b $CPT10 > bed/deDSBs/late_deDSB.bed

bedtools intersect -wa -u -a $CPT10 -b $CPT20 > bed/deDSBs/temp.bed
bedtools intersect -wa -u -b $CPT10 -a $CPT20 >> bed/deDSBs/temp.bed
bedtools sort -i temp.bed | bedtools merge -i - > bed/deDSBs/persistent_deDSB.bed
rm bed/deDSBs/temp.bed

# Print the number of transient, late and persistent deDSBs
wc -l bed/deDSBs/transient_deDSB.bed \
    bed/deDSBs/late_deDSB.bed \
    bed/deDSBs/persistent_deDSB.bed

## Find transient, late and persistent seDSBs
for strand in ${strands[@]};do
CPT10=bed/seDSBs/CPT10_"$strand"_seDSB_merged.bed
CPT20=bed/seDSBs/CPT20_"$strand"_seDSB_merged.bed
# -v Only report those entries in A that have no overlap in B. 
bedtools intersect -v -a $CPT10 -b $CPT20 > bed/seDSBs/transient_seDSBs_"$strand".bed

bedtools intersect -v -a $CPT20 -b $CPT10 > bed/seDSBs/late_seDSBs_"$strand".bed

bedtools intersect -wa -u -a $CPT10 -b $CPT20 > bed/seDSBs/temp.bed
bedtools intersect -wa -u -b $CPT10 -a $CPT20 >> bed/seDSBs/temp.bed
bedtools sort -i temp.bed | bedtools merge -i - > bed/seDSBs/persistent_seDSBs_"$strand".bed
rm bed/seDSBs/temp.bed
done

# Print the number of transient, late and persistent seDSBs
for strand in ${strands[@]};do
wc -l bed/seDSBs/transient_seDSBs_"$strand".bed \
    bed/seDSBs/late_seDSBs_"$strand".bed \
    bed/seDSBs/persistent_seDSBs_"$strand".bed
done


# Validazione picchi deDSBs ------------------------------------------------
conda activate deeptools
mkdir matrices

for condition in ${conditions[@]};do
    bw_fwd="bigwig/"$condition"_mean_fwd_spikein.bw"
    bw_rev="bigwig/"$condition"_mean_rev_spikein.bw"

    bw_NT_fwd=bigwig/NT_mean_fwd_spikein.bw
    bw_NT_rev=bigwig/NT_mean_rev_spikein.bw

    bed="bed/deDSBs/"$condition"_deDSB_merged.bed"

    matrix="matrices/matrix_deDSB_"$condition".gz"

    computeMatrix reference-point \
    -S $bw_fwd $bw_rev $bw_NT_fwd $bw_NT_rev \
    -R $bed \
    -o $matrix \
    -a 3000 -b 3000 \
    --missingDataAsZero \
    --skipZeros \
    --referencePoint center \
    -p $THR

    heatmap="heatmap_deDSB_"$condition".pdf"

    plotHeatmap -m $matrix -out $heatmap \
            --perGroup \
            --plotType se \
            --colorList white,blue \
            --heatmapWidth 9 --heatmapHeight 9 \
            --regionsLabel ""$condition" deDSBs" \
            --samplesLabel ""$condition" fwd" ""$condition" rev" "NT fwd" "NT rev"
done

# Validazione picchi seDSBs
for condition in ${conditions[@]};do
for strand in ${strands[@]};do
    bw_fwd="bigwig/"$condition"_mean_fwd_spikein.bw"
    bw_rev="bigwig/"$condition"_mean_rev_spikein.bw"

    bw_NT_fwd=bigwig/NT_mean_fwd_spikein.bw
    bw_NT_rev=bigwig/NT_mean_rev_spikein.bw

    bed="bed/seDSBs/"$condition"_"$strand"_seDSB_merged.bed"

    matrix="matrices/matrix_seDSB_"$condition"_"$strand".gz"

    computeMatrix reference-point \
    -S $bw_fwd $bw_rev $bw_NT_fwd $bw_NT_rev \
    -R $bed \
    -o $output \
    -a 3000 -b 3000 \
    --missingDataAsZero \
    --skipZeros \
    --referencePoint center \
    -p $THR

    heatmap="heatmap_seDSB_"$condition"_"$strand".pdf"

    plotHeatmap -m $matrix -out $HEAT_MAP \
            --perGroup \
            --plotType se \
            --colorList white,blue \
            --heatmapWidth 9 --heatmapHeight 5 \
            --regionsLabel ""$condition" seDSBs" \
            --samplesLabel ""$condition"_fwd" ""$condition"_rev" "NT_fwd" "NT_rev"
done
done


# Validazione transient, late and persistent deDSBs

BW_CPT10_fwd="bigwig/CPT10_mean_fwd_spikein.bw"
BW_CPT10_rev="bigwig/CPT10_mean_rev_spikein.bw"

BW_CPT20_fwd="bigwig/CPT20_mean_fwd_spikein.bw"
BW_CPT20_rev="bigwig/CPT20_mean_rev_spikein.bw"

BW_NT_fwd="bigwig/NT_mean_fwd_spikein.bw"
BW_NT_rev="bigwig/NT_mean_rev_spikein.bw"

BED_transient="bed/deDSBs/transient_deDSB.bed"
BED_late="bed/deDSBs/late_deDSB.bed"
BED_persistent="bed/deDSBs/persistent_deDSB.bed"

matrix="matrices/matrix_trans_late_pers_deDSBs.gz"

computeMatrix reference-point \
    -S $BW_CPT10_fwd $BW_CPT10_rev $BW_CPT20_fwd $BW_CPT20_rev $BW_NT_fwd $BW_NT_rev \
    -R $BED_trans $BED_late $BED_pers \
    -o $matrix \
    -a 3000 -b 3000 \
    --missingDataAsZero \
    --skipZeros \
    --referencePoint center \
    -p $THR

heatmap="heatmap_trans_late_pers_deDSBs.pdf"

plotHeatmap -m $matrix -out $heatmap \
    --perGroup \
    --plotType se \
    --colorList white,blue \
    --heatmapWidth 10 --heatmapHeight 10 \
    --regionsLabel "transient deDSBs" "late deDSBs" "persistent deDSBs" \
    --samplesLabel "10'CPT fwd" "10'CPT rev" "20'CPT fwd" "20'CPT rev" "NT fwd" "NT rev"