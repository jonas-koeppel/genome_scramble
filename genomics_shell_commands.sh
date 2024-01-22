# Basecalling with Guppy for each sample
while read SAMPLE; do
    bsub -G team154-vwork -o ./pod5/logs/guppy_${SAMPLE}.log -gpu - -q "gpu-normal" -R'select[mem>20000] rusage[mem=20000]' -M 20000 -J ${SAMPLE}_gpy \
        /path/to/guppy_basecaller --input_path ./pod5/${SAMPLE}/ --save_path ./fastq/${SAMPLE} -x cuda:0 -c dna_r10.4.1_e8.2_400bps_sup.cfg --num_callers 12
done < ./pod5/samples.txt

# Mapping reads to a reference genome using minimap
while read SAMPLE; do
    bsub -G team154-vwork -q normal -n 8 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M 32000 -o ./bam/minimap_${SAMPLE}.log -J minimap2 \
        '~/path/to/minimap2 --MD -ax map-ont /path/to/genome.fa -t 8 ./fastq/'${SAMPLE}'.fastq.gz > ./bam/'${SAMPLE}'.sam'
done < ./fastq/samples.txt

# Sorting mapped reads
while read SAMPLE; do
    bsub -G team154-vwork -q normal -n 8 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M 32000 -o ./bam/minimap_${SAMPLE}.log -J minimap2 \
    ~/path/to/samtools sort -@ 16 -o ./bam/${SAMPLE}.bam ./bam/${SAMPLE}.sam
done < ./fastq/samples.txt

# Indexing sorted reads
while read SAMPLE; do
    bsub -G team154-vwork -q normal -n 8 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M 32000 -o ./bam/minimap_${SAMPLE}.log -J minimap2 \
    ~/path/to/samtools index -@ 16 ./bam/${SAMPLE}.bam
done < ./fastq/samples.txt

# Calculating coverages using bedtools
while read SAMPLE; do
    bsub -G team154-vwork -q normal -n 8 -R'select[mem>10000] rusage[mem=10000] span[hosts=1]' -M 10000 -o ./${SAMPLE}/depth/${SAMPLE}.log -J ${SAMPLE}_depth \
        "samtools bedcov -Q 10 /path/to/hg38_windows_10k.bed ./${SAMPLE}/bam/${SAMPLE}.bam > ./${SAMPLE}/depth/${SAMPLE}.bedgraph"
done < ./samples.txt

# Calculating coverages using mosdepth
while read SAMPLE; do
    bsub -G team154-vwork -q normal -n 8 -R'select[mem>10000] rusage[mem=10000] span[hosts=1]' -M 10000 -o ./${SAMPLE}/depth/mosdepth_${SAMPLE}.log -J ${SAMPLE}_mosdepth \
    bash -c "cd ./${SAMPLE}/depth &&
    /path/to/mosdepth -n --by 50000 cov_50k ../bam/${SAMPLE}.bam &&
    cd ../.."
done < ./samples.txt

# Calculating B allele frequencies with amber
while read SAMPLE; do
     bsub -G team154-vwork -o ./${SAMPLE}/amber/${SAMPLE}.log  -q normal -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -n4 -J ${SAMPLE}_amber \
          java -Xmx32G -cp /path/to/amber-3.5.jar com.hartwig.hmftools.amber.AmberApplication \
               -tumor_only \
               -min_het_af_percent 0.1 \
               -max_het_af_percent 0.9 \
               -tumor ${SAMPLE} \
               -tumor_bam ./${SAMPLE}/bam/${SAMPLE}.bam \
               -output_dir ./${SAMPLE}/amber \
               -threads 16 \
               -loci /path/to/Analysis_supplements/GermlineHetPon.38.vcf.gz
done < ./samples.txt

# SV detection using Sniffles
while read SAMPLE; do
    /path/to/sniffles -i ./${SAMPLE}/bam/${SAMPLE}.bam \
    --reference /path/to/genome.fa -t 8 --mosaic --output-rnames --tandem-repeats /path/to/Analysis_supplements/tandem_repeat_annotation.bed --minsupport 1 --minsvlen 25 \
    --qc-output-all --mapq 5 -v ./${SAMPLE}/sniffles/${SAMPLE}.vcf.gz
done < ./samples.txt

# Structural variant analysis using NanoMonsv
# Setting up the environment for NanoMonsv
export PATH=/path/to/nanomonsv:$PATH
conda activate nanomonsv

# Parse
while read SAMPLE; do
     bsub -G team154-vwork -o ./${SAMPLE}/nanomonsv/${SAMPLE}.log  -q normal -R'select[mem>20000] rusage[mem=20000] span[hosts=1]' -M20000 -n8 -J nano_${SAMPLE} \
          nanomonsv parse ./${SAMPLE}/bam/${SAMPLE}.bam ./${SAMPLE}/nanomonsv/${SAMPLE} 
done < ./samples.txt

# Get SVs for HAP1
while read SAMPLE; do
    bsub -G team154-vwork -o ./nanomonsv/${SAMPLE}.log  -q week -R'select[mem>20000] rusage[mem=20000] span[hosts=1]' -M20000 -n8 -J nanomonsv_${SAMPLE} \
        nanomonsv get ./${SAMPLE}/nanomonsv/${SAMPLE} ./${SAMPLE}/bam/${SAMPLE}.bam /path/to/Analysis_supplements/genome.fa --control_prefix ./C516/nanomonsv/C516 --control_bam ./C516/bam/C516.bam --min_tumor_variant_read_num 1 --min_indel_size 1000 --processes 32 --min_tumor_VAF 0.02 --max_control_variant_read_num 0 --check_read_max_num 200 --control_panel_prefix /path/to/Analysis_supplements/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control
done < ./samples.txt

# Get SVs for HEK293T
while read SAMPLE; do
    bsub -G team154-vwork -o ./nanomonsv/${SAMPLE}.log  -q week -R'select[mem>20000] rusage[mem=20000] span[hosts=1]' -M20000 -n8 -J nanomonsv_${SAMPLE} \
        nanomonsv get ./${SAMPLE}/nanomonsv/${SAMPLE} ./${SAMPLE}/bam/${SAMPLE}.bam /path/to/Analysis_supplements/genome.fa --control_prefix ./F3R/nanomonsv/F3R --control_bam ./F3R/bam/F3R.bam --min_tumor_variant_read_num 1 --min_indel_size 1000 --processes 32 --min_tumor_VAF 0.02 --max_control_variant_read_num 0 --check_read_max_num 200 --control_panel_prefix /path/to/Analysis_supplements/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control
done < ./samples.txt

# RNA sequencing
# downloading the salmon index for hg38
bsub -G team154-vwork -q normal -n 1 -R'select[mem>10000] rusage[mem=10000] span[hosts=1]' -M 10000 -o wget.log -J wget \
    wget http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/salmon_sa_index__default.tgz

# running salmon quant
while read SAMPLE; do
    bsub -G team154-vwork -q normal -n 8 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M 32000 -o plex${SAMPLE}/salmon.log -J salmon${SAMPLE}\
        /path/to/salmon-latest_linux_x86_64/bin/salmon quant -p 8 -i /path/to/salmon-latest_linux_x86_64/index/default -l IU -1 ${SAMPLE}/${SAMPLE}_1.fq -2 ${SAMPLE}/${SAMPLE}_2.fq --validateMappings -o ${SAMPLE}/transcripts_quant
done < ./samples.txt
