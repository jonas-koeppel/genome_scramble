library("tximport")
library("tidyverse")
library("DESeq2")
library("purrr")
library("biomaRt")

setwd("/path/to/wd")

# Set up gene list
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 110)
gene_coordinates_RNA <- getBM(attributes=c('hgnc_symbol', 'chromosome_name', 'transcript_start', 'transcript_end', "ensembl_transcript_id", "transcript_biotype", "strand", "transcript_is_canonical"), mart=ensembl) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
  filter(chromosome_name %in% chr_list, transcript_biotype %in% c("protein_coding", "lncRNA", "snRNA", "snoRNA", "miRNA", "rRNA", "ribozyme"), hgnc_symbol != "", transcript_is_canonical == 1) %>%
  mutate(start = ifelse(strand == 1, transcript_start, transcript_end), end = ifelse(strand == 1, transcript_end, transcript_start)) %>%
  dplyr::select("chr" = "chromosome_name", "start", "end", "gene" = "hgnc_symbol", transcript_biotype)

# Get chromosome information
chr_list <- sprintf("chr%s",c(seq(1,22,1), "X", "Y"))
chr_sizes <- read_tsv("./input_data/GRCh38.chrom_sizes.txt", col_names = c("tmp", "chr", "end")) %>% mutate(start = 0, supp_reads = 0) %>% filter(chr != "chrY") %>%
  pivot_longer(c(start, end), names_to = "label", values_to = "start") %>% dplyr::select(chr, start, supp_reads) 


# HAP1 RNAseq Analysis
HAP1_RNAseq <- list()
clones_HAP1 <- c("C516R_c7", "C516R_c8", "C516R_c10", "C516R_c11")
clones_HAP1_G1 <- c("C516R_G1_A31", "C516R_G1_C31", "C516R_G1_D61", "C516R_G1_D21", "C516R_G1_D101", "C516R_G1_E12", "C516R_G1_G32" )
clones_HAP1_WT <- c(outer(c("C516", "C516R_c6", "C516R_c9"), c("_R1", "_R2"), paste0))
clones_HAP1_WT_G1 <- c(outer(c("C516R_G1_A71", "C516R_G1_C91", "C516R_G1_C22"), c("_R1", "_R2"), paste0))

for(i in 1:length(clones_HAP1)) {
  selection = c(clones_HAP1_WT, paste0(clones_HAP1[i], c("_R1", "_R2")))
  files_selected <- paste0("./prc_data/salmon_output/", selection, ".sf")
  txi_selected <- tximport(files_selected, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "scaledTPM")
  key_clone = tibble(name = selection, condition = c(rep("B", 6), rep("A", 2)))
  print(key_clone)
  dds <- DESeqDataSetFromTximport(txi_selected, key_clone, ~ condition)
  dds <- dds[rowSums(counts(dds)) > 9,] # filter out reads with at least 10 counts
  rld <- rlog(dds, blind = FALSE)
  dds <- estimateSizeFactors(dds)
  dds_diff <- DESeq(dds)
  
  
  dds_results <- results(dds_diff, contrast=c("condition","A","B"))
  dds_results <- as_tibble(dds_results) %>% mutate(gene = dds_results@rownames, sample = clones_HAP1[i])
  HAP1_RNAseq[[i]] <- dds_results
  write_tsv(dds_results, paste0("./prc_data/deseq2/", clones_HAP1[i], ".tsv"))
}

for(i in 1:length(clones_HAP1_G1)) {
  selection = c(clones_HAP1_WT_G1, paste0(clones_HAP1_G1[i], c("_R1", "_R2")))
  files_selected <- paste0("./prc_data/salmon_output/", selection, ".sf")
  txi_selected <- tximport(files_selected, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "scaledTPM")
  key_clone = tibble(name = selection, condition = c(rep("B", 6), rep("A", 2)))
  print(key_clone)
  dds <- DESeqDataSetFromTximport(txi_selected, key_clone, ~ condition)
  dds <- dds[rowSums(counts(dds)) > 9,] # filter out reads with at least 10 counts
  rld <- rlog(dds, blind = FALSE)
  dds <- estimateSizeFactors(dds)
  dds_diff <- DESeq(dds)
  
  
  dds_results <- results(dds_diff, contrast=c("condition","A","B"))
  dds_results <- as_tibble(dds_results) %>% mutate(gene = dds_results@rownames, sample = clones_HAP1_G1[i])
  HAP1_RNAseq[[i]] <- dds_results
  write_tsv(dds_results, paste0("./prc_data/deseq2/", clones_HAP1_G1[i], ".tsv"))
}

# HEK293T RNAseq Analysis
HEK293T_RNAseq <- list()
clones_HEK293T <- c("F3R_c6")
clones_HEK293T_WT <- c("F3RC_ctl_R1", "F3RC_ctl_R2")

for(i in 1:length(clones_HEK293T)) {
  selection = c(clones_HEK293T_WT, paste0(clones_HEK293T[i], c("_R1", "_R2")))
  files_selected <- paste0("./prc_data/salmon_output/", selection, ".sf")
  txi_selected <- tximport(files_selected, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "scaledTPM")
  key_clone = tibble(name = selection, condition = c(rep("B", 2), rep("A", 2)))
  print(key_clone)
  dds <- DESeqDataSetFromTximport(txi_selected, key_clone, ~ condition)
  dds <- dds[rowSums(counts(dds)) > 9,] # filter out reads with at least 10 counts
  rld <- rlog(dds, blind = FALSE)
  dds <- estimateSizeFactors(dds)
  dds_diff <- DESeq(dds)
  
  
  dds_results <- results(dds_diff, contrast=c("condition","A","B"))
  dds_results <- as_tibble(dds_results) %>% mutate(gene = dds_results@rownames, sample = clones_HEK293T[i])
  HAP1_RNAseq[[i]] <- dds_results
  write_tsv(dds_results, paste0("./prc_data/deseq2/", clones_HEK293T[i], ".tsv"))
}

# Annotate and filter RNA-seq data
# plot results across all clones
file_list <- list.files(pattern = "*.tsv", path = "./prc_data/deseq2", full.names = TRUE)
variant_genes <- read_tsv("./prc_data/structural_variation/filtered_variants/genes_in_variants.tsv")
RNAseq <- file_list %>% map_df(~ read_tsv(.x))
RNAseq_annotated <- RNAseq %>% left_join(gene_coordinates_RNA, by = "gene") %>% filter(!is.na(transcript_biotype), baseMean > 50) %>%
  left_join(dplyr::select(filter(variant_genes, category == "Deletion", location == "Within variant"), sample, gene, location), by = c("sample", "gene")) %>%
  mutate(location = ifelse(is.na(location), "Outside variant", "Within variant"))

