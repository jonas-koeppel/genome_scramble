# June, 13th, 2024. Jonas Koeppel
# This script takes raw variant calls from nanomonsv and performs the filtering steps to derive at the condensed set of variants used in the manuscript

library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(Repitools)
library(tidyverse)
library(plyranges)
library(spgs)

setwd("/Users/jk24/Library/CloudStorage/OneDrive-Personal/PhD/Scramble/submission/")

rc <- function(x) {toupper(spgs::reverseComplement(x))} # define a quick function to reverse compliment sequences
chr_sizes <- read_tsv("./input_data/GRCh38.chrom_sizes.txt", col_names = c("tmp", "chr", "end")) %>% mutate(start = 0, supp_reads = 0) %>% filter(chr != "chrY") %>%
  pivot_longer(c(start, end), names_to = "label", values_to = "start") %>% dplyr::select(chr, start, supp_reads) 
chr_list <- sprintf("chr%s",c(seq(1,22,1), "X", "Y"))
centromeres <- read_tsv("./input_data/hg38_centromeres.bed", col_names = c("chr", "start", "end", "name")) %>% group_by(chr) %>% summarise(start = min(start), end = min(end))
missing_variants <- read_tsv("/Users/jk24/Library/CloudStorage/OneDrive-Personal/PhD/Scramble/submission/prc_data/structural_variation/missing_variants.txt")

theme_sv <-   theme_bw(base_size = 7, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        text = element_text(family = "Helvetica"),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.25, "cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black"))

# read in all required files
LINE1_gr <- read_tsv("./input_data/LINE1_nick_3mm.tsv") %>% mutate(start = nick_coordinates - 300, end = nick_coordinates + 300) %>% dplyr::select(chr, start, end) %>% GRanges()
loxPsym_insertions <- read_tsv("./prc_data/insertions/loxPsym_insertions.tsv")
sample_annotation <- read_tsv("./prc_data/sample_annotation.txt")

C516_clonal <- filter(loxPsym_insertions, clone == "HAP1 loxPsym(301)", supp_reads > 5, AF > 0.5, chr %in% chr_list) %>% dplyr::select(chr, start, end) %>% mutate(start = start - 300, end = end + 300)
F3_insertions <- filter(loxPsym_insertions, clone == "HEK293T loxPsym(638)", supp_reads > 3, AF  > 0.1, chr %in% chr_list) %>% dplyr::select(chr, start, end) %>% mutate(start = start - 300, end = end + 300)

# find the exact start positions
liftover_coordinates <- function(x, insertion_sites) {
  print("Harmonizing coordinates") 
  x <- filter(x, SVTYPE != "INS", chr %in% chr_list, chr_2 %in% chr_list) %>% mutate(id = 1:nrow(.))
  # liftover start
  sv_liftover_start <- dplyr::select(x, chr, start) %>% mutate(end = start, sv_start = start) %>% GRanges()
  sv_liftover_start <- join_overlap_left(GRanges(insertion_sites), sv_liftover_start) %>% filter(!is.na(sv_start)) %>% annoGR2DF() %>% distinct()
  sv_liftover_start <- left_join(dplyr::select(sv_liftover_start, chr, start, sv_start), dplyr::select(x, chr,"sv_start" = "start", chr_2, end, SVLEN, supp_reads, SVTYPE, category, id), by = c("chr", "sv_start"), relationship = "many-to-many") %>% distinct()
  # liftover end
  sv_liftover_end <- dplyr::select(x, "chr" = "chr_2", "start" = "end") %>% mutate(end = start, sv_end = start) %>% GRanges()
  sv_liftover_end <- join_overlap_left(GRanges(insertion_sites), sv_liftover_end) %>% filter(!is.na(sv_end)) %>% annoGR2DF() %>% distinct()
  sv_liftover <- left_join(dplyr::select(sv_liftover_end, "chr_2" = "chr", "end" = "start", sv_end), dplyr::select(sv_liftover_start, chr, start, "sv_end" = "end", chr_2, end, supp_reads, SVTYPE, category, id), by = c("chr_2", "sv_end"), relationship = "many-to-many") %>% filter(!is.na(supp_reads)) %>% distinct()
  # final cleaning to combine reads that are now the same
  passing_rearrangements <- sv_liftover %>% group_by(chr, start, chr_2, end, SVTYPE, category) %>% summarise(supp_reads = sum(supp_reads), id = id[1]) %>% mutate(start = start + 300, end = end + 300) %>% ungroup()
}

# ==== Part1 read in SVs from nanomonsv and filter for Cre induced variants ====
sv_annotations <- tibble(
  orientation = c("+-c", "++c", "--c", "-+c", "+-t", "++t", "--t", "-+t", "+-l", "++l", "--l"),
  SVTYPE = c("DEL", "INVh2h", "INVt2t", "CIRC", "TRAh2t", "TRAh2h", "TRAt2t", "TRAt2h", "INS", "FBh2h", "FBt2t"),
  category = c("Deletion", "Inversion", "Inversion", "Circle", "Translocation", "Translocation", "Translocation", "Translocation", "Insertion", "Fold back", "Fold back"))

read_nanomonsv <- function(path) {
  read_tsv(path) %>%
    mutate(type = ifelse(Chr_1 != Chr_2, "t", ifelse(abs(Pos_2 - Pos_1) < 50, "l", "c")), orientation = paste0(Dir_1, Dir_2, type)) %>% left_join(sv_annotations, by = "orientation") %>%
    mutate(SVLEN = ifelse(!str_detect(SVTYPE , "TRA"), Pos_2 - Pos_1, 0)) %>% filter(Chr_1 %in% chr_list, Chr_2 %in% chr_list, type %in% c("l", "t") | SVLEN > 3000) %>%
    dplyr::select("chr" = "Chr_1", "start" = "Pos_1", "chr_2" = "Chr_2", "end" = "Pos_2", SVLEN, category, SVTYPE, orientation, "ins_seq" = "Inserted_Seq", "depth" = "Checked_Read_Num_Tumor", "supp_reads" = "Supporting_Read_Num_Tumor")
}

find_rearrangements_nano = function(x, insertion_sites = LINE1_gr, liftover_sites = C516_clonal) {
  # find start and end sites that are overlapping with loxPsym site insertions
  print("find start and end sites that are overlapping with loxPsym site insertions")
  start_sites <- x %>% filter(str_detect(ins_seq, "TAACTTCGTAT|ATACGAAGTTA|AATGTACATTAT|GTATAATGTAC")) %>% dplyr::select(chr, start) %>% mutate(end = start)
  end_sites <- x %>% filter(str_detect(ins_seq, "TAACTTCGTAT|ATACGAAGTTA|AATGTACATTAT|GTATAATGTAC")) %>% dplyr::select("chr" = "chr_2", "start" = "end") %>% mutate(end = start)
  gr_start <- makeGRangesFromDataFrame(start_sites, keep.extra.columns=T)
  gr_end <- makeGRangesFromDataFrame(end_sites, keep.extra.columns=T)
  overlaps_start <- subsetByOverlaps(gr_start, insertion_sites) %>% annoGR2DF() %>% mutate(region = paste0(chr, ":", start)) %>% `$`(region)
  overlaps_end <- subsetByOverlaps(gr_end, insertion_sites) %>% annoGR2DF() %>% mutate(region = paste0(chr, ":", start)) %>% `$`(region)
  
  print("filtering variants")
  passing_rearrangements <- x %>% mutate(region_start = paste0(chr, ":", start), region_end = paste0(chr_2, ":", end)) %>% filter(region_start %in% overlaps_start & region_end %in% overlaps_end, chr %in% chr_list, chr_2 %in% chr_list, SVTYPE != "INS",  chr %in% chr_list, chr_2 %in% chr_list)
  
  # deal with zero cases
  if(nrow(passing_rearrangements) == 0) {
    print("No rearrangements")
    write_tsv(passing_rearrangements, paste0("./prc_data/structural_variation/cre_induced//", x$sample[1], "_nano.tsv"))
  } else {
    passing_rearrangements <- mutate(passing_rearrangements, id = 1:nrow(passing_rearrangements))
    passing_rearrangements_lifted <- liftover_coordinates(passing_rearrangements, insertion_sites = liftover_sites)
    print(paste(nrow(passing_rearrangements_lifted), "rearrangements harmonized"))
    passing_rearrangements <- passing_rearrangements %>% left_join(dplyr::select(passing_rearrangements_lifted, id, "start_lift" = "start", "end_lift" = "end"), by = "id") %>%
      mutate(start = ifelse(!is.na(start_lift), start_lift, start), end = ifelse(!is.na(end_lift), end_lift, end)) %>% dplyr::select(-start_lift, -end_lift)
    write_tsv(passing_rearrangements, paste0("./prc_data/structural_variation/cre_induced/", x$sample[1], "_nano.tsv"))
  }
}

# parentals
C516_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516_par.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516", sample = "C516") %>% find_rearrangements_nano()
F3R_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R.nanomonsv.result.txt") %>% mutate(cell_line = "HEK293T", clone = "F3R", sample = "F3R") %>% find_rearrangements_nano(liftover_sites = F3_insertions)

# pools early
C516_d1_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516_5Cre_d1.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516", sample = "C516_d1") %>% find_rearrangements_nano()
C516_d1b_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516_5Cre_d1_as.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516", sample = "C516_d1as") %>% find_rearrangements_nano()
C516R_PB_d1_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_PB_d1.nanomonsv.result.txt") %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_PB_d1") %>% find_rearrangements_nano()
F3R_d1_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_d1.nanomonsv.result.txt") %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_d1") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_d1_mn_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_d1_mn.nanomonsv.result.txt") %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_d1_mn") %>% find_rearrangements_nano(liftover_sites = F3_insertions)

# pools late
C516_d14_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516_d14.nanomonsv.result.txt") %>% mutate(cell_line = "HAP1", clone = "C516", sample = "C516_d14") %>% find_rearrangements_nano()
C516_d14_mn_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516_5Cre_d14.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516", sample = "C516_d14_mn") %>% find_rearrangements_nano()
C516_d14b_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516_d14b.nanomonsv.result.txt") %>% mutate(cell_line = "HAP1", clone = "C516", sample = "C516_d14b") %>% find_rearrangements_nano()
C516R_d13_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_d13.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_d13") %>% find_rearrangements_nano()
C516R_d13_mn_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_d13_mn.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_d13_mn") %>% find_rearrangements_nano()
C516R_d13b_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_d13b.nanomonsv.result.txt") %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_d13b") %>% find_rearrangements_nano()
C516R_d13c_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_d13c.nanomonsv.result.txt") %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_d13c") %>% find_rearrangements_nano()
C516R_PB_d17_BFP_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_PB_d17_BFP.nanomonsv.result.txt") %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_PB_d17_BFP") %>% find_rearrangements_nano()
F3R_d15_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_d15.nanomonsv.result.txt") %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_d15") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_d15b_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_d15b.nanomonsv.result.txt") %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_d15b") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_d15c_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_d15c.nanomonsv.result.txt") %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_d15c") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_d15d_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_d15d.nanomonsv.result.txt") %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_d15d") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_d15_mn_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3_d15.nanomonsv.result.txt") %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_d15_mn") %>% find_rearrangements_nano(liftover_sites = F3_insertions)

# clones
C516R_c1t6_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_c1t6.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_c1t6") %>% find_rearrangements_nano()
C516R_c7t12_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_c7t12.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_c7t12") %>% find_rearrangements_nano()
C516R_c11_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_c11.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_c11") %>% find_rearrangements_nano()
F3R_c6_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3_c6.nanomonsv.result.txt") %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_c6") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_c1t4_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_c1t4.nanomonsv.result.txt") %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_c1t4") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_c5t8_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_c5t8.nanomonsv.result.txt") %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_c5t8") %>% find_rearrangements_nano(liftover_sites = F3_insertions)

C516R_c6_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_c6.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_c6") %>% find_rearrangements_nano()
C516R_c7_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_c7.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_c7") %>% find_rearrangements_nano()
C516R_c8_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_c8.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_c8") %>% find_rearrangements_nano()
C516R_c9_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_c9.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_c9") %>% find_rearrangements_nano()
C516R_c10_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_c10.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_c10_mn") %>% find_rearrangements_nano()
C516R_c11_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_c11.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_c11_mn") %>% find_rearrangements_nano()

C516R_G1_C22_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_C22.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_C22") %>% find_rearrangements_nano()
C516R_G1_C91_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_C91.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_C91") %>% find_rearrangements_nano()
C516R_G1_D61_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_D61.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_D61") %>% find_rearrangements_nano()
C516R_G1_D21_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_D21.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_D21") %>% find_rearrangements_nano()

C516R_G1_A71_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_A71.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_A71") %>% find_rearrangements_nano()
C516R_G1_A31_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_A31.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_A31") %>% find_rearrangements_nano()
C516R_G1_B41_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_B41.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_B41") %>% find_rearrangements_nano()
C516R_G1_C31_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_C31.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_C31") %>% find_rearrangements_nano()
C516R_G1_C32_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_C32.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_C32") %>% find_rearrangements_nano()

C516R_G1_C72_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_C72.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_C72") %>% find_rearrangements_nano()
C516R_G1_D22_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_D22.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_D22") %>% find_rearrangements_nano()
C516R_G1_D101_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_D101.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_D101") %>% find_rearrangements_nano()
C516R_G1_E12_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_E12.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_E12") %>% find_rearrangements_nano()
C516R_G1_G32_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/C516R_G1_G32.nanomonsv.result.txt")  %>% mutate(cell_line = "HAP1", clone = "C516R", sample = "C516R_G1_G32") %>% find_rearrangements_nano()

F3R_c1_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_c1.nanomonsv.result.txt")  %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_c1") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_c4_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_c4.nanomonsv.result.txt")  %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_c4") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_c5_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_c5.nanomonsv.result.txt")  %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_c5") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_c10_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_c10.nanomonsv.result.txt")  %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_c15") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_c11_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_c11.nanomonsv.result.txt")  %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_c11") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_c12_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_c12.nanomonsv.result.txt")  %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_c12") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_c13_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_c13.nanomonsv.result.txt")  %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_c13") %>% find_rearrangements_nano(liftover_sites = F3_insertions)
F3R_c15_nano <- read_nanomonsv("./prc_data/structural_variation/raw_variants/F3R_c15.nanomonsv.result.txt")  %>% mutate(cell_line = "HEK293T", clone = "HEK293T_F3R", sample = "F3R_c15") %>% find_rearrangements_nano(liftover_sites = F3_insertions)

# ==== Part2 read in filtered variants, merge and perform final filtering ====

# read in files
depth_files <- read_tsv("./prc_data/depth_files/samples.txt", col_names = "file")
rearrangement_files <- read_tsv("./prc_data/structural_variation/cre_induced/samples.txt", col_names = "file")

read_depth <- function(files, path) {
  files <- mutate(files, path = paste0(path, file))
  data = list()
  
  for(i in 1:nrow(files)) {
    print(files$path[i])
    data[[i]] = read_tsv(files$path[i], col_names = c("chr", "start", "end", "depth")) %>% mutate(filename = str_remove(files$file[i], ".bedgraph"))
  }
  
  names <- vector("character", length = length(data))
  for(i in 1:length(data)) {
    names[i] = str_remove(files$file[i], ".bedgraph")
  }
  
  names(data) = names
  
  return(data)
}
read_rearrangements <- function(files, path) {
  files <- mutate(files, path = paste0(path, file))
  data = list()
  
  for(i in 1:nrow(files)) {
    print(files$path[i])
    data[[i]] = read_tsv(files$path[i], col_types = c("c", "n", "c", "n", "n", "c", "c", "c", "c", "n", "n", "c", "c", "c", "c", "c", "c")) %>% mutate(filename = str_remove(files$file[i], "_nano.tsv"))
  }
  
  names <- vector("character", length = length(data))
  for(i in 1:length(data)) {
    names[i] = str_remove(files$file[i], ".tsv")
  }
  
  names(data) = names
  
  return(data)
}

depth <- read_depth(depth_files, path = "./prc_data/depth_files/") %>% bind_rows() %>% left_join(sample_annotation, by = "filename")
rearrangements <- read_rearrangements(rearrangement_files, path = "./prc_data/structural_variation/cre_induced/") %>% bind_rows() %>% dplyr::select(-sample) %>% left_join(sample_annotation, by = "filename")

# generate summary statistics and combine with coverages
coverages_filename <- depth %>%
  group_by(filename) %>%
  summarise(coverage = sum(depth)/3300000000)

coverages <- depth %>%
  group_by(sample) %>%
  summarise(coverage = sum(depth)/3300000000)

coverages_old <- depth_old %>%
  group_by(filename) %>%
  summarise(coverage = sum(depth)/3300000000)

coverages_replicates <- depth %>%
  group_by(sample, group) %>%
  summarise(coverage = sum(depth)/3300000000)

# filter rearrangements by depth
sv_filtered <- rearrangements %>% 
  mutate(sample = str_remove(sample, "_nano"), depth = ifelse(category == "Fold back", depth*2, depth)) %>%
  left_join(coverages_filename, by = "filename") %>%
  filter(depth < coverage * 5 | (chr == "chr15" & depth < coverage * 10)) %>% 
  group_by(chr, start, chr_2, end, category, SVTYPE, cell_line, arrest, timepoint, clone, sample) %>% # combine svs that are the same
  summarise(supp_reads = sum(supp_reads)) %>%
  left_join(coverages, by = "sample") %>%
  mutate(SVLEN = ifelse(category %in% c("Translocation", "Fold back"), 0, end-start)) %>%
  filter(!(timepoint == "clones" & supp_reads == 1 & category == "Fold back")) %>%
  bind_rows(dplyr::select(missing_variants, -orientation))

sv_filtered_replicate <- rearrangements %>% 
  mutate(depth = ifelse(category == "Fold back", depth*2, depth)) %>%
  left_join(coverages_filename, by = "filename") %>%
  filter(depth < coverage * 5 | (chr == "chr15" & depth < coverage * 10)) %>% 
  group_by(chr, start, chr_2, end, category, SVTYPE, cell_line, clone, sample, group, timepoint, replicate, arrest) %>% # combine svs that are the same
  summarise(supp_reads = sum(supp_reads)) %>%
  left_join(coverages_replicates,  by = c("sample", "group")) %>%
  mutate(SVLEN = ifelse(category %in% c("Translocation", "Fold back"), 0, end-start))


write_tsv(sv_filtered, "./prc_data/structural_variation/filtered_variants/sv_filtered.tsv")

# generate early and late sets
rearrangements_early <- filter(sv_filtered, timepoint == "early")
rearrangements_late_pool <- filter(sv_filtered, timepoint == "late")
rearrangements_clones <- filter(sv_filtered, timepoint == "clones")
rearrangements_late <- bind_rows(rearrangements_late_pool, mutate(rearrangements_clones, supp_reads = 1))

write_tsv(rearrangements_late_pool, "./prc_data/structural_variation/filtered_variants/late_pool_rearrangements.tsv")
write_tsv(rearrangements_early, "./prc_data/structural_variation/filtered_variants/early_rearrangements.tsv")
write_tsv(rearrangements_clones, "./prc_data/structural_variation/filtered_variants/clones_rearrangements.tsv")
write_tsv(rearrangements_late, "./prc_data/structural_variation/filtered_variants/late_rearrangements.tsv")

ungroup(rearrangements_early) %>% filter(clone == "C516", category %in% c("Deletion", "Inversion")) %>% dplyr::select(chr, start, chr_2, end, SVLEN, category, SVTYPE, supp_reads, clone, sample) %>% write_tsv("./prc_data/structural_variation/filtered_variants/early_rearrangements_thomas.tsv")
ungroup(rearrangements_late) %>% filter(cell_line == "HAP1", category %in% c("Deletion", "Inversion")) %>%  dplyr::select(chr, start, chr_2, end, SVLEN, category, SVTYPE, supp_reads, clone, sample) %>% write_tsv("./prc_data/structural_variation/filtered_variants/late_rearrangements_thomas.tsv")

