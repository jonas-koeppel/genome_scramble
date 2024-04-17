library(ShortRead)
library(tidyverse)

theme_sv <-   theme_bw(base_size = 7, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        text = element_text(family = "Helvetica"),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.25, "cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black"))


count_insertions <- function(path) {
  filename = sub("_S.*$", "", basename(path))
  print(filename)
  
  reads <- readFastq(path)
  reads <- as.character(reads@sread)
  matching_guide = str_subset(reads, "TTCTACCAGAGGTACA(.*?)AGG.GG") %>%
    str_extract("TTCTACCAGAGGTACA(.*?)AGG.GG")
  
  
  # first summing up all identical reads to save computation time later when pattern matching
  reads_sum <- tibble(sequence = matching_guide) %>%
    group_by(sequence) %>%
    summarise(count = dplyr::n()) %>%
    mutate(length = nchar(sequence))
  
  n_insert <- sum(filter(reads_sum, length == 56)$count)
  n_wt <- sum(filter(reads_sum, length == 22)$count)
  n_other <- sum(filter(reads_sum, !length %in% c(22, 56))$count)
  perc_guide <- sum(length(matching_guide)/length(reads)*100)
  perc_Ins <- n_insert/(n_insert+n_wt)*100
  return(tibble(n_insert, n_wt, n_other, perc_guide, perc_Ins, sample = filename))
}

setwd("/Users/jk24/Library/CloudStorage/OneDrive-Personal/PhD/Scramble/submission")

file_names_1 <- read_tsv("./raw_data/202211_LINE1_insertions/samples1/samples.txt", col_names = "file")
file_names_2 <- read_tsv("./raw_data/202211_LINE1_insertions/samples2/samples.txt", col_names = "file")
file_names_3 <- read_tsv("./raw_data/202211_LINE1_insertions/samples3/samples.txt", col_names = "file")
file_names_4 <- read_tsv("./raw_data/202211_LINE1_insertions/transfections/samples.txt", col_names = "file")
file_names_5 <- read_tsv("./raw_data/202211_LINE1_insertions/clones/samples.txt", col_names = "file")

# ==== Loading files ====
read_files <- function(files, path) {
  files <- mutate(files, path = paste0(path, file), sample = sub("_S.*$", "", file))
  data = list()
  
  for(i in 1:nrow(files)) {
    data[[i]] = count_insertions(files$path[i])
  }
  
  names(data) = files$sample
  
  return(data)
}


samples1 <- read_files(file_names_1, "./raw_data/202211_LINE1_insertions/samples1/") %>% bind_rows()
samples2 <- read_files(file_names_2, "./raw_data/202211_LINE1_insertions/samples2/") %>% bind_rows()
samples3 <- read_files(file_names_3, "./raw_data/202211_LINE1_insertions/samples3/") %>% bind_rows()
transfections <- read_files(file_names_4, "./raw_data/202211_LINE1_insertions/transfections/") %>% bind_rows()
clones <- read_files(file_names_5, "./raw_data/202211_LINE1_insertions/clones/") %>% bind_rows()

samples1 <- separate(samples1, col = sample, into = c("clone", "day", "replicate"), sep = "_", remove = F) %>%
  mutate(day = as.numeric(str_remove(day, "day")))

samples2 <- separate(samples2, col = sample, into = c("clone", "day", "replicate"), sep = "-", remove = F) %>%
  mutate(day = as.numeric(str_remove(day, "day")))

samples3 <- separate(samples3, col = sample, into = c("clone", "day", "replicate"), sep = "_", remove = F) %>%
  mutate(day = as.numeric(str_remove(day, "day")), replicate = str_remove(replicate, "R"))

transfections <- separate(transfections, col = sample, into = c("clone", "day", "replicate"), sep = "_", remove = F) %>%
  mutate(day = as.numeric(str_remove(day, "day")), replicate = str_remove(replicate, "R"), clone = "HEK293T")

clones <- separate(clones, col = sample, into = c("clone", "cell_line", "ID"), sep = "_", remove = F) %>%
  mutate(n = n_insert + n_wt + n_other, perc_other = n_other/(n_insert + n_wt)) %>% filter(n > 300)

timecourse <- bind_rows(samples1, samples2, samples3) %>%
  mutate(n = n_insert + n_wt + n_other, perc_other = n_other/(n_insert + n_wt), clone = ifelse(str_detect(clone, "c12"), "HAP1", "HEK293T")) %>% filter(replicate %in% c(1,2), n > 300)

write_tsv(timecourse, "./prc_data/amplicon_seq/editing_timecourse.tsv")
write_tsv(transfections, "./prc_data/amplicon_seq/transfection_timecourse.tsv")
write_tsv(clones, "./prc_data/amplicon_seq/editing_clones.tsv")
