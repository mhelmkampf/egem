### ===========================================================================================
### R code to prepare genotyping metadata file for
### "A plea for a conservative approach to the report of new species records from eDNA surveys"
### by Oscar Puebla and Martin Helmkampf (2022)
### ===========================================================================================

library("tidyverse")

setwd("/Users/martin/Documents/Projects/2_Other/egem/0_metadata")


# read in files
files <- read_lines("raw_data.fof")
labels <- read_delim("egem_labels.tsv", delim = "\t")


# convert to tibbles
r1r2 <- files %>%
  as.data.frame(x = matrix(., ncol = 2, byrow = TRUE)) %>%
  as_tibble() %>%
  rename("R1" = "V1", "R2" = "V2") %>%
  mutate(FileID = str_replace_all(R1, "([A-Z-]{1,5}[0-9]{5})[-_].*", "\\1")) %>%
  mutate(LaneF = str_replace_all(R1, ".*_L0{0,2}([1-9]).*", "\\1")) %>%
  select(FileID, R1, R2, LaneF) %>%
  arrange(FileID)

names <- labels %>%
  mutate(Name = paste(SampleID, Species, tolower(Site), sep = "")) %>%
  select(FileID, Name) %>%
  arrange(FileID)


# bind columns
meta <- left_join(r1r2, names) %>%
  select(Name, R1, R2, LaneF)


# write to file
write_delim(meta, "egem_meta.csv", col_names = FALSE, delim = ";")
