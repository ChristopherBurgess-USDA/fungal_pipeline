library(dada2)
library(tidyverse)
library(phyloseq)


asv_to_otu_phyloseq <- function(mf){
  phylo_obj <- read_rds(
    paste0("output_data/asv/its2_asv_", mf, "_phyloseq.rds")
  )

  asv_to_otu_map <- read_tsv(
    paste0("temp/04_vsearch/01_asv_otu_", mf, "_map.uc"),
    col_names = FALSE
  )

  asv_to_otu_hash <- setNames(asv_to_otu_map$X10, asv_to_otu_map$X9)

  otu_data <- phylo_obj %>%
    otu_table() %>%
    as("matrix") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "asv_id") %>%
    mutate(otu_id = asv_to_otu_hash[asv_id]) %>%
    select(-asv_id) %>%
    group_by(otu_id) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames(var = "otu_id") %>%
    as.matrix() %>%
    otu_table(taxa_are_rows = TRUE)

  taxa_data <- read_tsv(
      paste0("temp/04_vsearch/05_its2_otu_", mf, "_unite-04-04-2024_taxanomy.tsv"),
      col_names = c("otu_id", "Taxon", "Confidence"),
      skip = 1
    ) %>%
    mutate(Taxon = str_remove_all(Taxon, "[:alpha:]{1}__")) %>%
    separate_wider_delim(
      Taxon, ";",
      names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "sh_id"),
      too_few = "align_start"
    ) %>%
    select(otu_id, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    column_to_rownames(var = "otu_id") %>%
    as.matrix() %>%
    tax_table()
  taxa_data %>%
    as("matrix") %>%
    as.data.frame() %>%
    rownames_to_column(var = "otu_id") %>%
    write_tsv(paste0("output_data/its2_otu_", mf, "_taxonomy_unite_04_04_2024.tsv"))

  temp <- taxa_data %>%
        as.matrix() %>%
        as.data.frame() %>%
        mutate(across(Kingdom:Species, ~str_remove(.x, "[:alpha:]__"))) |>
        rownames_to_column(var = "otu_id") %>%
        unite("taxonomy", Kingdom:Species, sep = ";")
    
  otu_data %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column(var = "otu_id") %>%
      select(otu_id, last_col()) %>%
      left_join(temp) %>%
      write_tsv(paste0("temp/05_funguild/its2_otu_", mf, "_funguild_input.tsv"))

  otu_phylo <- phyloseq(sample_data(phylo_obj), taxa_data, otu_data)
  write_rds(
      otu_phylo, paste0("output_data/its2_otu_", mf, "_phyloseq.rds")
    )

}

ky_metadata_merged_phylo <- read_csv("output_data/ky_metadata_composite.csv") |>
  mutate(row_id = sample_id) |>
  column_to_rownames(var = "row_id") %>%
  sample_data()


asv_to_otu_phyloseq("forward")

asv_to_otu_phyloseq("merged")
