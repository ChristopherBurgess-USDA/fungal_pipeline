library(dada2)
library(tidyverse)
library(phyloseq)


saving_final_results <- function(seqtab_obj, taxa_obj, mf){

    phylo_obj <- phyloseq(
        ky_metadata_phylo,
        tax_table(taxa_obj),
        otu_table(seqtab_obj, taxa_are_rows = F)
    ) %>%
        prune_taxa(x = ., taxa_sums(x = .) >= 1)

    seqs <- Biostrings::DNAStringSet(taxa_names(phylo_obj))
    names(seqs) <- taxa_names(phylo_obj)
    phylo_obj <- merge_phyloseq(phylo_obj, seqs)
    taxa_names(phylo_obj) <- paste0("asv_", seq(1:length(seqs)))

    write_rds(
        phylo_obj,paste0("output_data/asv/its2_asv_", mf, "_phyloseq.rds")
    )

    refseq(phylo_obj) |>
        Biostrings::writeXStringSet(
            paste0("output_data/asv/its2_asv_", mf, "_rep_set.fasta"),
            append = F, compress = F, format = "fasta"
        )

    phylo_obj |>
        tax_table() |>
        as.matrix() |>
        as.data.frame() |>
        rownames_to_column(var = "asv_id") |>
        select(asv_id, everything()) |>
        write_tsv(paste0("output_data/asv/its2_asv_", mf, "_taxonomy_unite_04_04_2024.tsv"))

    phylo_obj |>
        otu_table() |>
        as.matrix() |>
        as.data.frame() |>
        rownames_to_column(var = "sample_seq_id") |>
        select(sample_seq_id, everything()) |>
        write_tsv(paste0("output_data/asv/its2_asv_", mf, "_table.tsv"))
    
    temp <- phylo_obj |>
        tax_table() |>
        as.matrix() |>
        as.data.frame() |>
        mutate(across(Kingdom:Species, ~str_remove(.x, "[:alpha:]__"))) |>
        rownames_to_column(var = "asv_id") |>
        unite("taxonomy", Kingdom:Species, sep = ";")

    phylo_obj |>
        otu_table() |>
        as.matrix() |>
        t() |>
        as.data.frame() |>
        rownames_to_column(var = "asv_id") |>
        select(asv_id, last_col()) |>
        left_join(temp) |>
        write_tsv(paste0("temp/05_funguild/its2_asv_", mf, "_funguild_input.tsv"))
    
}

seqtab_nochim_merged <- read_rds("temp/03_dada/05_merged_seqtab.rds")
seqtab_nochim_forward <- read_rds("temp/03_dada/05_forward_seqtab.rds")
taxa <- read_rds("temp/03_dada/06_merged_unite.rds")
taxa_forward <- read_rds("temp/03_dada/06_forward_unite.rds")

reads_tracking <- read_csv("temp/03_dada/reads_pipeline_tracking.csv")
mlsh_run_metadata <- read_csv("mlsh_import_data/ky21_mlsh_seq_run_meta_data.csv") |>
     mutate(
        seq_id = paste0("mlsh-", seq_id),
        sample_seq_id = sample_id,
        sample_id = str_remove_all(sample_id, "(KY21|MLSH)_ITS2_")
    ) |>
    select(sample_seq_id, sample_id, seq_id, sample_qc = sample_type, everything())

###############################################################################
#*****************************************************************************#
#                            Saving Final results                             #
#*****************************************************************************#
###############################################################################

#*****************************************************************************#
#                            formating metadata                               #
#*****************************************************************************#


ky_metadata <- reads_tracking |>
    mutate(
        sample_qc = case_when(
            str_detect(sample_id, "Mock|PCR|Blank|Balnk") ~ "QC",
            T ~ "KY21"
        ),
        seq_id = paste0("ky-", seq_id),
        sample_seq_id = sample_id
    ) |>
    select(sample_seq_id, sample_id, seq_id, sample_qc, everything()) |>
    bind_rows(mlsh_run_metadata) |>
    separate_wider_delim(
        sample_id,
        delim = "-",
        names = c("temp_id", "plate_loc"),
        too_few = "align_start",
        cols_remove = FALSE
    ) |>
    separate_wider_delim(
        temp_id,
        delim = "_",
        names = c("drop_col", "plot_id", "depth", "sampling_point", "sample_type"),
        too_few = "align_start"
    ) |>
    mutate(
        across(
            c("plot_id", "depth", "sampling_point", "sample_type"),
            ~if_else(sample_qc == "QC", NA, .x)
        ),
        sample_id = str_remove_all(sample_id, "run[:digit:]_|-.+")
    ) |>
    select(
        sample_seq_id, sample_id, sample_qc, seq_id,
        plot_id, depth, sampling_point, sample_type,
        everything(), -drop_col
    )

ky_metadata_merged <- ky_metadata |>
    group_by(
        sample_id, sample_qc, plot_id, depth, sampling_point, sample_type
    ) |>
    summarise(
        across(where(is.numeric), sum)
    ) |>
    ungroup()

write_csv(ky_metadata, "output_data/ky_metadata.csv")
write_csv(ky_metadata_merged, "output_data/ky_metadata_composite.csv")

ky_metadata_phylo <- ky_metadata %>%
    column_to_rownames(var = "sample_seq_id") %>%
    sample_data()

ky_metadata_merged_phylo <- ky_metadata_merged |>
    mutate(row_id = sample_id) |>
    column_to_rownames(var = "row_id") %>%
    sample_data()



# #*****************************************************************************#
# #                            Calling saving_final_results                     #
# #*****************************************************************************#

saving_final_results(seqtab_nochim_merged, taxa, "merged")

saving_final_results(seqtab_nochim_forward, taxa_forward, "forward")