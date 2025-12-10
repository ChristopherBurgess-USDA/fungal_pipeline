library(fs)
library(tidyverse)
library(ShortRead)
library(dada2)

###############################################################################
#*****************************************************************************#
#                            Defining Functions                               #
#*****************************************************************************#
###############################################################################

get_n <- function(x) sum(getUniques(x))

make_df <- function(input_data, c_name){
###############################################################################
# takes a named list of sequencing reads and returns them as a data frame.

# Args:
#    input_data:  named list of sequencing reads
#        c_name: a string with to name the sequencing read column by
# Returns:
#    A data frame with 2 columns sample names and the reads kept named c_name
###############################################################################
    input_data %>%
         tibble(sample_names = names(.), !!c_name := .) %>%
         select(sample_names, !!c_name)
}



its_pipeline <- function(seq_run){
###############################################################################
# This function returns denoised its reads for a given seq_run by 
# first filtering/trimming them, then denoising them and removing chimeras.

# Args:
#    seq_run: The seq_run id which is used to references the correct reads for a given sequencing run.

# Returns:
#    A list of 3 objects:
#          i) merged forward and reverse ITS2 as a dada2 sequence table obj
#          ii) only forward ITS2 reads as a dada2 sequence table obj
#          iii) a data frame of the number of reads kept at each step
###############################################################################


#*****************************************************************************#
#                            Data Import                                      #
#*****************************************************************************#


    filter_path <- "temp/02_qc_filtered"

    itsxpress_forward <- dir_ls("temp/00_itsxpress/reads", glob = paste0("*",seq_run, "_*R1.fastq.gz"))
    itsxpress_reverse <- dir_ls("temp/00_itsxpress/reads", glob = paste0("*",seq_run, "_*R2.fastq.gz"))
    sample_names <- str_remove_all(itsxpress_forward, "temp/00_itsxpress/reads/|_R1.fastq.gz")

    filter_forward <- paste0(sample_names, "_R1") %>%
        fs::path(filter_path, ., ext = "fastq.gz")
    
    filter_reverse <- paste0(sample_names, "_R2") %>%
        fs::path(filter_path, ., ext = "fastq.gz")
    
    names(filter_forward) <- sample_names
    names(filter_reverse) <- sample_names


#*****************************************************************************#
#                            Filter and Trimming                              #
#*****************************************************************************#

    out <- filterAndTrim(
        fwd = itsxpress_forward, filt = filter_forward,
        rev = itsxpress_reverse, filt.rev = filter_reverse,
        minLen = 50,
        maxN = 0, maxEE = c(2, 2), truncQ = 8,
        rm.phix = T, compress = T, multithread = T, verbose = T
    )


#*****************************************************************************#
#                            Tracking reads                                   #
#*****************************************************************************#

    out <- out %>%
        as_tibble() %>%
        mutate(sample_names = sample_names) %>%
        dplyr::rename(input = reads.in, filtered = reads.out)
    
    write_rds(out, paste0("temp/03_dada/reads_", seq_run, ".rds"))


#*****************************************************************************#
#                            Dada2 denoising                                  #
#*****************************************************************************#

    filter_forward <- dir_ls(filter_path, glob = paste0("*",seq_run, "_*R1.fastq.gz"))
    filter_reverse <- dir_ls(filter_path, glob = paste0("*",seq_run, "_*R2.fastq.gz"))
    sample_names <- str_remove_all(filter_forward, "temp/02_qc_filtered/|_R1.fastq.gz")

    names(filter_forward) <- sample_names
    names(filter_reverse) <- sample_names


    err_forward <- learnErrors(filter_forward, multithread = T)
    err_reverse <- learnErrors(filter_reverse, multithread = T)

    dada_forward <- dada(filter_forward, err = err_forward, multithread = TRUE)
    dada_reverse <- dada(filter_reverse, err = err_reverse, multithread = TRUE)

#*****************************************************************************#
#                            saving intermediate files                        #
#*****************************************************************************#
    write_rds(err_forward, paste0("temp/03_dada/01_err_forward_", seq_run, ".rds"))
    write_rds(err_reverse, paste0("temp/03_dada/01_err_reverse_", seq_run, ".rds"))

    write_rds(dada_forward, paste0("temp/03_dada/02_dada_forward_", seq_run, ".rds"))
    write_rds(dada_reverse, paste0("temp/03_dada/02_dada_reverse_", seq_run, ".rds"))



#*****************************************************************************#
#                            removing Chimeras                                #
#*****************************************************************************#

    merge_seq <- mergePairs(dada_forward, filter_forward, dada_reverse, filter_reverse, verbose=TRUE)

    seqtab <- makeSequenceTable(merge_seq)

    seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = T, verbose = T)
    
    seqtab_forward <- makeSequenceTable(dada_forward)

    seqtab_forward_nochim <- removeBimeraDenovo(seqtab_forward, method = "consensus", multithread = T, verbose = T)

#*****************************************************************************#
#                            saving intermediate files                        #
#*****************************************************************************#
    write_rds(seqtab, paste0("temp/03_dada/03_seqtab_", seq_run, ".rds"))
    write_rds(seqtab_nochim, paste0("temp/03_dada/04_seqtab_nochim_", seq_run, ".rds"))
    write_rds(seqtab_forward_nochim, paste0("temp/03_dada/04_seqtab_forward_nochim_", seq_run, ".rds"))


#*****************************************************************************#
#                            combining read counts at each step               #
#*****************************************************************************#

    reads_tracking <- out %>%
        left_join(make_df(map_int(dada_forward, get_n), "denoised_forward")) %>%
        left_join(make_df(map_int(dada_reverse, get_n), "denoised_reverse")) %>%
        left_join(make_df(rowSums(seqtab), "merged")) %>%
        left_join(make_df(rowSums(seqtab_nochim), "no_chimera")) %>%
        left_join(make_df(rowSums(seqtab_forward_nochim), "no_forward_chim")) %>%
        mutate(seq_id = seq_run) %>%
        select(
            sample_id = sample_names, seq_id,
            input, filtered,
            denoised_forward, denoised_reverse,
            merged, no_chimera, no_forward_chim
        )
    
#*****************************************************************************#
#                            Returns all needed data                          #
#*****************************************************************************#
    return(list(seqtab_nochim, seqtab_forward_nochim, reads_tracking))
}

###############################################################################
#*****************************************************************************#
#                            calling its_pipeline function                    #
#*****************************************************************************#
###############################################################################

run1_seqtabs <- its_pipeline("run1")

run2_seqtabs <- its_pipeline("run2")


#*****************************************************************************#
#                            Merging plates 1,2 mlsh runs                     #
#*****************************************************************************#
mlsh_run2_seqtabs_nochim <- read_rds("mlsh_import_data/04_seqtab_nochim_run2.rds")
mlsh_run2_forward_nochim <- read_rds("mlsh_import_data/04_seqtab_forward_nochim_run2.rds")

mlsh_run3_seqtabs_nochim <- read_rds("mlsh_import_data/04_seqtab_nochim_run3.rds")
mlsh_run3_forward_nochim <- read_rds("mlsh_import_data/04_seqtab_forward_nochim_run3.rds")


seqtab_nochim_merged <- mergeSequenceTables(
    run1_seqtabs[[1]],
    run2_seqtabs[[1]],
    mlsh_run2_seqtabs_nochim,
    mlsh_run3_seqtabs_nochim,
    repeats = "sum"
)

seqtab_nochim_forward <- mergeSequenceTables(
    run1_seqtabs[[2]],
    run2_seqtabs[[2]],
    mlsh_run2_forward_nochim,
    mlsh_run3_forward_nochim,
    repeats = "sum"
)

write_rds(seqtab_nochim_merged, "temp/03_dada/05_merged_seqtab.rds")
write_rds(seqtab_nochim_forward, "temp/03_dada/05_forward_seqtab.rds")

#*****************************************************************************#
#          Merging read tracking info with itsxpress results and saving       #
#*****************************************************************************#

reads_tracking <- bind_rows(run1_seqtabs[[3]], run2_seqtabs[[3]])

read_csv("temp/00_itsxpress/reads_track.csv") %>%
    left_join(reads_tracking) %>%
    write_csv("temp/03_dada/reads_pipeline_tracking.csv")

#*****************************************************************************#
#                            Assigning taxanomy                               #
#*****************************************************************************#

unite_ref <- "../../databases/sh_general_release_dynamic_04.04.2024.fasta"

taxa <- assignTaxonomy(seqtab_nochim_merged, unite_ref, multithread = T, tryRC = T)
taxa_forward <- assignTaxonomy(seqtab_nochim_forward, unite_ref, multithread = T, tryRC = T)

#*****************************************************************************#
#              Saving final sequencing data ready for phyloseq import         #
#*****************************************************************************#

write_rds(taxa, "temp/03_dada/06_merged_unite.rds")
write_rds(taxa_forward, "temp/03_dada/06_forward_unite.rds")