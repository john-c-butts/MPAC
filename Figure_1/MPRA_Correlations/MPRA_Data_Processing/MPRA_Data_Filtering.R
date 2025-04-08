# import packages
library(dplyr)
library(tidyr)

# Define emVar function
def_emVars <- function(df) {
  df <- df %>%
    # if there isn't a skew_logPadj value, replace it with zero
    tidyr::replace_na(list(Skew_logPadj = 0)) %>%
    # within that, emVars (lenients) are variants within active elements whose 
    dplyr::mutate(active = ifelse(logPadj_BF >= -log10(0.01) & abs(log2FC) >= 1, T, F),
                  emVar = ifelse(active & Skew_logPadj >= -log10(0.1) & !is.na(Skew_logPadj) & abs(Log2Skew) >= 0, T, F)) %>%
    dplyr::group_by(variant) %>%
    dplyr::mutate(active_any = ifelse(any(active), T, F),
                  emVar_any = ifelse(any(emVar), T, F),
                  emVar_all = ifelse(all(emVar),T,F))
  return(df) }

# Open GTEx RDS
gtex_mpra_paired_final20230117 <- readRDS("~/Dropbox (JAX)/for_john/empirical_ukbb_gtex_data/gtex_mpra_paired_final20230117.rds")
# Filter for best library, and Malinois Cell Types
gtex_mpra_paired_filter <- gtex_mpra_paired_final20230117 %>% dplyr::filter(best_library == 1,
                                                                                      type != "other_test")
# Add emVar and other columns
gtex_mpra_paired_filter_emvars <- def_emVars(gtex_mpra_paired_filter)

gtex_paired_no_filter_emvars <- def_emVars(gtex_mpra_paired_final20230117)
# Define PRCs
# Filtered
eqtl_mpra_prc_df <- gtex_mpra_paired_filter_emvars %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(pip = max(pip, na.rm = T)) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  filter(row_number() == 1) %>%
  ungroup() %>% 
  dplyr::mutate(causal = case_when((pip > 0.9 & abs(z) > qnorm(1 - 5*10^-8)) & (type !='3tissue_CS') ~ TRUE,
                                   (type == "3tissue_CS" & pip < 0.02) ~ FALSE,
                                   TRUE ~ NA)) %>% 
  dplyr::filter(!is.na(causal)) %>% 
  dplyr::filter(!consequence %in% c("synonymous","missense","LoF"))

sampsize <- min(table(eqtl_mpra_prc_df$causal))
eqtl_mpra_prc_df <- eqtl_mpra_prc_df %>%
  group_by(causal) %>%
  sample_n(sampsize)
# Calculate Recall
eqtl_total_causal <- length(which(eqtl_mpra_prc_df$causal, TRUE))
eqtl_mpra_pos <- length(which(eqtl_mpra_prc_df$emVar_any, TRUE))

eqtl_recall <- (eqtl_mpra_pos / eqtl_total_causal)
# No Filter
eqtl_mpra_prc_nf_df <- gtex_paired_no_filter_emvars %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(pip = max(pip, na.rm = T)) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  filter(row_number() == 1) %>%
  ungroup() %>% 
  dplyr::mutate(causal = case_when((pip > 0.9 & abs(z) > qnorm(1 - 5*10^-8)) & (type !='3tissue_CS') ~ TRUE,
                                   (type == "3tissue_CS" & pip < 0.02) ~ FALSE,
                                   TRUE ~ NA)) %>% 
  dplyr::filter(!is.na(causal)) %>% 
  dplyr::filter(!consequence %in% c("synonymous","missense","LoF"))

sampsize <- min(table(eqtl_mpra_prc_nf_df$causal))
eqtl_mpra_prc_nf_df <- eqtl_mpra_prc_nf_df %>%
  group_by(causal) %>%
  sample_n(sampsize)
# Calculate REcall
eqtl_nf_tota_causal <- length(which(eqtl_mpra_prc_nf_df$causal, TRUE))
eqtl_nf_mpra_pos <- length(which(eqtl_mpra_prc_nf_df$emVar_any, TRUE))

# Write emVar annotated GTEx DF to file
write.table(gtex_mpra_paired_filter_emvars, file='~/Dropbox (JAX)/Variant_Effects/ukbb_gtex_mpra/final_datasets/gtex_mpra_paired_filtered_emvar.txt',
            sep = '\t', row.names = FALSE, col.names = TRUE)
# Write GTEx PRC DF to file for checking in Python
write.table(eqtl_mpra_prc_df_2, file='~/Dropbox (JAX)/Variant_Effects/ukbb_gtex_mpra/final_datasets/gtex.paired.emvar.prc.txt',
            sep = '\t', row.names = FALSE, col.names = TRUE)
# Open Traits RDS
traits_mpra_paired_final20230117 <- readRDS("~/Dropbox (JAX)/for_john/empirical_ukbb_gtex_data/traits_mpra_paired_final20230117.rds")
# Filter for only best library
traits_mpra_paired_filter <- traits_mpra_paired_final20230117 %>% dplyr::filter(best_library == 1,
                                                                                type != "other_test")
# Add emVar and other columns
traits_mpra_paired_filter_emvars <- def_emVars(traits_mpra_paired_filter)
# Complex traits
set.seed(123)
traits_mpra_prc_df <- traits_mpra_paired_filter_emvars %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(pip = max(pip, na.rm = T)) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  filter(row_number() == 1) %>%
  ungroup() %>% 
  dplyr::mutate(causal = case_when(pip > 0.9 & (pchisq(chisq_marginal, 1, log.p = TRUE, lower.tail = F) / -log(10) > -log10(5 * 10^-8)) & (type %in% c("CS", "PIP10")) ~ TRUE,
                                   type == "CS" & pip < 0.01 ~ FALSE,
                                   TRUE ~ NA)) %>% 
  dplyr::filter(!is.na(causal)) %>% 
  dplyr::filter(! consequence %in% c("synonymous","missense","LoF")) 

sampsize <- min(table(traits_mpra_prc_df$causal))
traits_mpra_prc_df <- traits_mpra_prc_df %>%
  group_by(causal) %>%
  sample_n(sampsize)

# Write emVar annotated GTEx DF to file
write.table(traits_mpra_paired_filter_emvars, file='~/Dropbox (JAX)/Variant_Effects/ukbb_gtex_mpra/final_datasets/traits_mpra_paired_filtered_emvar.txt',
            sep = '\t', row.names = FALSE, col.names = TRUE)
# 

