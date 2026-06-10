# Data preparation utilities for INFUSE
# Source this file at the top of Data_prepare.Rmd:
#   source(here::here("code/data_prep_utils.R"))

library(data.table)
library(dplyr)

# ── Step 1: GWAS Harmonization ────────────────────────────────────────────────

#' Read, rename columns, QC-filter, and intersect one ancestry's summary stats
#' with the 1KG reference panel BIM file.
#' Expects GWAS files with columns BP/A1/A2/FRQ (regenie / SAIGE output).
read_sumstats <- function(file, ref, maf_thresh = 0.01) {
  dt <- fread(file)
  colnames(dt)[colnames(dt) == "BP"]  <- "POS"
  colnames(dt)[colnames(dt) == "A1"]  <- "ALT"
  colnames(dt)[colnames(dt) == "A2"]  <- "REF"
  colnames(dt)[colnames(dt) == "FRQ"] <- "MAF"
  dt <- GWAS_QC(dt, maf_thresh)
  dt %>%
    filter(CHR %in% 1:22, CHR_POS %in% ref$CHR_POS) %>%
    arrange(CHR, POS) %>%
    select(-CHISQ)
}

# ── Step 2: Align GWAS SNP list with 1KG LD reference ────────────────────────

#' Run plink2 to subset and flip a 1KG bfile to match the GWAS SNP list
#' for one ancestry (EUR or AFR).
align_ld_ref <- function(ancestry, trait,
                          kg_base  = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall",
                          ss_dir   = "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summstats_to_finemap",
                          out_base = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8",
                          plink2   = "~/gwas_software/plink2.0/plink2") {
  anc      <- tolower(ancestry)
  bfile    <- file.path(kg_base,  paste0("1kg_", anc, "_maf_unrelated_biallelic_nodup_extracted"))
  snplist  <- file.path(ss_dir,   paste0(trait, "_GRCh38_combined_format_noMHC_to_finemap_snplist"))
  ref_flip <- file.path(ss_dir,   paste0(trait, "_GRCh38_combined_format_noMHC_to_finemap_ref_flip"))
  outfile  <- file.path(out_base, paste0("1kg_", anc, "_maf_unrelated_biallelic_nodup_extracted_", trait, "_flipped"))
  system(paste0(plink2, " --bfile ", bfile, " --extract ", snplist,
                " --ref-allele ", ref_flip, " 2 1 --make-bed --out ", outfile),
         intern = TRUE)
  message("Aligned LD reference for ", ancestry, " (", trait, ")")
}

# ── Step 3: LD Matrix Construction ───────────────────────────────────────────

#' Impute missing values and scale each SNP column to zero mean / unit variance
preprocess_data <- function(data_matrix) {
  apply(data_matrix, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    scale(x)
  })
}

#' Compute symmetric LD correlation matrix and write to file
compute_cov <- function(data_matrix, file_name) {
  cov_mat <- cov2cor(crossprod(data_matrix))
  cov_mat <- (cov_mat + t(cov_mat)) / 2
  fwrite(cov_mat, file_name, sep = " ")
}

#' Subset a genomic region from a 1KG BIM file, run plink2 to extract genotypes,
#' identify monomorphic SNPs, and write the LD correlation matrix.
build_ld_matrix <- function(bim_ref, chr_num, start_pos, end_pos,
                             bfile, out_dir, locus_name = "loci_example",
                             plink2 = "~/gwas_software/plink2.0/plink2") {
  loci_snp     <- bim_ref[bim_ref$V4 >= start_pos & bim_ref$V4 <= end_pos &
                           bim_ref$V1 == chr_num, ]$V2
  dir.create(out_dir, showWarnings = FALSE)
  snplist_file <- file.path(out_dir, locus_name)
  write.table(loci_snp, snplist_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

  system(paste(plink2, "--bfile", bfile, "--extract", snplist_file,
               "--allow-no-sex --make-bed --out", snplist_file), intern = TRUE)

  plink_geno <- as(snpStats::read.plink(paste0(snplist_file, ".bed"))$genotypes, "numeric")
  mono_snps  <- names(which(apply(plink_geno, 2, function(x) length(unique(x)) == 1)))
  write.table(mono_snps, paste0(snplist_file, "_snp_to_exclude"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  mat <- preprocess_data(as.matrix(plink_geno))
  compute_cov(mat, paste0(snplist_file, ".ld"))
  invisible(snplist_file)
}

# ── Step 4: LD Mismatch Detection ────────────────────────────────────────────

#' Iteratively remove SNPs with LD-inconsistent Z-scores using kriging_rss.
#' Returns a list with the cleaned sumstats, final LD matrix, and SNP list.
filter_ld_mismatch <- function(geno_dir, locus_name, sumstats) {
  snplist_file <- file.path(geno_dir, locus_name)
  exclude_file <- paste0(snplist_file, "_snp_to_exclude")

  snplist <- fread(snplist_file, header = FALSE) %>% pull(V1)
  if (file.info(exclude_file)$size > 0)
    snplist <- setdiff(snplist, fread(exclude_file, header = FALSE) %>% pull(V1))

  idx_to_remove <- c()
  while (!identical(idx_to_remove, integer(0))) {
    sumstats_sub <- sumstats %>% filter(SNP %in% snplist)
    plink_geno   <- as(snpStats::read.plink(paste0(snplist_file, ".bed"))$genotypes, "numeric")
    mat          <- preprocess_data(as.matrix(plink_geno[, colnames(plink_geno) %in% snplist]))
    cov_mat      <- (cov2cor(crossprod(mat)) + t(cov2cor(crossprod(mat)))) / 2
    diag         <- susieR::kriging_rss(sumstats_sub$Z, cov_mat, n = median(sumstats_sub$N))
    print(diag$plot)
    idx_to_remove <- diag$plot$plot_env$idx
    if (!identical(idx_to_remove, integer(0))) snplist <- snplist[-idx_to_remove]
  }
  list(sumstats = sumstats_sub, cov = cov_mat, snplist = snplist)
}
