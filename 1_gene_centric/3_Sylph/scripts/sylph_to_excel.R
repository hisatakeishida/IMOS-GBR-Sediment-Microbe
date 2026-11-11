#!/usr/bin/env Rscript
# =============================================================================
# sylph_to_excel.R
# 
# Originally written by Julian Zaugg 
#
# Converts sylph-tax merged output tables (relative_abundance or
# sequence_abundance) into a multi-sheet Excel workbook with one sheet per
# taxonomy level (Phylum, Class, Order, Family, Genus, Species).
#
# Input files produced by:
#   sylph-tax merge $OUT_DIR/*.sylphmpa --column relative_abundance \
#     -o $OUT_DIR/merged_relative_abundance.tsv
#   sylph-tax merge $OUT_DIR/*.sylphmpa --column sequence_abundance \
#     -o $OUT_DIR/merged_sequence_abundance.tsv
#
# Usage (from R or Rscript):
#   source("sylph_to_excel.R")          # uses defaults below
#   # or
#   Rscript sylph_to_excel.R [relative_abundance.tsv] [sequence_abundance.tsv] [output.xlsx]
#
# Dependencies:
#   install.packages(c("dplyr", "tibble", "readr", "stringr", "purrr", "tidyr", "openxlsx"))
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(openxlsx)
})

setwd("<path to sylph results>")
# ---------------------------------------------------------------------------
# 0.  Configuration — edit these paths or pass them as command-line arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

rel_abund_file  <- if (length(args) >= 1) args[1] else "merged_relative_abundance.tsv"
seq_abund_file  <- if (length(args) >= 2) args[2] else "merged_sequence_abundance.tsv"
output_xlsx     <- if (length(args) >= 3) args[3] else "sylph_taxonomy_abundances.xlsx"

# Which input file(s) to process. Set to NA to skip one.
# Both are processed when both paths are supplied.
process_relative  <- !is.na(rel_abund_file)  && file.exists(rel_abund_file)
process_sequence  <- !is.na(seq_abund_file)  && file.exists(seq_abund_file)

if (!process_relative && !process_sequence) {
  stop("No input files found. Check rel_abund_file / seq_abund_file paths.")
}

# ---------------------------------------------------------------------------
# 1.  Helper: clean sample column names
#
#     sylph-tax merge can preserve full input read paths as sample headers, e.g.
#       /path/to/processed_reads/barcode10.fastq.gz
#     For output tables, reduce these to bare sample IDs, e.g.
#       barcode10
# ---------------------------------------------------------------------------
clean_sample_name <- function(x) {
  x <- basename(x)
  x <- stringr::str_remove(x, "\\.(fastq|fq|fasta|fa|fna)(\\.gz|\\.bz2|\\.xz)?$")
  x <- stringr::str_remove(x, "\\.sylphmpa$")
  x
}

clean_sample_column_names <- function(df) {
  sample_cols <- setdiff(names(df), "clade_name")
  cleaned_sample_cols <- clean_sample_name(sample_cols)
  
  if (anyDuplicated(cleaned_sample_cols)) {
    duplicates <- unique(cleaned_sample_cols[duplicated(cleaned_sample_cols)])
    stop(
      "Cleaning sample names created duplicate column names: ",
      paste(duplicates, collapse = ", "),
      ". Please make input filenames unique before removing extensions."
    )
  }
  
  names(df)[match(sample_cols, names(df))] <- cleaned_sample_cols
  df
}

# ---------------------------------------------------------------------------
# 1.  Taxonomy level definitions
#     prefix   : the letter prefix used in clade_name (e.g. "p__")
#     col_label: sheet name and first-column header in the output
# ---------------------------------------------------------------------------
TAX_LEVELS <- tibble::tibble(
  rank   = c("phylum", "class",  "order",  "family", "genus",  "species"),
  prefix = c("p__",    "c__",    "o__",    "f__",    "g__",    "s__"),
  sheet  = c("taxonomy_phylum", "taxonomy_class",  "taxonomy_order",
             "taxonomy_family", "taxonomy_genus",  "taxonomy_species")
)

# ---------------------------------------------------------------------------
# 2.  Helper: read a merged sylph-tax TSV
#
#     The merged file looks like:
#       clade_name  sample1  sample2  ...
#     Rows that contain higher-level clades are already separated by pipe (|)
#     or semicolon (;) depending on the sylph-tax version / database.
#     We handle both separators.
# ---------------------------------------------------------------------------
read_sylph_merged <- function(path) {
  message("Reading: ", path)
  df <- readr::read_tsv(path, comment = "#", show_col_types = FALSE)
  
  # Normalise separator: replace "|" with ";" so we can split uniformly
  df <- df %>%
    dplyr::mutate(clade_name = stringr::str_replace_all(.data$clade_name, "\\|", ";"))
  
  # Clean full FASTQ/FASTA paths in sample column headers to bare sample IDs
  # before downstream sheet generation.
  df <- clean_sample_column_names(df)
  
  df
}

# ---------------------------------------------------------------------------
# 3.  Helper: extract one taxonomy level from the merged table
#
#     Strategy:
#       - Keep only rows whose deepest token is at exactly this rank. This
#         avoids double-counting rows that belong to a deeper level (those
#         also contain the prefix string but their last token is deeper).
#       - Rows whose clade_name does NOT contain the rank prefix at all
#         represent lineages unresolved at this level. Their abundance is
#         aggregated under "Unassigned" for each sample.
#       - Build the full lineage string for the label column.
#       - Sum across duplicate lineages (safe dedup).
#       - Replace missing/NA abundances with 0.
# ---------------------------------------------------------------------------
extract_level <- function(df, rank_prefix, col_label) {
  sample_cols <- setdiff(names(df), "clade_name")
  rank_regex <- paste0("^", rank_prefix)
  empty_rank_regex <- paste0("^", rank_prefix, "$")
  
  # --- rows resolved to exactly this rank -----------------------------------
  df_resolved <- df %>%
    dplyr::mutate(
      tokens     = stringr::str_split(.data$clade_name, ";"),
      last_token = purrr::map_chr(.data$tokens, dplyr::last)
    ) %>%
    dplyr::filter(stringr::str_detect(.data$last_token, rank_regex))
  df_resolved <- df_resolved[, setdiff(names(df_resolved), c("tokens", "last_token")), drop = FALSE]
  
  # Rename clade_name to the dynamic taxonomy-level column.
  names(df_resolved)[names(df_resolved) == "clade_name"] <- col_label
  
  # --- rows that have NO representation at this rank (unassigned) -----------
  # These are rows whose deepest token is above this rank (e.g. a species row
  # contributes to genus even if the genus token is "g__" with no name given).
  # More practically: sylph-tax emits one row per strain/species, and higher
  # aggregates; if a genus token is literally "g__" (empty label after prefix)
  # we label it Unassigned.
  unassigned_rows <- df %>%
    dplyr::mutate(
      tokens     = stringr::str_split(.data$clade_name, ";"),
      last_token = purrr::map_chr(.data$tokens, dplyr::last),
      # Rows that ARE at this rank but have an empty/unknown label
      is_this_rank_empty = stringr::str_detect(.data$last_token, empty_rank_regex)
    ) %>%
    dplyr::filter(.data$is_this_rank_empty)
  unassigned_rows <- unassigned_rows[, setdiff(names(unassigned_rows), c("tokens", "last_token", "is_this_rank_empty")), drop = FALSE]
  
  if (nrow(unassigned_rows) > 0) {
    unassigned_summed <- unassigned_rows %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(sample_cols), function(x) sum(x, na.rm = TRUE)))
    
    # Add the dynamic taxonomy-level column.
    unassigned_summed[[col_label]] <- "Unassigned"
    unassigned_summed <- unassigned_summed %>%
      dplyr::select(dplyr::all_of(col_label), dplyr::everything())
    
    df_resolved <- dplyr::bind_rows(df_resolved, unassigned_summed)
  }
  
  # --- deduplicate and clean ------------------------------------------------
  df_level <- df_resolved %>%
    dplyr::group_by(.data[[col_label]]) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(sample_cols), function(x) sum(x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(sample_cols), function(x) tidyr::replace_na(x, 0)))
  
  # Sort: Unassigned last, everything else alphabetically
  is_unassigned <- df_level[[col_label]] == "Unassigned"
  df_level <- dplyr::bind_rows(
    df_level[!is_unassigned, ] %>% dplyr::arrange(.data[[col_label]]),
    df_level[is_unassigned, ]
  )
  
  df_level
}

# ---------------------------------------------------------------------------
# 4.  Helper: build and style one sheet in the workbook
# ---------------------------------------------------------------------------
add_sheet <- function(wb, sheet_name, df) {
  openxlsx::addWorksheet(wb, sheet_name)
  
  # Header style: bold, light blue background
  header_style <- openxlsx::createStyle(
    fontName       = "Arial",
    fontSize       = 11,
    textDecoration = "bold",
    fgFill         = "#DAEEF3",
    border         = "Bottom",
    borderColour   = "#4BACC6",
    wrapText       = FALSE,
    halign         = "center",
    valign         = "center"
  )
  
  # Data style: normal font
  data_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 10,
    halign   = "center",
    valign   = "center"
  )
  
  # Number style for abundance values
  number_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 10,
    numFmt   = "0.0000",
    halign   = "center",
    valign   = "center"
  )
  
  openxlsx::writeData(wb, sheet_name, df, startRow = 1, startCol = 1,
                      headerStyle = header_style, keepNA = FALSE, na.string = "0")
  
  n_rows <- nrow(df)
  n_cols <- ncol(df)
  
  if (n_rows > 0) {
    # Apply data style to all data cells
    openxlsx::addStyle(wb, sheet_name, data_style,
                       rows = 2:(n_rows + 1), cols = 1, gridExpand = TRUE)
    if (n_cols > 1) {
      openxlsx::addStyle(wb, sheet_name, number_style,
                         rows = 2:(n_rows + 1), cols = 2:n_cols, gridExpand = TRUE)
      
      # Colour-scale abundance values from white (0) to red (100).
      openxlsx::conditionalFormatting(
        wb,
        sheet = sheet_name,
        cols  = 2:n_cols,
        rows  = 2:(n_rows + 1),
        type  = "colourScale",
        style = c("#FFFFFF", "#FF0000"),
        rule  = c(0, 100)
      )
    }
  }
  
  # Auto-size the taxonomy label column; fixed width for sample columns
  openxlsx::setColWidths(wb, sheet_name, cols = 1,       widths = "auto")
  if (n_cols > 1) {
    openxlsx::setColWidths(wb, sheet_name, cols = 2:n_cols, widths = 14)
  }
  
  # Freeze the header row and first column
  openxlsx::freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
}

# ---------------------------------------------------------------------------
# 5.  Main processing function
# ---------------------------------------------------------------------------
process_file <- function(tsv_path, abundance_type) {
  df_raw <- read_sylph_merged(tsv_path)
  
  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 10)
  
  for (i in seq_len(nrow(TAX_LEVELS))) {
    rank_prefix <- TAX_LEVELS$prefix[i]
    sheet_name  <- TAX_LEVELS$sheet[i]
    
    message("  Processing level: ", sheet_name)
    
    df_level <- extract_level(df_raw, rank_prefix, sheet_name)
    
    if (nrow(df_level) == 0) {
      message("    No rows found for this level — adding empty sheet.")
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, data.frame(sheet_name = character(0)))
    } else {
      add_sheet(wb, sheet_name, df_level)
      message("    ", nrow(df_level), " lineages written.")
    }
  }
  
  wb
}

# ---------------------------------------------------------------------------
# 6.  Run
# ---------------------------------------------------------------------------
if (process_relative && process_sequence) {
  # Process both; write separate workbooks named after abundance type
  base   <- tools::file_path_sans_ext(output_xlsx)
  ext    <- tools::file_ext(output_xlsx)
  
  out_rel <- paste0(base, "_relative_abundance.", ext)
  out_seq <- paste0(base, "_sequence_abundance.", ext)
  
  wb_rel <- process_file(rel_abund_file, "relative_abundance")
  dir.create(dirname(out_rel), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb_rel, out_rel, overwrite = TRUE)
  message("Saved: ", out_rel)
  
  wb_seq <- process_file(seq_abund_file, "sequence_abundance")
  dir.create(dirname(out_seq), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb_seq, out_seq, overwrite = TRUE)
  message("Saved: ", out_seq)
  
} else if (process_relative) {
  wb <- process_file(rel_abund_file, "relative_abundance")
  dir.create(dirname(output_xlsx), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  message("Saved: ", output_xlsx)
  
} else {
  wb <- process_file(seq_abund_file, "sequence_abundance")
  dir.create(dirname(output_xlsx), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  message("Saved: ", output_xlsx)
}

message("Done.")
