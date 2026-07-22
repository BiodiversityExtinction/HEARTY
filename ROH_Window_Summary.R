#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript ROH_Window_Summary.R --mode <single|batch> [options]",
  "",
  "Modes:",
  "  single: one window-summary txt -> one PDF",
  "  batch : sample list -> summary TSV + combined PDF",
  "",
  "Common options:",
  "  --min-length-mb 20",
  "  --max-scaffolds 40",
  "  --fill-cap 0.0005",
  "  --fill-cap-scale 1.0",
  "  --bin-size-bp 500000",
  "  --fallback-resolution 100",
  "  --roh-threshold 0.0001",
  "  --roh-threshold-mode absolute",
  "",
  "Single-mode required:",
  "  --input <path/to/sample.windows.txt.gz>",
  "  --output-pdf <path/to/output.pdf>",
  "",
  "Batch-mode required:",
  "  --sample-list <path/to/sample_list.tsv>",
  "  --summary-out <path/to/summary.tsv>",
  "  --combined-pdf <path/to/combined.pdf>",
  sep = "\n"
)

parse_named_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s\n\n%s", key, usage), call. = FALSE)
    if (i == length(x)) stop(sprintf("Missing value for %s\n\n%s", key, usage), call. = FALSE)
    val <- x[[i + 1L]]
    if (startsWith(val, "--")) stop(sprintf("Missing value for %s\n\n%s", key, usage), call. = FALSE)
    out[[sub("^--", "", key)]] <- val
    i <- i + 2L
  }
  out
}

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  cat(usage, "\n")
  quit(status = 0)
}

opt <- parse_named_args(args)
get_opt <- function(name, default = NULL, required = FALSE) {
  if (!is.null(opt[[name]])) return(opt[[name]])
  if (required) stop(sprintf("Missing required option --%s\n\n%s", name, usage), call. = FALSE)
  default
}

as_num <- function(x, nm) {
  v <- suppressWarnings(as.numeric(x))
  if (!is.finite(v)) stop(sprintf("--%s must be numeric", nm), call. = FALSE)
  v
}
as_int <- function(x, nm) {
  v <- suppressWarnings(as.integer(x))
  if (!is.finite(v)) stop(sprintf("--%s must be integer", nm), call. = FALSE)
  v
}

parse_roh_txt <- function(path) {
  if (!file.exists(path)) stop(sprintf("Input window-summary file not found: %s", path), call. = FALSE)
  lines <- readr::read_lines(path)
  df <- tibble(line = lines) %>%
    tidyr::extract(
      line,
      into = c("chrom", "start", "end", "total_sites", "het_sites", "het_label", "prop"),
      regex = "^\\s*(\\S+)\\s+(\\d+)-(\\d+):\\s+(\\d+)\\s+sites,\\s+(\\d+)\\s+sites\\s+with\\s+'(HET(?:_trv)?)',\\s+Proportion:\\s+([0-9.]+)\\s*$",
      convert = TRUE
    ) %>%
    filter(!is.na(chrom), !is.na(start), !is.na(end), !is.na(het_label), !is.na(prop))
  if (nrow(df) == 0) stop(sprintf("No parsable rows in %s", path), call. = FALSE)
  df
}

describe_dataset_type <- function(df, path) {
  labels <- unique(stats::na.omit(df$het_label))
  if (length(labels) == 1) {
    if (identical(labels[[1]], "HET_trv")) return("transversion-only")
    if (identical(labels[[1]], "HET")) return("all-sites")
  }

  lower_path <- tolower(path)
  if (grepl("windows_trv", lower_path, fixed = TRUE)) return("transversion-only")
  if (grepl("windows", lower_path, fixed = TRUE)) return("all-sites")
  "unknown"
}

resolve_fill_cap <- function(df_list, fill_cap_raw, fill_cap_scale) {
  if (tolower(fill_cap_raw) != "auto") {
    return(as_num(fill_cap_raw, "fill-cap"))
  }
  all_props <- unlist(lapply(df_list, function(d) d$prop), use.names = FALSE)
  cap <- mean(all_props, na.rm = TRUE) * fill_cap_scale
  if (!is.finite(cap) || cap <= 0) {
    stop("Automatic fill-cap could not be calculated from the input windows", call. = FALSE)
  }
  cap
}

apply_fill_cap <- function(df, fill_cap) {
  df %>% mutate(prop_cap = pmin(prop, fill_cap))
}

parse_sample_list <- function(path) {
  if (!file.exists(path)) stop(sprintf("Sample list not found: %s", path), call. = FALSE)
  lines <- readr::read_lines(path)
  out <- list()
  for (i in seq_along(lines)) {
    ln <- stringr::str_trim(lines[[i]])
    if (ln == "" || stringr::str_starts(ln, "#")) next
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2) stop(sprintf("%s:%d expected at least 2 tab-separated columns", path, i), call. = FALSE)
    s <- stringr::str_trim(parts[[1]])
    f <- stringr::str_trim(parts[[2]])
    if (i == 1 && tolower(s) %in% c("sample", "sample_id", "sample_name") && grepl("path|file|txt", tolower(f))) next
    if (s == "" || f == "") stop(sprintf("%s:%d has empty sample or path", path, i), call. = FALSE)
    out[[length(out) + 1]] <- tibble(sample = s, file = f)
  }
  if (length(out) == 0) stop(sprintf("No samples parsed from %s", path), call. = FALSE)
  bind_rows(out)
}

calc_roh_segments <- function(df_windows, flag_col) {
  d <- df_windows %>%
    arrange(chrom, start, end) %>%
    filter(.data[[flag_col]]) %>%
    mutate(
      new_seg = dplyr::row_number() == 1L | chrom != dplyr::lag(chrom) | start > (dplyr::lag(end) + 1L),
      seg_id = cumsum(new_seg)
    )
  if (nrow(d) == 0) {
    return(tibble(chrom = character(), seg_id = integer(), seg_start = integer(), seg_end = integer(), len_bp = integer()))
  }
  d %>%
    group_by(chrom, seg_id) %>%
    summarise(seg_start = min(start), seg_end = max(end), len_bp = seg_end - seg_start + 1L, .groups = "drop")
}

calc_froh_bins <- function(seg_df, genome_size_bp) {
  bins <- c(1e5, 1e6, 5e6)
  labels <- c("gt100kb", "gt1Mb", "gt5Mb")
  roh_bp <- vapply(bins, function(b) sum(seg_df$len_bp[seg_df$len_bp > b], na.rm = TRUE), numeric(1))
  tibble(metric = labels, roh_bp = roh_bp, froh = roh_bp / genome_size_bp)
}

calc_het_stats <- function(df) {
  values <- df$prop[is.finite(df$prop)]
  if (length(values) == 0) {
    return(tibble(mean_het = NA_real_, median_het = NA_real_))
  }
  tibble(
    mean_het = mean(values),
    median_het = median(values)
  )
}

calc_non_roh_stats <- function(df_windows, flag_col) {
  non_roh <- df_windows %>% filter(!.data[[flag_col]])
  stats <- calc_het_stats(non_roh)
  tibble(
    n_windows = nrow(non_roh),
    mean_het = stats$mean_het,
    median_het = stats$median_het
  )
}

build_window_flags <- function(df, roh_threshold) {
  df %>%
    arrange(chrom, start, end) %>%
    group_by(chrom) %>%
    mutate(
      roh_raw = prop <= roh_threshold,
      roh_bridged = roh_raw | (!roh_raw & lag(roh_raw, default = FALSE) & lead(roh_raw, default = FALSE))
    ) %>%
    ungroup()
}

resolve_roh_threshold <- function(df, threshold_value, threshold_mode) {
  if (threshold_mode == "absolute") {
    return(threshold_value)
  }
  basis <- if (threshold_mode == "mean") {
    mean(df$prop, na.rm = TRUE)
  } else {
    median(df$prop, na.rm = TRUE)
  }
  if (!is.finite(basis)) {
    stop(sprintf("Could not calculate %s heterozygosity for relative ROH threshold", threshold_mode), call. = FALSE)
  }
  threshold_value * basis
}

format_roh_caption <- function(threshold_input, threshold_mode, threshold_used = NULL) {
  if (threshold_mode == "absolute") {
    return(sprintf("Black ticks: raw bins <= %.6f | Orange ticks: bridged ROH (single-bin gap filled if ROH on both sides)", threshold_input))
  }
  if (!is.null(threshold_used)) {
    return(sprintf(
      "Black ticks: raw bins <= %.6f (%s x %s window heterozygosity) | Orange ticks: bridged ROH (single-bin gap filled if ROH on both sides)",
      threshold_used,
      threshold_input,
      threshold_mode
    ))
  }
  sprintf(
    "Black ticks: raw bins <= sample-specific threshold (%s x %s window heterozygosity) | Orange ticks: bridged ROH (single-bin gap filled if ROH on both sides)",
    threshold_input,
    threshold_mode
  )
}

build_single_plot_df <- function(df, min_length_mb, max_scaffolds, bin_size_bp, roh_threshold) {
  sizes_all <- df %>% group_by(chrom) %>% summarise(size_bp = max(end, na.rm = TRUE), .groups = "drop")
  sizes_sel <- sizes_all %>% filter(size_bp >= min_length_mb * 1e6) %>% arrange(desc(size_bp)) %>% slice_head(n = max_scaffolds)
  if (nrow(sizes_sel) == 0) stop("No scaffolds pass min-length threshold", call. = FALSE)
  d <- df %>% filter(chrom %in% sizes_sel$chrom)
  chrom_levels <- sizes_sel %>% arrange(desc(size_bp)) %>% pull(chrom)
  bin_size <- if (!is.na(bin_size_bp)) bin_size_bp else {
    b <- median(d$end - d$start, na.rm = TRUE)
    if (!is.finite(b) || b <= 0) 1e5 else b
  }

  d_bins <- d %>%
    mutate(bin_id = as.integer(floor(start / bin_size))) %>%
    group_by(chrom, bin_id) %>%
    summarise(prop = mean(prop, na.rm = TRUE), prop_cap = mean(prop_cap, na.rm = TRUE), roh_hit = any(prop <= roh_threshold, na.rm = TRUE), .groups = "drop")

  full_grid <- sizes_sel %>%
    rowwise() %>%
    mutate(bins = list(tibble(bin_id = 0:floor(size_bp / bin_size)))) %>%
    unnest(bins) %>%
    mutate(start_Mb = bin_id * bin_size / 1e6, end_Mb = (bin_id + 1) * bin_size / 1e6) %>%
    ungroup() %>%
    mutate(chrom = factor(chrom, levels = chrom_levels, ordered = TRUE))

  pdat <- full_grid %>%
    left_join(d_bins, by = c("chrom", "bin_id")) %>%
    group_by(chrom) %>%
    arrange(bin_id, .by_group = TRUE) %>%
    mutate(
      roh_hit_raw = if_else(is.na(roh_hit), FALSE, roh_hit),
      roh_hit_bridged = roh_hit_raw | (!roh_hit_raw & lag(roh_hit_raw, default = FALSE) & lead(roh_hit_raw, default = FALSE))
    ) %>%
    ungroup()

  list(plot_df = pdat, scaffolds_plotted = nrow(sizes_sel), rows_used = nrow(d), bin_size = bin_size)
}

plot_single <- function(plot_df, sample_label, fill_cap, caption_text, fallback_resolution, out_pdf) {
  n_panels <- n_distinct(plot_df$chrom)
  p <- ggplot(plot_df) +
    geom_rect(aes(xmin = start_Mb, xmax = end_Mb, ymin = 0, ymax = 0.82, fill = prop_cap)) +
    geom_segment(data = dplyr::filter(plot_df, roh_hit_raw), aes(x = start_Mb, xend = end_Mb, y = 0.92, yend = 0.92), inherit.aes = FALSE, linewidth = 2.0, color = "black") +
    geom_segment(data = dplyr::filter(plot_df, roh_hit_bridged), aes(x = start_Mb, xend = end_Mb, y = 1.04, yend = 1.04), inherit.aes = FALSE, linewidth = 1.8, color = "#E69F00") +
    facet_grid(rows = vars(chrom), switch = "y") +
    scale_fill_gradient(name = "HET proportion", limits = c(0, fill_cap), low = "red", high = "blue", na.value = "grey90") +
    labs(x = "Position (Mb)", y = NULL, title = sample_label,
         caption = caption_text) +
    coord_cartesian(ylim = c(0, 1.10), expand = FALSE) +
    theme_minimal(base_size = 11) +
    theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid = element_blank(), strip.placement = "outside", strip.text.y = element_blank(),
          legend.title = element_text(size = 13), legend.text = element_text(size = 11), legend.key.height = grid::unit(0.9, "cm"),
          plot.caption = element_text(size = 10, hjust = 0), plot.title = element_text(face = "bold", size = 14))

  pdf_h <- 0.55 * max(1, n_panels)
  ggsave(out_pdf, p, width = 12, height = pdf_h, units = "in", device = cairo_pdf,
         fallback_resolution = fallback_resolution, limitsize = FALSE)
}

build_batch_plot_df <- function(df, sample_id, min_length_mb, max_scaffolds, bin_size_bp, roh_threshold) {
  sizes_all <- df %>% group_by(chrom) %>% summarise(size_bp = max(end, na.rm = TRUE), .groups = "drop")
  sizes_sel <- sizes_all %>% filter(size_bp >= min_length_mb * 1e6) %>% arrange(desc(size_bp)) %>% slice_head(n = max_scaffolds)
  if (nrow(sizes_sel) == 0) return(list(plot_df = tibble(), n_chrom = 0, rows_used = 0, bin_size = NA_real_))
  d <- df %>% filter(chrom %in% sizes_sel$chrom)
  chrom_levels <- sizes_sel %>% arrange(desc(size_bp)) %>% pull(chrom)
  chrom_rank <- tibble(chrom = chrom_levels, chrom_idx = rev(seq_along(chrom_levels)))
  bin_size <- if (!is.na(bin_size_bp)) bin_size_bp else {
    b <- median(d$end - d$start, na.rm = TRUE)
    if (!is.finite(b) || b <= 0) 1e5 else b
  }

  d_bins <- d %>%
    mutate(bin_id = as.integer(floor(start / bin_size))) %>%
    group_by(chrom, bin_id) %>%
    summarise(prop = mean(prop, na.rm = TRUE), prop_cap = mean(prop_cap, na.rm = TRUE), roh_hit = any(prop <= roh_threshold, na.rm = TRUE), .groups = "drop")

  full_grid <- sizes_sel %>%
    rowwise() %>%
    mutate(bins = list(tibble(bin_id = 0:floor(size_bp / bin_size)))) %>%
    unnest(bins) %>%
    mutate(start_Mb = bin_id * bin_size / 1e6, end_Mb = (bin_id + 1) * bin_size / 1e6) %>%
    ungroup()

  pdat <- full_grid %>%
    left_join(d_bins, by = c("chrom", "bin_id")) %>%
    left_join(chrom_rank, by = "chrom") %>%
    mutate(sample = sample_id, chrom = as.character(chrom)) %>%
    group_by(sample, chrom) %>%
    arrange(bin_id, .by_group = TRUE) %>%
    mutate(
      roh_hit_raw = if_else(is.na(roh_hit), FALSE, roh_hit),
      roh_hit_bridged = roh_hit_raw | (!roh_hit_raw & lag(roh_hit_raw, default = FALSE) & lead(roh_hit_raw, default = FALSE)),
      y_rect_min = chrom_idx - 0.35,
      y_rect_max = chrom_idx + 0.25,
      y_raw_tick = chrom_idx + 0.37,
      y_bridged_tick = chrom_idx + 0.52
    ) %>%
    ungroup()

  list(plot_df = pdat, n_chrom = nrow(sizes_sel), rows_used = nrow(d), bin_size = bin_size)
}

plot_batch <- function(plot_df_all, caption_text, fallback_resolution, out_pdf, sample_order) {
  plot_df_all$sample <- factor(plot_df_all$sample, levels = sample_order, ordered = TRUE)
  sample_counts <- plot_df_all %>% distinct(sample, chrom) %>% count(sample, name = "n_panels") %>%
    mutate(sample_order = as.integer(sample), row_id = ((sample_order - 1) %/% 2) + 1)

  p <- ggplot(plot_df_all) +
    geom_rect(aes(xmin = start_Mb, xmax = end_Mb, ymin = y_rect_min, ymax = y_rect_max, fill = prop_cap)) +
    geom_segment(data = dplyr::filter(plot_df_all, roh_hit_raw), aes(x = start_Mb, xend = end_Mb, y = y_raw_tick, yend = y_raw_tick), inherit.aes = FALSE, linewidth = 2.0, color = "black") +
    geom_segment(data = dplyr::filter(plot_df_all, roh_hit_bridged), aes(x = start_Mb, xend = end_Mb, y = y_bridged_tick, yend = y_bridged_tick), inherit.aes = FALSE, linewidth = 1.8, color = "#E69F00") +
    facet_wrap(~sample, ncol = 2, scales = "free") +
    scale_fill_gradient(name = "HET proportion", low = "red", high = "blue", na.value = "grey90") +
    labs(x = "Position (Mb)", y = NULL, caption = caption_text) +
    coord_cartesian(expand = FALSE) +
    theme_minimal(base_size = 11) +
    theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 11), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank(), strip.placement = "outside", strip.text = element_text(face = "bold", size = 12),
          panel.spacing = grid::unit(10, "mm"), legend.title = element_text(size = 13), legend.text = element_text(size = 11),
          legend.key.height = grid::unit(0.9, "cm"), plot.caption = element_text(size = 10, hjust = 0))

  row_heights <- sample_counts %>% group_by(row_id) %>% summarise(row_h = max(2.8, 0.17 * n_panels), .groups = "drop")
  pdf_h <- sum(row_heights$row_h) + 1.5
  ggsave(out_pdf, p, width = 13, height = pdf_h, units = "in", device = cairo_pdf,
         fallback_resolution = fallback_resolution, limitsize = FALSE)
}

mode <- get_opt("mode", required = TRUE)
if (!mode %in% c("single", "batch")) stop("--mode must be either 'single' or 'batch'", call. = FALSE)

min_length_mb <- as_num(get_opt("min-length-mb", "20"), "min-length-mb")
max_scaffolds <- as_int(get_opt("max-scaffolds", "40"), "max-scaffolds")
fill_cap_raw <- get_opt("fill-cap", "0.0005")
fill_cap_scale <- as_num(get_opt("fill-cap-scale", "1.0"), "fill-cap-scale")
bin_size_bp <- as_num(get_opt("bin-size-bp", "500000"), "bin-size-bp")
fallback_resolution <- as_int(get_opt("fallback-resolution", "100"), "fallback-resolution")
roh_threshold <- as_num(get_opt("roh-threshold", "0.0001"), "roh-threshold")
roh_threshold_mode <- tolower(get_opt("roh-threshold-mode", "absolute"))

if (min_length_mb <= 0 || max_scaffolds <= 0 || fill_cap_scale <= 0 || bin_size_bp <= 0 || fallback_resolution <= 0 || roh_threshold < 0) {
  stop("Invalid parameter values", call. = FALSE)
}
if (!roh_threshold_mode %in% c("absolute", "mean", "median")) {
  stop("--roh-threshold-mode must be one of: absolute, mean, median", call. = FALSE)
}

WINDOW_BP_FOR_FROH <- 100000

if (identical(mode, "single")) {
  input <- get_opt("input", required = TRUE)
  out_pdf <- get_opt("output-pdf", required = TRUE)
  sample_label <- basename(input)

  df <- parse_roh_txt(input)
  dataset_type <- describe_dataset_type(df, input)
  fill_cap <- resolve_fill_cap(list(df), fill_cap_raw, fill_cap_scale)
  df <- apply_fill_cap(df, fill_cap)
  het_stats_all <- calc_het_stats(df)
  roh_threshold_used <- resolve_roh_threshold(df, roh_threshold, roh_threshold_mode)
  roh_caption <- format_roh_caption(roh_threshold, roh_threshold_mode, roh_threshold_used)
  genome_bp <- nrow(df) * WINDOW_BP_FOR_FROH
  flags <- build_window_flags(df, roh_threshold_used)
  non_roh_raw <- calc_non_roh_stats(flags, "roh_raw")
  non_roh_bridged <- calc_non_roh_stats(flags, "roh_bridged")
  seg_raw <- calc_roh_segments(flags, "roh_raw")
  seg_br <- calc_roh_segments(flags, "roh_bridged")
  froh_raw <- calc_froh_bins(seg_raw, genome_bp)
  froh_br <- calc_froh_bins(seg_br, genome_bp)

  built <- build_single_plot_df(df, min_length_mb, max_scaffolds, bin_size_bp, roh_threshold_used)
  plotted_df <- df %>% filter(chrom %in% unique(as.character(built$plot_df$chrom)))
  het_stats_plotted <- calc_het_stats(plotted_df)
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  plot_single(built$plot_df, sprintf("%s [%s]", sample_label, dataset_type), fill_cap, roh_caption, fallback_resolution, out_pdf)

  cat(sprintf("Wrote PDF: %s\n", out_pdf))
  cat(sprintf("Dataset type: %s\n", dataset_type))
  cat(sprintf("HET label parsed: %s\n", paste(unique(df$het_label), collapse = ", ")))
  cat(sprintf("Input rows parsed: %d\n", nrow(df)))
  cat(sprintf("Mean heterozygosity (all parsed windows): %.7f\n", het_stats_all$mean_het))
  cat(sprintf("Median heterozygosity (all parsed windows): %.7f\n", het_stats_all$median_het))
  cat(sprintf("Rows used for plotted scaffolds: %d\n", built$rows_used))
  cat(sprintf("Mean heterozygosity (plotted scaffolds): %.7f\n", het_stats_plotted$mean_het))
  cat(sprintf("Median heterozygosity (plotted scaffolds): %.7f\n", het_stats_plotted$median_het))
  cat(sprintf("Scaffolds plotted: %d (min length %.2f Mb; max panels %d)\n", built$scaffolds_plotted, min_length_mb, max_scaffolds))
  cat(sprintf("Bin size used: %.0f bp\n", built$bin_size))
  cat(sprintf("PDF fallback resolution: %d dpi\n", fallback_resolution))
  cat(sprintf("ROH threshold mode: %s\n", roh_threshold_mode))
  cat(sprintf("ROH threshold input: %.6f\n", roh_threshold))
  cat(sprintf("ROH threshold used: %.6f\n", roh_threshold_used))
  cat(sprintf("Non-ROH windows (raw definition): %d\n", non_roh_raw$n_windows))
  cat(sprintf("Mean heterozygosity (non-ROH, raw definition): %.7f\n", non_roh_raw$mean_het))
  cat(sprintf("Median heterozygosity (non-ROH, raw definition): %.7f\n", non_roh_raw$median_het))
  cat(sprintf("Non-ROH windows (bridged definition): %d\n", non_roh_bridged$n_windows))
  cat(sprintf("Mean heterozygosity (non-ROH, bridged definition): %.7f\n", non_roh_bridged$mean_het))
  cat(sprintf("Median heterozygosity (non-ROH, bridged definition): %.7f\n", non_roh_bridged$median_het))
  cat(sprintf("Raw ROH bins flagged: %d\n", sum(built$plot_df$roh_hit_raw, na.rm = TRUE)))
  cat(sprintf("Bridged ROH bins flagged: %d\n", sum(built$plot_df$roh_hit_bridged, na.rm = TRUE)))
  cat(sprintf("Genome size used for FROH: %.0f bp (n_windows=%d x %d bp)\n", genome_bp, nrow(df), WINDOW_BP_FOR_FROH))

  cat("\nFROH (unbridged/raw; ROH segments longer than cutoff)\n")
  for (i in seq_len(nrow(froh_raw))) cat(sprintf("  >%s: %.6f (ROH bp=%.0f)\n", sub("gt", "", froh_raw$metric[i]), froh_raw$froh[i], froh_raw$roh_bp[i]))
  cat("FROH (bridged; ROH segments longer than cutoff)\n")
  for (i in seq_len(nrow(froh_br))) cat(sprintf("  >%s: %.6f (ROH bp=%.0f)\n", sub("gt", "", froh_br$metric[i]), froh_br$froh[i], froh_br$roh_bp[i]))

} else {
  sample_list <- get_opt("sample-list", required = TRUE)
  summary_out <- get_opt("summary-out", required = TRUE)
  combined_pdf <- get_opt("combined-pdf", required = TRUE)

  samples <- parse_sample_list(sample_list)
  samples$sample <- as.character(samples$sample)

  dir.create(dirname(summary_out), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(combined_pdf), recursive = TRUE, showWarnings = FALSE)

  parsed_rows <- list()
  summary_rows <- list()
  plot_rows <- list()

  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample[[i]]
    path <- samples$file[[i]]
    df <- parse_roh_txt(path)
    parsed_rows[[length(parsed_rows) + 1]] <- df
  }

  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample[[i]]
    path <- samples$file[[i]]
    fill_cap <- resolve_fill_cap(list(parsed_rows[[i]]), fill_cap_raw, fill_cap_scale)
    df <- apply_fill_cap(parsed_rows[[i]], fill_cap)
    dataset_type <- describe_dataset_type(df, path)
    het_stats_all <- calc_het_stats(df)
    roh_threshold_used <- resolve_roh_threshold(df, roh_threshold, roh_threshold_mode)

    genome_bp <- nrow(df) * WINDOW_BP_FOR_FROH
    flags <- build_window_flags(df, roh_threshold_used)
    non_roh_raw <- calc_non_roh_stats(flags, "roh_raw")
    non_roh_bridged <- calc_non_roh_stats(flags, "roh_bridged")
    seg_raw <- calc_roh_segments(flags, "roh_raw")
    seg_br <- calc_roh_segments(flags, "roh_bridged")
    froh_raw <- calc_froh_bins(seg_raw, genome_bp)
    froh_br <- calc_froh_bins(seg_br, genome_bp)

    raw_map <- setNames(froh_raw$froh, froh_raw$metric)
    raw_bp_map <- setNames(froh_raw$roh_bp, froh_raw$metric)
    br_map <- setNames(froh_br$froh, froh_br$metric)
    br_bp_map <- setNames(froh_br$roh_bp, froh_br$metric)

    built <- build_batch_plot_df(df, sid, min_length_mb, max_scaffolds, bin_size_bp, roh_threshold_used)
    plotted_df <- df %>% filter(chrom %in% unique(built$plot_df$chrom))
    het_stats_plotted <- if (nrow(plotted_df) > 0) calc_het_stats(plotted_df) else tibble(mean_het = NA_real_, median_het = NA_real_)
    if (nrow(built$plot_df) > 0) plot_rows[[length(plot_rows) + 1]] <- built$plot_df

    summary_rows[[length(summary_rows) + 1]] <- tibble(
      sample = sid,
      roh_file = path,
      dataset_type = dataset_type,
      het_label = paste(unique(df$het_label), collapse = ","),
      fill_cap_used = fill_cap,
      n_windows = nrow(df),
      mean_heterozygosity_all_windows = het_stats_all$mean_het,
      median_heterozygosity_all_windows = het_stats_all$median_het,
      n_non_roh_raw_windows = non_roh_raw$n_windows,
      mean_heterozygosity_non_roh_raw_windows = non_roh_raw$mean_het,
      median_heterozygosity_non_roh_raw_windows = non_roh_raw$median_het,
      n_non_roh_bridged_windows = non_roh_bridged$n_windows,
      mean_heterozygosity_non_roh_bridged_windows = non_roh_bridged$mean_het,
      median_heterozygosity_non_roh_bridged_windows = non_roh_bridged$median_het,
      genome_bp_for_froh = genome_bp,
      roh_threshold_mode = roh_threshold_mode,
      roh_threshold_input = roh_threshold,
      roh_threshold = roh_threshold_used,
      n_raw_segments = nrow(seg_raw),
      n_bridged_segments = nrow(seg_br),
      raw_froh_gt100kb = raw_map[["gt100kb"]],
      raw_froh_gt1Mb = raw_map[["gt1Mb"]],
      raw_froh_gt5Mb = raw_map[["gt5Mb"]],
      raw_roh_bp_gt100kb = raw_bp_map[["gt100kb"]],
      raw_roh_bp_gt1Mb = raw_bp_map[["gt1Mb"]],
      raw_roh_bp_gt5Mb = raw_bp_map[["gt5Mb"]],
      bridged_froh_gt100kb = br_map[["gt100kb"]],
      bridged_froh_gt1Mb = br_map[["gt1Mb"]],
      bridged_froh_gt5Mb = br_map[["gt5Mb"]],
      bridged_roh_bp_gt100kb = br_bp_map[["gt100kb"]],
      bridged_roh_bp_gt1Mb = br_bp_map[["gt1Mb"]],
      bridged_roh_bp_gt5Mb = br_bp_map[["gt5Mb"]],
      scaffolds_plotted = built$n_chrom,
      rows_used_for_plot = built$rows_used,
      mean_heterozygosity_plotted_windows = het_stats_plotted$mean_het,
      median_heterozygosity_plotted_windows = het_stats_plotted$median_het,
      bin_size_bp_used = built$bin_size
    )
  }

  summary_df <- bind_rows(summary_rows)
  readr::write_tsv(summary_df, summary_out)
  if (length(plot_rows) == 0) stop("No samples had scaffolds passing min-length threshold", call. = FALSE)

  plot_df_all <- bind_rows(plot_rows)
  roh_caption <- format_roh_caption(roh_threshold, roh_threshold_mode)
  plot_batch(plot_df_all, roh_caption, fallback_resolution, combined_pdf, samples$sample)

  sample_counts <- plot_df_all %>% distinct(sample, chrom) %>% nrow()
  cat(sprintf("Wrote summary table: %s\n", summary_out))
  cat(sprintf("Wrote combined PDF: %s\n", combined_pdf))
  cat(sprintf("Samples processed: %d\n", nrow(summary_df)))
  cat(sprintf("Total plotted scaffold panels: %d\n", sample_counts))
  cat(sprintf("ROH threshold mode: %s\n", roh_threshold_mode))
  cat(sprintf("ROH threshold input: %.6f\n", roh_threshold))
  if (roh_threshold_mode != "absolute") {
    cat("ROH threshold used: sample-specific values written to the summary table\n")
  }
  cat(sprintf("Min scaffold length for plot: %.2f Mb\n", min_length_mb))
  cat(sprintf("Max scaffolds per sample: %d\n", max_scaffolds))
  if (tolower(fill_cap_raw) == "auto") {
    cat(sprintf("Fill cap mode: auto (sample-specific) with scale %.3f\n", fill_cap_scale))
  } else {
    cat(sprintf("Fill cap: %.6f\n", as_num(fill_cap_raw, "fill-cap")))
  }
  cat(sprintf("Bin size: %.0f bp\n", bin_size_bp))
}
