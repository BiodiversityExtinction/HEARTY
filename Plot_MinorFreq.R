#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

show_help <- function() {
  cat(
    paste(
      "Plot_MinorFreq.R - Plot minor allele frequency overlays from HEARTY output",
      "",
      "Single-sample usage:",
      "  Rscript Plot_MinorFreq.R --minorfreq-file <file> --out <output_prefix> [--format pdf|svg|png] [--log-y] [--y-cap-after-x FLOAT]",
      "",
      "Batch usage:",
      "  Rscript Plot_MinorFreq.R --sample-list <list.txt> [--batch-out <combined_output>] [--format pdf|svg|png] [--log-y] [--y-cap-after-x FLOAT]",
      "",
      "Single-sample arguments:",
      "  --minorfreq-file   Path to HEARTY consolidated .minorfreq.txt file",
      "  --out              Output prefix, e.g. results/plots/sample01",
      "",
      "Batch arguments:",
      "  --sample-list      Text file with: column 1 = minorfreq file, column 2 = output prefix",
      "  --batch-out        Optional combined batch plot path",
      "",
      "Optional flags:",
      "  --format           Output format: pdf, svg, or png [default: pdf]",
      "  --log-y            Draw the y axis on a log10 scale",
      "  --y-cap-after-x    Set y max to 2x the highest count at freq >= cutoff",
      "",
      "Outputs:",
      "  <output_prefix>_minorfreq_overlay.<format>",
      "  <output_prefix>_minorfreq_summary.txt",
      "  <batch_output>     Combined faceted batch plot when --sample-list is used",
      sep = "\n"
    ),
    "\n"
  )
}

print_run_summary <- function(output_prefixes, fmt, log_y, y_cap_after_x, batch_output = NULL) {
  cat("Plot_MinorFreq.R completed\n")
  cat(sprintf("Output format: %s\n", fmt))
  cat(sprintf("Y axis mode: %s\n", if (log_y) "log10" else "linear"))
  if (!is.na(y_cap_after_x)) {
    cat(sprintf("Y cap rule: 2x highest count at freq >= %s\n", y_cap_after_x))
  } else {
    cat("Y cap rule: none\n")
  }
  cat("Generated files:\n")
  if (!is.null(batch_output)) {
    cat(sprintf("  %s\n", batch_output))
  }
  for (prefix in output_prefixes) {
    cat(sprintf("  %s_minorfreq_summary.txt\n", prefix))
  }
}

parse_named_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    }
    if (key %in% c("--help", "--log-y")) {
      out[[sub("^--", "", key)]] <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(x)) {
      stop(sprintf("Missing value for %s", key), call. = FALSE)
    }
    val <- x[[i + 1L]]
    if (startsWith(val, "--")) {
      stop(sprintf("Missing value for %s", key), call. = FALSE)
    }
    out[[sub("^--", "", key)]] <- val
    i <- i + 2L
  }
  out
}

if (length(args) == 0) {
  show_help()
  quit(status = 0)
}

opt <- parse_named_args(args)
if (!is.null(opt$help)) {
  show_help()
  quit(status = 0)
}

open_device <- function(path, fmt) {
  if (fmt == "pdf") {
    pdf(path, width = 11, height = 8)
  } else if (fmt == "svg") {
    svg(path, width = 11, height = 8)
  } else {
    png(path, width = 1400, height = 1000, res = 150)
  }
}

strip_extension <- function(path) {
  sub("\\.[^.]+$", "", path)
}

pair_cols <- c("AC", "AG", "AT", "CG", "CT", "GT")
pair_types <- c(AC = "Transversion", AG = "Transition", AT = "Transversion", CG = "Transversion", CT = "Transition", GT = "Transversion")

read_minorfreq <- function(path) {
  d <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE)
  required <- c("freq", pair_cols)
  missing <- setdiff(required, colnames(d))
  if (length(missing) > 0) {
    stop(sprintf("Missing columns in %s: %s", path, paste(missing, collapse = ", ")), call. = FALSE)
  }
  long <- do.call(
    rbind,
    lapply(pair_cols, function(pair) {
      data.frame(
        freq = d$freq,
        pair = pair,
        count = d[[pair]],
        type = pair_types[[pair]],
        stringsAsFactors = FALSE
      )
    })
  )
  long <- long[long$count > 0, , drop = FALSE]
  if (nrow(long) == 0) {
    stop(sprintf("No non-zero rows found in %s", path), call. = FALSE)
  }
  long
}

write_summary <- function(d, output_prefix) {
  summary_df <- aggregate(count ~ pair + type, d, sum)
  summary_df$mean_freq <- sapply(summary_df$pair, function(p) {
    idx <- d$pair == p
    weighted.mean(d$freq[idx], d$count[idx])
  })
  write.table(
    summary_df,
    file = paste0(output_prefix, "_minorfreq_summary.txt"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}

plot_overlay_gg <- function(d, sample_label, output_path, fmt, log_y, y_max, pairs, cols) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  d$pair <- factor(d$pair, levels = pairs, ordered = TRUE)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = freq, y = count, color = pair)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_manual(values = cols[pairs], drop = FALSE) +
    ggplot2::labs(
      title = sample_label,
      x = "Minor Allele Frequency",
      y = if (log_y) "Site Count (log10)" else "Site Count",
      color = "State"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))

  if (log_y) {
    p <- p + ggplot2::scale_y_log10(limits = c(1, y_max))
  } else {
    p <- p + ggplot2::coord_cartesian(ylim = c(0, y_max))
  }

  ggplot2::ggsave(output_path, p, width = 11, height = 8, units = "in", device = fmt)
}

build_panel_data <- function(d, pairs) {
  freqs <- sort(unique(d$freq))
  panel_data <- lapply(pairs, function(p) {
    subset <- d[d$pair == p, ]
    out <- rep(0, length(freqs))
    names(out) <- freqs
    if (nrow(subset) > 0) {
      out[match(subset$freq, freqs)] <- subset$count
    }
    data.frame(freq = freqs, count = out, pair = p)
  })
  list(freqs = freqs, panel_data = do.call(rbind, panel_data))
}

calc_y_cap <- function(df, cutoff) {
  if (!is.finite(cutoff)) {
    return(NA_real_)
  }
  subset <- df[df$freq >= cutoff & df$count > 0, , drop = FALSE]
  if (nrow(subset) == 0) {
    return(NA_real_)
  }
  2 * max(subset$count, na.rm = TRUE)
}

calc_y_max <- function(panel_data, cutoff) {
  y_max <- max(panel_data$count, na.rm = TRUE)
  global_y_cap <- calc_y_cap(panel_data, cutoff)
  if (is.finite(global_y_cap)) {
    y_max <- min(y_max, global_y_cap)
  }
  if (!is.finite(y_max) || y_max <= 0) {
    y_max <- 1
  }
  y_max
}

plot_overlay <- function(d, output_prefix, fmt, log_y, y_max, pairs, cols) {
  built <- build_panel_data(d, pairs)
  count_matrix <- sapply(pairs, function(p) {
    subset <- built$panel_data[built$panel_data$pair == p, ]
    subset$count
  })

  overlay_file <- paste0(output_prefix, "_minorfreq_overlay.", fmt)
  open_device(overlay_file, fmt)
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

  matplot(
    built$freqs,
    count_matrix,
    type = "b",
    lty = 1,
    lwd = 2,
    pch = 16,
    col = cols[pairs],
    log = if (log_y) "y" else "",
    ylim = c(if (log_y) 1 else 0, y_max),
    xlab = "Minor Allele Frequency",
    ylab = if (log_y) "Site Count (log10)" else "Site Count",
    main = "Minor Allele Frequency Overlay"
  )
  legend("topright", legend = pairs, col = cols[pairs], lty = 1, pch = 16, cex = 0.8)

  dev.off()
  par(oldpar)
}

plot_batch_combined <- function(d, output_path, fmt, log_y, y_max, pairs, cols) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  d$pair <- factor(d$pair, levels = pairs, ordered = TRUE)
  d$sample <- factor(d$sample, levels = unique(d$sample), ordered = TRUE)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = freq, y = count, color = pair)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_wrap(~sample, scales = "fixed") +
    ggplot2::scale_color_manual(values = cols[pairs], drop = FALSE) +
    ggplot2::labs(
      title = "Minor Allele Frequency Overlay",
      x = "Minor Allele Frequency",
      y = if (log_y) "Site Count (log10)" else "Site Count",
      color = "State"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      strip.text = ggplot2::element_text(face = "bold")
    )

  if (log_y) {
    p <- p + ggplot2::scale_y_log10(limits = c(1, y_max))
  } else {
    p <- p + ggplot2::coord_cartesian(ylim = c(0, y_max))
  }

  n_samples <- length(unique(d$sample))
  ncol <- if (n_samples <= 4) 2 else 3
  nrow <- ceiling(n_samples / ncol)
  height <- max(8, 3.2 * nrow)
  width <- if (ncol == 2) 12 else 15

  ggplot2::ggsave(output_path, p, width = width, height = height, units = "in", device = fmt)
}

read_sample_list <- function(path) {
  lines <- readLines(path, warn = FALSE)
  rows <- list()
  for (line in lines) {
    line <- trimws(line)
    if (line == "" || startsWith(line, "#")) {
      next
    }
    parts <- strsplit(line, "[[:space:]]+")[[1]]
    if (length(parts) < 2) {
      stop("Sample list rows must contain at least 2 columns: minorfreq_file output_prefix", call. = FALSE)
    }
    rows[[length(rows) + 1]] <- data.frame(
      minorfreq_file = parts[[1]],
      output_prefix = parts[[2]],
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) {
    stop("No valid rows found in sample list", call. = FALSE)
  }
  do.call(rbind, rows)
}

fmt <- if (!is.null(opt[["format"]])) tolower(opt[["format"]]) else "pdf"
if (!fmt %in% c("pdf", "svg", "png")) {
  stop("--format must be one of: pdf, svg, png", call. = FALSE)
}

log_y <- isTRUE(opt[["log-y"]])
y_cap_after_x <- if (!is.null(opt[["y-cap-after-x"]])) suppressWarnings(as.numeric(opt[["y-cap-after-x"]])) else NA_real_
if (!is.na(y_cap_after_x) && !is.finite(y_cap_after_x)) {
  stop("--y-cap-after-x must be numeric", call. = FALSE)
}

pairs <- pair_cols
cols <- c(
  "AC" = "#1b9e77",
  "AG" = "#d95f02",
  "AT" = "#7570b3",
  "CG" = "#e7298a",
  "CT" = "#66a61e",
  "GT" = "#e6ab02"
)

if (!is.null(opt[["sample-list"]])) {
  sample_df <- read_sample_list(opt[["sample-list"]])
  datasets <- vector("list", nrow(sample_df))
  y_values <- numeric(nrow(sample_df))

  for (i in seq_len(nrow(sample_df))) {
    d <- read_minorfreq(sample_df$minorfreq_file[[i]])
    built <- build_panel_data(d, pairs)
    datasets[[i]] <- d
    y_values[[i]] <- calc_y_max(built$panel_data, y_cap_after_x)
  }

  shared_y_max <- max(y_values, na.rm = TRUE)
  if (!is.finite(shared_y_max) || shared_y_max <= 0) {
    shared_y_max <- 1
  }

  batch_rows <- vector("list", nrow(sample_df))
  for (i in seq_len(nrow(sample_df))) {
    d <- datasets[[i]]
    output_prefix <- sample_df$output_prefix[[i]]
    sample_label <- basename(output_prefix)
    write_summary(d, output_prefix)
    batch_rows[[i]] <- data.frame(
      freq = d$freq,
      pair = d$pair,
      count = d$count,
      type = d$type,
      sample = sample_label,
      stringsAsFactors = FALSE
    )
  }
  batch_output <- if (!is.null(opt[["batch-out"]])) {
    opt[["batch-out"]]
  } else {
    paste0(strip_extension(sample_df$minorfreq_file[[1]]), "_batch_overlay.", fmt)
  }
  plot_batch_combined(do.call(rbind, batch_rows), batch_output, fmt, log_y, shared_y_max, pairs, cols)
  print_run_summary(sample_df$output_prefix, fmt, log_y, y_cap_after_x, batch_output = batch_output)
} else {
  if (is.null(opt[["minorfreq-file"]]) || is.null(opt[["out"]])) {
    show_help()
    quit(status = 1)
  }

  d <- read_minorfreq(opt[["minorfreq-file"]])
  built <- build_panel_data(d, pairs)
  y_max <- calc_y_max(built$panel_data, y_cap_after_x)
  write_summary(d, opt[["out"]])
  plot_overlay_gg(d, basename(opt[["out"]]), paste0(opt[["out"]], "_minorfreq_overlay.", fmt), fmt, log_y, y_max, pairs, cols)
  print_run_summary(opt[["out"]], fmt, log_y, y_cap_after_x)
}
