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
      "  Rscript Plot_MinorFreq.R --sample-list <list.txt> [--format pdf|svg|png] [--log-y] [--y-cap-after-x FLOAT]",
      "",
      "Single-sample arguments:",
      "  --minorfreq-file   Path to HEARTY consolidated .minorfreq.txt file",
      "  --out              Output prefix, e.g. results/plots/sample01",
      "",
      "Batch arguments:",
      "  --sample-list      Text file with: column 1 = minorfreq file, column 2 = output prefix",
      "",
      "Optional flags:",
      "  --format           Output format: pdf, svg, or png [default: pdf]",
      "  --log-y            Draw the y axis on a log10 scale",
      "  --y-cap-after-x    Set y max to 2x the highest count at freq >= cutoff",
      "",
      "Outputs:",
      "  <output_prefix>_minorfreq_overlay.<format>",
      "  <output_prefix>_minorfreq_summary.txt",
      sep = "\n"
    ),
    "\n"
  )
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

  for (i in seq_len(nrow(sample_df))) {
    d <- datasets[[i]]
    output_prefix <- sample_df$output_prefix[[i]]
    write_summary(d, output_prefix)
    plot_overlay(d, output_prefix, fmt, log_y, shared_y_max, pairs, cols)
  }
} else {
  if (is.null(opt[["minorfreq-file"]]) || is.null(opt[["out"]])) {
    show_help()
    quit(status = 1)
  }

  d <- read_minorfreq(opt[["minorfreq-file"]])
  built <- build_panel_data(d, pairs)
  y_max <- calc_y_max(built$panel_data, y_cap_after_x)
  write_summary(d, opt[["out"]])
  plot_overlay(d, opt[["out"]], fmt, log_y, y_max, pairs, cols)
}
