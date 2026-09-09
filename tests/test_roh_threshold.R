# Run with: Rscript tests/test_roh_threshold.R (requires the plotter's R packages).
# Add --unit-only to check the threshold calculation without plotting dependencies.
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- file.path(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))), "ROH_Window_Summary.R")

run_tests <- function() {
  tmp <- tempfile("hearty-roh-threshold-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  # Load only the threshold function to check edge cases without running the CLI.
  scope <- new.env()
  for (expr in parse(script)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("resolve_roh_threshold"))) eval(expr, scope)
  }
  values <- c(0, 0.0002, 0.001, 0.002, 0.01)
  df <- data.frame(prop = values)
  stopifnot(isTRUE(all.equal(scope$resolve_roh_threshold(df, 0.2, "median"), 0.0002)))
  stopifnot(scope$resolve_roh_threshold(df, 0.2, "absolute") == 0.2)
  stopifnot(scope$resolve_roh_threshold(data.frame(prop = c(0, 0, 1)), 0.2, "median") == 0)
  stopifnot(inherits(try(scope$resolve_roh_threshold(df, 0.2, "mean"), silent = TRUE), "try-error"))
  cat("PASS: median calculation, absolute cutoff, zero median, and rejection of mean mode\n")
  if ("--unit-only" %in% commandArgs(trailingOnly = TRUE)) return(invisible(NULL))

  files <- file.path(tmp, c("sampleA.windows.txt", "sampleB.windows_trv.txt"))
  for (i in seq_along(files)) {
    label <- if (i == 1) "HET" else "HET_trv"
    props <- values / i
    writeLines(sprintf("chr1 %d-%d: 100000 sites, %d sites with '%s', Proportion: %.7f",
                       (0:4) * 100000 + 1, (1:5) * 100000,
                       as.integer(props * 100000), label, props), files[i])
  }
  samples <- file.path(tmp, "samples.tsv")
  writeLines(c("sample\tfile", paste(c("sampleA", "sampleB"), files, sep = "\t")), samples)

  cli <- function(args, success = TRUE) {
    log <- file.path(tmp, "command.log")
    status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(script, args)), stdout = log, stderr = log)
    output <- paste(readLines(log, warn = FALSE), collapse = "\n")
    if ((status == 0) != success) stop(output)
    output
  }
  tuning <- c("--min-length-mb", "0.1", "--bin-size-bp", "100000")
  relative <- c("--roh-threshold-mode", "median", "--roh-threshold", "0.2")
  pdf <- file.path(tmp, "single.pdf")
  output <- cli(c("--mode", "single", "--input", files[1], "--output-pdf", pdf, tuning, relative))
  stopifnot(grepl("0.000200", output, fixed = TRUE), file.info(pdf)$size > 100)
  summary <- file.path(tmp, "summary.tsv")
  batch <- c("--mode", "batch", "--sample-list", samples, "--summary-out", summary,
             "--combined-pdf", file.path(tmp, "batch.pdf"), tuning)
  cli(c(batch, relative))
  result <- read.delim(summary)
  stopifnot(all(result$roh_threshold_mode == "median"))
  stopifnot(isTRUE(all.equal(result$roh_threshold, c(0.0002, 0.0001))))
  stopifnot(all(result$raw_froh_gt100kb == 0.4))
  stopifnot(isTRUE(all.equal(result$mean_heterozygosity_all_windows, c(mean(values), mean(values) / 2))))

  # Omitting the mode retains the existing absolute default.
  cli(c(batch, "--roh-threshold", "0.0001"))
  result <- read.delim(summary)
  stopifnot(all(result$roh_threshold_mode == "absolute"), all(result$roh_threshold == 0.0001))
  for (mode in c("single", "batch")) {
    output <- cli(c("--mode", mode, "--roh-threshold-mode", "mean"), success = FALSE)
    stopifnot(grepl("must be one of: absolute, median", output, fixed = TRUE))
  }
  cat("PASS: median cutoffs in single/batch, sample-specific values, zero median, absolute default, unchanged mean summaries, and rejection of mean mode\n")
}

run_tests()
