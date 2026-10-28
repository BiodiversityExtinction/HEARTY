# Run with: Rscript tests/test_roh_threshold.R (requires the plotter's R packages).
# Add --unit-only to check the threshold calculation without plotting dependencies.
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- file.path(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))), "ROH_Window_Summary.R")

run_tests <- function() {
  tmp <- tempfile("hearty-roh-threshold-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  # Help must exit before loading plotting dependencies.
  for (help_args in list("--help", "-h", character())) {
    help_log <- file.path(tmp, "help.log")
    status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(script, help_args)),
                      stdout = help_log, stderr = help_log)
    stopifnot(status == 0, any(grepl("Usage:", readLines(help_log), fixed = TRUE)))
  }
  cat("PASS: --help, -h, and no-argument usage without plotting dependencies\n")

  # Load dependency-free helpers without running the CLI.
  scope <- new.env()
  for (expr in parse(script)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in%
        c("resolve_roh_threshold", "default_single_summary_path", "write_single_summary")) eval(expr, scope)
  }
  values <- c(0, 0.0002, 0.001, 0.002, 0.01)
  df <- data.frame(prop = values)
  stopifnot(isTRUE(all.equal(scope$resolve_roh_threshold(df, 0.2, "median"), 0.0002)))
  stopifnot(scope$resolve_roh_threshold(df, 0.2, "absolute") == 0.2)
  stopifnot(scope$resolve_roh_threshold(data.frame(prop = c(0, 0, 1)), 0.2, "median") == 0)
  stopifnot(inherits(try(scope$resolve_roh_threshold(df, 0.2, "mean"), silent = TRUE), "try-error"))
  cat("PASS: median calculation, absolute cutoff, zero median, and rejection of mean mode\n")
  stopifnot(scope$default_single_summary_path("plots/sample.pdf") == "plots/sample_summary.txt")
  stopifnot(scope$default_single_summary_path("sample.PDF") == "sample_summary.txt")
  stopifnot(scope$default_single_summary_path("sample") == "sample_summary.txt")
  text_file <- file.path(tmp, "nested", "summary.txt")
  report <- c("Mean heterozygosity: 0.001", "", "FROH: 0.2")
  screen <- capture.output(scope$write_single_summary(report, text_file))
  stopifnot(identical(readLines(text_file), report), identical(screen[seq_along(report)], report))
  cat("PASS: automatic summary filenames and identical saved/printed report text\n")
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
  report <- paste(readLines(file.path(tmp, "single_summary.txt")), collapse = "\n")
  stopifnot(grepl(report, output, fixed = TRUE), grepl("FROH (bridged", report, fixed = TRUE))
  stopifnot(grepl("Mean heterozygosity (non-ROH, raw definition)", report, fixed = TRUE))
  custom <- file.path(tmp, "reports", "transversion_summary.txt")
  output <- cli(c("--mode", "single", "--input", files[2], "--output-pdf", pdf,
                  "--summary-out", custom, "--froh-cutoffs", "100kb,2mb", tuning, relative))
  report <- paste(readLines(custom), collapse = "\n")
  stopifnot(grepl(report, output, fixed = TRUE), grepl("Dataset type: transversion-only", report, fixed = TRUE))
  stopifnot(grepl(">2 Mb:", report, fixed = TRUE))
  before <- readLines(files[1])
  cli(c("--mode", "single", "--input", files[1], "--output-pdf", pdf,
        "--summary-out", files[1], tuning, relative), success = FALSE)
  stopifnot(identical(before, readLines(files[1])))
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
