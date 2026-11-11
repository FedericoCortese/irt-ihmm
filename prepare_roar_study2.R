# Prepare ROAR Study 2 XYZ and describe the data. No model is fitted.
# Run source("prepare_roar_study2.R") from the project directory. Base R only.
input_dir <- file.path(path.expand("~"), "Library", "CloudStorage",
  "OneDrive-UniversitàdegliStudidiMilano", "Cortese-Mignemi", "roar_study2")
input_file <- file.path(input_dir, "roar_study2_xyz.csv")
output_dir <- file.path("output", "roar_study2")
standardize_covariates <- TRUE
# Input contract: y is already natural-log response time in seconds. Never log it again.
# Only the requested XYZ file is read; no joins to earlier datasets or *_std files.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Input not found; edit input_dir: ", input_file)
raw <- read.csv(input_file, na.strings = c("", "NA"), stringsAsFactors = FALSE)
required <- c("subject", "t", "item", "x", "y", "z1", "z2")
if (!setequal(names(raw), required) || any(!vapply(raw, is.numeric, logical(1))))
  stop("Expected numeric columns: ", paste(required, collapse = ", "))
if (any(!is.finite(as.matrix(raw[c("subject", "t", "item", "x", "z1", "z2")]))))
  stop("Keys, accuracy and covariates must be complete and finite.")
if (any(!raw$x %in% 0:1) || anyDuplicated(raw[c("subject", "t")]) ||
    anyDuplicated(raw[c("subject", "item")])) stop("Invalid accuracy or duplicate trial keys.")
subject_ids <- sort(unique(raw$subject)); item_ids <- sort(unique(raw$item))
P <- length(subject_ids); I <- length(item_ids)
if (nrow(raw) != P * I) stop("Expected a complete respondent by item panel.")
for (id in subject_ids) {
  d <- raw[raw$subject == id, ]
  if (!setequal(d$t, seq_len(I)) || !setequal(d$item, item_ids))
    stop("Incomplete panel or invalid presentation positions: ", id)
  for (name in c("z1", "z2")) if (length(unique(d[[name]])) != 1L)
    stop("Covariate must be constant within respondent: ", name)
}
trials <- raw[order(raw$subject, raw$t), ]; row.names(trials) <- NULL
subjects <- trials[match(subject_ids, trials$subject), c("subject", "z1", "z2")]
row.names(subjects) <- NULL
z_raw <- as.matrix(subjects[c("z1", "z2")]); rownames(z_raw) <- as.character(subject_ids)
# Estimate scaling across unique respondents, not repeated trial rows.
centers <- colMeans(z_raw); scales <- apply(z_raw, 2L, sd)
if (any(!is.finite(scales) | scales <= 0)) stop("Covariates must vary across respondents.")
z_standardized <- sweep(sweep(z_raw, 2L, centers, "-"), 2L, scales, "/")
z <- if (standardize_covariates) z_standardized else z_raw
scaling <- data.frame(covariate = colnames(z_raw), center = centers, sd = scales,
  standardized_for_z = standardize_covariates, row.names = NULL)
trials$valid_y <- is.finite(trials$y)
trials$rt_s <- NA_real_
trials$rt_s[trials$valid_y] <- exp(trials$y[trials$valid_y])
if (any(!is.finite(trials$rt_s[trials$valid_y]) | trials$rt_s[trials$valid_y] <= 0))
  stop("Finite log-times must imply finite positive response times.")
trials$outside_reference_bounds <- NA
trials$outside_reference_bounds[trials$valid_y] <-
  trials$rt_s[trials$valid_y] < .2 | trials$rt_s[trials$valid_y] > 5

make_matrix <- function(value, keys, columns) {
  out <- matrix(NA_real_, P, I, dimnames = list(as.character(subject_ids), as.character(columns)))
  out[cbind(match(trials$subject, subject_ids), match(keys, columns))] <- value
  out
}
chronological <- list(x = make_matrix(trials$x, trials$t, seq_len(I)),
  y = make_matrix(trials$y, trials$t, seq_len(I)),
  item_id = make_matrix(trials$item, trials$t, seq_len(I)), z = z)
item_aligned <- list(x = make_matrix(trials$x, trials$item, item_ids),
  y = make_matrix(trials$y, trials$item, item_ids),
  position = make_matrix(trials$t, trials$item, item_ids), z = z)
complete_times <- rowSums(!is.finite(chronological$y)) == 0L
common_order <- all(vapply(seq_len(P), function(p)
  identical(chronological$item_id[p, ], chronological$item_id[1, ]), logical(1)))
availability <- data.frame(subject = subject_ids, complete_log_times = complete_times,
  complete_covariates = complete.cases(z_raw))

metric_summary <- function(d) {
  rt <- d$rt_s[d$valid_y]
  q <- if (length(rt)) quantile(rt, c(.25, .5, .75, .95, 1), names = FALSE) else rep(NA_real_, 5)
  data.frame(n_trials = nrow(d), accuracy = mean(d$x), n_valid_times = length(rt),
    n_invalid_times = sum(!d$valid_y),
    n_outside_bounds = sum(d$outside_reference_bounds, na.rm = TRUE),
    outside_fraction_of_valid_times = if (length(rt)) mean(d$outside_reference_bounds[d$valid_y]) else NA_real_,
    mean_rt_s = if (length(rt)) mean(rt) else NA_real_,
    p25_rt_s = q[1], median_rt_s = q[2], p75_rt_s = q[3], p95_rt_s = q[4], max_rt_s = q[5])
}
by_group <- function(field) {
  out <- lapply(split(seq_len(nrow(trials)), trials[[field]]), function(idx)
    cbind(trials[idx[1], field, drop = FALSE], metric_summary(trials[idx, ])))
  out <- do.call(rbind, out); row.names(out) <- NULL
  out[order(out[[field]]), , drop = FALSE]
}
overall <- metric_summary(trials)
by_subject <- merge(subjects, by_group("subject"), by = "subject", sort = TRUE)
by_item <- by_group("item"); by_position <- by_group("t")
covariate_summary <- do.call(rbind, lapply(colnames(z_raw), function(name) {
  v <- z_raw[, name]
  data.frame(covariate = name, n_subjects = length(v), mean = mean(v), sd = sd(v),
    minimum = min(v), median = median(v), maximum = max(v))
}))
prepared <- list(source = list(csv = normalizePath(input_file),
  md5 = unname(tools::md5sum(input_file)), y_definition = "natural-log response time in seconds"),
  P = P, I = I, subject_ids = subject_ids, item_ids = item_ids,
  trials = trials, subjects = subjects, chronological = chronological,
  item_aligned = item_aligned, z_raw = z_raw, z_standardized = z_standardized,
  standardize_covariates = standardize_covariates, scaling = scaling,
  availability = availability, fit_ready = FALSE,
  fit_blockers = c(if (!common_order)
    "Different item orders: the sampler needs an item-index matrix for each respondent/position.",
    if (!all(complete_times)) "Missing/nonfinite log-times require an explicit handling policy."),
  summaries = list(overall = overall, by_subject = by_subject, by_item = by_item,
    by_position = by_position, covariate_summary = covariate_summary))
saveRDS(prepared, file.path(output_dir, "prepared_data.rds"))
tables <- c(prepared$summaries, list(covariate_scaling = scaling, subject_metadata = subjects,
  sample_availability = availability,
  flagged_times = trials[!trials$valid_y | (!is.na(trials$outside_reference_bounds) &
    trials$outside_reference_bounds), c(required, "rt_s", "outside_reference_bounds")]))
for (name in names(tables)) write.csv(tables[[name]],
  file.path(output_dir, paste0(name, ".csv")), row.names = FALSE, na = "")

png_plot <- function(name, draw) {
  png(file.path(output_dir, paste0(name, ".png")), width = 1800, height = 1200, res = 150)
  on.exit(dev.off()); draw()
}
png_plot("overview", function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  hist(trials$y[trials$valid_y], breaks = 60, col = "#648eae", border = "white",
    main = "Observed log-response times (all finite values)", xlab = "Natural log of response time (s)")
  abline(v = log(c(.2, 5)), col = "#c04b39", lty = 2)
  hist(by_subject$accuracy, breaks = 15, col = "#648eae", border = "white",
    main = "Respondent accuracy", xlab = "Fraction correct", xlim = c(0, 1))
  for (name in c("z1", "z2")) hist(z_raw[, name], breaks = 15,
    col = "#648eae", border = "white", main = paste("Covariate", name, "(original units)"),
    xlab = paste(name, "(one observation per respondent)"))
})
png_plot("respondents_and_items", function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  plot(by_subject$median_rt_s, by_subject$accuracy, pch = 19, col = "#244b78",
    main = "Respondents", xlab = "Median observed response time (s)", ylab = "Fraction correct")
  plot(by_item$median_rt_s, by_item$accuracy, pch = 19, col = "#244b78",
    main = "Items", xlab = "Median observed response time (s)", ylab = "Fraction correct")
  for (name in c("z1", "z2")) plot(by_subject[[name]], by_subject$accuracy,
    pch = 19, col = "#244b78", main = paste(name, "and accuracy"),
    xlab = paste(name, "(original units)"), ylab = "Fraction correct")
})
png_plot("presentation_order", function() {
  par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 1))
  plot(by_position$t, by_position$accuracy, type = "l", col = "#244b78",
    main = "Accuracy along the actual presentation sequence", xlab = "Trial position", ylab = "Fraction correct")
  plot(by_position$t, by_position$median_rt_s, type = "l", col = "#244b78",
    main = "Response time along the actual presentation sequence", xlab = "Trial position",
    ylab = "Median observed time (s)", ylim = range(by_position[c("p25_rt_s", "p75_rt_s")]))
  polygon(c(by_position$t, rev(by_position$t)), c(by_position$p25_rt_s, rev(by_position$p75_rt_s)),
    col = "#648eae40", border = NA)
  lines(by_position$t, by_position$median_rt_s, col = "#244b78")
  legend("topright", c("Median", "Interquartile range across respondents"),
    col = c("#244b78", "#648eae"), lty = c(1, NA), pch = c(NA, 15), bty = "n", cex = .8)
})
report <- c("ROAR Study 2 XYZ: preparation and descriptive analysis", "",
  paste("Source:", basename(input_file)),
  sprintf("Respondents: %d; items: %d; trials: %d.", P, I, nrow(trials)),
  sprintf("Accuracy: %.4f; median observed RT: %.4f s; maximum RT: %.4f s.",
    overall$accuracy, overall$median_rt_s, overall$max_rt_s),
  sprintf("Missing/nonfinite log-times: %d; affected respondents: %d; complete respondents: %d.",
    overall$n_invalid_times, sum(!complete_times), sum(complete_times)),
  sprintf("Times outside [0.2, 5] s: %d (%.2f%% of valid times). Not trimmed.",
    overall$n_outside_bounds, 100 * overall$outside_fraction_of_valid_times),
  paste("Covariates standardized for the saved z:", standardize_covariates), "",
  "z1 and z2 are retained under their supplied names; their meaning is not documented in this CSV.",
  "z_raw stores original units; z_standardized stores (z-center)/sample_SD across unique subjects.",
  "scaling records the center and SD to reuse when transforming future respondents.",
  "Standardizing is not mathematically required, but changes the effective slope prior.",
  "With beta~N(0,1), standardized continuous covariates make effects comparable per one SD",
  "and avoid an overly diffuse probit linear predictor at the raw covariate scale.",
  "x and y are unaltered: y is interpreted as already natural-log seconds, never logged again.",
  "Missing y values remain missing; their original raw times cannot be inferred from this CSV.",
  "No trimming, winsorization, imputation or respondent exclusion was performed.", "",
  "Load with: prepared <- readRDS('output/roar_study2/prepared_data.rds')",
  "chronological$x/y: respondents by actual presentation position.",
  "chronological$item_id: identity of the item at each respondent/position.",
  "chronological$z: chosen covariate matrix, without an intercept, aligned by subject.",
  "item_aligned$x/y: respondents by common item ID; this is not HMM temporal order.",
  "availability flags complete log-time profiles without dropping rows.",
  "No model has been fitted. fit_ready is FALSE. Remaining blockers:", prepared$fit_blockers,
  "Do not fit either view with the current common-order sampler before adapting item indexing.",
  "Only roar_study2_xyz.csv is used. Earlier long and pre-standardized files are not read.")
writeLines(report, file.path(output_dir, "preparation_report.txt"))
cat(paste(report[1:9], collapse = "\n"), "\n")
message("Prepared XYZ data, tables and figures saved in ", output_dir, "/.")
