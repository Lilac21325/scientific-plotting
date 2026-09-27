args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) stop("Arguments: skill_dir output_dir [R_library]")
skill <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
out <- args[2]
if (length(args) == 3) .libPaths(c(args[3], .libPaths()))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(skill, "scripts", "figure-tools.R"))
for (f in list.files(file.path(skill, "recipes", "r"), pattern = "[.]R$", full.names = TRUE)) source(f)
sp_require(c("ggplot2", "patchwork", "ggrepel", "ggridges", "sf", "circlize", "ggseqlogo", "ragg", "svglite", "jsonlite"))
set.seed(20260927)
checks <- list()
check <- function(name, value) {
  if (!isTRUE(value)) stop("FAILED: ", name)
  checks[[name]] <<- TRUE
}
reject <- function(name, expr) check(name, inherits(tryCatch({force(expr); NULL}, error = identity), "error"))
save_plot <- function(p, name, w = 150, h = 100, formats = c("png", "pdf")) {
  ggplot2::ggplot_build(p)
  export_figure(p, file.path(out, name), w, h, formats = formats, overwrite = TRUE)
  check(paste0(name, "_exported"), all(file.info(paste0(file.path(out, name), ".", formats))$size > 1000))
}

d <- data.frame(group = factor(rep(c("Control", "Treatment", "Reference"), each = 45),
                               levels = c("Control", "Treatment", "Reference")),
                x = rep(seq(1, 9, length.out = 45), 3))
d$response <- 2 + 0.65 * d$x + rep(c(0, 1.6, -0.8), each = 45) + rnorm(135, sd = 1)
d$mass <- runif(135, 1, 30)
d$id <- paste0("S", seq_len(nrow(d)))
single <- distribution_plot(d, "response", binwidth = 0.8)
save_plot(single, "01-distribution", formats = c("png", "pdf", "svg"))
groups <- group_distribution(d, "group", "response")
save_plot(groups, "02-group-distribution")
check("all_raw_points_retained", nrow(ggplot2::ggplot_build(groups)$data[[3]]) == nrow(d))
save_plot(ridge_plot(d, "group", "response", bandwidth = 0.7), "03-ridges")
fit <- scatter_fit(d, "x", "response", "group")
save_plot(fit$plot, "04-regression")
check("model_n_and_grouping", nrow(fit$statistics) == 3 && all(fit$statistics$n == 45))
check("CI_bounds_ordered", all(fit$predictions$lwr <= fit$predictions$fit & fit$predictions$upr >= fit$predictions$fit))
check("R2_from_actual_fit", abs(fit$statistics$r_squared[1] - summary(fit$models[[1]])$r.squared) < 1e-12)
overall <- scatter_fit(d, "x", "response", "group", fit_scope = "overall")
save_plot(overall$plot, "05-overall-fit")
check("overall_uses_all_rows", overall$statistics$n == nrow(d))
save_plot(bubble_plot(d[seq(1, 135, 3), ], "x", "response", "mass", "group"), "06-bubbles")
save_plot(density2d_plot(d, "x", "response"), "07-density2d")
save_plot(highlight_scatter(d, "x", "response", "id", d$response > quantile(d$response, 0.97)), "08-highlight")

tr <- expand.grid(time = 1:6, group = factor(c("Control", "Treatment"), levels = c("Control", "Treatment")))
tr$response <- 5 + tr$time * ifelse(tr$group == "Control", 0.3, 0.6) + rnorm(nrow(tr), sd = 0.15)
save_plot(trend_plot(tr, "time", "response", "group"), "09-trends")
paired <- expand.grid(id = paste0("P", 1:18), occasion = c("Before", "After"), stringsAsFactors = FALSE)
baseline <- rnorm(18, 5, 1)
paired$value <- c(baseline, baseline + rnorm(18, 0.8, 0.45))
pair_plot <- paired_slope(paired, "id", "occasion", "value", c("Before", "After"))
save_plot(pair_plot, "10-paired")
check("pair_values_unchanged", identical(sort(ggplot2::ggplot_build(pair_plot)$data[[2]]$y), sort(paired$value)))
check("one_line_group_per_subject", length(unique(ggplot2::ggplot_build(pair_plot)$data[[1]]$group)) == 18)
reject("incomplete_pair_rejected", paired_slope(paired[-1, ], "id", "occasion", "value", c("Before", "After")))
reject("duplicate_pair_rejected", paired_slope(rbind(paired, paired[1, ]), "id", "occasion", "value", c("Before", "After")))

ranked <- data.frame(category = c("Long outcome A", "Outcome B", "Outcome C", "Outcome D", "Outcome E"), value = c(18, 12, 24, 5, 9))
save_plot(ordered_values(ranked, "category", "value"), "11-lollipop")
save_plot(ordered_values(ranked, "category", "value", kind = "bar"), "12-ordered-bars")
bars <- expand.grid(category = c("A", "B", "C"), group = c("Control", "Treatment"), stringsAsFactors = FALSE)
bars$value <- c(12, 18, 15, 21, 20, 16)
save_plot(grouped_bars(bars, "category", "group", "value"), "13-grouped-bars")
save_plot(grouped_bars(bars, "category", "group", "value", position = "fill"), "14-proportions")

samples <- expand.grid(category = c("Outcome A", "Outcome B", "Outcome C"), group = c("Control", "Treatment"), replicate = 1:24)
samples$value <- 5 + as.numeric(samples$category) * 0.5 + ifelse(samples$group == "Treatment", 0.8, 0) + rnorm(nrow(samples))
keys <- unique(samples[c("category", "group")])
est <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  z <- subset(samples, category == keys$category[i] & group == keys$group[i])$value
  margin <- qt(0.975, length(z) - 1) * sd(z) / sqrt(length(z))
  data.frame(keys[i, ], estimate = mean(z), lower = mean(z) - margin, upper = mean(z) + margin, n = length(z))
}))
ci <- estimate_plot(est, "category", "estimate", "lower", "upper", "95% t CI of independent simulated means; n=24 per cell", group = "group")
save_plot(ci, "15-estimates-ci")
ci_layers <- ggplot2::ggplot_build(ci)$data
check("CI_point_dodge_alignment", identical(ci_layers[[1]]$x, ci_layers[[2]]$x))
check("CI_values_preserved", identical(sort(ci_layers[[1]]$ymin), sort(est$lower)))
bad <- est; bad$lower[1] <- bad$estimate[1] + 1
reject("bad_interval_rejected", estimate_plot(bad, "category", "estimate", "lower", "upper", "95% CI", group = "group"))

hm <- expand.grid(gene = paste0("Gene ", 1:5), sample = paste0("S", 1:7), stringsAsFactors = FALSE)
hm$expression <- rnorm(nrow(hm)); hm$expression[3] <- NA
heat <- heatmap_plot(hm, "gene", "sample", "expression", midpoint = 0)
save_plot(heat, "16-heatmap")
check("missing_heatmap_cell_retained", nrow(ggplot2::ggplot_build(heat)$data[[1]]) == nrow(hm))
reject("incomplete_matrix_rejected", heatmap_plot(hm[-1, ], "gene", "sample", "expression"))
correl <- correlation_plot(d, c("x", "response", "mass"))
save_plot(correl$plot, "17-correlation")
check("correlation_matches_R", isTRUE(all.equal(correl$matrix, cor(d[c("x", "response", "mass")]))))

multi <- compose_panels(list(single + ggplot2::labs(title = "Distribution"),
                              groups + ggplot2::labs(title = "Group comparison"),
                              fit$plot + ggplot2::labs(title = "Association"),
                              ci + ggplot2::labs(title = "Estimates and uncertainty")), ncol = 2)
export_figure(multi, file.path(out, "18-multipanel"), 180, 160, overwrite = TRUE)
check("multipanel_exported", file.info(file.path(out, "18-multipanel.png"))$size > 1000)
reject("unsafe_guide_collection_rejected", compose_panels(list(single, groups), collect = TRUE))

seqs <- replicate(80, paste0(vapply(seq_len(10), function(i) sample(c("A", "C", "G", "T"), 1,
                          prob = if (i %% 3 == 0) c(0.05, 0.1, 0.8, 0.05) else c(0.55, 0.2, 0.15, 0.1)), character(1)), collapse = ""))
save_plot(sequence_logo(seqs, "bits"), "19-sequence-logo", w = 155, h = 80)
save_plot(sequence_logo(seqs, "prob"), "20-logo-probability", w = 155, h = 80)
reject("unaligned_sequences_rejected", sequence_logo(c("ACGT", "AG")))
links <- matrix(c(12, 5, 8, 9, 11, 4, 3, 8, 15), nrow = 3,
                dimnames = list(c("Source A", "Source B", "Source C"), c("Target X", "Target Y", "Target Z")))
ragg::agg_png(file.path(out, "21-chord.png"), width = 140, height = 140, units = "mm", res = 300, background = "white")
chord_result <- draw_chord(links, directional = TRUE)
grDevices::dev.off()
grDevices::cairo_pdf(file.path(out, "21-chord.pdf"), width = 140/25.4, height = 140/25.4)
draw_chord(links, directional = TRUE)
grDevices::dev.off()
check("chord_relationships_preserved", nrow(chord_result) == length(links))
check("chord_state_cleared", length(circlize::get.all.sector.index()) == 0)

make_square <- function(x, y) sf::st_polygon(list(matrix(c(x, y, x+1, y, x+1, y+1, x, y+1, x, y), ncol = 2, byrow = TRUE)))
geo <- sf::st_sf(region = LETTERS[1:6], geometry = sf::st_sfc(lapply(0:5, function(i) make_square(i %% 3, i %/% 3)), crs = 3857))
spatial_data <- data.frame(region = LETTERS[1:5], rate = c(2, 4, 3, 7, NA))
reject("unmatched_map_requires_opt_in", choropleth_plot(geo, spatial_data, "region", "region", "rate", "Rate / 1000"))
map <- choropleth_plot(geo, spatial_data, "region", "region", "rate", "Rate / 1000", allow_unmatched = TRUE)
save_plot(map$plot + ggplot2::coord_sf(datum = NA) + ggplot2::labs(title = "Synthetic regions: missing values remain visible"), "22-map", w = 140, h = 100)
check("all_regions_retained", nrow(map$data) == 6 && map$missing_values == 2 && identical(map$unmatched_geometry, "F"))
reject("duplicate_map_join_rejected", choropleth_plot(geo, rbind(spatial_data, spatial_data[1, ]), "region", "region", "rate", "Rate", allow_unmatched = TRUE))

reject("missing_observation_rejected", distribution_plot(data.frame(value = c(1, NA, 3)), "value", 1))
reject("negative_bubble_size_rejected", bubble_plot(data.frame(x = 1:3, y = 1:3, m = c(2,-1,3)), "x", "y", "m"))
reject("palette_not_silently_recycled", scientific_palette(letters[1:8]))
reject("export_overwrite_guard", export_figure(single, file.path(out, "01-distribution"), 150, 100))
numeric_groups <- d; numeric_groups$group <- as.integer(numeric_groups$group)
check("numeric_group_bubbles", nrow(ggplot2::ggplot_build(bubble_plot(numeric_groups, "x", "response", "mass", "group"))$data[[1]]) == nrow(d))
numeric_trend <- tr; numeric_trend$group <- as.integer(numeric_trend$group)
check("numeric_group_trends", nrow(ggplot2::ggplot_build(trend_plot(numeric_trend, "time", "response", "group"))$data[[2]]) == nrow(tr))
numeric_bars <- bars; numeric_bars$group <- rep(1:2, each = 3)
check("numeric_group_bars", nrow(ggplot2::ggplot_build(grouped_bars(numeric_bars, "category", "group", "value"))$data[[1]]) == nrow(bars))
connected <- estimate_plot(est, "category", "estimate", "lower", "upper", "95% CI", "group", connect = TRUE)
check("connected_CI_alignment", identical(ggplot2::ggplot_build(connected)$data[[1]]$x, ggplot2::ggplot_build(connected)$data[[3]]$x))
check("no_group_OLS", scatter_fit(d, "x", "response")$statistics$n == nrow(d))
check("nonviolin_raw_points", nrow(ggplot2::ggplot_build(group_distribution(d, "group", "response", violin = FALSE))$data[[2]]) == nrow(d))
zeros <- bars; zeros$value <- 0
reject("zero_denominator_rejected", grouped_bars(zeros, "category", "group", "value", "fill"))
sp_require("systemfonts")
fonts <- systemfonts::system_fonts()
cjk <- unique(fonts$family[grepl("Microsoft YaHei", fonts$family)])[1]
if (!is.na(cjk)) {
  cjk_plot <- ci + ggplot2::labs(title = "\u7ec4\u95f4\u6bd4\u8f83", y = "\u4f30\u8ba1\u503c", x = NULL) +
    theme_scientific(base_family = cjk)
  save_plot(cjk_plot, "23-cjk-font", w = 150, h = 100)
  check("CJK_font_selected", nzchar(cjk))
}
jsonlite::write_json(list(status = "passed", checks = checks, seed = 20260927,
  packages = lapply(c("ggplot2", "patchwork", "ggseqlogo", "circlize", "sf", "ragg", "ggridges", "ggrepel", "svglite"),
                    function(p) list(name = p, version = as.character(packageVersion(p))))),
  file.path(out, "validation.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
cat("PASS", length(checks), "checks\n")
