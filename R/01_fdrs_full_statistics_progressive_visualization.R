############################################################
# FDRS integrated script
# Full statistics + progressive subsampling + visualization
############################################################

# ==========================================================
# 0. Working directory and user configuration
# ==========================================================

setwd("F:/Data Model/Figure/Aov")

set.seed(20260524)

# 输入文件，需放在 WORK_DIR 目录下
INPUT_FILES <- c("RawData.txt", "ErrData.txt")
INPUT_NAMES <- c("RawData", "ErrData")

# 小数位设置
SINGLE_DIGIT_POS <- 3
JOINT_START_POS  <- 3
JOINT_LENGTH     <- 2
DIGITS_KEEP      <- 10

# 渐进式抽样设置
B_PROGRESSIVE <- 1000
N_PROGRESSIVE_POINTS <- 12
MIN_PROGRESSIVE_N <- 30

# 输出文件夹
OUTPUT_DIR <- "FDRS_integrated_results"

# ==========================================================
# 1. Load packages
# ==========================================================

pkg_needed <- c("ggplot2", "dplyr", "tidyr", "stringr", "scales")

for (pkg in pkg_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==========================================================
# 2. Color palette and figure theme
# ==========================================================

# 与前面图中 Con、0.25、25 一致的配色
color_con <- "#A9D1E8"
color_025 <- "#F6C1C3"
color_25  <- "#E68589"

theme_fdrs_style <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text.x = element_text(
        color = "black",
        size = base_size - 2,
        face = "bold"
      ),
      axis.text.y = element_text(
        color = "black",
        size = base_size - 2
      ),
      axis.title.x = element_text(
        color = "black",
        size = base_size,
        face = "bold"
      ),
      axis.title.y = element_text(
        color = "black",
        size = base_size,
        face = "bold"
      ),
      axis.line = element_line(
        color = "black",
        linewidth = 0.8
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.9
      ),
      panel.grid.major = element_line(
        color = "grey90",
        linewidth = 0.4
      ),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.8
      ),
      strip.text = element_text(
        color = "black",
        face = "bold",
        size = base_size - 1
      ),
      legend.title = element_text(
        color = "black",
        size = base_size - 1,
        face = "bold"
      ),
      legend.text = element_text(
        color = "black",
        size = base_size - 2
      ),
      plot.title = element_text(
        color = "black",
        size = base_size + 1,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        color = "black",
        size = base_size - 1,
        hjust = 0.5
      ),
      plot.background = element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = element_rect(
        fill = "white",
        color = NA
      )
    )
}

save_plot <- function(plot_obj, filename, width = 7, height = 5.5) {
  ggsave(
    filename = file.path(OUTPUT_DIR, paste0(filename, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(OUTPUT_DIR, paste0(filename, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    bg = "white"
  )
}

format_metric <- function(x) {
  sapply(x, function(z) {
    if (is.na(z)) return("")
    if (z != 0 && abs(z) < 0.0001) {
      return(format(z, scientific = TRUE, digits = 3))
    }
    sprintf("%.4f", z)
  })
}

scale01 <- function(v) {
  if (all(is.na(v))) {
    return(rep(NA_real_, length(v)))
  }
  
  rng <- range(v, na.rm = TRUE)
  
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(0.5, length(v)))
  }
  
  scales::rescale(v, to = c(0, 1), from = rng)
}

# ==========================================================
# 3. Data reading and decimal digit extraction
# ==========================================================

read_single_column_txt <- function(file) {
  raw <- readLines(file, warn = FALSE, encoding = "UTF-8")
  raw <- trimws(raw)
  raw <- raw[raw != ""]
  
  raw <- sapply(strsplit(raw, "[,\t ]+"), `[`, 1)
  raw <- trimws(raw)
  
  num_check <- suppressWarnings(as.numeric(raw))
  raw <- raw[!is.na(num_check)]
  
  return(as.character(raw))
}

get_decimal_part <- function(x, digits_keep = 10) {
  x <- trimws(as.character(x))
  x <- gsub(",", "", x)
  out <- rep(NA_character_, length(x))
  
  for (i in seq_along(x)) {
    z <- x[i]
    num <- suppressWarnings(as.numeric(z))
    
    if (is.na(num)) {
      out[i] <- NA_character_
      next
    }
    
    if (grepl("[eE]", z) || !grepl("\\.", z)) {
      z2 <- formatC(num, format = "f", digits = digits_keep)
    } else {
      z2 <- z
    }
    
    z2 <- gsub("^[-+]", "", z2)
    
    if (!grepl("\\.", z2)) {
      z2 <- paste0(z2, ".")
    }
    
    parts <- strsplit(z2, "\\.")[[1]]
    dec <- ifelse(length(parts) >= 2, parts[2], "")
    dec <- gsub("[^0-9]", "", dec)
    dec <- stringr::str_pad(
      dec,
      width = digits_keep,
      side = "right",
      pad = "0"
    )
    
    out[i] <- dec
  }
  
  return(out)
}

extract_single_digit <- function(x, pos = 3, digits_keep = 10) {
  dec <- get_decimal_part(x, digits_keep = digits_keep)
  d <- substr(dec, pos, pos)
  d[d == ""] <- NA_character_
  return(as.integer(d))
}

extract_digit_group <- function(
    x,
    start_pos = 3,
    group_length = 2,
    digits_keep = 10
) {
  dec <- get_decimal_part(x, digits_keep = digits_keep)
  g <- substr(dec, start_pos, start_pos + group_length - 1)
  g[g == ""] <- NA_character_
  return(g)
}

# ==========================================================
# 4. Single decimal-digit statistics
# ==========================================================

calculate_single_digit_stats <- function(
    x,
    dataset_name,
    digit_pos = 3,
    digits_keep = 10
) {
  digits <- extract_single_digit(
    x,
    pos = digit_pos,
    digits_keep = digits_keep
  )
  
  digits <- digits[!is.na(digits)]
  
  cats <- as.character(0:9)
  obs <- as.numeric(table(factor(digits, levels = 0:9)))
  names(obs) <- cats
  
  n <- sum(obs)
  k <- 10
  expected <- n / k
  
  chi_square <- sum((obs - expected)^2 / expected)
  df <- k - 1
  p_value <- pchisq(chi_square, df = df, lower.tail = FALSE)
  cramers_v <- sqrt(chi_square / (n * (k - 1)))
  
  prop <- obs / n
  prop_nonzero <- prop[prop > 0]
  
  shannon_entropy <- -sum(prop_nonzero * log(prop_nonzero))
  h_star <- shannon_entropy / log(k)
  kl <- sum(prop_nonzero * log(prop_nonzero / (1 / k)))
  
  residual <- (obs - expected) / sqrt(expected)
  
  prop_named <- setNames(prop, cats)
  
  pref_05 <- sum(prop_named[c("0", "5")], na.rm = TRUE)
  pref_258 <- sum(prop_named[c("2", "5", "8")], na.rm = TRUE)
  pref_odd <- sum(prop_named[c("1", "3", "5", "7", "9")], na.rm = TRUE)
  pref_extreme <- sum(prop_named[c("0", "1", "8", "9")], na.rm = TRUE)
  pref_middle <- sum(prop_named[c("3", "4", "5", "6")], na.rm = TRUE)
  
  freq_df <- data.frame(
    Dataset = dataset_name,
    Test = paste0("Single decimal digit ", digit_pos),
    Category = cats,
    Observed = obs,
    Expected = expected,
    Proportion = prop,
    Residual = residual,
    stringsAsFactors = FALSE
  )
  
  summary_df <- data.frame(
    Dataset = dataset_name,
    Test = paste0("Single decimal digit ", digit_pos),
    N = n,
    Categories = k,
    Expected_per_category = expected,
    ChiSquare = chi_square,
    DF = df,
    Pvalue = p_value,
    CramersV = cramers_v,
    ShannonEntropy = shannon_entropy,
    H_star = h_star,
    KL = kl,
    Pref_05 = pref_05,
    Pref_258 = pref_258,
    Pref_odd = pref_odd,
    Pref_extreme = pref_extreme,
    Pref_middle = pref_middle,
    Pref_same_pair = NA_real_,
    Pref_neat_combo = NA_real_,
    Pref_sequential_combo = NA_real_,
    stringsAsFactors = FALSE
  )
  
  return(list(summary = summary_df, frequency = freq_df))
}

# ==========================================================
# 5. Joint decimal-digit statistics
# ==========================================================

calculate_joint_digit_stats <- function(
    x,
    dataset_name,
    start_pos = 3,
    group_length = 2,
    digits_keep = 10
) {
  groups <- extract_digit_group(
    x,
    start_pos = start_pos,
    group_length = group_length,
    digits_keep = digits_keep
  )
  
  groups <- groups[!is.na(groups)]
  
  k <- 10^group_length
  cats <- sprintf(paste0("%0", group_length, "d"), 0:(k - 1))
  
  obs <- as.numeric(table(factor(groups, levels = cats)))
  names(obs) <- cats
  
  n <- sum(obs)
  expected <- n / k
  
  chi_square <- sum((obs - expected)^2 / expected)
  df <- k - 1
  p_value <- pchisq(chi_square, df = df, lower.tail = FALSE)
  cramers_v <- sqrt(chi_square / (n * (k - 1)))
  
  prop <- obs / n
  prop_nonzero <- prop[prop > 0]
  
  shannon_entropy <- -sum(prop_nonzero * log(prop_nonzero))
  h_star <- shannon_entropy / log(k)
  kl <- sum(prop_nonzero * log(prop_nonzero / (1 / k)))
  
  residual <- (obs - expected) / sqrt(expected)
  
  freq_df <- data.frame(
    Dataset = dataset_name,
    Test = paste0(
      "Joint decimal digits ",
      start_pos,
      "-",
      start_pos + group_length - 1
    ),
    Category = cats,
    Observed = obs,
    Expected = expected,
    Proportion = prop,
    Residual = residual,
    stringsAsFactors = FALSE
  )
  
  if (group_length == 2) {
    freq_df$Digit1 <- substr(freq_df$Category, 1, 1)
    freq_df$Digit2 <- substr(freq_df$Category, 2, 2)
    
    same_pair <- sprintf("%d%d", 0:9, 0:9)
    neat_combo <- c("00", "25", "50", "75")
    sequential_combo <- paste0(0:8, 1:9)
    
    prop_named <- setNames(prop, cats)
    
    pref_same_pair <- sum(
      prop_named[intersect(same_pair, cats)],
      na.rm = TRUE
    )
    
    pref_neat_combo <- sum(
      prop_named[intersect(neat_combo, cats)],
      na.rm = TRUE
    )
    
    pref_sequential_combo <- sum(
      prop_named[intersect(sequential_combo, cats)],
      na.rm = TRUE
    )
  } else {
    freq_df$Digit1 <- NA_character_
    freq_df$Digit2 <- NA_character_
    pref_same_pair <- NA_real_
    pref_neat_combo <- NA_real_
    pref_sequential_combo <- NA_real_
  }
  
  summary_df <- data.frame(
    Dataset = dataset_name,
    Test = paste0(
      "Joint decimal digits ",
      start_pos,
      "-",
      start_pos + group_length - 1
    ),
    N = n,
    Categories = k,
    Expected_per_category = expected,
    ChiSquare = chi_square,
    DF = df,
    Pvalue = p_value,
    CramersV = cramers_v,
    ShannonEntropy = shannon_entropy,
    H_star = h_star,
    KL = kl,
    Pref_05 = NA_real_,
    Pref_258 = NA_real_,
    Pref_odd = NA_real_,
    Pref_extreme = NA_real_,
    Pref_middle = NA_real_,
    Pref_same_pair = pref_same_pair,
    Pref_neat_combo = pref_neat_combo,
    Pref_sequential_combo = pref_sequential_combo,
    stringsAsFactors = FALSE
  )
  
  return(list(summary = summary_df, frequency = freq_df))
}

# ==========================================================
# 6. Progressive subsampling analysis
# ==========================================================

generate_sample_sizes <- function(N, min_n = 30, n_points = 12) {
  if (N <= min_n) {
    return(unique(round(seq(max(5, floor(N / 3)), N, length.out = min(n_points, N)))))
  }
  
  sample_sizes <- unique(round(seq(min_n, N, length.out = n_points)))
  sample_sizes <- sample_sizes[sample_sizes <= N]
  
  if (!N %in% sample_sizes) {
    sample_sizes <- c(sample_sizes, N)
  }
  
  return(unique(sample_sizes))
}

progressive_subsampling <- function(
    x,
    dataset_name,
    start_pos = 3,
    group_length = 2,
    digits_keep = 10,
    B = 1000,
    min_n = 30,
    n_points = 12
) {
  N <- length(x)
  sample_sizes <- generate_sample_sizes(
    N,
    min_n = min_n,
    n_points = n_points
  )
  
  result_list <- list()
  counter <- 1
  
  for (n_i in sample_sizes) {
    B_use <- ifelse(n_i >= N, 1, B)
    
    for (b in seq_len(B_use)) {
      idx <- if (n_i >= N) {
        seq_len(N)
      } else {
        sample(seq_len(N), size = n_i, replace = FALSE)
      }
      
      tmp_x <- x[idx]
      
      tmp_stat <- calculate_joint_digit_stats(
        tmp_x,
        dataset_name = dataset_name,
        start_pos = start_pos,
        group_length = group_length,
        digits_keep = digits_keep
      )$summary
      
      result_list[[counter]] <- data.frame(
        Dataset = dataset_name,
        SampleSize = n_i,
        Iteration = b,
        ChiSquare = tmp_stat$ChiSquare,
        DF = tmp_stat$DF,
        Pvalue = tmp_stat$Pvalue,
        CramersV = tmp_stat$CramersV,
        ShannonEntropy = tmp_stat$ShannonEntropy,
        H_star = tmp_stat$H_star,
        KL = tmp_stat$KL,
        stringsAsFactors = FALSE
      )
      
      counter <- counter + 1
    }
  }
  
  return(dplyr::bind_rows(result_list))
}

get_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 2) {
    return(NA_real_)
  }
  
  unname(coef(lm(y[ok] ~ log(x[ok])))[2])
}

# ==========================================================
# 7. Run analysis
# ==========================================================

missing_files <- INPUT_FILES[!file.exists(INPUT_FILES)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "These input files were not found in the working directory:",
      paste(missing_files, collapse = ", ")
    )
  )
}

raw_data_list <- list()

for (i in seq_along(INPUT_FILES)) {
  raw_data_list[[INPUT_NAMES[i]]] <- read_single_column_txt(INPUT_FILES[i])
}

single_summary_list <- list()
single_freq_list <- list()
joint_summary_list <- list()
joint_freq_list <- list()
progressive_raw_list <- list()

for (dataset_name in names(raw_data_list)) {
  x <- raw_data_list[[dataset_name]]
  
  message("Processing dataset: ", dataset_name, " | n = ", length(x))
  
  single_res <- calculate_single_digit_stats(
    x,
    dataset_name = dataset_name,
    digit_pos = SINGLE_DIGIT_POS,
    digits_keep = DIGITS_KEEP
  )
  
  joint_res <- calculate_joint_digit_stats(
    x,
    dataset_name = dataset_name,
    start_pos = JOINT_START_POS,
    group_length = JOINT_LENGTH,
    digits_keep = DIGITS_KEEP
  )
  
  prog_res <- progressive_subsampling(
    x,
    dataset_name = dataset_name,
    start_pos = JOINT_START_POS,
    group_length = JOINT_LENGTH,
    digits_keep = DIGITS_KEEP,
    B = B_PROGRESSIVE,
    min_n = MIN_PROGRESSIVE_N,
    n_points = N_PROGRESSIVE_POINTS
  )
  
  single_summary_list[[dataset_name]] <- single_res$summary
  single_freq_list[[dataset_name]] <- single_res$frequency
  
  joint_summary_list[[dataset_name]] <- joint_res$summary
  joint_freq_list[[dataset_name]] <- joint_res$frequency
  
  progressive_raw_list[[dataset_name]] <- prog_res
}

single_summary_df <- dplyr::bind_rows(single_summary_list)
single_freq_df <- dplyr::bind_rows(single_freq_list)

joint_summary_df <- dplyr::bind_rows(joint_summary_list)
joint_freq_df <- dplyr::bind_rows(joint_freq_list)

full_summary_df <- dplyr::bind_rows(single_summary_df, joint_summary_df)
progressive_raw_df <- dplyr::bind_rows(progressive_raw_list)

progressive_summary_df <- progressive_raw_df %>%
  group_by(Dataset, SampleSize) %>%
  summarise(
    across(
      c(ChiSquare, Pvalue, CramersV, ShannonEntropy, H_star, KL),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        q025 = ~quantile(.x, 0.025, na.rm = TRUE),
        q975 = ~quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    Iterations = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Pvalue_median_safe = pmax(Pvalue_median, .Machine$double.xmin),
    Pvalue_q025_safe = pmax(Pvalue_q025, .Machine$double.xmin),
    Pvalue_q975_safe = pmax(Pvalue_q975, .Machine$double.xmin),
    minus_log10P_median = -log10(Pvalue_median_safe),
    minus_log10P_q025 = -log10(Pvalue_q975_safe),
    minus_log10P_q975 = -log10(Pvalue_q025_safe)
  )

progressive_trend_df <- progressive_summary_df %>%
  group_by(Dataset) %>%
  summarise(
    CramersV_slope = get_slope(SampleSize, CramersV_mean),
    H_star_slope = get_slope(SampleSize, H_star_mean),
    KL_slope = get_slope(SampleSize, KL_mean),
    Pvalue_slope = get_slope(SampleSize, Pvalue_median),
    minus_log10P_slope = get_slope(SampleSize, minus_log10P_median),
    .groups = "drop"
  )

single_top_residuals <- single_freq_df %>%
  group_by(Dataset, Test) %>%
  arrange(desc(abs(Residual)), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

joint_top_residuals <- joint_freq_df %>%
  group_by(Dataset, Test) %>%
  arrange(desc(abs(Residual)), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

# ==========================================================
# 8. Export result tables
# ==========================================================

write.csv(
  full_summary_df,
  file.path(OUTPUT_DIR, "FDRS_full_statistics_summary.csv"),
  row.names = FALSE
)

write.csv(
  single_freq_df,
  file.path(OUTPUT_DIR, "FDRS_single_digit_frequency.csv"),
  row.names = FALSE
)

write.csv(
  joint_freq_df,
  file.path(OUTPUT_DIR, "FDRS_joint_digit_frequency.csv"),
  row.names = FALSE
)

write.csv(
  single_top_residuals,
  file.path(OUTPUT_DIR, "FDRS_single_digit_top_residuals.csv"),
  row.names = FALSE
)

write.csv(
  joint_top_residuals,
  file.path(OUTPUT_DIR, "FDRS_joint_digit_top_residuals.csv"),
  row.names = FALSE
)

write.csv(
  progressive_raw_df,
  file.path(OUTPUT_DIR, "FDRS_progressive_raw_results.csv"),
  row.names = FALSE
)

write.csv(
  progressive_summary_df,
  file.path(OUTPUT_DIR, "FDRS_progressive_summary.csv"),
  row.names = FALSE
)

write.csv(
  progressive_trend_df,
  file.path(OUTPUT_DIR, "FDRS_progressive_trend_slopes.csv"),
  row.names = FALSE
)

# ==========================================================
# 9. Dataset colors
# ==========================================================

dataset_levels <- unique(full_summary_df$Dataset)

dataset_colors <- setNames(
  colorRampPalette(c(color_con, color_025, color_25))(length(dataset_levels)),
  dataset_levels
)

if ("RawData" %in% dataset_levels) {
  dataset_colors["RawData"] <- color_con
}

if ("ErrData" %in% dataset_levels) {
  dataset_colors["ErrData"] <- color_25
}

# ==========================================================
# 10. Figure 1: Single digit frequency
# ==========================================================

single_freq_df$Category <- factor(
  single_freq_df$Category,
  levels = as.character(0:9)
)

p_single <- ggplot(
  single_freq_df,
  aes(x = Category, y = Observed, fill = Dataset)
) +
  geom_col(
    color = "black",
    linewidth = 0.25,
    width = 0.75
  ) +
  geom_hline(
    aes(yintercept = Expected),
    linetype = "dashed",
    linewidth = 0.7
  ) +
  facet_wrap(~ Dataset, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Single decimal-digit distribution",
    subtitle = paste0("Decimal position: ", SINGLE_DIGIT_POS),
    x = "Digit",
    y = "Observed count",
    fill = "Dataset"
  ) +
  theme_fdrs_style()

print(p_single)

save_plot(
  p_single,
  "FDRS_figure_01_single_digit_frequency",
  width = 7,
  height = 6.5
)

# ==========================================================
# 11. Figure 2: Single digit standardized residuals
# ==========================================================

p_single_res <- ggplot(
  single_freq_df,
  aes(x = Category, y = Residual, fill = Residual)
) +
  geom_col(
    color = "black",
    linewidth = 0.25,
    width = 0.75
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = c(-2, 2),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  facet_wrap(~ Dataset, ncol = 1) +
  scale_fill_gradient2(
    low = color_con,
    mid = color_025,
    high = color_25,
    midpoint = 0,
    name = "Residual"
  ) +
  labs(
    title = "Single decimal-digit standardized residuals",
    subtitle = paste0("Decimal position: ", SINGLE_DIGIT_POS),
    x = "Digit",
    y = "Standardized residual"
  ) +
  theme_fdrs_style()

print(p_single_res)

save_plot(
  p_single_res,
  "FDRS_figure_02_single_digit_residuals",
  width = 7,
  height = 6.5
)

# ==========================================================
# 12. Figure 3: Joint digit heatmap
# ==========================================================

if (JOINT_LENGTH == 2) {
  
  joint_heat_df <- joint_freq_df %>%
    mutate(
      Digit1 = factor(Digit1, levels = as.character(0:9)),
      Digit2 = factor(Digit2, levels = as.character(0:9))
    )
  
  p_joint_heat <- ggplot(
    joint_heat_df,
    aes(x = Digit2, y = Digit1, fill = Residual)
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.35
    ) +
    facet_wrap(~ Dataset, ncol = length(dataset_levels)) +
    scale_fill_gradient2(
      low = color_con,
      mid = color_025,
      high = color_25,
      midpoint = 0,
      name = "Residual"
    ) +
    labs(
      title = "Joint distribution heatmap of decimal digits",
      subtitle = paste0(
        "Decimal positions: ",
        JOINT_START_POS,
        "-",
        JOINT_START_POS + JOINT_LENGTH - 1
      ),
      x = "Second digit",
      y = "First digit"
    ) +
    theme_fdrs_style()
  
  print(p_joint_heat)
  
  save_plot(
    p_joint_heat,
    "FDRS_figure_03_joint_digit_heatmap_residuals",
    width = 9,
    height = 5.5
  )
}

# ==========================================================
# 13. Figure 4: Full statistics bubble plot
# ==========================================================

full_long_df <- full_summary_df %>%
  select(Dataset, Test, ChiSquare, Pvalue, CramersV, H_star, KL) %>%
  pivot_longer(
    cols = c(ChiSquare, Pvalue, CramersV, H_star, KL),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  group_by(Metric) %>%
  mutate(Value_scaled = scale01(Value)) %>%
  ungroup() %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("ChiSquare", "Pvalue", "CramersV", "H_star", "KL")
    ),
    Label = format_metric(Value)
  )

p_full_bubble <- ggplot(
  full_long_df,
  aes(x = Metric, y = Dataset)
) +
  geom_point(
    aes(size = Value_scaled, color = Value_scaled),
    alpha = 0.9
  ) +
  geom_text(
    aes(label = Label),
    size = 3.2,
    color = "black",
    vjust = -1.2
  ) +
  facet_wrap(~ Test, ncol = 1) +
  scale_size_continuous(
    range = c(4, 12),
    limits = c(0, 1),
    name = "Relative value"
  ) +
  scale_color_gradientn(
    colors = c(color_con, color_025, color_25),
    limits = c(0, 1),
    name = "Relative value"
  ) +
  labs(
    title = "Full-sample statistical summary",
    x = "Metric",
    y = "Dataset"
  ) +
  theme_fdrs_style() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

print(p_full_bubble)

save_plot(
  p_full_bubble,
  "FDRS_figure_04_full_statistics_bubble",
  width = 8.5,
  height = 7
)

# ==========================================================
# 14. Figure 5: Full statistics faceted point plot
# ==========================================================

p_full_point <- ggplot(
  full_long_df,
  aes(x = Dataset, y = Value, color = Dataset)
) +
  geom_point(size = 4) +
  geom_line(
    aes(group = Test),
    linewidth = 0.7,
    alpha = 0.6
  ) +
  facet_grid(Metric ~ Test, scales = "free_y") +
  scale_color_manual(values = dataset_colors) +
  labs(
    title = "Full-sample statistical metrics",
    x = "",
    y = "Value",
    color = "Dataset"
  ) +
  theme_fdrs_style()

print(p_full_point)

save_plot(
  p_full_point,
  "FDRS_figure_05_full_statistics_facet_point",
  width = 10,
  height = 8
)

# ==========================================================
# 15. Figure 6: Progressive Cramér's V
# ==========================================================

p_prog_v <- ggplot(
  progressive_summary_df,
  aes(
    x = SampleSize,
    y = CramersV_mean,
    color = Dataset,
    fill = Dataset
  )
) +
  geom_ribbon(
    aes(
      ymin = CramersV_q025,
      ymax = CramersV_q975
    ),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  scale_color_manual(values = dataset_colors) +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Progressive subsampling analysis",
    subtitle = "Cramér's V across increasing sample sizes",
    x = "Sample size (n)",
    y = "Cramér's V",
    color = "Dataset",
    fill = "Dataset"
  ) +
  theme_fdrs_style()

print(p_prog_v)

save_plot(
  p_prog_v,
  "FDRS_figure_06_progressive_CramersV",
  width = 7.5,
  height = 5.8
)

# ==========================================================
# 16. Figure 7: Progressive normalized entropy
# ==========================================================

p_prog_h <- ggplot(
  progressive_summary_df,
  aes(
    x = SampleSize,
    y = H_star_mean,
    color = Dataset,
    fill = Dataset
  )
) +
  geom_ribbon(
    aes(
      ymin = H_star_q025,
      ymax = H_star_q975
    ),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  scale_color_manual(values = dataset_colors) +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Progressive subsampling analysis",
    subtitle = "Normalized entropy across increasing sample sizes",
    x = "Sample size (n)",
    y = "Normalized entropy (H*)",
    color = "Dataset",
    fill = "Dataset"
  ) +
  theme_fdrs_style()

print(p_prog_h)

save_plot(
  p_prog_h,
  "FDRS_figure_07_progressive_Hstar",
  width = 7.5,
  height = 5.8
)

# ==========================================================
# 17. Figure 8: Progressive KL divergence
# ==========================================================

p_prog_kl <- ggplot(
  progressive_summary_df,
  aes(
    x = SampleSize,
    y = KL_mean,
    color = Dataset,
    fill = Dataset
  )
) +
  geom_ribbon(
    aes(
      ymin = KL_q025,
      ymax = KL_q975
    ),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  scale_color_manual(values = dataset_colors) +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Progressive subsampling analysis",
    subtitle = "KL divergence across increasing sample sizes",
    x = "Sample size (n)",
    y = "KL divergence",
    color = "Dataset",
    fill = "Dataset"
  ) +
  theme_fdrs_style()

print(p_prog_kl)

save_plot(
  p_prog_kl,
  "FDRS_figure_08_progressive_KL",
  width = 7.5,
  height = 5.8
)

# ==========================================================
# 18. Figure 9: Progressive -log10(P)
# ==========================================================

p_prog_p <- ggplot(
  progressive_summary_df,
  aes(
    x = SampleSize,
    y = minus_log10P_median,
    color = Dataset,
    fill = Dataset
  )
) +
  geom_ribbon(
    aes(
      ymin = minus_log10P_q025,
      ymax = minus_log10P_q975
    ),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.8,
    color = "black"
  ) +
  scale_color_manual(values = dataset_colors) +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Progressive subsampling analysis",
    subtitle = expression(paste("-log"[10], "(P) across increasing sample sizes")),
    x = "Sample size (n)",
    y = expression(paste("-log"[10], "(P)")),
    color = "Dataset",
    fill = "Dataset"
  ) +
  theme_fdrs_style()

print(p_prog_p)

save_plot(
  p_prog_p,
  "FDRS_figure_09_progressive_minus_log10P",
  width = 7.5,
  height = 5.8
)

# ==========================================================
# 19. Figure 10: Combined progressive metrics
# ==========================================================

progressive_long_df <- progressive_summary_df %>%
  transmute(
    Dataset = Dataset,
    SampleSize = SampleSize,
    `Cramér's V` = CramersV_mean,
    `Normalized entropy H*` = H_star_mean,
    `KL divergence` = KL_mean,
    `-log10(P)` = minus_log10P_median
  ) %>%
  pivot_longer(
    cols = c(
      `Cramér's V`,
      `Normalized entropy H*`,
      `KL divergence`,
      `-log10(P)`
    ),
    names_to = "Metric",
    values_to = "Value"
  )

p_prog_combined <- ggplot(
  progressive_long_df,
  aes(x = SampleSize, y = Value, color = Dataset)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.4) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  scale_color_manual(values = dataset_colors) +
  labs(
    title = "Combined progressive subsampling metrics",
    x = "Sample size (n)",
    y = "Metric value",
    color = "Dataset"
  ) +
  theme_fdrs_style()

print(p_prog_combined)

save_plot(
  p_prog_combined,
  "FDRS_figure_10_progressive_combined_metrics",
  width = 9,
  height = 7
)

# ==========================================================
# 20. Console summary
# ==========================================================

cat("\n============================================================\n")
cat("FDRS integrated analysis finished.\n")
cat("Working directory: ", getwd(), "\n", sep = "")
cat("Results saved in folder: ", OUTPUT_DIR, "\n", sep = "")
cat("\nMain output tables:\n")
cat("  - FDRS_full_statistics_summary.csv\n")
cat("  - FDRS_single_digit_frequency.csv\n")
cat("  - FDRS_joint_digit_frequency.csv\n")
cat("  - FDRS_progressive_raw_results.csv\n")
cat("  - FDRS_progressive_summary.csv\n")
cat("  - FDRS_progressive_trend_slopes.csv\n")
cat("\nMain figures:\n")
cat("  - FDRS_figure_01_single_digit_frequency\n")
cat("  - FDRS_figure_02_single_digit_residuals\n")
cat("  - FDRS_figure_03_joint_digit_heatmap_residuals\n")
cat("  - FDRS_figure_04_full_statistics_bubble\n")
cat("  - FDRS_figure_05_full_statistics_facet_point\n")
cat("  - FDRS_figure_06_progressive_CramersV\n")
cat("  - FDRS_figure_07_progressive_Hstar\n")
cat("  - FDRS_figure_08_progressive_KL\n")
cat("  - FDRS_figure_09_progressive_minus_log10P\n")
cat("  - FDRS_figure_10_progressive_combined_metrics\n")
cat("============================================================\n")