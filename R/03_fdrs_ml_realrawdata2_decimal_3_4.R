############################################################
# FDRS integrated machine-learning framework
# Full statistics + progressive subsampling + semi-supervised ML
# Author: Zhuohua Cao
############################################################

# ==========================================================
# 0. Working directory and user configuration
# ==========================================================

# 修改为你的工作目录
setwd("F:/Data Model/Figure/Aov")

set.seed(20260524)

# 待检测数据文件
TARGET_FILES <- c("RealRawData2.txt")
TARGET_NAMES <- c("RealRawData2")

# 输出文件夹
OUTPUT_DIR <- "FDRS_ML_integrated_RealRawData2_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 小数位设置
DIGITS_KEEP <- 10
SINGLE_DIGIT_POSITIONS <- 1:4
PRIMARY_SINGLE_POS <- 3
JOINT_START_POS <- 3
JOINT_LENGTH <- 2

# 半监督模拟训练集设置
N_PER_DATASET <- 255

# 正式论文运行可设为 500、500、500
# 调试时可先设为 100、100、100
N_NORMAL_SETS <- 500
N_IRREGULAR_SETS <- 500
N_MIXED_SETS <- 500

# 机器学习特征中的渐进式抽样设置
# 注意：这个值不宜太大，否则 1500 个训练集会非常慢
B_PROGRESSIVE_ML <- 20
N_PROGRESSIVE_POINTS_ML <- 5
MIN_PROGRESSIVE_N_ML <- 30

# 对 RawData / ErrData 单独输出详细渐进式抽样结果
RUN_TARGET_PROGRESSIVE_REPORT <- TRUE

# 正式论文可设为 1000；如果运行慢，先设为 300
B_PROGRESSIVE_REPORT <- 1000
N_PROGRESSIVE_POINTS_REPORT <- 12
MIN_PROGRESSIVE_N_REPORT <- 30

# 训练集划分比例
TRAIN_RATIO <- 0.8

# 集成模型权重
# supervised risk = RF / Logistic / SVM 平均
# final risk = LAMBDA_SUPERVISED * supervised risk + (1 - LAMBDA_SUPERVISED) * isolation risk
LAMBDA_SUPERVISED <- 0.75

# 风险分级阈值
risk_grade_fun <- function(score) {
  ifelse(
    score < 0.20, "Grade 0: no apparent irregularity",
    ifelse(
      score < 0.40, "Grade 1: mild irregularity",
      ifelse(
        score < 0.60, "Grade 2: moderate irregularity",
        ifelse(
          score < 0.80, "Grade 3: high irregularity",
          "Grade 4: very high irregularity"
        )
      )
    )
  )
}

# ==========================================================
# 1. Load packages
# ==========================================================

pkg_needed <- c(
  "ggplot2", "dplyr", "tidyr", "stringr", "scales",
  "ranger", "glmnet", "e1071", "isotree", "pROC"
)

options(repos = c(CRAN = "https://cloud.r-project.org"))

for (pkg in pkg_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# ==========================================================
# 2. Color palette and figure theme
# ==========================================================

color_con <- "#A9D1E8"
color_025 <- "#F6C1C3"
color_25  <- "#E68589"

theme_fdrs_style <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text.x = element_text(color = "black", size = base_size - 2, face = "bold"),
      axis.text.y = element_text(color = "black", size = base_size - 2),
      axis.title.x = element_text(color = "black", size = base_size, face = "bold"),
      axis.title.y = element_text(color = "black", size = base_size, face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
      strip.text = element_text(color = "black", face = "bold", size = base_size - 1),
      legend.title = element_text(color = "black", size = base_size - 1, face = "bold"),
      legend.text = element_text(color = "black", size = base_size - 2),
      plot.title = element_text(color = "black", size = base_size + 1, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(color = "black", size = base_size - 1, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
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

# ==========================================================
# 3. Basic helper functions
# ==========================================================

safe_p <- function(p) {
  pmax(p, .Machine$double.xmin)
}

scale01 <- function(v) {
  if (all(is.na(v))) return(rep(NA_real_, length(v)))
  rng <- range(v, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(0.5, length(v)))
  }
  scales::rescale(v, to = c(0, 1), from = rng)
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

get_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2) return(NA_real_)
  unname(coef(lm(y[ok] ~ log(x[ok])))[2])
}

# ==========================================================
# 4. Data reading and decimal extraction
# ==========================================================

read_single_column_txt <- function(file) {
  raw <- readLines(file, warn = FALSE, encoding = "UTF-8")
  raw <- trimws(raw)
  raw <- raw[raw != ""]
  
  raw <- sapply(strsplit(raw, "[,\t ]+"), `[`, 1)
  raw <- trimws(raw)
  
  num_check <- suppressWarnings(as.numeric(raw))
  raw <- raw[!is.na(num_check)]
  
  as.character(raw)
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
    dec <- stringr::str_pad(dec, width = digits_keep, side = "right", pad = "0")
    out[i] <- dec
  }
  
  out
}

extract_single_digit <- function(x, pos = 3, digits_keep = 10) {
  dec <- get_decimal_part(x, digits_keep = digits_keep)
  d <- substr(dec, pos, pos)
  d[d == ""] <- NA_character_
  as.integer(d)
}

extract_digit_group <- function(x, start_pos = 3, group_length = 2, digits_keep = 10) {
  dec <- get_decimal_part(x, digits_keep = digits_keep)
  g <- substr(dec, start_pos, start_pos + group_length - 1)
  g[g == ""] <- NA_character_
  g
}

# ==========================================================
# 5. Distribution statistics
# ==========================================================

calc_stats_from_counts <- function(obs, category_names) {
  obs <- as.numeric(obs)
  names(obs) <- category_names
  
  n <- sum(obs)
  k <- length(obs)
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
  
  list(
    n = n,
    k = k,
    expected = expected,
    chi_square = chi_square,
    df = df,
    p_value = p_value,
    cramers_v = cramers_v,
    prop = prop,
    shannon_entropy = shannon_entropy,
    h_star = h_star,
    kl = kl,
    residual = residual
  )
}

calculate_single_digit_stats <- function(x, dataset_name, digit_pos = 3, digits_keep = 10) {
  digits <- extract_single_digit(x, pos = digit_pos, digits_keep = digits_keep)
  digits <- digits[!is.na(digits)]
  
  cats <- as.character(0:9)
  obs <- as.numeric(table(factor(digits, levels = 0:9)))
  names(obs) <- cats
  
  s <- calc_stats_from_counts(obs, cats)
  prop_named <- setNames(s$prop, cats)
  
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
    Expected = s$expected,
    Proportion = s$prop,
    Residual = s$residual,
    stringsAsFactors = FALSE
  )
  
  summary_df <- data.frame(
    Dataset = dataset_name,
    Test = paste0("Single decimal digit ", digit_pos),
    N = s$n,
    Categories = s$k,
    Expected_per_category = s$expected,
    ChiSquare = s$chi_square,
    DF = s$df,
    Pvalue = s$p_value,
    minus_log10P = -log10(safe_p(s$p_value)),
    CramersV = s$cramers_v,
    ShannonEntropy = s$shannon_entropy,
    H_star = s$h_star,
    KL = s$kl,
    Pref_05 = pref_05,
    Pref_258 = pref_258,
    Pref_odd = pref_odd,
    Pref_extreme = pref_extreme,
    Pref_middle = pref_middle,
    Pref_same_pair = NA_real_,
    Pref_neat_combo = NA_real_,
    Pref_sequential_combo = NA_real_,
    MaxAbsResidual = max(abs(s$residual), na.rm = TRUE),
    NonzeroCategories = sum(obs > 0),
    stringsAsFactors = FALSE
  )
  
  list(summary = summary_df, frequency = freq_df)
}

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
  
  s <- calc_stats_from_counts(obs, cats)
  
  freq_df <- data.frame(
    Dataset = dataset_name,
    Test = paste0("Joint decimal digits ", start_pos, "-", start_pos + group_length - 1),
    Category = cats,
    Observed = obs,
    Expected = s$expected,
    Proportion = s$prop,
    Residual = s$residual,
    stringsAsFactors = FALSE
  )
  
  if (group_length == 2) {
    freq_df$Digit1 <- substr(freq_df$Category, 1, 1)
    freq_df$Digit2 <- substr(freq_df$Category, 2, 2)
    
    same_pair <- sprintf("%d%d", 0:9, 0:9)
    neat_combo <- c("00", "25", "50", "75")
    sequential_combo <- paste0(0:8, 1:9)
    
    prop_named <- setNames(s$prop, cats)
    
    pref_same_pair <- sum(prop_named[intersect(same_pair, cats)], na.rm = TRUE)
    pref_neat_combo <- sum(prop_named[intersect(neat_combo, cats)], na.rm = TRUE)
    pref_sequential_combo <- sum(prop_named[intersect(sequential_combo, cats)], na.rm = TRUE)
  } else {
    freq_df$Digit1 <- NA_character_
    freq_df$Digit2 <- NA_character_
    pref_same_pair <- NA_real_
    pref_neat_combo <- NA_real_
    pref_sequential_combo <- NA_real_
  }
  
  summary_df <- data.frame(
    Dataset = dataset_name,
    Test = paste0("Joint decimal digits ", start_pos, "-", start_pos + group_length - 1),
    N = s$n,
    Categories = s$k,
    Expected_per_category = s$expected,
    ChiSquare = s$chi_square,
    DF = s$df,
    Pvalue = s$p_value,
    minus_log10P = -log10(safe_p(s$p_value)),
    CramersV = s$cramers_v,
    ShannonEntropy = s$shannon_entropy,
    H_star = s$h_star,
    KL = s$kl,
    Pref_05 = NA_real_,
    Pref_258 = NA_real_,
    Pref_odd = NA_real_,
    Pref_extreme = NA_real_,
    Pref_middle = NA_real_,
    Pref_same_pair = pref_same_pair,
    Pref_neat_combo = pref_neat_combo,
    Pref_sequential_combo = pref_sequential_combo,
    MaxAbsResidual = max(abs(s$residual), na.rm = TRUE),
    NonzeroCategories = sum(obs > 0),
    stringsAsFactors = FALSE
  )
  
  list(summary = summary_df, frequency = freq_df)
}

# ==========================================================
# 6. Progressive subsampling
# ==========================================================

generate_sample_sizes <- function(N, min_n = 30, n_points = 5) {
  if (N <= min_n) {
    return(unique(round(seq(max(5, floor(N / 3)), N, length.out = min(n_points, N)))))
  }
  
  sample_sizes <- unique(round(seq(min_n, N, length.out = n_points)))
  sample_sizes <- sample_sizes[sample_sizes <= N]
  
  if (!N %in% sample_sizes) {
    sample_sizes <- c(sample_sizes, N)
  }
  
  unique(sample_sizes)
}

calc_joint_stats_from_groups <- function(groups, group_length = 2) {
  k <- 10^group_length
  cats <- sprintf(paste0("%0", group_length, "d"), 0:(k - 1))
  obs <- as.numeric(table(factor(groups, levels = cats)))
  calc_stats_from_counts(obs, cats)
}

progressive_subsampling <- function(
    x,
    dataset_name,
    start_pos = 3,
    group_length = 2,
    digits_keep = 10,
    B = 20,
    min_n = 30,
    n_points = 5
) {
  groups_all <- extract_digit_group(
    x,
    start_pos = start_pos,
    group_length = group_length,
    digits_keep = digits_keep
  )
  
  groups_all <- groups_all[!is.na(groups_all)]
  N <- length(groups_all)
  sample_sizes <- generate_sample_sizes(N, min_n = min_n, n_points = n_points)
  
  result_list <- vector("list", length = length(sample_sizes) * B)
  counter <- 1
  
  for (n_i in sample_sizes) {
    B_use <- ifelse(n_i >= N, 1, B)
    
    for (b in seq_len(B_use)) {
      idx <- if (n_i >= N) {
        seq_len(N)
      } else {
        sample(seq_len(N), size = n_i, replace = FALSE)
      }
      
      s <- calc_joint_stats_from_groups(groups_all[idx], group_length = group_length)
      
      result_list[[counter]] <- data.frame(
        Dataset = dataset_name,
        SampleSize = n_i,
        Iteration = b,
        ChiSquare = s$chi_square,
        DF = s$df,
        Pvalue = s$p_value,
        minus_log10P = -log10(safe_p(s$p_value)),
        CramersV = s$cramers_v,
        ShannonEntropy = s$shannon_entropy,
        H_star = s$h_star,
        KL = s$kl,
        stringsAsFactors = FALSE
      )
      
      counter <- counter + 1
    }
  }
  
  dplyr::bind_rows(result_list[seq_len(counter - 1)])
}

summarize_progressive <- function(progressive_raw_df) {
  progressive_raw_df %>%
    group_by(Dataset, SampleSize) %>%
    summarise(
      ChiSquare_mean = mean(ChiSquare, na.rm = TRUE),
      ChiSquare_median = median(ChiSquare, na.rm = TRUE),
      ChiSquare_q025 = quantile(ChiSquare, 0.025, na.rm = TRUE),
      ChiSquare_q975 = quantile(ChiSquare, 0.975, na.rm = TRUE),
      
      Pvalue_mean = mean(Pvalue, na.rm = TRUE),
      Pvalue_median = median(Pvalue, na.rm = TRUE),
      Pvalue_q025 = quantile(Pvalue, 0.025, na.rm = TRUE),
      Pvalue_q975 = quantile(Pvalue, 0.975, na.rm = TRUE),
      
      minus_log10P_mean = mean(minus_log10P, na.rm = TRUE),
      minus_log10P_median = median(minus_log10P, na.rm = TRUE),
      minus_log10P_q025 = quantile(minus_log10P, 0.025, na.rm = TRUE),
      minus_log10P_q975 = quantile(minus_log10P, 0.975, na.rm = TRUE),
      
      CramersV_mean = mean(CramersV, na.rm = TRUE),
      CramersV_median = median(CramersV, na.rm = TRUE),
      CramersV_q025 = quantile(CramersV, 0.025, na.rm = TRUE),
      CramersV_q975 = quantile(CramersV, 0.975, na.rm = TRUE),
      
      H_star_mean = mean(H_star, na.rm = TRUE),
      H_star_median = median(H_star, na.rm = TRUE),
      H_star_q025 = quantile(H_star, 0.025, na.rm = TRUE),
      H_star_q975 = quantile(H_star, 0.975, na.rm = TRUE),
      
      KL_mean = mean(KL, na.rm = TRUE),
      KL_median = median(KL, na.rm = TRUE),
      KL_q025 = quantile(KL, 0.025, na.rm = TRUE),
      KL_q975 = quantile(KL, 0.975, na.rm = TRUE),
      
      Iterations = n(),
      .groups = "drop"
    )
}

progressive_features <- function(
    x,
    prefix = "prog",
    start_pos = 3,
    group_length = 2,
    digits_keep = 10,
    B = 20,
    min_n = 30,
    n_points = 5
) {
  raw <- progressive_subsampling(
    x,
    dataset_name = "tmp",
    start_pos = start_pos,
    group_length = group_length,
    digits_keep = digits_keep,
    B = B,
    min_n = min_n,
    n_points = n_points
  )
  
  sm <- summarize_progressive(raw)
  
  last_row <- sm[which.max(sm$SampleSize), ]
  first_row <- sm[which.min(sm$SampleSize), ]
  
  out <- c(
    first_ChiSquare = first_row$ChiSquare_mean,
    final_ChiSquare = last_row$ChiSquare_mean,
    mean_ChiSquare = mean(sm$ChiSquare_mean, na.rm = TRUE),
    slope_ChiSquare = get_slope(sm$SampleSize, sm$ChiSquare_mean),
    
    first_Pvalue = first_row$Pvalue_median,
    final_Pvalue = last_row$Pvalue_median,
    mean_Pvalue = mean(sm$Pvalue_median, na.rm = TRUE),
    slope_Pvalue = get_slope(sm$SampleSize, sm$Pvalue_median),
    
    first_minus_log10P = first_row$minus_log10P_median,
    final_minus_log10P = last_row$minus_log10P_median,
    mean_minus_log10P = mean(sm$minus_log10P_median, na.rm = TRUE),
    slope_minus_log10P = get_slope(sm$SampleSize, sm$minus_log10P_median),
    
    first_CramersV = first_row$CramersV_mean,
    final_CramersV = last_row$CramersV_mean,
    mean_CramersV = mean(sm$CramersV_mean, na.rm = TRUE),
    max_CramersV = max(sm$CramersV_mean, na.rm = TRUE),
    slope_CramersV = get_slope(sm$SampleSize, sm$CramersV_mean),
    
    first_H_star = first_row$H_star_mean,
    final_H_star = last_row$H_star_mean,
    mean_H_star = mean(sm$H_star_mean, na.rm = TRUE),
    min_H_star = min(sm$H_star_mean, na.rm = TRUE),
    slope_H_star = get_slope(sm$SampleSize, sm$H_star_mean),
    
    first_KL = first_row$KL_mean,
    final_KL = last_row$KL_mean,
    mean_KL = mean(sm$KL_mean, na.rm = TRUE),
    max_KL = max(sm$KL_mean, na.rm = TRUE),
    slope_KL = get_slope(sm$SampleSize, sm$KL_mean)
  )
  
  names(out) <- paste0(prefix, "_", names(out))
  out
}

# ==========================================================
# 7. Dataset-level feature extraction for machine learning
# ==========================================================

extract_dataset_features <- function(
    x,
    dataset_name = "Dataset",
    digits_keep = 10,
    single_positions = 1:4,
    primary_single_pos = 3,
    joint_start_pos = 3,
    joint_length = 2,
    B_progressive = 20,
    n_progressive_points = 5,
    min_progressive_n = 30
) {
  feature_vec <- c()
  
  # Basic sample information
  feature_vec["N"] <- length(x)
  
  # Single-digit full statistics and digit frequencies
  for (pos in single_positions) {
    res <- calculate_single_digit_stats(
      x,
      dataset_name = dataset_name,
      digit_pos = pos,
      digits_keep = digits_keep
    )
    
    s <- res$summary
    
    feature_vec[paste0("single_pos", pos, "_ChiSquare")] <- s$ChiSquare
    feature_vec[paste0("single_pos", pos, "_Pvalue")] <- s$Pvalue
    feature_vec[paste0("single_pos", pos, "_minus_log10P")] <- s$minus_log10P
    feature_vec[paste0("single_pos", pos, "_CramersV")] <- s$CramersV
    feature_vec[paste0("single_pos", pos, "_H_star")] <- s$H_star
    feature_vec[paste0("single_pos", pos, "_KL")] <- s$KL
    feature_vec[paste0("single_pos", pos, "_MaxAbsResidual")] <- s$MaxAbsResidual
    feature_vec[paste0("single_pos", pos, "_Pref_05")] <- s$Pref_05
    feature_vec[paste0("single_pos", pos, "_Pref_258")] <- s$Pref_258
    feature_vec[paste0("single_pos", pos, "_Pref_odd")] <- s$Pref_odd
    feature_vec[paste0("single_pos", pos, "_Pref_extreme")] <- s$Pref_extreme
    feature_vec[paste0("single_pos", pos, "_Pref_middle")] <- s$Pref_middle
    
    freq <- res$frequency
    for (d in 0:9) {
      val <- freq$Proportion[freq$Category == as.character(d)]
      feature_vec[paste0("single_pos", pos, "_freq_", d)] <- val
    }
  }
  
  # Primary joint full statistics
  joint_res <- calculate_joint_digit_stats(
    x,
    dataset_name = dataset_name,
    start_pos = joint_start_pos,
    group_length = joint_length,
    digits_keep = digits_keep
  )
  
  js <- joint_res$summary
  
  feature_vec["joint_34_ChiSquare"] <- js$ChiSquare
  feature_vec["joint_34_Pvalue"] <- js$Pvalue
  feature_vec["joint_34_minus_log10P"] <- js$minus_log10P
  feature_vec["joint_34_CramersV"] <- js$CramersV
  feature_vec["joint_34_H_star"] <- js$H_star
  feature_vec["joint_34_KL"] <- js$KL
  feature_vec["joint_34_MaxAbsResidual"] <- js$MaxAbsResidual
  feature_vec["joint_34_NonzeroCategories"] <- js$NonzeroCategories
  feature_vec["joint_34_Pref_same_pair"] <- js$Pref_same_pair
  feature_vec["joint_34_Pref_neat_combo"] <- js$Pref_neat_combo
  feature_vec["joint_34_Pref_sequential_combo"] <- js$Pref_sequential_combo
  
  joint_freq <- joint_res$frequency
  for (cat_i in joint_freq$Category) {
    feature_vec[paste0("joint_34_freq_", cat_i)] <- joint_freq$Proportion[joint_freq$Category == cat_i]
  }
  
  # Progressive subsampling features
  prog_feat <- progressive_features(
    x,
    prefix = "prog_34",
    start_pos = joint_start_pos,
    group_length = joint_length,
    digits_keep = digits_keep,
    B = B_progressive,
    min_n = min_progressive_n,
    n_points = n_progressive_points
  )
  
  feature_vec <- c(feature_vec, prog_feat)
  
  feature_vec[!is.finite(feature_vec)] <- NA_real_
  feature_vec
}

# ==========================================================
# 8. Simulated dataset generation
# ==========================================================

format_values <- function(v, digits = 4) {
  v <- pmax(pmin(v, 9.9999), 0)
  formatC(v, format = "f", digits = digits)
}

generate_normal_dataset <- function(n = 255) {
  type <- sample(c("uniform", "normal", "lognormal", "gamma"), 1)
  
  if (type == "uniform") {
    v <- runif(n, 0, 9.9999)
  } else if (type == "normal") {
    v <- rnorm(n, mean = runif(1, 3, 7), sd = runif(1, 0.8, 2.0))
    v <- pmax(pmin(v, 9.9999), 0)
  } else if (type == "lognormal") {
    v <- rlnorm(n, meanlog = runif(1, 0.5, 1.4), sdlog = runif(1, 0.2, 0.7))
    v <- pmax(pmin(v, 9.9999), 0)
  } else {
    v <- rgamma(n, shape = runif(1, 1.5, 5), rate = runif(1, 0.4, 1.2))
    v <- pmax(pmin(v, 9.9999), 0)
  }
  
  format_values(v, digits = 4)
}

generate_irregular_dataset <- function(n = 255) {
  int_part <- sample(0:9, n, replace = TRUE)
  d1 <- sample(0:9, n, replace = TRUE)
  d2 <- sample(0:9, n, replace = TRUE)
  
  mode <- sample(c("258", "05", "neat_pair", "same_pair", "mixed_pref"), 1)
  
  if (mode == "258") {
    prob <- rep(0.06, 10)
    prob[c(3, 6, 9)] <- c(0.18, 0.22, 0.18)
    prob <- prob / sum(prob)
    d3 <- sample(0:9, n, replace = TRUE, prob = prob)
    d4 <- sample(0:9, n, replace = TRUE, prob = prob)
  } else if (mode == "05") {
    prob <- rep(0.07, 10)
    prob[c(1, 6)] <- c(0.23, 0.21)
    prob <- prob / sum(prob)
    d3 <- sample(0:9, n, replace = TRUE, prob = prob)
    d4 <- sample(0:9, n, replace = TRUE, prob = prob)
  } else if (mode == "neat_pair") {
    pair_pool <- c(
      "00", "25", "50", "75", "55", "88",
      sprintf("%02d", 0:99)
    )
    pair_prob <- c(rep(0.06, 6), rep((1 - 0.06 * 6) / 100, 100))
    pair <- sample(pair_pool, n, replace = TRUE, prob = pair_prob)
    d3 <- as.integer(substr(pair, 1, 1))
    d4 <- as.integer(substr(pair, 2, 2))
  } else if (mode == "same_pair") {
    same_pair <- sprintf("%d%d", 0:9, 0:9)
    pair_pool <- c(same_pair, sprintf("%02d", 0:99))
    pair_prob <- c(rep(0.035, 10), rep((1 - 0.035 * 10) / 100, 100))
    pair <- sample(pair_pool, n, replace = TRUE, prob = pair_prob)
    d3 <- as.integer(substr(pair, 1, 1))
    d4 <- as.integer(substr(pair, 2, 2))
  } else {
    prob3 <- rep(0.075, 10)
    prob3[c(3, 6, 9)] <- c(0.15, 0.18, 0.14)
    prob3 <- prob3 / sum(prob3)
    
    prob4 <- rep(0.08, 10)
    prob4[c(1, 6)] <- c(0.18, 0.18)
    prob4 <- prob4 / sum(prob4)
    
    d3 <- sample(0:9, n, replace = TRUE, prob = prob3)
    d4 <- sample(0:9, n, replace = TRUE, prob = prob4)
  }
  
  paste0(int_part, ".", d1, d2, d3, d4)
}

generate_mixed_dataset <- function(n = 255) {
  normal <- generate_normal_dataset(n)
  irregular <- generate_irregular_dataset(n)
  rho <- runif(1, 0.15, 0.75)
  replace_flag <- rbinom(n, size = 1, prob = rho) == 1
  out <- normal
  out[replace_flag] <- irregular[replace_flag]
  out
}

# ==========================================================
# 9. Build semi-supervised training feature matrix
# ==========================================================

make_training_features <- function() {
  total_sets <- N_NORMAL_SETS + N_IRREGULAR_SETS + N_MIXED_SETS
  
  feature_list <- vector("list", total_sets)
  meta_list <- vector("list", total_sets)
  
  pb <- txtProgressBar(min = 0, max = total_sets, style = 3)
  counter <- 1
  
  for (i in seq_len(N_NORMAL_SETS)) {
    x <- generate_normal_dataset(N_PER_DATASET)
    feature_list[[counter]] <- extract_dataset_features(
      x,
      dataset_name = paste0("Normal_", i),
      B_progressive = B_PROGRESSIVE_ML,
      n_progressive_points = N_PROGRESSIVE_POINTS_ML,
      min_progressive_n = MIN_PROGRESSIVE_N_ML
    )
    meta_list[[counter]] <- data.frame(
      DatasetID = paste0("Normal_", i),
      DataType = "Simulated normal",
      Label = 0,
      stringsAsFactors = FALSE
    )
    setTxtProgressBar(pb, counter)
    counter <- counter + 1
  }
  
  for (i in seq_len(N_IRREGULAR_SETS)) {
    x <- generate_irregular_dataset(N_PER_DATASET)
    feature_list[[counter]] <- extract_dataset_features(
      x,
      dataset_name = paste0("Irregular_", i),
      B_progressive = B_PROGRESSIVE_ML,
      n_progressive_points = N_PROGRESSIVE_POINTS_ML,
      min_progressive_n = MIN_PROGRESSIVE_N_ML
    )
    meta_list[[counter]] <- data.frame(
      DatasetID = paste0("Irregular_", i),
      DataType = "Simulated digit-preference irregular",
      Label = 1,
      stringsAsFactors = FALSE
    )
    setTxtProgressBar(pb, counter)
    counter <- counter + 1
  }
  
  for (i in seq_len(N_MIXED_SETS)) {
    x <- generate_mixed_dataset(N_PER_DATASET)
    feature_list[[counter]] <- extract_dataset_features(
      x,
      dataset_name = paste0("Mixed_", i),
      B_progressive = B_PROGRESSIVE_ML,
      n_progressive_points = N_PROGRESSIVE_POINTS_ML,
      min_progressive_n = MIN_PROGRESSIVE_N_ML
    )
    meta_list[[counter]] <- data.frame(
      DatasetID = paste0("Mixed_", i),
      DataType = "Mixed-contamination",
      Label = 1,
      stringsAsFactors = FALSE
    )
    setTxtProgressBar(pb, counter)
    counter <- counter + 1
  }
  
  close(pb)
  
  all_names <- unique(unlist(lapply(feature_list, names)))
  
  feature_mat <- do.call(
    rbind,
    lapply(feature_list, function(v) {
      out <- rep(NA_real_, length(all_names))
      names(out) <- all_names
      out[names(v)] <- v
      out
    })
  )
  
  feature_df <- as.data.frame(feature_mat)
  meta_df <- dplyr::bind_rows(meta_list)
  
  cbind(meta_df, feature_df)
}

message("Building semi-supervised training feature matrix...")
training_df <- make_training_features()

write.csv(
  training_df,
  file.path(OUTPUT_DIR, "FDRS_ML_training_feature_matrix_raw.csv"),
  row.names = FALSE
)

# ==========================================================
# 10. Feature preprocessing
# ==========================================================

feature_cols <- setdiff(colnames(training_df), c("DatasetID", "DataType", "Label"))

X_raw <- training_df[, feature_cols, drop = FALSE]
X_raw <- as.data.frame(lapply(X_raw, as.numeric))

# Remove columns with all NA
all_na_cols <- names(X_raw)[sapply(X_raw, function(z) all(is.na(z)))]
if (length(all_na_cols) > 0) {
  X_raw <- X_raw[, setdiff(names(X_raw), all_na_cols), drop = FALSE]
}

# Replace infinite with NA
X_raw[] <- lapply(X_raw, function(z) {
  z[!is.finite(z)] <- NA_real_
  z
})

# Median imputation
feature_medians <- sapply(X_raw, function(z) median(z, na.rm = TRUE))
feature_medians[!is.finite(feature_medians)] <- 0

for (nm in names(X_raw)) {
  X_raw[[nm]][is.na(X_raw[[nm]])] <- feature_medians[nm]
}

# Remove zero-variance columns
feature_sd <- sapply(X_raw, sd, na.rm = TRUE)
keep_cols <- names(feature_sd)[is.finite(feature_sd) & feature_sd > 0]
X_raw <- X_raw[, keep_cols, drop = FALSE]

# Standardization
feature_means <- sapply(X_raw, mean, na.rm = TRUE)
feature_sds <- sapply(X_raw, sd, na.rm = TRUE)
feature_sds[!is.finite(feature_sds) | feature_sds == 0] <- 1

X_scaled <- as.data.frame(scale(X_raw, center = feature_means, scale = feature_sds))

y <- training_df$Label

processed_training_df <- cbind(
  training_df[, c("DatasetID", "DataType", "Label")],
  X_scaled
)

write.csv(
  processed_training_df,
  file.path(OUTPUT_DIR, "FDRS_ML_training_feature_matrix_processed.csv"),
  row.names = FALSE
)

# ==========================================================
# 11. Train-test split
# ==========================================================

idx0 <- which(y == 0)
idx1 <- which(y == 1)

train0 <- sample(idx0, size = floor(length(idx0) * TRAIN_RATIO))
train1 <- sample(idx1, size = floor(length(idx1) * TRAIN_RATIO))

train_idx <- sort(c(train0, train1))
test_idx <- setdiff(seq_along(y), train_idx)

X_train <- X_scaled[train_idx, , drop = FALSE]
X_test  <- X_scaled[test_idx, , drop = FALSE]

y_train <- y[train_idx]
y_test  <- y[test_idx]

train_model_df <- data.frame(Label = factor(y_train, levels = c(0, 1)), X_train)
test_model_df  <- data.frame(Label = factor(y_test, levels = c(0, 1)), X_test)

# ==========================================================
# 12. Model training
# ==========================================================

message("Training Random Forest...")
rf_model <- ranger(
  Label ~ .,
  data = train_model_df,
  probability = TRUE,
  importance = "impurity",
  num.trees = 500,
  seed = 20260524
)

message("Training Elastic-net Logistic Regression...")
x_train_mat <- as.matrix(X_train)
x_test_mat <- as.matrix(X_test)

elastic_model <- cv.glmnet(
  x = x_train_mat,
  y = y_train,
  family = "binomial",
  alpha = 0.5,
  nfolds = 5,
  type.measure = "auc"
)

message("Training SVM radial...")
svm_model <- e1071::svm(
  Label ~ .,
  data = train_model_df,
  kernel = "radial",
  probability = TRUE,
  scale = FALSE
)

message("Training Isolation Forest...")
normal_train_x <- X_train[y_train == 0, , drop = FALSE]

if_model <- isotree::isolation.forest(
  normal_train_x,
  ntrees = 500,
  sample_size = min(256, nrow(normal_train_x)),
  ndim = 1,
  seed = 20260524
)

# ==========================================================
# 13. Prediction helper functions
# ==========================================================

as_prob_vector <- function(x) {
  if (is.null(x)) {
    return(numeric(0))
  }
  
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  
  if (is.matrix(x)) {
    if (ncol(x) == 1) {
      x <- x[, 1]
    } else {
      x <- x[, ncol(x)]
    }
  }
  
  if (is.list(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  return(x)
}

predict_rf <- function(model, newx_df) {
  pred <- predict(model, data = newx_df)$predictions
  
  if (is.matrix(pred) || is.data.frame(pred)) {
    if ("1" %in% colnames(pred)) {
      prob <- pred[, "1"]
    } else {
      prob <- pred[, ncol(pred)]
    }
  } else {
    prob <- pred
  }
  
  as_prob_vector(prob)
}

predict_elastic <- function(model, newx_df) {
  prob <- predict(
    model,
    newx = as.matrix(newx_df),
    s = "lambda.min",
    type = "response"
  )
  
  as_prob_vector(prob)
}

predict_svm <- function(model, newx_df) {
  pred <- predict(
    model,
    newdata = newx_df,
    probability = TRUE
  )
  
  prob <- attr(pred, "probabilities")
  
  if (is.null(prob)) {
    stop("SVM probability output is NULL. Please make sure probability = TRUE was used in svm().")
  }
  
  if ("1" %in% colnames(prob)) {
    prob <- prob[, "1"]
  } else {
    prob <- prob[, ncol(prob)]
  }
  
  as_prob_vector(prob)
}

predict_if_score <- function(model, newx_df, reference_scores = NULL) {
  score <- predict(
    model,
    newdata = newx_df,
    type = "score"
  )
  
  score <- as_prob_vector(score)
  
  if (is.null(reference_scores)) {
    return(scale01(score))
  }
  
  reference_scores <- as_prob_vector(reference_scores)
  
  q01 <- quantile(reference_scores, 0.01, na.rm = TRUE)
  q99 <- quantile(reference_scores, 0.99, na.rm = TRUE)
  
  if (!is.finite(q01) || !is.finite(q99) || q99 == q01) {
    return(rep(0.5, length(score)))
  }
  
  score_clip <- pmin(pmax(score, q01), q99)
  risk <- (score_clip - q01) / (q99 - q01)
  risk[!is.finite(risk)] <- 0.5
  
  return(risk)
}

# Test predictions
rf_test_prob <- predict_rf(rf_model, X_test)
elastic_test_prob <- predict_elastic(elastic_model, X_test)
svm_test_prob <- predict_svm(svm_model, X_test)
if_test_prob <- predict_if_score(if_model, X_test, reference_scores = if_train_raw)

supervised_test_prob <- rowMeans(
  cbind(rf_test_prob, elastic_test_prob, svm_test_prob),
  na.rm = TRUE
)

ensemble_test_prob <- LAMBDA_SUPERVISED * supervised_test_prob +
  (1 - LAMBDA_SUPERVISED) * if_test_prob

# ==========================================================
# 14. Model evaluation
# ==========================================================
get_best_threshold <- function(y_true, prob) {
  y_true <- as_prob_vector(y_true)
  prob <- as_prob_vector(prob)
  
  ok <- is.finite(y_true) & is.finite(prob)
  y_true <- y_true[ok]
  prob <- prob[ok]
  
  thresholds <- sort(unique(prob))
  
  if (length(thresholds) == 0) {
    return(0.5)
  }
  
  youden_values <- sapply(thresholds, function(thr) {
    pred <- ifelse(prob >= thr, 1, 0)
    
    TP <- sum(pred == 1 & y_true == 1)
    TN <- sum(pred == 0 & y_true == 0)
    FP <- sum(pred == 1 & y_true == 0)
    FN <- sum(pred == 0 & y_true == 1)
    
    sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
    specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
    
    sensitivity + specificity - 1
  })
  
  if (all(is.na(youden_values))) {
    return(0.5)
  }
  
  thresholds[which.max(youden_values)]
}

evaluate_model <- function(y_true, prob, threshold = 0.5, model_name = "Model") {
  y_true <- as_prob_vector(y_true)
  prob <- as_prob_vector(prob)
  
  ok <- is.finite(y_true) & is.finite(prob)
  y_true <- y_true[ok]
  prob <- prob[ok]
  
  pred <- ifelse(prob >= threshold, 1, 0)
  
  TP <- sum(pred == 1 & y_true == 1)
  TN <- sum(pred == 0 & y_true == 0)
  FP <- sum(pred == 1 & y_true == 0)
  FN <- sum(pred == 0 & y_true == 1)
  
  accuracy <- (TP + TN) / length(y_true)
  sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  npv <- ifelse((TN + FN) == 0, NA, TN / (TN + FN))
  
  f1 <- ifelse(
    is.na(precision) | is.na(sensitivity) | (precision + sensitivity) == 0,
    NA,
    2 * precision * sensitivity / (precision + sensitivity)
  )
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  auc_value <- auc_manual(y_true, prob)
  brier <- mean((prob - y_true)^2)
  
  data.frame(
    Model = model_name,
    Threshold = threshold,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    NPV = npv,
    F1_score = f1,
    Balanced_accuracy = balanced_accuracy,
    AUC = auc_value,
    Brier_score = brier,
    stringsAsFactors = FALSE
  )
}

# ==========================================================
# 15. Feature importance
# ==========================================================

rf_importance <- data.frame(
  Feature = names(rf_model$variable.importance),
  RF_importance = as.numeric(rf_model$variable.importance),
  stringsAsFactors = FALSE
)

elastic_coef <- coef(elastic_model, s = "lambda.min")
elastic_importance <- data.frame(
  Feature = rownames(elastic_coef),
  Elastic_coef = as.numeric(elastic_coef),
  stringsAsFactors = FALSE
) %>%
  filter(Feature != "(Intercept)") %>%
  mutate(Elastic_abs_coef = abs(Elastic_coef))

feature_importance <- full_join(rf_importance, elastic_importance, by = "Feature") %>%
  mutate(
    RF_importance = ifelse(is.na(RF_importance), 0, RF_importance),
    Elastic_abs_coef = ifelse(is.na(Elastic_abs_coef), 0, Elastic_abs_coef),
    RF_scaled = scale01(RF_importance),
    Elastic_scaled = scale01(Elastic_abs_coef),
    Integrated_importance = rowMeans(cbind(RF_scaled, Elastic_scaled), na.rm = TRUE)
  ) %>%
  arrange(desc(Integrated_importance))

write.csv(
  feature_importance,
  file.path(OUTPUT_DIR, "FDRS_ML_feature_importance.csv"),
  row.names = FALSE
)

# ==========================================================
# 16. Apply model to target datasets
# ==========================================================

missing_files <- TARGET_FILES[!file.exists(TARGET_FILES)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "These target files were not found in working directory:",
      paste(missing_files, collapse = ", ")
    )
  )
}

target_data_list <- list()

for (i in seq_along(TARGET_FILES)) {
  target_data_list[[TARGET_NAMES[i]]] <- read_single_column_txt(TARGET_FILES[i])
}

target_feature_list <- list()

message("Extracting features for target datasets...")

for (nm in names(target_data_list)) {
  message("Target dataset: ", nm, " | n = ", length(target_data_list[[nm]]))
  
  target_feature_list[[nm]] <- extract_dataset_features(
    target_data_list[[nm]],
    dataset_name = nm,
    digits_keep = DIGITS_KEEP,
    single_positions = SINGLE_DIGIT_POSITIONS,
    primary_single_pos = PRIMARY_SINGLE_POS,
    joint_start_pos = JOINT_START_POS,
    joint_length = JOINT_LENGTH,
    B_progressive = B_PROGRESSIVE_ML,
    n_progressive_points = N_PROGRESSIVE_POINTS_ML,
    min_progressive_n = MIN_PROGRESSIVE_N_ML
  )
}

target_all_names <- names(X_raw)

target_feature_mat <- do.call(
  rbind,
  lapply(target_feature_list, function(v) {
    out <- rep(NA_real_, length(target_all_names))
    names(out) <- target_all_names
    common <- intersect(names(v), target_all_names)
    out[common] <- v[common]
    out
  })
)

target_feature_df_raw <- as.data.frame(target_feature_mat)
target_feature_df_raw[] <- lapply(target_feature_df_raw, as.numeric)

for (nm in names(target_feature_df_raw)) {
  target_feature_df_raw[[nm]][!is.finite(target_feature_df_raw[[nm]])] <- NA_real_
  target_feature_df_raw[[nm]][is.na(target_feature_df_raw[[nm]])] <- feature_medians[nm]
}

target_feature_df_raw <- target_feature_df_raw[, keep_cols, drop = FALSE]

target_feature_scaled <- as.data.frame(
  scale(target_feature_df_raw, center = feature_means[keep_cols], scale = feature_sds[keep_cols])
)

target_rf <- predict_rf(rf_model, target_feature_scaled)
target_elastic <- predict_elastic(elastic_model, target_feature_scaled)
target_svm <- predict_svm(svm_model, target_feature_scaled)
target_if <- predict_if_score(if_model, target_feature_scaled, reference_scores = if_train_raw)

target_supervised <- rowMeans(
  cbind(target_rf, target_elastic, target_svm),
  na.rm = TRUE
)

target_ensemble <- LAMBDA_SUPERVISED * target_supervised +
  (1 - LAMBDA_SUPERVISED) * target_if

target_risk_df <- data.frame(
  Dataset = rownames(target_feature_scaled),
  N = sapply(target_data_list, length),
  Random_Forest = target_rf,
  Elastic_net_Logistic = target_elastic,
  SVM_radial = target_svm,
  Isolation_Forest = target_if,
  Supervised_mean = target_supervised,
  Ensemble_risk = target_ensemble,
  Risk_grade = risk_grade_fun(target_ensemble),
  stringsAsFactors = FALSE
)

write.csv(
  target_risk_df,
  file.path(OUTPUT_DIR, "FDRS_target_dataset_risk_scores.csv"),
  row.names = FALSE
)

target_feature_export <- cbind(
  Dataset = rownames(target_feature_scaled),
  target_feature_df_raw
)

write.csv(
  target_feature_export,
  file.path(OUTPUT_DIR, "FDRS_target_dataset_features_raw.csv"),
  row.names = FALSE
)

# ==========================================================
# 17. Full statistics and detailed progressive report for targets
# ==========================================================

single_summary_list <- list()
single_freq_list <- list()
joint_summary_list <- list()
joint_freq_list <- list()
target_progressive_raw_list <- list()

for (nm in names(target_data_list)) {
  x <- target_data_list[[nm]]
  
  single_res <- calculate_single_digit_stats(
    x,
    dataset_name = nm,
    digit_pos = PRIMARY_SINGLE_POS,
    digits_keep = DIGITS_KEEP
  )
  
  joint_res <- calculate_joint_digit_stats(
    x,
    dataset_name = nm,
    start_pos = JOINT_START_POS,
    group_length = JOINT_LENGTH,
    digits_keep = DIGITS_KEEP
  )
  
  single_summary_list[[nm]] <- single_res$summary
  single_freq_list[[nm]] <- single_res$frequency
  
  joint_summary_list[[nm]] <- joint_res$summary
  joint_freq_list[[nm]] <- joint_res$frequency
  
  if (RUN_TARGET_PROGRESSIVE_REPORT) {
    message("Running detailed progressive subsampling for target: ", nm)
    
    target_progressive_raw_list[[nm]] <- progressive_subsampling(
      x,
      dataset_name = nm,
      start_pos = JOINT_START_POS,
      group_length = JOINT_LENGTH,
      digits_keep = DIGITS_KEEP,
      B = B_PROGRESSIVE_REPORT,
      min_n = MIN_PROGRESSIVE_N_REPORT,
      n_points = N_PROGRESSIVE_POINTS_REPORT
    )
  }
}

target_single_summary_df <- bind_rows(single_summary_list)
target_single_freq_df <- bind_rows(single_freq_list)
target_joint_summary_df <- bind_rows(joint_summary_list)
target_joint_freq_df <- bind_rows(joint_freq_list)

target_full_summary_df <- bind_rows(
  target_single_summary_df,
  target_joint_summary_df
)

write.csv(
  target_full_summary_df,
  file.path(OUTPUT_DIR, "FDRS_target_full_statistics_summary.csv"),
  row.names = FALSE
)

write.csv(
  target_single_freq_df,
  file.path(OUTPUT_DIR, "FDRS_target_single_digit_frequency.csv"),
  row.names = FALSE
)

write.csv(
  target_joint_freq_df,
  file.path(OUTPUT_DIR, "FDRS_target_joint_digit_frequency.csv"),
  row.names = FALSE
)

if (RUN_TARGET_PROGRESSIVE_REPORT) {
  target_progressive_raw_df <- bind_rows(target_progressive_raw_list)
  target_progressive_summary_df <- summarize_progressive(target_progressive_raw_df)
  
  target_progressive_trend_df <- target_progressive_summary_df %>%
    group_by(Dataset) %>%
    summarise(
      ChiSquare_slope = get_slope(SampleSize, ChiSquare_mean),
      Pvalue_slope = get_slope(SampleSize, Pvalue_median),
      minus_log10P_slope = get_slope(SampleSize, minus_log10P_median),
      CramersV_slope = get_slope(SampleSize, CramersV_mean),
      H_star_slope = get_slope(SampleSize, H_star_mean),
      KL_slope = get_slope(SampleSize, KL_mean),
      .groups = "drop"
    )
  
  write.csv(
    target_progressive_raw_df,
    file.path(OUTPUT_DIR, "FDRS_target_progressive_raw_results.csv"),
    row.names = FALSE
  )
  
  write.csv(
    target_progressive_summary_df,
    file.path(OUTPUT_DIR, "FDRS_target_progressive_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(
    target_progressive_trend_df,
    file.path(OUTPUT_DIR, "FDRS_target_progressive_trend_slopes.csv"),
    row.names = FALSE
  )
}

# ==========================================================
# 18. Visualizations
# ==========================================================

dataset_levels <- unique(c(training_df$DataType, target_risk_df$Dataset))

target_colors <- setNames(
  colorRampPalette(c(color_con, color_025, color_25))(length(unique(target_risk_df$Dataset))),
  unique(target_risk_df$Dataset)
)

if ("RawData" %in% names(target_colors)) target_colors["RawData"] <- color_con
if ("ErrData" %in% names(target_colors)) target_colors["ErrData"] <- color_25

# ==========================================================
# Rebuild model_performance before visualization
# ==========================================================

as_prob_vector <- function(x) {
  if (is.null(x)) return(numeric(0))
  
  if (is.data.frame(x)) x <- as.matrix(x)
  
  if (is.matrix(x)) {
    if (ncol(x) == 1) {
      x <- x[, 1]
    } else {
      x <- x[, ncol(x)]
    }
  }
  
  if (is.list(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  return(x)
}

auc_manual <- function(y_true, prob) {
  y_true <- as_prob_vector(y_true)
  prob <- as_prob_vector(prob)
  
  ok <- is.finite(y_true) & is.finite(prob)
  y_true <- y_true[ok]
  prob <- prob[ok]
  
  n_pos <- sum(y_true == 1)
  n_neg <- sum(y_true == 0)
  
  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }
  
  ranks <- rank(prob, ties.method = "average")
  auc <- (sum(ranks[y_true == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  return(as.numeric(auc))
}

get_best_threshold <- function(y_true, prob) {
  y_true <- as_prob_vector(y_true)
  prob <- as_prob_vector(prob)
  
  ok <- is.finite(y_true) & is.finite(prob)
  y_true <- y_true[ok]
  prob <- prob[ok]
  
  thresholds <- sort(unique(prob))
  
  if (length(thresholds) == 0) {
    return(0.5)
  }
  
  youden_values <- sapply(thresholds, function(thr) {
    pred <- ifelse(prob >= thr, 1, 0)
    
    TP <- sum(pred == 1 & y_true == 1)
    TN <- sum(pred == 0 & y_true == 0)
    FP <- sum(pred == 1 & y_true == 0)
    FN <- sum(pred == 0 & y_true == 1)
    
    sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
    specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
    
    sensitivity + specificity - 1
  })
  
  if (all(is.na(youden_values))) {
    return(0.5)
  }
  
  thresholds[which.max(youden_values)]
}

evaluate_model <- function(y_true, prob, threshold = 0.5, model_name = "Model") {
  y_true <- as_prob_vector(y_true)
  prob <- as_prob_vector(prob)
  
  ok <- is.finite(y_true) & is.finite(prob)
  y_true <- y_true[ok]
  prob <- prob[ok]
  
  pred <- ifelse(prob >= threshold, 1, 0)
  
  TP <- sum(pred == 1 & y_true == 1)
  TN <- sum(pred == 0 & y_true == 0)
  FP <- sum(pred == 1 & y_true == 0)
  FN <- sum(pred == 0 & y_true == 1)
  
  accuracy <- (TP + TN) / length(y_true)
  sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  npv <- ifelse((TN + FN) == 0, NA, TN / (TN + FN))
  
  f1 <- ifelse(
    is.na(precision) | is.na(sensitivity) | (precision + sensitivity) == 0,
    NA,
    2 * precision * sensitivity / (precision + sensitivity)
  )
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  auc_value <- auc_manual(y_true, prob)
  brier <- mean((prob - y_true)^2)
  
  data.frame(
    Model = model_name,
    Threshold = threshold,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    NPV = npv,
    F1_score = f1,
    Balanced_accuracy = balanced_accuracy,
    AUC = auc_value,
    Brier_score = brier,
    stringsAsFactors = FALSE
  )
}

# 重新确保所有预测概率都是 numeric vector
rf_train_prob <- as_prob_vector(rf_train_prob)
elastic_train_prob <- as_prob_vector(elastic_train_prob)
svm_train_prob <- as_prob_vector(svm_train_prob)
if_train_prob <- as_prob_vector(if_train_prob)
ensemble_train_prob <- as_prob_vector(ensemble_train_prob)

rf_test_prob <- as_prob_vector(rf_test_prob)
elastic_test_prob <- as_prob_vector(elastic_test_prob)
svm_test_prob <- as_prob_vector(svm_test_prob)
if_test_prob <- as_prob_vector(if_test_prob)
ensemble_test_prob <- as_prob_vector(ensemble_test_prob)

# 计算最佳阈值
thr_rf <- get_best_threshold(y_train, rf_train_prob)
thr_elastic <- get_best_threshold(y_train, elastic_train_prob)
thr_svm <- get_best_threshold(y_train, svm_train_prob)
thr_if <- get_best_threshold(y_train, if_train_prob)
thr_ensemble <- get_best_threshold(y_train, ensemble_train_prob)

# 生成模型评价表
model_performance <- dplyr::bind_rows(
  evaluate_model(y_test, rf_test_prob, thr_rf, "Random Forest"),
  evaluate_model(y_test, elastic_test_prob, thr_elastic, "Elastic-net Logistic"),
  evaluate_model(y_test, svm_test_prob, thr_svm, "SVM radial"),
  evaluate_model(y_test, if_test_prob, thr_if, "Isolation Forest"),
  evaluate_model(y_test, ensemble_test_prob, thr_ensemble, "Ensemble")
)

print(model_performance)

write.csv(
  model_performance,
  file.path(OUTPUT_DIR, "FDRS_ML_model_performance.csv"),
  row.names = FALSE
)
############################################################
# 18. Integrated visualization module
# Red-blue palette consistent with previous figures
############################################################

# ==========================================================
# 18.0 Packages
# ==========================================================

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggrepel)

# ==========================================================
# 18.1 Red-blue palette
# ==========================================================

color_con <- "#A9D1E8"   # Con-like light blue
color_025 <- "#F6C1C3"   # 0.25-like light pink
color_25  <- "#E68589"   # 25-like red

model_colors <- c(
  "Random Forest" = "#A9D1E8",
  "Elastic-net Logistic" = "#6FAFD3",
  "SVM radial" = "#F6C1C3",
  "Isolation Forest" = "#E68589",
  "Ensemble" = "#B94C5A"
)

model_linetypes <- c(
  "Random Forest" = "solid",
  "Elastic-net Logistic" = "longdash",
  "SVM radial" = "solid",
  "Isolation Forest" = "dotdash",
  "Ensemble" = "solid"
)

class_colors <- c(
  "Single digit" = "#A9D1E8",
  "Joint digit" = "#6FAFD3",
  "Progressive" = "#F6C1C3",
  "Information theory" = "#E68589",
  "Preference index" = "#B94C5A",
  "Other" = "#BDBDBD"
)

training_colors <- c(
  "Simulated normal" = "#A9D1E8",
  "Mixed-contamination" = "#F6C1C3",
  "Simulated digit-preference irregular" = "#E68589"
)

target_colors <- c(
  "RawData" = "#A9D1E8",
  "ErrData" = "#E68589"
)

# ==========================================================
# 18.2 Helper functions
# ==========================================================

as_prob_vector <- function(x) {
  if (is.null(x)) return(numeric(0))
  
  if (is.data.frame(x)) x <- as.matrix(x)
  
  if (is.matrix(x)) {
    if (ncol(x) == 1) {
      x <- x[, 1]
    } else {
      x <- x[, ncol(x)]
    }
  }
  
  if (is.list(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  return(x)
}

auc_manual <- function(y_true, prob) {
  y_true <- as_prob_vector(y_true)
  prob <- as_prob_vector(prob)
  
  ok <- is.finite(y_true) & is.finite(prob)
  y_true <- y_true[ok]
  prob <- prob[ok]
  
  n_pos <- sum(y_true == 1)
  n_neg <- sum(y_true == 0)
  
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  
  ranks <- rank(prob, ties.method = "average")
  auc <- (sum(ranks[y_true == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  
  return(as.numeric(auc))
}

brier_manual <- function(y_true, prob) {
  y_true <- as_prob_vector(y_true)
  prob <- as_prob_vector(prob)
  
  ok <- is.finite(y_true) & is.finite(prob)
  
  mean((prob[ok] - y_true[ok])^2)
}

get_best_threshold <- function(y_true, prob) {
  y_true <- as_prob_vector(y_true)
  prob <- as_prob_vector(prob)
  
  ok <- is.finite(y_true) & is.finite(prob)
  y_true <- y_true[ok]
  prob <- prob[ok]
  
  thresholds <- sort(unique(prob))
  
  if (length(thresholds) == 0) return(0.5)
  
  youden_values <- sapply(thresholds, function(thr) {
    pred <- ifelse(prob >= thr, 1, 0)
    
    TP <- sum(pred == 1 & y_true == 1)
    TN <- sum(pred == 0 & y_true == 0)
    FP <- sum(pred == 1 & y_true == 0)
    FN <- sum(pred == 0 & y_true == 1)
    
    sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
    specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
    
    sensitivity + specificity - 1
  })
  
  if (all(is.na(youden_values))) return(0.5)
  
  thresholds[which.max(youden_values)]
}

evaluate_model <- function(y_true, prob, threshold = 0.5, model_name = "Model") {
  y_true <- as_prob_vector(y_true)
  prob <- as_prob_vector(prob)
  
  ok <- is.finite(y_true) & is.finite(prob)
  y_true <- y_true[ok]
  prob <- prob[ok]
  
  pred <- ifelse(prob >= threshold, 1, 0)
  
  TP <- sum(pred == 1 & y_true == 1)
  TN <- sum(pred == 0 & y_true == 0)
  FP <- sum(pred == 1 & y_true == 0)
  FN <- sum(pred == 0 & y_true == 1)
  
  accuracy <- (TP + TN) / length(y_true)
  sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  npv <- ifelse((TN + FN) == 0, NA, TN / (TN + FN))
  
  f1 <- ifelse(
    is.na(precision) | is.na(sensitivity) | (precision + sensitivity) == 0,
    NA,
    2 * precision * sensitivity / (precision + sensitivity)
  )
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  auc_value <- auc_manual(y_true, prob)
  brier <- brier_manual(y_true, prob)
  
  data.frame(
    Model = model_name,
    Threshold = threshold,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    NPV = npv,
    F1_score = f1,
    Balanced_accuracy = balanced_accuracy,
    AUC = auc_value,
    Brier_score = brier,
    stringsAsFactors = FALSE
  )
}

make_roc_df <- function(y_true, prob, model_name) {
  y_true <- as_prob_vector(y_true)
  prob <- as_prob_vector(prob)
  
  ok <- is.finite(y_true) & is.finite(prob)
  y_true <- y_true[ok]
  prob <- prob[ok]
  
  thresholds <- sort(unique(prob), decreasing = TRUE)
  thresholds <- c(Inf, thresholds, -Inf)
  
  roc_list <- lapply(thresholds, function(thr) {
    pred <- ifelse(prob >= thr, 1, 0)
    
    TP <- sum(pred == 1 & y_true == 1)
    TN <- sum(pred == 0 & y_true == 0)
    FP <- sum(pred == 1 & y_true == 0)
    FN <- sum(pred == 0 & y_true == 1)
    
    sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
    specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
    
    data.frame(
      Model = model_name,
      Threshold = thr,
      FPR = 1 - specificity,
      Sensitivity = sensitivity,
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(roc_list)
}

make_calibration_df <- function(y_true, prob, model_name, n_bins = 10) {
  df <- data.frame(
    y_true = as_prob_vector(y_true),
    Probability = as_prob_vector(prob)
  )
  
  df <- df[is.finite(df$y_true) & is.finite(df$Probability), ]
  
  df <- df %>%
    dplyr::mutate(
      Bin = dplyr::ntile(Probability, n_bins)
    ) %>%
    dplyr::group_by(Bin) %>%
    dplyr::summarise(
      Mean_predicted = mean(Probability, na.rm = TRUE),
      Observed_rate = mean(y_true, na.rm = TRUE),
      N = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Model = model_name
    )
  
  return(df)
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

# ==========================================================
# 18.3 Rebuild model performance if needed
# ==========================================================

rf_train_prob <- as_prob_vector(rf_train_prob)
elastic_train_prob <- as_prob_vector(elastic_train_prob)
svm_train_prob <- as_prob_vector(svm_train_prob)
if_train_prob <- as_prob_vector(if_train_prob)
ensemble_train_prob <- as_prob_vector(ensemble_train_prob)

rf_test_prob <- as_prob_vector(rf_test_prob)
elastic_test_prob <- as_prob_vector(elastic_test_prob)
svm_test_prob <- as_prob_vector(svm_test_prob)
if_test_prob <- as_prob_vector(if_test_prob)
ensemble_test_prob <- as_prob_vector(ensemble_test_prob)

if (!exists("thr_rf")) thr_rf <- get_best_threshold(y_train, rf_train_prob)
if (!exists("thr_elastic")) thr_elastic <- get_best_threshold(y_train, elastic_train_prob)
if (!exists("thr_svm")) thr_svm <- get_best_threshold(y_train, svm_train_prob)
if (!exists("thr_if")) thr_if <- get_best_threshold(y_train, if_train_prob)
if (!exists("thr_ensemble")) thr_ensemble <- get_best_threshold(y_train, ensemble_train_prob)

model_performance <- dplyr::bind_rows(
  evaluate_model(y_test, rf_test_prob, thr_rf, "Random Forest"),
  evaluate_model(y_test, elastic_test_prob, thr_elastic, "Elastic-net Logistic"),
  evaluate_model(y_test, svm_test_prob, thr_svm, "SVM radial"),
  evaluate_model(y_test, if_test_prob, thr_if, "Isolation Forest"),
  evaluate_model(y_test, ensemble_test_prob, thr_ensemble, "Ensemble")
)

write.csv(
  model_performance,
  file.path(OUTPUT_DIR, "FDRS_ML_model_performance.csv"),
  row.names = FALSE
)

# ==========================================================
# 18.4 Model performance bubble plot
# ==========================================================

performance_long <- model_performance %>%
  tidyr::pivot_longer(
    cols = c(
      Accuracy,
      Sensitivity,
      Specificity,
      Precision,
      F1_score,
      Balanced_accuracy,
      AUC,
      Brier_score
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::group_by(Metric) %>%
  dplyr::mutate(Value_scaled = scale01(Value)) %>%
  dplyr::ungroup()

p_perf <- ggplot(performance_long, aes(x = Metric, y = Model)) +
  geom_point(aes(size = Value_scaled, color = Value_scaled), alpha = 0.9) +
  geom_text(aes(label = sprintf("%.3f", Value)), size = 3, vjust = -1.1) +
  scale_size_continuous(range = c(3, 11), limits = c(0, 1), name = "Relative value") +
  scale_color_gradientn(
    colors = c(color_con, color_025, color_25),
    limits = c(0, 1),
    name = "Relative value"
  ) +
  labs(
    title = "Model performance summary",
    x = "Metric",
    y = "Model"
  ) +
  theme_fdrs_style() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

print(p_perf)

save_plot(
  p_perf,
  "FDRS_ML_figure_01_model_performance",
  width = 10,
  height = 6
)

# ==========================================================
# 18.5 Target risk score plot
# ==========================================================

risk_long <- target_risk_df %>%
  dplyr::select(
    Dataset,
    Random_Forest,
    Elastic_net_Logistic,
    SVM_radial,
    Isolation_Forest,
    Ensemble_risk
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      Random_Forest,
      Elastic_net_Logistic,
      SVM_radial,
      Isolation_Forest,
      Ensemble_risk
    ),
    names_to = "Model",
    values_to = "Risk"
  )

p_risk <- ggplot(risk_long, aes(x = Model, y = Risk, fill = Dataset)) +
  geom_col(
    position = position_dodge(width = 0.75),
    color = "black",
    linewidth = 0.25,
    width = 0.68
  ) +
  geom_hline(
    yintercept = c(0.2, 0.4, 0.6, 0.8),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = target_colors) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Risk scores of target datasets",
    x = "Model",
    y = "Risk score",
    fill = "Dataset"
  ) +
  theme_fdrs_style() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p_risk)

save_plot(
  p_risk,
  "FDRS_ML_figure_02_target_risk_scores",
  width = 9,
  height = 5.8
)

# ==========================================================
# 18.6 ROC curves with AUC labels
# ==========================================================

roc_df <- dplyr::bind_rows(
  make_roc_df(y_test, rf_test_prob, "Random Forest"),
  make_roc_df(y_test, elastic_test_prob, "Elastic-net Logistic"),
  make_roc_df(y_test, svm_test_prob, "SVM radial"),
  make_roc_df(y_test, if_test_prob, "Isolation Forest"),
  make_roc_df(y_test, ensemble_test_prob, "Ensemble")
)

auc_df <- data.frame(
  Model = c(
    "Random Forest",
    "Elastic-net Logistic",
    "SVM radial",
    "Isolation Forest",
    "Ensemble"
  ),
  AUC = c(
    auc_manual(y_test, rf_test_prob),
    auc_manual(y_test, elastic_test_prob),
    auc_manual(y_test, svm_test_prob),
    auc_manual(y_test, if_test_prob),
    auc_manual(y_test, ensemble_test_prob)
  ),
  stringsAsFactors = FALSE
)

auc_df <- auc_df %>%
  dplyr::mutate(
    AUC_label = paste0(Model, ": AUC = ", sprintf("%.3f", AUC)),
    x = 0.58,
    y = seq(0.30, 0.30 - 0.075 * (dplyr::n() - 1), by = -0.075)
  )

legend_labels <- setNames(
  paste0(auc_df$Model, " (AUC=", sprintf("%.3f", auc_df$AUC), ")"),
  auc_df$Model
)

p_roc <- ggplot(
  roc_df,
  aes(x = FPR, y = Sensitivity, color = Model, linetype = Model)
) +
  geom_line(linewidth = 1.15) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey55"
  ) +
  geom_label(
    data = auc_df,
    aes(x = x, y = y, label = AUC_label, color = Model),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.7,
    fill = "white",
    label.size = 0.25,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = model_colors,
    labels = legend_labels
  ) +
  scale_linetype_manual(
    values = model_linetypes,
    labels = legend_labels
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "ROC curves",
    subtitle = "AUC values are shown on the plot",
    x = "False positive rate",
    y = "Sensitivity",
    color = "Model",
    linetype = "Model"
  ) +
  theme_fdrs_style() +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

print(p_roc)

save_plot(
  p_roc,
  "FDRS_ML_figure_03_ROC_curves_with_AUC",
  width = 8.5,
  height = 6.5
)

# ==========================================================
# 18.7 PCA feature space
# ==========================================================

pca_train <- prcomp(X_scaled, center = FALSE, scale. = FALSE)

pca_train_df <- data.frame(
  PC1 = pca_train$x[, 1],
  PC2 = pca_train$x[, 2],
  DataType = training_df$DataType,
  Label = factor(training_df$Label),
  stringsAsFactors = FALSE
)

target_pca <- predict(pca_train, newdata = target_feature_scaled)

target_pca_df <- data.frame(
  PC1 = target_pca[, 1],
  PC2 = target_pca[, 2],
  Dataset = rownames(target_feature_scaled),
  stringsAsFactors = FALSE
)

pc_var <- summary(pca_train)$importance[2, ]
pc1_var <- round(pc_var[1] * 100, 1)
pc2_var <- round(pc_var[2] * 100, 1)

x_lim <- quantile(pca_train_df$PC1, probs = c(0.01, 0.99), na.rm = TRUE)
y_lim <- quantile(pca_train_df$PC2, probs = c(0.01, 0.99), na.rm = TRUE)

x_pad <- 0.08 * diff(range(c(x_lim, target_pca_df$PC1), na.rm = TRUE))
y_pad <- 0.08 * diff(range(c(y_lim, target_pca_df$PC2), na.rm = TRUE))

x_min <- min(c(x_lim[1], target_pca_df$PC1), na.rm = TRUE) - x_pad
x_max <- max(c(x_lim[2], target_pca_df$PC1), na.rm = TRUE) + x_pad
y_min <- min(c(y_lim[1], target_pca_df$PC2), na.rm = TRUE) - y_pad
y_max <- max(c(y_lim[2], target_pca_df$PC2), na.rm = TRUE) + y_pad

target_shapes <- c(
  "RawData" = 24,
  "ErrData" = 21
)

missing_shapes <- setdiff(target_pca_df$Dataset, names(target_shapes))

if (length(missing_shapes) > 0) {
  extra_shapes <- c(22, 23, 25, 8, 4, 15)
  target_shapes[missing_shapes] <- extra_shapes[seq_along(missing_shapes)]
}

p_pca <- ggplot() +
  geom_point(
    data = pca_train_df,
    aes(x = PC1, y = PC2, color = DataType),
    alpha = 0.30,
    size = 1.6
  ) +
  stat_ellipse(
    data = pca_train_df,
    aes(x = PC1, y = PC2, color = DataType),
    type = "norm",
    level = 0.95,
    linewidth = 0.9,
    alpha = 0.95
  ) +
  geom_point(
    data = target_pca_df,
    aes(x = PC1, y = PC2, shape = Dataset),
    size = 4.8,
    stroke = 1.2,
    fill = "black",
    color = "black"
  ) +
  ggrepel::geom_label_repel(
    data = target_pca_df,
    aes(x = PC1, y = PC2, label = Dataset),
    size = 4,
    fontface = "bold",
    box.padding = 0.45,
    point.padding = 0.35,
    label.size = 0.3,
    fill = "white",
    color = "black",
    min.segment.length = 0
  ) +
  scale_color_manual(values = training_colors) +
  scale_shape_manual(values = target_shapes) +
  coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
  labs(
    title = "PCA of digit-pattern feature space",
    subtitle = "Training clusters and target datasets",
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    color = "Training data",
    shape = "Target dataset"
  ) +
  theme_fdrs_style()

print(p_pca)

save_plot(
  p_pca,
  "FDRS_ML_figure_04_PCA_feature_space",
  width = 8.5,
  height = 6.8
)

# ==========================================================
# 18.8 Target full statistics bubble plot
# ==========================================================

target_full_long <- target_full_summary_df %>%
  dplyr::select(
    Dataset,
    Test,
    ChiSquare,
    Pvalue,
    CramersV,
    H_star,
    KL,
    Pref_same_pair,
    Pref_neat_combo,
    Pref_sequential_combo
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      ChiSquare,
      Pvalue,
      CramersV,
      H_star,
      KL,
      Pref_same_pair,
      Pref_neat_combo,
      Pref_sequential_combo
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::group_by(Metric) %>%
  dplyr::mutate(Value_scaled = scale01(Value)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Label = format_metric(Value))

p_full <- ggplot(target_full_long, aes(x = Metric, y = Dataset)) +
  geom_point(aes(size = Value_scaled, color = Value_scaled), alpha = 0.9) +
  geom_text(aes(label = Label), size = 3, vjust = -1.15) +
  facet_wrap(~ Test, ncol = 1) +
  scale_size_continuous(
    range = c(3, 10),
    limits = c(0, 1),
    name = "Relative value"
  ) +
  scale_color_gradientn(
    colors = c(color_con, color_025, color_25),
    limits = c(0, 1),
    name = "Relative value"
  ) +
  labs(
    title = "Full-sample statistical summary of target datasets",
    x = "Metric",
    y = "Dataset"
  ) +
  theme_fdrs_style() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

print(p_full)

save_plot(
  p_full,
  "FDRS_ML_figure_05_target_full_statistics",
  width = 10,
  height = 7
)

# ==========================================================
# 18.9 Progressive subsampling curves
# ==========================================================

if (exists("target_progressive_summary_df")) {
  
  progressive_long <- target_progressive_summary_df %>%
    dplyr::transmute(
      Dataset = Dataset,
      SampleSize = SampleSize,
      `Chi-square` = ChiSquare_mean,
      `Median P value` = Pvalue_median,
      `Cramér's V` = CramersV_mean,
      `Normalized entropy H*` = H_star_mean,
      `KL divergence` = KL_mean,
      `-log10(P)` = minus_log10P_median
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        `Chi-square`,
        `Median P value`,
        `Cramér's V`,
        `Normalized entropy H*`,
        `KL divergence`,
        `-log10(P)`
      ),
      names_to = "Metric",
      values_to = "Value"
    )
  
  p_prog <- ggplot(progressive_long, aes(x = SampleSize, y = Value, color = Dataset)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.4) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
    scale_color_manual(values = target_colors) +
    labs(
      title = "Progressive subsampling metrics of target datasets",
      x = "Sample size (n)",
      y = "Metric value",
      color = "Dataset"
    ) +
    theme_fdrs_style()
  
  print(p_prog)
  
  save_plot(
    p_prog,
    "FDRS_ML_figure_06_target_progressive_metrics",
    width = 10,
    height = 8
  )
}

# ==========================================================
# 18.10 Confusion matrices
# ==========================================================

model_prob_df <- dplyr::bind_rows(
  data.frame(
    Model = "Random Forest",
    y_true = as_prob_vector(y_test),
    Probability = as_prob_vector(rf_test_prob),
    Threshold = thr_rf
  ),
  data.frame(
    Model = "Elastic-net Logistic",
    y_true = as_prob_vector(y_test),
    Probability = as_prob_vector(elastic_test_prob),
    Threshold = thr_elastic
  ),
  data.frame(
    Model = "SVM radial",
    y_true = as_prob_vector(y_test),
    Probability = as_prob_vector(svm_test_prob),
    Threshold = thr_svm
  ),
  data.frame(
    Model = "Isolation Forest",
    y_true = as_prob_vector(y_test),
    Probability = as_prob_vector(if_test_prob),
    Threshold = thr_if
  ),
  data.frame(
    Model = "Ensemble",
    y_true = as_prob_vector(y_test),
    Probability = as_prob_vector(ensemble_test_prob),
    Threshold = thr_ensemble
  )
)

model_prob_df <- model_prob_df %>%
  dplyr::mutate(
    Predicted = ifelse(Probability >= Threshold, 1, 0),
    Actual_class = ifelse(y_true == 1, "Irregular", "Normal"),
    Predicted_class = ifelse(Predicted == 1, "Irregular", "Normal"),
    Actual_class = factor(Actual_class, levels = c("Irregular", "Normal")),
    Predicted_class = factor(Predicted_class, levels = c("Normal", "Irregular"))
  )

confusion_df <- model_prob_df %>%
  dplyr::group_by(Model, Actual_class, Predicted_class) %>%
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
  tidyr::complete(
    Model,
    Actual_class,
    Predicted_class,
    fill = list(Count = 0)
  ) %>%
  dplyr::group_by(Model) %>%
  dplyr::mutate(
    Percent = Count / sum(Count) * 100,
    Label = paste0(Count, "\n", sprintf("%.1f", Percent), "%")
  ) %>%
  dplyr::ungroup()

p_confusion <- ggplot(
  confusion_df,
  aes(x = Predicted_class, y = Actual_class, fill = Count)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = Label),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(~ Model, ncol = 3) +
  scale_fill_gradientn(
    colors = c("#F7FBFF", color_con, color_025, color_25),
    name = "Count"
  ) +
  labs(
    title = "Confusion matrices of machine-learning models",
    x = "Predicted class",
    y = "Actual class"
  ) +
  theme_fdrs_style() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank()
  )

print(p_confusion)

save_plot(
  p_confusion,
  "FDRS_ML_figure_07_confusion_matrices",
  width = 10,
  height = 6.8
)

write.csv(
  confusion_df,
  file.path(OUTPUT_DIR, "FDRS_ML_confusion_matrix_summary.csv"),
  row.names = FALSE
)

# ==========================================================
# 18.11 Calibration curves
# ==========================================================

calibration_df <- dplyr::bind_rows(
  make_calibration_df(y_test, rf_test_prob, "Random Forest"),
  make_calibration_df(y_test, elastic_test_prob, "Elastic-net Logistic"),
  make_calibration_df(y_test, svm_test_prob, "SVM radial"),
  make_calibration_df(y_test, if_test_prob, "Isolation Forest"),
  make_calibration_df(y_test, ensemble_test_prob, "Ensemble")
)

calibration_label_df <- data.frame(
  Model = c(
    "Random Forest",
    "Elastic-net Logistic",
    "SVM radial",
    "Isolation Forest",
    "Ensemble"
  ),
  Brier = c(
    brier_manual(y_test, rf_test_prob),
    brier_manual(y_test, elastic_test_prob),
    brier_manual(y_test, svm_test_prob),
    brier_manual(y_test, if_test_prob),
    brier_manual(y_test, ensemble_test_prob)
  ),
  stringsAsFactors = FALSE
)

calibration_legend_labels <- setNames(
  paste0(
    calibration_label_df$Model,
    " (Brier=",
    sprintf("%.3f", calibration_label_df$Brier),
    ")"
  ),
  calibration_label_df$Model
)

p_calibration <- ggplot(
  calibration_df,
  aes(
    x = Mean_predicted,
    y = Observed_rate,
    color = Model,
    linetype = Model
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.8
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(aes(size = N), alpha = 0.9) +
  scale_color_manual(
    values = model_colors,
    labels = calibration_legend_labels
  ) +
  scale_linetype_manual(
    values = model_linetypes,
    labels = calibration_legend_labels
  ) +
  scale_size_continuous(range = c(2.5, 6), name = "Bin size") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Calibration curves",
    subtitle = "Observed irregularity rate versus predicted risk",
    x = "Mean predicted risk",
    y = "Observed irregularity rate",
    color = "Model",
    linetype = "Model"
  ) +
  theme_fdrs_style()

print(p_calibration)

save_plot(
  p_calibration,
  "FDRS_ML_figure_08_calibration_curves",
  width = 8.5,
  height = 6.5
)

write.csv(
  calibration_df,
  file.path(OUTPUT_DIR, "FDRS_ML_calibration_curve_summary.csv"),
  row.names = FALSE
)

# ==========================================================
# 18.12 Feature importance plot
# ==========================================================

if (!"Integrated_importance" %in% colnames(feature_importance)) {
  
  if (!"RF_scaled" %in% colnames(feature_importance) &&
      "RF_importance" %in% colnames(feature_importance)) {
    feature_importance$RF_scaled <- scale01(feature_importance$RF_importance)
  }
  
  if (!"Elastic_scaled" %in% colnames(feature_importance) &&
      "Elastic_abs_coef" %in% colnames(feature_importance)) {
    feature_importance$Elastic_scaled <- scale01(feature_importance$Elastic_abs_coef)
  }
  
  feature_importance$Integrated_importance <- rowMeans(
    cbind(
      feature_importance$RF_scaled,
      feature_importance$Elastic_scaled
    ),
    na.rm = TRUE
  )
}

feature_importance_plot_df <- feature_importance %>%
  dplyr::mutate(
    Feature_class = dplyr::case_when(
      grepl("^single_pos", Feature) ~ "Single digit",
      grepl("^joint_34", Feature) ~ "Joint digit",
      grepl("^prog_34", Feature) ~ "Progressive",
      grepl("H_star|KL|Entropy", Feature) ~ "Information theory",
      grepl("Pref|pair|combo|odd|middle|extreme", Feature) ~ "Preference index",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::arrange(desc(Integrated_importance)) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::mutate(
    Feature = factor(Feature, levels = rev(Feature))
  )

p_importance <- ggplot(
  feature_importance_plot_df,
  aes(x = Integrated_importance, y = Feature)
) +
  geom_segment(
    aes(
      x = 0,
      xend = Integrated_importance,
      y = Feature,
      yend = Feature,
      color = Feature_class
    ),
    linewidth = 0.9,
    alpha = 0.85
  ) +
  geom_point(
    aes(
      size = Integrated_importance,
      color = Feature_class
    ),
    alpha = 0.95
  ) +
  scale_color_manual(values = class_colors, name = "Feature class") +
  scale_size_continuous(range = c(3, 8), name = "Importance") +
  labs(
    title = "Top feature importance",
    subtitle = "Integrated importance from machine-learning feature ranking",
    x = "Integrated importance",
    y = "Feature"
  ) +
  theme_fdrs_style() +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

print(p_importance)

save_plot(
  p_importance,
  "FDRS_ML_figure_09_feature_importance_top30",
  width = 10,
  height = 8
)

write.csv(
  feature_importance_plot_df,
  file.path(OUTPUT_DIR, "FDRS_ML_top30_feature_importance.csv"),
  row.names = FALSE
)

# ==========================================================
# 18.13 Isolation Forest anomaly score
# ==========================================================

if_all_raw <- as_prob_vector(
  predict(
    if_model,
    newdata = X_scaled,
    type = "score"
  )
)

if_reference_raw <- as_prob_vector(
  predict(
    if_model,
    newdata = X_scaled[training_df$Label == 0, , drop = FALSE],
    type = "score"
  )
)

q01 <- quantile(if_reference_raw, 0.01, na.rm = TRUE)
q99 <- quantile(if_reference_raw, 0.99, na.rm = TRUE)

if_all_clip <- pmin(pmax(if_all_raw, q01), q99)
if_all_risk <- (if_all_clip - q01) / (q99 - q01)
if_all_risk[!is.finite(if_all_risk)] <- 0.5

if_training_df <- data.frame(
  DatasetID = training_df$DatasetID,
  DataType = training_df$DataType,
  Label = training_df$Label,
  AnomalyRisk = if_all_risk,
  stringsAsFactors = FALSE
)

if_target_raw <- as_prob_vector(
  predict(
    if_model,
    newdata = target_feature_scaled,
    type = "score"
  )
)

if_target_clip <- pmin(pmax(if_target_raw, q01), q99)
if_target_risk <- (if_target_clip - q01) / (q99 - q01)
if_target_risk[!is.finite(if_target_risk)] <- 0.5

if_target_df <- data.frame(
  Dataset = rownames(target_feature_scaled),
  AnomalyRisk = if_target_risk,
  stringsAsFactors = FALSE
)

if ("Isolation_Forest" %in% colnames(target_risk_df)) {
  if_target_df <- target_risk_df %>%
    dplyr::select(Dataset, Isolation_Forest) %>%
    dplyr::rename(AnomalyRisk = Isolation_Forest)
}

p_if_density <- ggplot() +
  geom_density(
    data = if_training_df,
    aes(x = AnomalyRisk, fill = DataType),
    alpha = 0.35,
    linewidth = 0.6,
    color = "grey30"
  ) +
  geom_vline(
    data = if_target_df,
    aes(xintercept = AnomalyRisk, color = Dataset),
    linewidth = 1.3,
    linetype = "solid"
  ) +
  geom_label(
    data = if_target_df,
    aes(
      x = AnomalyRisk,
      y = Inf,
      label = paste0(Dataset, "\n", sprintf("%.3f", AnomalyRisk)),
      color = Dataset
    ),
    vjust = 1.15,
    hjust = 0.5,
    size = 3.8,
    fill = "white",
    label.size = 0.25,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = training_colors, name = "Training data") +
  scale_color_manual(values = target_colors, name = "Target dataset") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "Isolation Forest anomaly score distribution",
    subtitle = "Target datasets are marked by vertical lines",
    x = "Isolation Forest anomaly risk",
    y = "Density"
  ) +
  theme_fdrs_style()

print(p_if_density)

save_plot(
  p_if_density,
  "FDRS_ML_figure_10_isolation_forest_anomaly_density",
  width = 8.5,
  height = 6.5
)

p_if_violin <- ggplot(
  if_training_df,
  aes(x = DataType, y = AnomalyRisk, fill = DataType)
) +
  geom_violin(
    alpha = 0.45,
    color = "black",
    linewidth = 0.3,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    alpha = 0.75,
    color = "black",
    linewidth = 0.35
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.12,
    size = 0.8,
    color = "black"
  ) +
  geom_point(
    data = if_target_df,
    aes(x = Dataset, y = AnomalyRisk, color = Dataset),
    inherit.aes = FALSE,
    size = 4.5,
    stroke = 1.2
  ) +
  geom_label(
    data = if_target_df,
    aes(
      x = Dataset,
      y = AnomalyRisk,
      label = paste0(Dataset, "\n", sprintf("%.3f", AnomalyRisk)),
      color = Dataset
    ),
    inherit.aes = FALSE,
    size = 3.6,
    fill = "white",
    label.size = 0.25,
    vjust = -0.8,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = training_colors, name = "Training data") +
  scale_color_manual(values = target_colors, name = "Target dataset") +
  scale_y_continuous(limits = c(0, 1.08), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Isolation Forest anomaly risk",
    subtitle = "Training distributions and target datasets",
    x = "",
    y = "Anomaly risk"
  ) +
  theme_fdrs_style() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

print(p_if_violin)

save_plot(
  p_if_violin,
  "FDRS_ML_figure_11_isolation_forest_anomaly_violin",
  width = 9,
  height = 6.5
)

write.csv(
  if_training_df,
  file.path(OUTPUT_DIR, "FDRS_ML_isolation_forest_training_scores.csv"),
  row.names = FALSE
)

write.csv(
  if_target_df,
  file.path(OUTPUT_DIR, "FDRS_ML_isolation_forest_target_scores.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("Integrated visualization module finished.\n")
cat("Figures saved in folder: ", OUTPUT_DIR, "\n", sep = "")
cat("Generated figures:\n")
cat("  - FDRS_ML_figure_01_model_performance\n")
cat("  - FDRS_ML_figure_02_target_risk_scores\n")
cat("  - FDRS_ML_figure_03_ROC_curves_with_AUC\n")
cat("  - FDRS_ML_figure_04_PCA_feature_space\n")
cat("  - FDRS_ML_figure_05_target_full_statistics\n")
cat("  - FDRS_ML_figure_06_target_progressive_metrics\n")
cat("  - FDRS_ML_figure_07_confusion_matrices\n")
cat("  - FDRS_ML_figure_08_calibration_curves\n")
cat("  - FDRS_ML_figure_09_feature_importance_top30\n")
cat("  - FDRS_ML_figure_10_isolation_forest_anomaly_density\n")
cat("  - FDRS_ML_figure_11_isolation_forest_anomaly_violin\n")
cat("============================================================\n")