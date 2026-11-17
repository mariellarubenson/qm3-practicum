### Quant 3 Practicum ###

library(ggplot2)
library(estimatr)
library(tidyverse)

set.seed(123)

## Randomization

fat_levels  <- c("3.25", "2", "1", "0")
temp_levels <- c("fridge", "warm")
temp_target <- c(fridge = 4, warm = 40)
blocks <- 3

conditions <- expand.grid(fat = fat_levels, temp = temp_levels,
                          stringsAsFactors = FALSE)

design <- do.call(rbind, lapply(seq_len(blocks), function(b) {
  d <- conditions[sample(nrow(conditions)), ]
  d$block <- b
  d$order_in_block <- seq_len(nrow(d))
  d
}))
design$run_id <- seq_len(nrow(design))
rownames(design) <- NULL
design <- design[, c("run_id", "block", "order_in_block", "fat", "temp")]

# check: each condition appears exactly once per block
stopifnot(all(table(design$fat, design$temp, design$block) == 1))
design

## Data collection sheet

sheet <- design
sheet$start_temp_c  <- NA
sheet$total_ml_t0   <- NA
sheet$liquid_ml_t0  <- NA
sheet$total_ml_t10  <- NA
sheet$liquid_ml_t10 <- NA
sheet$rater1 <- NA
sheet$rater2 <- NA
sheet$rater3 <- NA
sheet$rerun  <- 0

write.csv(sheet, "/Users/mariellarubenson/Library/CloudStorage/Dropbox/Quant 3/practicum/randomization_sheet.csv", row.names = FALSE)

## Cleaning and derived variables

dat <- read.csv("/Users/mariellarubenson/Library/CloudStorage/Dropbox/Quant 3/practicum/practicum_data.csv",
                colClasses = c(fat = "character", temp = "character"))
dat$fat <- factor(dat$fat, levels = c("0", "1", "2", "3.25"))
dat$temp <- factor(dat$temp, levels = temp_levels)

dat$foam_t0 <- dat$total_ml_t0  - dat$liquid_ml_t0
dat$foam_t10 <- dat$total_ml_t10 - dat$liquid_ml_t10
dat$stability <- dat$foam_t10 / dat$foam_t0
dat$foam_share <- dat$foam_t0 / dat$total_ml_t0
dat$texture   <- rowMeans(dat[, c("rater1", "rater2", "rater3")], na.rm = TRUE)

# manipulation check: starting temperature by condition
dat$temp_compliant <- abs(dat$start_temp_c - temp_target[as.character(dat$temp)]) <= 2
aggregate(start_temp_c ~ temp, data = dat, FUN = mean)
table(dat$temp, dat$temp_compliant)

## Distribution of outcomes by condition

aggregate(foam_t0 ~ fat + temp, data = dat,
          FUN = function(x) round(c(mean = mean(x), sd = sd(x)), 1))

dist <- dat %>%
  group_by(fat, temp) %>%
  summarise(n = n(),
            across(c(foam_t0, stability, texture),
                   list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))),
            .groups = "drop")
knitr::kable(dist, digits = 2, format = "latex", booktabs = TRUE)

## Analysis

ate <- function(y, var, treat, control, sims = 10000) {
  sub   <- dat[dat[[var]] %in% c(treat, control), ]
  sub$Z <- as.numeric(sub[[var]] == treat)
  fit   <- lm_robust(as.formula(paste(y, "~ Z")), data = sub)
  
  other  <- if (var == "fat") "temp" else "fat"
  strata <- interaction(dat$block, dat[[other]])
  lab    <- as.character(dat[[var]])
  dim    <- function(l) mean(dat[[y]][l == treat], na.rm = TRUE) -
    mean(dat[[y]][l == control], na.rm = TRUE)
  obs  <- dim(lab)
  null <- replicate(sims, dim(ave(lab, strata, FUN = sample)))
  
  data.frame(outcome  = y,
             contrast = paste(treat, "vs", control),
             ate      = unname(coef(fit)["Z"]),
             se       = unname(fit$std.error["Z"]),
             ci_low   = unname(fit$conf.low["Z"]),
             ci_high  = unname(fit$conf.high["Z"]),
             p_ri     = mean(abs(null) >= abs(obs) - 1e-9))
}

set.seed(123)
results <- rbind(
  # H1: foam volume by fat
  ate("foam_t0", "fat", "1", "0"),
  ate("foam_t0", "fat", "2", "0"),
  ate("foam_t0", "fat", "3.25", "0"),
  # H2: stability by fat
  ate("stability", "fat", "1", "0"),
  ate("stability", "fat", "2", "0"),
  ate("stability", "fat", "3.25", "0"),
  # H3: foam volume by starting temperature
  ate("foam_t0", "temp", "warm", "fridge"),
  # H4: texture by fat
  ate("texture", "fat", "1", "0"),
  ate("texture", "fat", "2", "0"),
  ate("texture", "fat", "3.25", "0")
)
results

knitr::kable(results, digits = 2)

rbind(ate("foam_share", "fat", "3.25", "0"),
      ate("foam_share", "temp", "warm", "fridge"))

# main effects in a regression
summary(lm_robust(foam_t0 ~ fat + temp, data = dat))

## Figures

temp_labs <- c(fridge = "Fridge", warm = "Pre-warmed")

ci95 <- function(y) {
  y <- na.omit(y); m <- mean(y); se <- sd(y) / sqrt(length(y))
  data.frame(y = m, ymin = m - qt(0.975, length(y) - 1) * se,
             ymax = m + qt(0.975, length(y) - 1) * se)}

plot_outcome <- function(y, ylab) {
  ggplot(dat, aes(x = fat, y = .data[[y]])) +
    geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1),
               colour = "grey60", alpha = 0.4, size = 2.5) +
    stat_summary(fun.data = ci95, geom = "pointrange",
                 colour = "black", size = 0.8, linewidth = 0.8) +
    facet_wrap(~ temp, labeller = as_labeller(temp_labs)) +
    scale_x_discrete(labels = function(x) paste0(x, "%")) +
    labs(x = "Milk fat", y = ylab) +
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", linewidth = 0.8),
          strip.background = element_blank())}

plot_outcome("foam_t0", "Foam volume at t0 (ml)")
plot_outcome("stability", "Share of foam remaining after 10 min")
plot_outcome("foam_share", "Share of foam")
plot_outcome("texture", "Texture score (1-5)")
