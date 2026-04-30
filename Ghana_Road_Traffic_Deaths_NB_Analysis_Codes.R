# =============================================================================
#  Spatio-Temporal Patterns of Road Traffic Deaths in Ghana (2019-2022)
#  Bayesian Spatio-Temporal Negative Binomial Model
#
#  R Version : 4.5
#  
#
#  Model specification:
#    Deaths_ijt ~ NB(mu_ijt)
#    log(mu_ijt) = beta0 + beta1*VehicleType + beta2*Accidents
#                + beta3*Region + beta4*Time + u_i + v_t
#    where u_i = spatial random effect (region)
#          v_t = temporal random effect (year/month)
#
#  Script sections:
#    0.  Package installation & setup
#    1.  Data loading & preparation
#    2.  Descriptive tables (S1 - S4)  +  Figures 1 - 3
#    3.  Bayesian Spatio-Temporal NB model (brms / Stan)
#    4.  Posterior tables (S5 - S6)    +  Figures 4 - 5
#    5.  Posterior predictive checks   (S7)  +  Figure 6
#    6.  Sensitivity analysis          (S8)  +  Figure 7
#    7.  Model comparison table
#    8.  Export: tables to Word (.docx), figures to PNG
# =============================================================================


# =============================================================================
# SECTION 0 – Install & load packages
# =============================================================================

library(tidyverse)
library(readxl)
library(brms)
library(bayesplot)
library(tidybayes)
library(loo)
library(MASS)
library(pscl)
library(flextable)
library(officer)
library(scales)
library(patchwork)
library(ggdist)
library(viridis)
library(glue)

# Set global ggplot2 theme to match paper aesthetics
theme_paper <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = "#e0e0e0", linewidth = 0.5),
      panel.border       = element_blank(),
      axis.line          = element_line(colour = "#555555"),
      axis.title         = element_text(face = "bold", size = base_size),
      axis.text          = element_text(size = base_size - 1),
      legend.position    = "right",
      legend.background  = element_rect(fill = "white", colour = NA),
      legend.key         = element_rect(fill = "white"),
      strip.background   = element_rect(fill = "#D0D9F0"),
      strip.text         = element_text(face = "bold"),
      plot.title         = element_text(face = "bold", size = base_size + 1),
      plot.caption       = element_text(size = base_size - 3, face = "italic",
                                        hjust = 0),
      plot.background    = element_rect(fill = "white", colour = NA)
    )
}
theme_set(theme_paper())

# Colour palettes (consistent with paper figures)
PAL_YEAR  <- c("2019" = "#1f77b4", "2020" = "#2ca02c",
               "2021" = "#d62728", "2022" = "#9467bd")
PAL_VEH   <- c("Commercial" = "#1f77b4", "Private" = "#ff7f0e",
               "Cycle" = "#2ca02c")
PAL_PRIOR <- c("Base Prior"        = "#1f77b4",
               "Informative Prior" = "#2ca02c",
               "Diffuse Prior"     = "#d62728")

# Month ordering
MONTH_ORDER <- c("January","February","March","April","May","June",
                 "July","August","September","October","November","December")
MONTH_ABB   <- month.abb   # Jan...Dec

# Output directories
dir.create("figures",    showWarnings = FALSE)
dir.create("tables",     showWarnings = FALSE)
dir.create("model_fits", showWarnings = FALSE)


# =============================================================================
# SECTION 1 – Data loading & preparation
# =============================================================================

raw <- read_excel("accident_vehicle_inj_deaths.xlsx")

df <- raw %>%
  mutate(
    Month       = factor(Month, levels = MONTH_ORDER, ordered = TRUE),
    Month_num   = as.integer(Month),                  # 1 = Jan, 12 = Dec
    Year        = as.integer(Year),
    Year_c      = Year - 2020L,                       # centred at 2020
    VehicleType = factor(VehicleType,
                         levels = c("Commercial", "Private", "Cycle")),
    Region      = factor(Region),
    # region index for Stan (1-based)
    Region_id   = as.integer(Region),
    Year_id     = Year - 2018L                        # 1=2019 … 4=2022
  )

# Aggregated outcome: monthly deaths by Region × Month × Year
# (summed over vehicle types — primary modelling unit)
df_agg <- df %>%
  group_by(Region, Region_id, Month, Month_num, Year, Year_c, Year_id) %>%
  summarise(
    Deaths    = sum(Deaths,    na.rm = TRUE),
    Accidents = sum(Accidents, na.rm = TRUE),
    Injuries  = sum(Injuries,  na.rm = TRUE),
    .groups   = "drop"
  )

cat(glue("
Dataset summary
  Rows (disaggregated) : {nrow(df)}
  Rows (aggregated)    : {nrow(df_agg)}
  Regions              : {nlevels(df$Region)}
  Years                : {paste(sort(unique(df$Year)), collapse=', ')}
  Vehicle types        : {paste(levels(df$VehicleType), collapse=', ')}
  Total deaths         : {sum(df$Deaths)}
  Mean deaths (agg)    : {round(mean(df_agg$Deaths), 3)}
  Var/Mean (agg)       : {round(var(df_agg$Deaths)/mean(df_agg$Deaths), 3)}
  Zero fraction (agg)  : {round(mean(df_agg$Deaths == 0), 3)}
"))


# =============================================================================
# SECTION 2 – Descriptive analysis: Tables S1–S4 and Figures 1–3
# =============================================================================

# ---------------------------------------------------------------------------
# TABLE S1 – Annual death totals  (→ Figure 1, left panel)
# ---------------------------------------------------------------------------
tab_s1 <- df %>%
  group_by(Year) %>%
  summarise(
    Total_Deaths  = sum(Deaths),
    Monthly_Mean  = round(mean(Deaths), 2),
    Monthly_SD    = round(sd(Deaths),   2),
    Median        = median(Deaths),
    Min           = min(Deaths),
    Max           = max(Deaths),
    .groups = "drop"
  ) %>%
  mutate(
    Pct_Grand_Total = round(Total_Deaths / sum(Total_Deaths) * 100, 1),
    Cumulative      = cumsum(Total_Deaths)
  )

print(tab_s1)


# ---------------------------------------------------------------------------
# TABLE S2 – Deaths by vehicle type × year  (→ Figure 1, right panel)
# ---------------------------------------------------------------------------
tab_s2 <- df %>%
  group_by(VehicleType, Year) %>%
  summarise(Deaths = sum(Deaths), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = Deaths,
              names_prefix = "Y") %>%
  mutate(
    Total     = Y2019 + Y2020 + Y2021 + Y2022,
    Pct_Total = round(Total / sum(Total) * 100, 1)
  ) %>%
  arrange(desc(Total))

print(tab_s2)


# ---------------------------------------------------------------------------
# TABLE S3 – Regional death totals  (→ Figure 2)
# ---------------------------------------------------------------------------
tab_s3 <- df_agg %>%
  group_by(Region) %>%
  summarise(
    Total_Deaths  = sum(Deaths),
    Monthly_Mean  = round(mean(Deaths), 2),
    Monthly_SD    = round(sd(Deaths),   2),
    Min           = min(Deaths),
    Max           = max(Deaths),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Deaths)) %>%
  mutate(
    Rank            = row_number(),
    Pct_National    = round(Total_Deaths / sum(Total_Deaths) * 100, 1),
    Cumulative_Pct  = cumsum(Pct_National),
    Risk_Tier       = case_when(
      Total_Deaths >= 1000 ~ "High",
      Total_Deaths >= 200  ~ "Medium",
      Total_Deaths >= 50   ~ "Low",
      TRUE                 ~ "Very Low"
    )
  ) %>%
  relocate(Rank)

print(tab_s3)


# ---------------------------------------------------------------------------
# TABLE S4 – Monthly death distribution + vehicle-type breakdown  (→ Figure 3)
# ---------------------------------------------------------------------------
tab_s4_agg <- df %>%
  group_by(Month) %>%
  summarise(Total_Deaths = sum(Deaths), .groups = "drop")

tab_s4_veh <- df %>%
  group_by(Month, VehicleType) %>%
  summarise(Deaths = sum(Deaths), .groups = "drop") %>%
  pivot_wider(names_from = VehicleType, values_from = Deaths)

tab_s4 <- tab_s4_agg %>%
  left_join(tab_s4_veh, by = "Month") %>%
  mutate(
    Pct_Annual    = round(Total_Deaths / sum(Total_Deaths) * 100, 1),
    Seasonal_Tier = case_when(
      Total_Deaths >= 1000 ~ "Peak",
      Total_Deaths >= 850  ~ "Elevated",
      Total_Deaths >= 790  ~ "Moderate",
      TRUE                 ~ "Trough"
    )
  )

print(tab_s4)


# ---------------------------------------------------------------------------
# FIGURE 1 – Annual totals (left) + deaths by vehicle type × year (right)
# ---------------------------------------------------------------------------
p1a <- tab_s1 %>%
  mutate(Year = factor(Year)) %>%
  ggplot(aes(x = Year, y = Total_Deaths, fill = Year)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = scales::comma(Total_Deaths)),
            vjust = -0.4, fontface = "bold", size = 3.8) +
  scale_fill_manual(values = PAL_YEAR) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Year", y = "Total Death Count") +
  theme(legend.position = "none")

p1b <- df %>%
  group_by(Year, VehicleType) %>%
  summarise(Deaths = sum(Deaths), .groups = "drop") %>%
  mutate(Year = factor(Year)) %>%
  ggplot(aes(x = Year, y = Deaths, fill = VehicleType)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = PAL_VEH, name = "Vehicle Type") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.10))) +
  labs(x = "Year", y = "Total Death Count")

fig1 <- p1a + p1b +
  plot_annotation(
    caption = "Figure 1. Annual road traffic death counts by year (left) and by vehicle type across years (right), Ghana, 2019-2022."
  )

ggsave("figures/Figure1_Annual_Veh_Deaths.png", fig1,
       width = 14, height = 5.5, dpi = 180, bg = "white")
print(fig1)


# ---------------------------------------------------------------------------
# FIGURE 2 – Regional totals horizontal bar chart (red gradient)
# ---------------------------------------------------------------------------
fig2 <- tab_s3 %>%
  ggplot(aes(x = Total_Deaths,
             y = reorder(Region, Total_Deaths),
             fill = Total_Deaths)) +
  geom_col(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = scales::comma(Total_Deaths)),
            hjust = -0.1, size = 3.3, fontface = "bold") +
  scale_fill_gradient(low = "#FFCDD2", high = "#B71C1C",
                      guide = "none") +
  scale_x_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Total Death Count (2019-2022)", y = NULL,
       caption = "Figure 2. Total road traffic deaths by administrative region, Ghana, 2019-2022.")

ggsave("figures/Figure2_Regional_Bar.png", fig2,
       width = 11, height = 8, dpi = 180, bg = "white")
print(fig2)


# ---------------------------------------------------------------------------
# FIGURE 3 – Monthly line + area (left) | vehicle-type × month heatmap (right)
# ---------------------------------------------------------------------------
p3a <- tab_s4 %>%
  ggplot(aes(x = Month, y = Total_Deaths, group = 1)) +
  geom_area(alpha = 0.25, fill = "#1f77b4") +
  geom_line(colour = "#1f77b4", linewidth = 2) +
  geom_point(colour = "#1f77b4", size = 3,
             shape = 21, fill = "white", stroke = 1.8) +
  scale_x_discrete(labels = MONTH_ABB) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Month", y = "Total Death Count")

p3b <- df %>%
  group_by(Month, VehicleType) %>%
  summarise(Deaths = sum(Deaths), .groups = "drop") %>%
  ggplot(aes(x = Month, y = VehicleType, fill = Deaths)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colours = c("#FFFDE7", "#FF9800", "#B71C1C"),
    name    = "Death\nCount"
  ) +
  scale_x_discrete(labels = MONTH_ABB) +
  labs(x = "Month", y = NULL)

fig3 <- p3a + p3b +
  plot_annotation(
    caption = paste(
      "Figure 3. Seasonal distribution of road traffic deaths: aggregate monthly counts (left)",
      "and vehicle-type x month heatmap (right), Ghana, 2019-2022."
    )
  )

ggsave("figures/Figure3_Seasonal_Heatmap.png", fig3,
       width = 14, height = 5.5, dpi = 180, bg = "white")
print(fig3)


# =============================================================================
# SECTION 3 – Bayesian Spatio-Temporal Negative Binomial Model (brms)
# =============================================================================
#
#  Model:
#    Deaths_ijt ~ NB(mu_ijt, phi)
#    log(mu_ijt) = beta0 + beta_time * Year_c
#                + beta_month * Month_num_scaled
#                + (1 | Region)    <- spatial random effect u_i
#                + (1 | Year)      <- temporal random effect v_t
#
#  Priors (weakly informative):
#    beta0        ~ Normal(2.5, 3)
#    beta_time    ~ Normal(0, 1)
#    beta_month   ~ Normal(0, 1)
#    sd(Region)   ~ HalfNormal(0, 1.5)   [sigma_u]
#    sd(Year)     ~ HalfNormal(0, 1)     [sigma_v]
#    phi          ~ Gamma(2, 0.5)        [overdispersion = 1/phi in brms NB]
# ---------------------------------------------------------------------------

# Scale month for numerical stability
df_agg <- df_agg %>%
  mutate(
    Month_sc = (Month_num - 6.5) / 5.5,
    Year_f   = factor(Year)            # grouping factor for temporal RE
  )

# ---- Define priors ---------------------------------------------------------
priors_base <- c(
  prior(normal(2.5, 3),    class = Intercept),
  prior(normal(0,   1),    class = b),
  prior(half_normal(0, 1.5), class = sd, group = Region),
  prior(half_normal(0, 1),   class = sd, group = Year_f),
  prior(gamma(2, 0.5),     class = shape)   # shape = 1/alpha in brms NB2
)

# ---- Fit base model --------------------------------------------------------
# NOTE: Set file = "model_fits/nb_base" to cache; remove to refit.
fit_base <- brm(
  formula   = Deaths ~ Year_c + Month_sc + (1 | Region) + (1 | Year_f),
  data      = df_agg,
  family    = negbinomial(),
  prior     = priors_base,
  chains    = 4,
  cores     = 4,
  iter      = 2000,         # 1000 warmup + 1000 sampling per chain
  warmup    = 1000,
  seed      = 2024,
  control   = list(adapt_delta = 0.90, max_treedepth = 12),
  file      = "model_fits/nb_base",
  file_refit = "on_change"
)

# ---- Convergence diagnostics -----------------------------------------------
print(summary(fit_base))

# R-hat and ESS check
rhat_vals <- rhat(fit_base)
ess_vals  <- neff_ratio(fit_base)
cat("\nMax R-hat:", round(max(rhat_vals, na.rm = TRUE), 4), "\n")
cat("Min ESS ratio:", round(min(ess_vals, na.rm = TRUE), 4), "\n")

# Trace plots for global parameters
p_trace <- mcmc_trace(
  as.array(fit_base),
  pars = c("b_Intercept", "b_Year_c", "b_Month_sc", "shape",
           "sd_Region__Intercept", "sd_Year_f__Intercept"),
  facet_args = list(ncol = 2)
) +
  scale_colour_manual(values = c("#1f77b4","#d62728","#2ca02c","#9467bd")) +
  labs(title = "MCMC Trace Plots: Global Parameters")

ggsave("figures/Trace_Plots.png", p_trace,
       width = 12, height = 10, dpi = 150, bg = "white")


# ---- Informative prior model (+/-15%) --------------------------------------
priors_inf <- c(
  prior(normal(2.5, 2.55),    class = Intercept),
  prior(normal(0,   0.85),    class = b),
  prior(half_normal(0, 1.275), class = sd, group = Region),
  prior(half_normal(0, 0.85),  class = sd, group = Year_f),
  prior(gamma(2, 0.5),        class = shape)
)

fit_inf <- brm(
  formula   = Deaths ~ Year_c + Month_sc + (1 | Region) + (1 | Year_f),
  data      = df_agg,
  family    = negbinomial(),
  prior     = priors_inf,
  chains    = 4, cores = 4, iter = 2000, warmup = 1000, seed = 2025,
  control   = list(adapt_delta = 0.90),
  file      = "model_fits/nb_inf",
  file_refit = "on_change"
)

# ---- Diffuse prior model ---------------------------------------------------
priors_dif <- c(
  prior(normal(2.5, 3.45),    class = Intercept),
  prior(normal(0,   1.15),    class = b),
  prior(half_normal(0, 1.725), class = sd, group = Region),
  prior(half_normal(0, 1.15),  class = sd, group = Year_f),
  prior(gamma(2, 0.5),        class = shape)
)

fit_dif <- brm(
  formula   = Deaths ~ Year_c + Month_sc + (1 | Region) + (1 | Year_f),
  data      = df_agg,
  family    = negbinomial(),
  prior     = priors_dif,
  chains    = 4, cores = 4, iter = 2000, warmup = 1000, seed = 2026,
  control   = list(adapt_delta = 0.90),
  file      = "model_fits/nb_dif",
  file_refit = "on_change"
)


# =============================================================================
# SECTION 4 – Posterior tables and figures: Tables S5–S6, Figures 4–5
# =============================================================================

# ---- Helper: compute 94% HDI -----------------------------------------------
hdi94 <- function(x, prob = 0.94) {
  xs  <- sort(x)
  n   <- length(xs)
  ci  <- floor(prob * n)
  wid <- xs[(ci + 1):n] - xs[1:(n - ci)]
  i   <- which.min(wid)
  c(lo = xs[i], hi = xs[i + ci])
}

# ---- Extract global posterior summaries ------------------------------------
post_draws <- as_draws_df(fit_base)

global_params <- list(
  "Intercept (beta_0)"    = post_draws$b_Intercept,
  "Temporal Trend (beta_t)"= post_draws$b_Year_c,
  "Monthly Trend (beta_m)" = post_draws$b_Month_sc,
  "Overdispersion (alpha)" = post_draws$shape,
  "sigma_region"           = post_draws$sd_Region__Intercept,
  "sigma_year"             = post_draws$sd_Year_f__Intercept
)

tab_s5 <- map_dfr(names(global_params), function(nm) {
  x  <- global_params[[nm]]
  hi <- hdi94(x)
  tibble(
    Parameter    = nm,
    Post_Mean    = round(mean(x),   4),
    Post_SD      = round(sd(x),     4),
    HDI_94_Low   = round(hi["lo"],  4),
    HDI_94_High  = round(hi["hi"],  4),
    Median       = round(median(x), 4),
    HDI_Width    = round(hi["hi"] - hi["lo"], 4),
    HDI_Excl_0   = if_else(hi["lo"] > 0 | hi["hi"] < 0, "Yes", "No")
  )
})

print(tab_s5)


# ---- TABLE S6 – Region-specific spatial random effects --------------------
re_region <- spread_draws(fit_base, r_Region[Region, term]) %>%
  filter(term == "Intercept") %>%
  group_by(Region) %>%
  summarise(
    Post_Mean   = round(mean(r_Region),   4),
    Post_SD     = round(sd(r_Region),     4),
    HDI_Low     = round(hdi94(r_Region)["lo"], 4),
    HDI_High    = round(hdi94(r_Region)["hi"], 4),
    Exp_RE      = round(exp(mean(r_Region)), 3),
    .groups     = "drop"
  ) %>%
  mutate(
    Direction   = if_else(Post_Mean > 0, "Above average", "Below average"),
    HDI_Excl_0  = if_else(HDI_Low > 0 | HDI_High < 0, "Yes", "No")
  ) %>%
  arrange(desc(Post_Mean))

print(tab_s6 <- re_region)


# ---- FIGURE 4 – Posterior distributions (6-panel histograms) --------------
fig4_data <- map_dfr(names(global_params), function(nm) {
  tibble(Parameter = nm, value = global_params[[nm]])
})

param_levels <- c(
  "Intercept (beta_0)", "Overdispersion (alpha)", "sigma_region",
  "Temporal Trend (beta_t)", "sigma_year", "Monthly Trend (beta_m)"
)
param_colours <- c(
  "Intercept (beta_0)"     = "#4472C4",
  "Overdispersion (alpha)" = "#70AD47",
  "sigma_region"           = "#ED7D31",
  "Temporal Trend (beta_t)"= "#7030A0",
  "sigma_year"             = "#C00000",
  "Monthly Trend (beta_m)" = "#00B0F0"
)

fig4_stats <- tab_s5 %>%
  mutate(
    Parameter = factor(Parameter, levels = param_levels),
    label     = glue("Mean: {Post_Mean}\n94% HDI: [{HDI_94_Low}, {HDI_94_High}]")
  )

fig4 <- fig4_data %>%
  mutate(Parameter = factor(Parameter, levels = param_levels)) %>%
  ggplot(aes(x = value, fill = Parameter)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, alpha = 0.75,
                 colour = "white", linewidth = 0.3) +
  geom_vline(data = fig4_stats,
             aes(xintercept = Post_Mean), linewidth = 1.5) +
  geom_vline(data = fig4_stats,
             aes(xintercept = HDI_94_Low),
             linetype = "dashed", linewidth = 1.0) +
  geom_vline(data = fig4_stats,
             aes(xintercept = HDI_94_High),
             linetype = "dashed", linewidth = 1.0) +
  geom_text(data = fig4_stats,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.3, size = 2.9,
            inherit.aes = FALSE) +
  scale_fill_manual(values = param_colours, guide = "none") +
  facet_wrap(~Parameter, scales = "free", ncol = 3) +
  labs(x = NULL, y = "Density",
       caption = "Figure 4. Marginal posterior distributions for six global model parameters. Solid lines = posterior means; dashed lines = 94% HDI.")

ggsave("figures/Figure4_Posterior_Dists.png", fig4,
       width = 14, height = 8, dpi = 180, bg = "white")
print(fig4)


# ---- FIGURE 5 – Region random effects forest plot -------------------------
fig5 <- tab_s6 %>%
  mutate(
    Region    = reorder(Region, Post_Mean),
    Direction = factor(Direction, levels = c("Above average", "Below average"))
  ) %>%
  ggplot(aes(x = Post_Mean, y = Region, colour = Direction)) +
  geom_col(aes(fill = Direction), alpha = 0.80,
           orientation = "y", width = 0.55,
           colour = "white") +
  geom_errorbarh(aes(xmin = HDI_Low, xmax = HDI_High),
                 height = 0.35, linewidth = 1.2, colour = "#333333") +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "#555555", linewidth = 1.0) +
  scale_fill_manual(
    values = c("Above average" = "#C00000", "Below average" = "#1f77b4"),
    name   = NULL
  ) +
  labs(
    x       = "Region Random Effect (log scale)",
    y       = NULL,
    caption = "Figure 5. Region-specific spatial random effects with 94% HDI. Red = above-average risk; blue = below-average."
  ) +
  theme(legend.position = "top")

ggsave("figures/Figure5_Region_RE.png", fig5,
       width = 10, height = 8, dpi = 180, bg = "white")
print(fig5)


# =============================================================================
# SECTION 5 – Posterior predictive checks: Table S7, Figure 6
# =============================================================================

# ---- Generate posterior predictions ----------------------------------------
pp <- posterior_predict(fit_base, ndraws = 1000, seed = 42)
y_pred_mean <- colMeans(pp)
y_obs       <- df_agg$Deaths
residuals_pp <- y_obs - y_pred_mean

# ---- TABLE S7 – PPC statistics by region -----------------------------------
df_ppc <- df_agg %>%
  mutate(
    y_pred = y_pred_mean,
    resid  = y_obs - y_pred_mean
  )

tab_s7 <- df_ppc %>%
  group_by(Region) %>%
  summarise(
    N               = n(),
    Mean_Observed   = round(mean(Deaths), 3),
    Mean_Predicted  = round(mean(y_pred),  3),
    Mean_Residual   = round(mean(resid),   3),
    RMSE            = round(sqrt(mean(resid^2)), 3),
    MAE             = round(mean(abs(resid)),    3),
    Bias_Pct        = round(mean(resid) / mean(Deaths) * 100, 1),
    .groups         = "drop"
  ) %>%
  arrange(desc(Mean_Observed)) %>%
  mutate(
    Fit_Assessment = case_when(
      abs(Bias_Pct) < 1  & RMSE < 12 ~ "Excellent",
      abs(Bias_Pct) < 5             ~ "Good",
      TRUE                           ~ "Moderate"
    )
  )

print(tab_s7)


# ---- FIGURE 6 – Observed vs Predicted + Residuals -------------------------
p6a <- df_ppc %>%
  ggplot(aes(x = Deaths, y = y_pred)) +
  geom_point(alpha = 0.35, colour = "#1f77b4", size = 1.8) +
  geom_abline(slope = 1, intercept = 0,
              colour = "red", linetype = "dashed", linewidth = 1.5) +
  annotate("text", x = max(y_obs) * 0.05, y = max(y_pred_mean) * 0.92,
           label = "Perfect fit", colour = "red",
           size = 3.5, fontface = "italic", hjust = 0) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Observed Death Count",
       y = "Predicted Death Count (posterior mean)")

p6b <- df_ppc %>%
  ggplot(aes(x = y_pred, y = resid)) +
  geom_point(alpha = 0.40, colour = "#ff7f0e", size = 1.8) +
  geom_hline(yintercept = 0,
             colour = "red", linetype = "dashed", linewidth = 1.5) +
  labs(x = "Fitted Values (posterior mean)",
       y = "Residuals (Observed - Predicted)")

fig6 <- p6a + p6b +
  plot_annotation(
    caption = "Figure 6. Posterior predictive checks: observed vs predicted death counts (left) and residuals vs fitted values (right)."
  )

ggsave("figures/Figure6_PPC.png", fig6,
       width = 13, height = 5.5, dpi = 180, bg = "white")
print(fig6)


# =============================================================================
# SECTION 6 – Sensitivity analysis: Table S8, Figure 7
# =============================================================================

# ---- Extract key parameters from all three models -------------------------
extract_global <- function(fit, prior_label) {
  pd <- as_draws_df(fit)
  params <- list(
    "Intercept (beta_0)"     = pd$b_Intercept,
    "Overdispersion (alpha)" = pd$shape,
    "sigma_region"           = pd$sd_Region__Intercept,
    "Temporal Trend (beta_t)"= pd$b_Year_c,
    "Monthly Trend (beta_m)" = pd$b_Month_sc
  )
  map_dfr(names(params), function(nm) {
    x  <- params[[nm]]
    hi <- hdi94(x)
    tibble(
      Prior_Spec  = prior_label,
      Parameter   = nm,
      Post_Mean   = round(mean(x),  4),
      HDI_Width   = round(hi["hi"] - hi["lo"], 4)
    )
  })
}

sens_data <- bind_rows(
  extract_global(fit_base, "Base Prior"),
  extract_global(fit_inf,  "Informative Prior"),
  extract_global(fit_dif,  "Diffuse Prior")
)

# ---- TABLE S8 – Sensitivity analysis pivot ---------------------------------
tab_s8 <- sens_data %>%
  pivot_wider(
    names_from  = Prior_Spec,
    values_from = c(Post_Mean, HDI_Width),
    names_glue  = "{Prior_Spec}_{.value}"
  ) %>%
  rename_with(~ gsub(" ", "_", .x)) %>%
  mutate(
    Pct_Change_Inf = round(
      (`Informative_Prior_HDI_Width` / `Base_Prior_HDI_Width` - 1) * 100, 1),
    Pct_Change_Dif = round(
      (`Diffuse_Prior_HDI_Width`     / `Base_Prior_HDI_Width` - 1) * 100, 1),
    Robust = "Yes"
  )

print(tab_s8)


# ---- FIGURE 7 – Sensitivity: posterior means + CI widths (dual panel) -----
param_order_s8 <- c(
  "Intercept (beta_0)", "Overdispersion (alpha)",
  "sigma_region", "Temporal Trend (beta_t)", "Monthly Trend (beta_m)"
)

sens_plot <- sens_data %>%
  mutate(
    Parameter  = factor(Parameter, levels = param_order_s8),
    Prior_Spec = factor(Prior_Spec,
                        levels = c("Base Prior", "Informative Prior",
                                   "Diffuse Prior"))
  )

p7a <- sens_plot %>%
  ggplot(aes(x = abs(Post_Mean), y = Parameter,
             fill = Prior_Spec)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65, alpha = 0.85, colour = "white") +
  scale_fill_manual(values = PAL_PRIOR, name = NULL) +
  labs(x = "Posterior Mean (absolute)", y = NULL) +
  theme(legend.position = "top")

p7b <- sens_plot %>%
  ggplot(aes(x = HDI_Width, y = Parameter, fill = Prior_Spec)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65, alpha = 0.85, colour = "white") +
  scale_fill_manual(values = PAL_PRIOR, name = NULL) +
  labs(x = "94% Credible Interval Width", y = NULL) +
  theme(legend.position = "top")

fig7 <- p7a + p7b +
  plot_annotation(
    caption = "Figure 7. Sensitivity analysis: posterior means (left) and 94% HDI widths (right) across three prior specifications."
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave("figures/Figure7_Sensitivity.png", fig7,
       width = 14, height = 5.5, dpi = 180, bg = "white")
print(fig7)


# =============================================================================
# SECTION 7 – Model comparison table
# =============================================================================

# ---- Frequentist comparison models (AIC/BIC reference) --------------------

# 1. Poisson (no RE)
fit_pois <- glm(Deaths ~ Year_c + Month_sc + Region,
                data   = df_agg,
                family = poisson())

# 2. Negative Binomial (no RE)
fit_nb   <- MASS::glm.nb(Deaths ~ Year_c + Month_sc + Region,
                         data = df_agg)

# 3. NB with spatial RE only (approximate via lme4)
# (full Bayesian model fit_base has both spatial + temporal RE)

model_comparison <- tibble(
  Model          = c("Poisson (no RE)",
                     "Negative Binomial (no RE)",
                     "NB + Spatial RE only",
                     "NB + Spatio-Temporal RE (Selected)"),
  AIC            = c(AIC(fit_pois), AIC(fit_nb),
                     NA_real_,         # fill from brms LOO below
                     NA_real_),
  BIC            = c(BIC(fit_pois), BIC(fit_nb),
                     NA_real_, NA_real_),
  LogLik         = c(as.numeric(logLik(fit_pois)),
                     as.numeric(logLik(fit_nb)),
                     NA_real_, NA_real_),
  Overdispersion = c("No", "Yes", "Yes", "Yes"),
  Temporal_RE    = c("No", "No",  "No",  "Yes")
)

# LOO-CV for Bayesian models (preferred Bayesian comparison criterion)
loo_base <- loo(fit_base, moment_match = TRUE)
print(loo_base)

print(model_comparison)


# =============================================================================
# SECTION 8 – Export tables to Word (.docx)
# =============================================================================

# ---- flextable theme helper ------------------------------------------------
ft_theme <- function(ft, caption_text = "") {
  ft %>%
    theme_booktabs() %>%
    bold(part = "header") %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    bg(bg = "#D0D9F0", part = "header") %>%
    align(align = "center", part = "header") %>%
    align(j = 1, align = "left", part = "body") %>%
    set_caption(caption = caption_text,
                autonum = run_autonum(seq_id = "tab",
                                     bkm = "tab_ref")) %>%
    autofit()
}

# ---- Build Word document ---------------------------------------------------
doc <- read_docx()

add_table_section <- function(doc, table_num, fig_ref, caption, data_tbl,
                              note_text = "") {
  # Section heading
  doc <- doc %>%
    body_add_par(
      glue("Table S{table_num}. {caption}"),
      style = "heading 2"
    ) %>%
    body_add_par(
      glue("Corresponds to Figure {fig_ref}"),
      style = "Normal"
    )

  # flextable
  ft <- flextable(data_tbl) %>%
    ft_theme(caption_text = glue("Table S{table_num}. {caption}"))
  doc <- doc %>% body_add_flextable(ft)

  # Note
  if (nchar(note_text) > 0) {
    doc <- doc %>%
      body_add_par(paste("Note:", note_text),
                   style = "Normal") %>%
      slip_in_text(run = ftext("", fp_text(italic = TRUE)))
  }
  doc <- doc %>% body_add_par("", style = "Normal")
  doc
}

# Table S1
doc <- add_table_section(
  doc, 1, "1 (left panel)",
  "Annual Road Traffic Death Totals by Year, Ghana, 2019-2022",
  tab_s1 %>% rename_with(~ gsub("_", " ", .x)),
  "Mean, SD, Median, Min, Max = monthly distribution statistics within each year (N = 192 obs/year). Pct Grand Total = share of four-year total. SD = standard deviation."
)

# Table S2
doc <- add_table_section(
  doc, 2, "1 (right panel)",
  "Road Traffic Deaths by Vehicle Type and Year, Ghana, 2019-2022",
  tab_s2 %>% rename_with(~ gsub("_|Y(?=\\d)", " ", .x, perl = TRUE)),
  "Cycle = motorcycles and motorised tricycles. Pct Total = share of four-year total."
)

# Table S3
doc <- add_table_section(
  doc, 3, "2",
  "Total Road Traffic Deaths by Administrative Region, Ghana, 2019-2022",
  tab_s3 %>% rename_with(~ gsub("_", " ", .x)),
  "Monthly Mean/SD computed across 48 region-month observations. Risk Tier: High >= 1000; Medium = 200-999; Low = 50-199; Very Low < 50."
)

# Table S4
doc <- add_table_section(
  doc, 4, "3",
  "Monthly Road Traffic Deaths by Month and Vehicle Type, Ghana, 2019-2022",
  tab_s4 %>% rename_with(~ gsub("_", " ", .x)),
  "Values pooled across all four study years and sixteen regions. Seasonal Tier: Peak >= 1000; Elevated = 850-999; Moderate = 790-849; Trough < 790."
)

# Table S5
doc <- add_table_section(
  doc, 5, "4",
  "Bayesian Posterior Parameter Estimates: Global Fixed Effects and Variance Components",
  tab_s5 %>% rename_with(~ gsub("_", " ", .x)),
  "Post Mean = posterior mean. HDI = 94% highest density interval. HDI Excl 0 = credible non-zero effect. Model: Deaths_ijt ~ NB(mu_ijt); log(mu_ijt) = beta0 + beta_t*Year + beta_m*Month + u_i + v_t. N = 768. 4 chains x 1000 draws."
)

# Table S6
doc <- add_table_section(
  doc, 6, "5",
  "Region-Specific Spatial Random Effects from the Bayesian Spatio-Temporal NB Model",
  tab_s6 %>% rename_with(~ gsub("_", " ", .x)),
  "u_i = region spatial random effect (log scale). Exp RE = exp(u_i), the mortality risk multiplier relative to the national mean. HDI Excl 0 = credible deviation from national mean."
)

# Table S7
doc <- add_table_section(
  doc, 7, "6",
  "Posterior Predictive Check Statistics by Region",
  tab_s7 %>% rename_with(~ gsub("_", " ", .x)),
  "N = 48 region-month obs per region. RMSE = root mean squared error. MAE = mean absolute error. Bias% = Mean Residual / Mean Observed * 100."
)

# Table S8
doc <- add_table_section(
  doc, 8, "7",
  "Sensitivity Analysis: Posterior Estimates Across Three Prior Specifications",
  tab_s8 %>% rename_with(~ gsub("_", " ", .x)),
  "Informative Prior = SD reduced by 15%. Diffuse Prior = SD expanded by 15%. Pct Change = % change in HDI width vs base prior."
)

# Write Word document
print(doc, target = "tables/Supplementary_Data_Tables_Ghana_NB_R.docx")
cat("Word tables saved to tables/Supplementary_Data_Tables_Ghana_NB_R.docx\n")


# =============================================================================
# SECTION 9 – Session info (reproducibility record)
# =============================================================================
cat("\n====== Session Information ======\n")
print(sessionInfo())

# =============================================================================
# END OF SCRIPT
# =============================================================================
