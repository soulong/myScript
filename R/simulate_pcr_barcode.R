# simulate_pcr_barcode.R
#
# Mimic a multiplexed barcode PCR + NGS quantification experiment.
#
# Model (recursive per-cycle, "continuous conversion" form):
#   - gDNA templates persist and keep producing amplicons every cycle with
#     efficiency alpha (scalar or per-barcode vector)
#   - amplicons amplify with efficiency beta (scalar or per-barcode vector)
#   - optional saturation: eff_cycle = base * (1 - A_prev / K)
#   - optional stochasticity: Poisson input draw + binomial amplification +
#     multinomial read draw
#
# Deterministic closed form (scalar alpha/beta, no saturation):
#   A(n) = G * alpha * ((1 + beta)^n - 1) / beta
#
# Usage (demo): Rscript simulate_pcr_barcode.R

script_dir <- if (interactive()) getwd() else {
  args <- commandArgs(FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file_arg) && nzchar(file_arg)) dirname(normalizePath(file_arg)) else getwd()
}
out <- function(f) file.path(script_dir, f)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Core simulation engine
# -----------------------------------------------------------------------------
simulate_pcr <- function(ratios       = c(100, 10, 1, 0.1, 0.01),
                         total_gDNA   = 1e6,
                         alpha        = 0.3,   # scalar or per-barcode vector
                         beta         = 0.9,   # scalar or per-barcode vector
                         cycles       = 35,
                         saturate     = FALSE,
                         K            = 1e13,  # plateau capacity (amplicons)
                         stochastic   = FALSE,
                         reps         = 20,
                         reads        = NULL,  # NGS reads per replicate (NULL = skip)
                         scenario     = "A") {

  n   <- length(ratios)
  barcode <- factor(paste0("B", seq_len(n), " (", ratios, ")"),
                    levels = paste0("B", seq_len(n), " (", ratios, ")"))

  alpha <- rep(alpha, length.out = n)
  beta  <- rep(beta,  length.out = n)

  expected <- total_gDNA * ratios / sum(ratios)
  input_prop <- expected / sum(expected)

  # --- one replicate, returns (A per cycle, final read counts) ---------------
  one_rep <- function(G0) {
    A      <- numeric(n)
    A_track <- matrix(0, nrow = cycles, ncol = n)
    for (cyc in seq_len(cycles)) {
      scale <- if (saturate) pmax(1 - A / K, 0) else rep(1, n)
      a_eff <- alpha * scale
      b_eff <- beta  * scale
      if (stochastic) {
        A <- A + rbinom(n, size = G0, prob = a_eff) + rbinom(n, size = A, prob = b_eff)
      } else {
        A <- A + G0 * a_eff + A * b_eff
      }
      A_track[cyc, ] <- A
    }
    rd <- NULL
    if (!is.null(reads)) {
      if (sum(A) > 0) rd <- as.integer(rmultinom(1, size = reads, prob = A)) else rd <- rep(0L, n)
    }
    list(A = A_track, reads = rd, G0 = G0)
  }

  # --- run replicates --------------------------------------------------------
  nreps <- if (stochastic) reps else 1
  out <- vector("list", nreps)
  for (r in seq_len(nreps)) {
    G0 <- if (stochastic) rpois(n, expected) else expected
    out[[r]] <- one_rep(G0)
  }

  # --- assemble trajectory ---------------------------------------------------
  traj <- lapply(seq_len(nreps), function(r) {
    A  <- out[[r]]$A
    df <- as.data.frame(A)
    names(df) <- paste0("b", seq_len(n))
    df <- df %>%
      mutate(cycle = seq_len(cycles)) %>%
      pivot_longer(-cycle, names_to = "b", values_to = "amount") %>%
      mutate(bi = as.integer(sub("^b", "", b)),
             barcode = barcode[bi],
             rep = if (stochastic) r else 1L,
             scenario = scenario) %>%
      group_by(scenario, barcode, cycle, rep) %>%
      mutate(proportion = amount / sum(amount)) %>%
      ungroup() %>%
      select(scenario, barcode, cycle, rep, amount, proportion)
  })
  traj <- bind_rows(traj)

  # --- assemble final (observed) table --------------------------------------
  final <- lapply(seq_len(nreps), function(r) {
    A  <- out[[r]]$A[cycles, ]
    rd <- if (!is.null(reads)) out[[r]]$reads else A
    tibble(scenario  = scenario,
           barcode   = barcode,
           rep       = if (stochastic) r else 1L,
           input_molecules   = out[[r]]$G0,
           input_proportion  = out[[r]]$G0 / sum(out[[r]]$G0),
           amplicon_amount   = A,
           observed          = rd,
           observed_prop     = if (sum(rd) > 0) rd / sum(rd) else rep(0, n))
  })
  final <- bind_rows(final)

  list(trajectory = traj, final = final)
}

# -----------------------------------------------------------------------------
# Validation against the closed-form solution
# -----------------------------------------------------------------------------
validate_pcr <- function(alpha = 0.3, beta = 0.9, ratios = c(100, 10, 1, 0.1, 0.01),
                         total_gDNA = 1e6, cycles = 20) {
  n  <- length(ratios)
  G  <- total_gDNA * ratios / sum(ratios)
  sim <- simulate_pcr(ratios = ratios, total_gDNA = total_gDNA,
                      alpha = alpha, beta = beta, cycles = cycles,
                      saturate = FALSE, stochastic = FALSE, scenario = "val")

  closed <- as.vector(outer(seq_len(cycles), G,
                            function(cyc, G) G * alpha * ((1 + beta)^cyc - 1) / beta))
  got <- sim$trajectory %>%
    arrange(barcode, cycle) %>%
    pull(amount)
  rel_err <- abs(got - rep(closed, times = n)) / rep(closed, times = n)
  list(passed = max(rel_err) < 1e-9, max_rel_error = max(rel_err))
}

# -----------------------------------------------------------------------------
# Demo (runs only when executed directly via Rscript, not when sourced)
# -----------------------------------------------------------------------------
if ((sys.nframe() == 0 || Sys.getenv("PCR_DEMO") == "1") && !interactive()) {
set.seed(42)

scenA <- simulate_pcr(saturate = FALSE, stochastic = FALSE, scenario = "A: exponential baseline")
scenB <- simulate_pcr(saturate = TRUE,  stochastic = FALSE, scenario = "B: saturation")
scenC <- simulate_pcr(alpha = c(0.45, 0.40, 0.30, 0.20, 0.15),
                      saturate = FALSE, stochastic = FALSE, scenario = "C: per-barcode alpha")
scenD <- simulate_pcr(saturate = TRUE, stochastic = TRUE, reps = 20, reads = 1e5,
                      scenario = "D: stochastic + saturation + reads")

traj  <- bind_rows(scenA$trajectory, scenB$trajectory, scenC$trajectory, scenD$trajectory)
final <- bind_rows(scenA$final,      scenB$final,      scenC$final,      scenD$final)

# stochastic trajectories -> median over replicates (for amount/proportion plots)
traj_med <- traj %>%
  group_by(scenario, barcode, cycle) %>%
  summarise(amount = median(amount), proportion = median(proportion), .groups = "drop")

# observed final proportion (deterministic = amplicon prop; stochastic = median read prop)
obs <- final %>%
  group_by(scenario, barcode) %>%
  summarise(observed_prop = median(observed_prop), .groups = "drop") %>%
  left_join(final %>% distinct(scenario, barcode, input_proportion), by = c("scenario", "barcode"))

# --- PDF 1: amplicon amount vs cycle (log scale) -----------------------------
p1 <- ggplot(traj_med, aes(cycle, amount, color = barcode)) +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  facet_wrap(~scenario, scales = "free_y") +
  labs(title = "Amplicon amount per barcode over PCR cycles",
       x = "Cycle", y = "Amplicon molecules (log10)") +
  theme_bw()
ggsave(out("pcr_amount_vs_cycle.pdf"), p1, device = "pdf", width = 11, height = 7)

# --- PDF 2: proportion vs cycle ----------------------------------------------
p2 <- ggplot(traj_med, aes(cycle, proportion, color = barcode)) +
  geom_line(linewidth = 1) +
  facet_wrap(~scenario, scales = "free_y") +
  labs(title = "Barcode proportion over PCR cycles",
       x = "Cycle", y = "Proportion of amplicons") +
  theme_bw()
ggsave(out("pcr_proportion_vs_cycle.pdf"), p2, device = "pdf", width = 11, height = 7)

# --- PDF 3: observed vs input proportion (log-log) ---------------------------
p3 <- ggplot(obs, aes(input_proportion, observed_prop, color = barcode)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  geom_point(size = 2.5) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~scenario) +
  labs(title = "Observed (post-PCR) vs input proportion",
       x = "Input proportion (gDNA molecules)", y = "Observed proportion (post-PCR)") +
  theme_bw()
ggsave(out("pcr_observed_vs_input.pdf"), p3, device = "pdf", width = 11, height = 7)

# --- PDF 4: replicate boxplot (scenario D) -----------------------------------
p4 <- ggplot(scenD$final, aes(barcode, observed_prop)) +
  geom_boxplot() +
  geom_point(aes(barcode, input_proportion), color = "red", size = 2.5, shape = 4) +
  labs(title = "Scenario D: observed read-proportion across 20 replicates",
       subtitle = "Red X = input proportion; zero-count barcodes = dropout",
       x = "Barcode", y = "Observed proportion (reads)") +
  theme_bw()
ggsave(out("pcr_replicate_boxplot.pdf"), p4, device = "pdf", width = 8, height = 5)

# --- Scenario D grid: total_gDNA x cycles x eff (alpha = beta = eff) ---------
demo_ratios <- c(100, 10, 1, 0.1, 0.01)
input_theory <- demo_ratios / sum(demo_ratios)
grid_df <- expand.grid(total_gDNA = c(1e3, 1e4, 1e5, 1e6),
                       cycles = c(10, 20, 30),
                       eff = c(0.3, 0.5, 0.9))
grid_final <- lapply(seq_len(nrow(grid_df)), function(i) {
  g <- grid_df[i, ]
  sim <- simulate_pcr(total_gDNA = g$total_gDNA, alpha = g$eff, beta = g$eff,
                      cycles = g$cycles, saturate = TRUE, stochastic = TRUE,
                      reps = 20, reads = 1e5,
                      scenario = paste0("input=", g$total_gDNA, ", cycles=", g$cycles,
                                        ", eff=", g$eff))
  sim$final %>% mutate(total_gDNA = g$total_gDNA, cycles = g$cycles, eff = g$eff,
                       input_theory = rep(input_theory, length.out = n()))
}) %>% bind_rows() %>%
  mutate(combo = paste0("input=", total_gDNA, ", cycles=", cycles, ", eff=", eff),
         combo = factor(combo, levels = unique(combo)))

# PDF: observed vs input (log-log) per grid combo, dropout floored at 1e-6
cap_model <- paste0(
  "Integrated model:  A_i(t) = A_i(t-1) + Bin(G_i, \u03b1\u00b7(1 - A_i(t-1)/K)) + Bin(A_i(t-1), \u03b2\u00b7(1 - A_i(t-1)/K))", "\n",
  "Input & reads:  G_i ~ Pois(\u03bb_i);  reads ~ Multinom(1e5, A_i(n)/\u03a3A_j(n));  observed p_hat_i = reads_i/\u03a3reads_j", "\n",
  "Variables:", "\n",
  "  A_i(t) = barcode i amplicons after t cycles;  G_i = input gDNA molecules;  \u03bb_i = N\u00b7r_i/\u03a3r_j = expected input molecules", "\n",
  "  \u03b1 = gDNA conversion efficiency/cycle;  \u03b2 = amplicon amplification efficiency/cycle;  K = plateau capacity;  n = total cycles;  Bin = stochastic amplification", "\n",
  "Here alpha = beta = eff (0.3/0.5/0.9), K = 1e13, 20 replicates, 1e5 reads")
p6 <- ggplot(grid_final,
             aes(pmax(input_proportion, 1e-6), pmax(observed_prop, 1e-6),
                 color = barcode, shape = factor(eff))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  geom_point(size = 1.3, alpha = 0.6) +
  scale_x_log10() + scale_y_log10() +
  facet_grid(total_gDNA ~ cycles, labeller = label_both) +
  labs(title = "Scenario D grid: observed vs input proportion (log-log)",
       subtitle = "alpha = beta = eff (0.3/0.5/0.9); stochastic PCR + saturation (K=1e13),\n
20 replicates, 1e5 reads; zero input/reads floored at 1e-6; dashed line = identity",
       caption = cap_model,
       x = "Input proportion (gDNA molecules)", y = "Observed proportion (reads)",
       shape = "alpha = beta") +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(plot.caption = element_text(size = 6))
ggsave(out("pcr_scenarioD_grid_observed_vs_input.pdf"), p6, device = cairo_pdf,
       width = 8, height = 8)

unlink(out("pcr_scenarioD_grid_boxplot.pdf"))

grid_summary <- grid_final %>%
  group_by(combo, barcode) %>%
  summarise(median_obs = median(observed_prop),
            input      = first(input_theory),
            dropout_pct = mean(observed == 0) * 100, .groups = "drop") %>%
  filter(barcode == "B5 (0.01)")

# --- low-input dropout check (console only) ----------------------------------
scenD_low <- simulate_pcr(total_gDNA = 1e4, saturate = TRUE, stochastic = TRUE,
                          reps = 200, reads = 1e5, scenario = "low")
dropout <- scenD_low$final %>%
  group_by(barcode) %>%
  summarise(dropout_rate = mean(observed == 0) * 100)

# --- console summary ---------------------------------------------------------
val <- validate_pcr()
cat("=== Validation vs closed form (deterministic, no saturation) ===\n")
cat(sprintf("passed = %s, max relative error = %g\n\n", val$passed, val$max_rel_error))

cat("=== Final observed vs input proportions (median over reps) ===\n")
print(obs %>% arrange(scenario, barcode), n = Inf)
cat("\n")

cat("=== Scenario D: dropout rate of the 0.01 barcode at low input (total_gDNA = 1e4) ===\n")
print(dropout, n = Inf)
cat("\n")

cat("=== Scenario D grid (alpha = beta): 0.01 barcode median observed & dropout per combo ===\n")
print(grid_summary %>% arrange(combo, barcode), n = Inf)
cat("\n")

cat("=== Notes ===\n")
cat(" - In this model, per-barcode alpha yields a CONSTANT multiplicative bias\n")
cat("   (proportions shifted but flat over cycles); cycle-dependent proportion drift\n")
cat("   is driven by saturation (scenario B) or differential beta (pass beta as a\n")
cat("   per-barcode vector to see exponential drift).\n")
cat(" - PDFs written to: R/pcr_amount_vs_cycle.pdf, R/pcr_proportion_vs_cycle.pdf,\n")
cat("   R/pcr_observed_vs_input.pdf, R/pcr_replicate_boxplot.pdf,\n")
cat("   R/pcr_scenarioD_grid_observed_vs_input.pdf (alpha = beta grid; boxplot file removed)\n")
if (length(warnings())) { cat("\n--- Warnings ---\n"); print(warnings()) }
} # end demo guard
