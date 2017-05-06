# Main analysis script for Marchand et al., "A stochastic movement model reproduces
#  "patterns of site fidelity and long-distance dispersal in a population of 
#  Fowler’s Toads (Anaxyrus fowleri)"


library(tidyverse)
library(abc)
library(cowplot)
library(gridExtra)

#### Load data and functions ####

source("prep_day_location_data.R")
source("toad_move_funcs.R")


#### Plot inter-refuge distance distribution ####

dt_vals <- c(1, 2, 4, 8)

day_pairs <- inner_join(radio, radio, by = "Toad") %>%
    dplyr::select(Toad, Day.x, Day.y, x.x, x.y) %>%
    filter((Day.y - Day.x) %in% dt_vals) %>%
    mutate(dt = Day.y - Day.x, dx = abs(x.y - x.x))

pdf("fig1_density_dx.pdf", width = 6, height = 6)
ggplot(data = day_pairs, aes(x = dx, color = as.factor(dt))) +
    labs(x = expression(abs(Delta*x) ~~ (m)), y = "Density",
         color = "lag (days)") +
    geom_density() +
    geom_vline(xintercept = 10, linetype = "dotted") +
    scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000),
                  limits = c(0.01, 1000), 
                  labels = c("0.01", "0.1", "1", "10", "100", "1000")) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_brewer(labels = c("1 (604 pairs)", "2 (487 pairs)",
                                  "4 (311 pairs)", "8 (170 pairs)"),
                       palette = "Dark2") +
    theme_cowplot() +
    theme(legend.position = c(0.3, 0.8))
dev.off()


#### Run simulations and process results ####

# Get observed toad-day pairs (to subset simulation results)
toad_days <- dplyr::select(radio, Toad, Day)


# Generate nsim sets of parameters from uniform prior distributions
nsim <- 10000
pars <- data.frame(alph = runif(nsim, 1, 2), scal = runif(nsim, 10, 100),
                   p0 = runif(nsim, 0, 1), d0 = runif(nsim, 20, 2000))

# Run the three model versions with same parameter sets
nmod <- 3 # Number of models
npars <- c(3, 3, 4) # Number of parameters for each model
model_res <- list()
for (m in 1:nmod) {
    model_res[[m]] <- lapply(1:nrow(pars), function(i) {
        sim_toad(unlist(pars[i, ]), model_version = m, 
                 ntoad = max(toad_days$Toad), nday = max(toad_days$Day),
                 obs_toad_days = toad_days)
    })
}

# saveRDS(pars, "sim_pars.RData")
# saveRDS(model_res, "sim_res_all.RData")

###

# pars <- readRDS("sim_pars.RData")
# model_res <- readRDS("sim_res_all.RData")

# Get the summary statistics for all model results into one data.frame
res_stat <- do.call(rbind, lapply(model_res, function(res_list) {
    do.call(rbind, lapply(res_list, get_sum_stats))
}))
res_stat <- as.data.frame(res_stat)

# Vector of model indices
model_id <- rep(1:nmod, each = nsim)

# Get observed summary statistics
obs_stat <- get_sum_stats(radio)


#### Approximate Bayesian computation ####

# Model selection

# Cross-validation
cv_sel <- cv4postpr(index = model_id, sumstat = res_stat, method = "mnlogistic",
                    nval = 100, tols = c(0.005, 0.01, 0.05, 0.1))
# Estimation
mod_sel <- postpr(target = obs_stat, index = model_id, sumstat = res_stat,
                  tol = 0.05, method = "mnlogistic")


# Parameter estimation

# Parameter transformations for local-linear regression
transf <- c("logit", "none", "none", "log")

cv_res <- list()
abc_rej <- list()
for (m in 1:nmod) {
    # Cross-validation
    cv_res[[m]] <- cv4abc(param = pars[, 1:npars[m]],
                          sumstat = res_stat[model_id == m, ],
                          nval = 100, tols = c(0.005, 0.01, 0.05, 0.1),
                          method = "loclinear", transf = transf[1:npars[m]],
                          logit.bounds = cbind(1, 2))
    
    # Estimation
    abc_rej[[m]] <- abc(target = obs_stat, param = pars[, 1:npars[m]],
                        sumstat = res_stat[model_id == m, ], 
                        nval = 100, tol = 0.05, method = "loclinear", 
                        transf = transf[1:npars[m]], logit.bounds = cbind(1, 2))
}


#### Cross-validation plot with 95% credible intervals ####

# Pick model
m <- 1

# Sample 100 simulations and run ABC using the chosen simulation as data
stat_df <- res_stat[model_id == m, ]
cv_ind <- sample.int(nsim, 100)
cv_out <- lapply(cv_ind, function(i) {
    cbind(truth = unlist(pars[i, 1:npars[m]]),
          abc_summary(target = unlist(stat_df[i, ]), param = pars[-i, 1:npars[m]],
                      sumstat = stat_df[-i, ], tol = 0.05, method = "loclinear",
                      transf = transf[1:npars[m]], logit.bounds = cbind(1, 2)))
})
cv_out <- bind_rows(cv_out)
levels(cv_out$param) <- list("alpha" = "alph", "gamma" = "scal", "p[0]" = "p0", 
                             "d[0]" = "d0")

pdf("fig_s1_crossval.pdf", width = 8, height = 8)
ggplot(cv_out, aes(x = truth, y = med)) +
    labs(x = "True value", y = "Estimated value") +
    geom_point() +
    geom_errorbar(aes(ymin = lo, ymax = hi), alpha = 0.5) +
    geom_abline() +
    facet_wrap(~param, scales = "free", ncol = 2, labeller = label_parsed) +
    theme_cowplot() +
    theme(strip.background = element_blank())
dev.off()


#### Plots of estimate variance vs. number of sims ####

nrep = 100
nvals = seq(4000, 16000, 2000)

# Pick model
m <- 1

# For each value of N, sample 100 sets of N simulations with replacement (bootstrap)
#  and run ABC 
stat_df <- res_stat[model_id == m, ]
abc_boot <- lapply(nvals, function(ns) {
    replicate(nrep, {
        inds <- sample.int(nrow(pars), ns, replace = TRUE)
        abc_summary(target = obs_stat, param = pars[inds, 1:npars[m]],
                    sumstat = stat_df[inds, ], tol = 0.05, method = "loclinear", 
                    transf = transf[1:npars[m]], logit.bounds = cbind(1, 2))     
    }, simplify = FALSE)
})
names(abc_boot) <- nvals

# Get the distibution of estimates for the 2.5%, 50% and 97.5% quantiles 
#  per parameter and per value of N
abc_boot <- bind_rows(lapply(abc_boot, bind_rows), .id = "nsim") %>%
    mutate(nsim = as.integer(nsim)) %>%
    gather(key = "quantile", value = value, lo:hi)

abc_boot <- group_by(abc_boot, param, nsim, quantile) %>%
    summarize(p025 = quantile(value, 0.025), median = quantile(value, 0.5),
              p975 = quantile(value, 0.975))
levels(abc_boot$param) <- list("alpha" = "alph", "gamma" = "scal", "p[0]" = "p0", 
                               "d[0]" = "d0")

pdf("fig_s2_estimates_nsim.pdf", width = 8, height = 8)
ggplot(abc_boot, aes(x = as.factor(nsim), y = median, color = quantile)) +
    labs(x = "Number of simulations", y = "Estimate", color = "Quantile") +
    geom_point() +
    geom_errorbar(aes(ymin = p025, ymax = p975), width = 0.1) +
    scale_color_brewer(breaks = c("hi", "med", "lo"), palette = "Dark2",
                       labels = c("97.5%", "50%", "2.5%")) +
    facet_wrap(~param, scales = "free_y", labeller = label_parsed, dir = "v",
               strip.position = "top") +
    theme_cowplot() +
    theme(strip.background = element_blank())
dev.off()


#### Posterior predictive checks ####

# Run simulations using the ABC posterior parameter distribution for each model
#  (500 parameter values with weights)
pars_post <- lapply(abc_rej, function(x) x$adj.values)
post_res <- list()
for (m in 1:nmod) {
    post_res[[m]] <- lapply(1:nrow(pars_post[[m]]), function(i) {
        sim_toad(unlist(pars_post[[m]][i, ]), model_version = m, 
                 ntoad = max(toad_days$Toad), nday = max(toad_days$Day),
                 obs_toad_days = toad_days)
    })    
}
post_weights <- lapply(abc_rej, function(x) x$weights / sum(x$weights))

# Get summary statistics for each simulation
post_stat <- do.call(rbind, lapply(post_res, function(res_list) {
    do.call(rbind, lapply(res_list, get_sum_stats))
}))
post_stat <- as.data.frame(post_stat)
# Add model ID and weight columns
post_stat$model_id <- rep(1:nmod, each = 0.05 * nsim)
post_stat$wgt <- do.call(c, post_weights)
# Convert stats data in "tall format" (statistic, time lag, value)
post_stat <- gather(post_stat, key = "var_dt", value = "value",
                    mld21:p_ret8) %>%
    extract(var_dt, c("var", "dt"), "(.*)([[:digit:]])") %>% 
    mutate(var = factor(var, c("mld2", "sdld2", "p_ret"),
                        labels = c("mean ~~ log ~ (Delta*x)^2", 
                                   "s.d. ~~ log ~~ (Delta*x)^2", 
                                   "p[ret]")), 
           dt = factor(dt, c("8", "4", "2", "1")))

# Convert vector of observed statistics into same format
obs_stat <- data.frame(var_dt = names(obs_stat), value = obs_stat)
obs_stat <- extract(obs_stat, var_dt, c("var", "dt"), "(.*)([[:digit:]])") %>%
    mutate(var = factor(var, c("mld2", "sdld2", "p_ret"), 
                        labels = levels(post_stat$var)), 
           dt = factor(dt, c("8", "4", "2", "1")))

# Facet plot for the posterior distribution of each summmary statistic
#  compared with the observed value in the dataset
pdf("fig2_post_sim.pdf", width = 9, height = 6)
ggplot() + 
    labs(x = "Statistic", y = "Time lag (days)", linetype = "Model") +
    geom_density(data = post_stat, aes(x = value, ..scaled.., weight = wgt, 
                                       linetype = as.factor(model_id))) +
    geom_vline(data = obs_stat, aes(xintercept = value), color = "red") +
    facet_grid(dt ~ var, scales = "free", switch = "both",
               labeller = label_parsed) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    theme_cowplot(font_size = 16, line_size = 1) +
    theme(strip.background = element_blank(), strip.placement = "outside", 
          strip.text.y = element_text(angle = 180), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), axis.text.y = element_blank())
dev.off()


# For each model, count the average number of distinct refuges per Toad
ref_sim <- list()
for (m in 1:nmod) {
    ref_dist <- vapply(post_res[[m]], function(df) {
        df2 <- group_by(df, Toad) %>%
            summarize(nref = count_refuges(x))
        df2$nref
    }, rep(0, max(radio$Toad)))
    ref_sim[[m]] <- data.frame(
        Toad = 1:max(radio$Toad),
        nref_mean = apply(ref_dist, 1, 
                          function(x) weighted.mean(x, post_weights[[m]]))
    )
}
names(ref_sim) <- paste('Model', 1:nmod)
ref_sim <- bind_rows(ref_sim, .id = "source")

# Count the observed number of refuges per toad, and add to dataframe
#  along with the number of observations per toad
ref_obs <- group_by(radio, Toad) %>%
    summarize(npts = n(), source = "Observed", nref_mean = count_refuges(x))
ref_obs$source <- "Observed"

ref_sim <- inner_join(dplyr::select(ref_obs, Toad, npts), ref_sim)
ref_sim <- bind_rows(ref_sim, ref_obs)

# Plot the observed and (mean) simulated number of refuges per Toad against 
#  the number of observations per toad, and add linear regression lines in log space
pdf("fig3_nrefuges.pdf", width = 6, height = 6)
ggplot(data = ref_sim, aes(x = npts, y = nref_mean, color = source, fill = source)) +
    labs(x = "Number of observations", y = "Number of refuges",
         color = "") +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, alpha = 0.2, show.legend = FALSE) +
    scale_x_log10(breaks = c(1, 2, 4, 8, 16, 32)) +
    scale_y_log10(breaks = c(1, 3, 6, 12)) +
    scale_fill_brewer(guide = "none", palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    theme_cowplot() +
    theme(legend.position = c(0.3, 0.9))
dev.off()

summary(lm(log(nref_mean) ~ log(npts), data = ref_sim))    

