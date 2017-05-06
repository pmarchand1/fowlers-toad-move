# Functions to simulate toad movement and analyze model output #
####


# Main function simulating daytime refuges of multiple toads over multiple days
#  along one dimension. 
#
#  Foraging at night is modelled as a random walk with final position 
#  following a stable distribution with stability parameter alph and
#  scale parameter scal. When done foraging, the toad either takes refuge
#  at its current location, or returns to a previous refuge. 
# 
# There are three model versions, differing in how returns are determined:
#
#  1. A toad has a constant probability p0 of returning to a previous refuge.
#    In that case, the refuge is picked at random from all previous days' refuges.
#  2. Still a constant probability of return p0, but the toad always returns
#    to the nearest previous refuge from its current position.
#  3. The probability of returning to each *distinct* previous refuge decays
#   with distance with a characteristic distance d0, i.e. pret = p0 * exp(-d/d0).

sim_toad <- function(params, model_version = 1, ntoad = 50, nday = 30,
                     obs_toad_days = NA) {
    # params is the vectors of parameters to be estimated via ABC,
    # model_version is the index for model selection,
    # ntoad and nday (number of toads and days to simulate) are fixed
    # obs_toad_days (optional) subsets the output to the specified Toad-Day pairs
    alph <- params[1]
    scal <- params[2]
    p0 <- params[3]
    d0 <- params[4]
    
    # xs[nday, ntoad] are the daytime locations
    # xn is the position of the toad after nighttime foraging
    xs <- matrix(0, nrow = nday, ncol = ntoad)
    for (j in 1:ntoad) {
        if (model_version == 3) {
            # We need to keep track of distinct refuges for this version
            # starting with one refuge at x = 0
            nref <- 1
            xref <- rep(0, nday)    
        }
        for (i in 1:(nday - 1)) {
            xn <- xs[i, j] + rstable(1, scale = scal, alpha = alph)
            if (model_version == 3) {
                # Distance vector to all previous refuges
                d <- abs(xn - xref[1:nref])
                # Distance-based probability factor
                pdist <- p0 * exp(-d/d0)
                # Probability of not returning is product of probability of not returning
                #  to each of the previous distinct refuges
                p_noret <- prod(1 - pdist) 
            } else {# Constant probability of return
                p_noret <- 1 - p0
            }

            rand <- runif(1)
            if (runif(1) < p_noret) {
                # If no return, take refuge here
                xs[i + 1, j] <- xn
                if (model_version == 3) {
                    nref <- nref + 1
                    xref[nref] <- xn
                }
            } else {# Return case
                if (model_version == 1) {
                    # Pick a previous refuge at random
                    ret_pt <- sample.int(i, 1)
                    xs[i + 1, j] <- xs[ret_pt, j]
                } else if (model_version == 2) {
                    # Return to nearest refuge
                    ret_pt <- which.min(abs(xn - xs[1:i, j]))
                    xs[i + 1, j] <- xs[ret_pt, j]
                } else {# model_version ==  3
                    # Pick a refuge with prob. proportional to pdist
                    p_cumul <- cumsum(pdist / sum(pdist))
                    rand <- runif(1)
                    ret_pt <- which(p_cumul > rand)[1]
                    xs[i + 1, j] <- xref[ret_pt]
                }
            }
        }
    }
    # Convert xs matrix to result data frame (Toad, Day, x)
    colnames(xs) <- 1:ntoad
    xs <- as.data.frame(cbind(Day = 1:nday, xs)) %>%
        gather(key = "Toad", value = x, -Day) %>%
        mutate(Toad = as.integer(Toad))
    if (is.data.frame(obs_toad_days)) {
        xs <- semi_join(xs, obs_toad_days, by = c("Toad", "Day"))
    }
    xs
}


# Draw n values from a symmetric, zero-centered stable distribution with
#  given scale and stability (alpha) parameters, using the CMS algorithm
rstable <- function(n, scale, alpha) {
    if (alpha > 2 | alpha <= 0) {
        stop("rstable is not defined for alpha outside the interval 0 < alpha <= 2\n")
    }
    if (alpha == 1) {
        return(scale * rcauchy(n, location = 0, scale = 1))
    }
    if (alpha == 2) {
        return(rnorm(n, mean = 0, sd = sqrt(2) * scale))
    }
    u <- pi * (runif(n, min = 0, max = 1) - 0.5)
    v <- rexp(n, rate = 1)
    t <- sin(alpha * u)/(cos(u)^(1/alpha))
    s <- (cos((1 - alpha) * u)/v)^((1 - alpha)/alpha)
    return(scale * t * s)
}


# Summary statistics from result data frame (output of sim_mrw)
#  optional parameters: 
#  - dt_vals: which day intervals to compute statistics for
#  - d2min: minimum distance-squared for non-returns (d2 < d2min is a return event)      
get_sum_stats <- function(result_df, dt_vals = c(1, 2, 4, 8), d2min = 100) {
    # Calculate statistics for all location pairs separated by dt_vals days
    day_pairs <- inner_join(result_df, result_df, by = "Toad") %>%
        dplyr::select(Toad, Day.x, Day.y, x.x, x.y) %>%
        filter((Day.y - Day.x) %in% dt_vals) %>%
        mutate(dt = Day.y - Day.x, d2x = (x.y - x.x)^2, ret = d2x < d2min)
    # Calculate the mean and s.d. of log distance-squared for non-returns,
    #  and the probability of return
    res_stats <- group_by(day_pairs, dt) %>%
                 summarise(mld2 = mean(log(d2x[!ret])), 
                           sdld2 = sd(log(d2x[!ret])), p_ret = sum(ret) / n())
    snames <- as.vector(t(outer(colnames(res_stats[, -1]), res_stats[[1]], paste0)))
    res_stats <- as.vector(as.matrix(res_stats[, -1]))
    names(res_stats) <- snames
    res_stats
}


# Calculate number of refuges with distance threshold "d"
#  from a vector of locations "x" using hierarchical clustering with 
#  "complete" method (so that all points in cluster are within "d" of each other)
count_refuges <- function(x, d = 10) {
    cl <- hclust(dist(x), method = "complete")
    max(cutree(cl, h = d))
}


# Wrapper around 'summary.abc' that outputs the medians and 95% credible intervals
abc_summary <- function(...) {
    stats_keep <- c("Weighted 2.5 % Perc.:", "Weighted Median:",
                    "Weighted 97.5 % Perc.:")
    abc_out <- summary(abc(...), print = FALSE)
    abc_out <- abc_out[stats_keep, ]
    rownames(abc_out) <- c("lo", "med", "hi")
    abc_out <- data.frame(abc_out)
    colnames(abc_out) <- c("stat", "param", "value")
    spread(abc_out, key = stat, value = value)
}


