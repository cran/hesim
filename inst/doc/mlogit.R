## ---- out.width = "500px", echo = FALSE---------------------------------------
knitr::include_graphics("illness-death.png")

## -----------------------------------------------------------------------------
tmat <- rbind(c(0, 1, 2),
              c(NA, 0, 1),
              c(NA, NA, NA))
colnames(tmat) <- rownames(tmat) <- c("Healthy", "Sick", "Dead")
print(tmat)

## ---- warning = FALSE, message = FALSE----------------------------------------
library("hesim")
library("data.table")
transitions_data <- data.table(multinom3_exdata$transitions)
head(transitions_data)

## -----------------------------------------------------------------------------
n_patients <- 100
patients <- transitions_data[year == 1, .(patient_id, age, female)][
  sample.int(nrow(transitions_data[year == 1]), n_patients)][
  , grp_id := 1:n_patients]
head(patients)

## -----------------------------------------------------------------------------
hesim_dat <- hesim_data(
  patients = patients,
  strategies = data.table(strategy_id = 1:2,
                          strategy_name = c("Reference", "Intervention")),
  states = data.table(state_id = c(1, 2),
                     state_name = c("Healthy", "Sick")) # Non-death health states
)

## -----------------------------------------------------------------------------
library("nnet")
library("splines")

# Transitions from healthy state
data_healthy <- transitions_data[state_from == "Healthy"]
fit_healthy <- multinom(state_to ~ strategy_name + female + 
                          ns(age, df = 5) + year_cat, 
                        data = data_healthy, trace = FALSE)

# Transitions from sick state
data_sick <- droplevels(transitions_data[state_from == "Sick"])
fit_sick <- multinom(state_to ~ strategy_name + female + 
                       ns(age, df = 5) + year_cat, 
                     data = data_sick, trace = FALSE)

## -----------------------------------------------------------------------------
transfits <- multinom_list(healthy = fit_healthy, sick = fit_sick)

## -----------------------------------------------------------------------------
utility_tbl <- stateval_tbl(multinom3_exdata$utility,
                            dist = "beta",
                            hesim_data = hesim_dat)
head(utility_tbl)

## -----------------------------------------------------------------------------
drugcost_tbl <- stateval_tbl(multinom3_exdata$costs$drugs,
                            dist = "fixed",
                            hesim_data = hesim_dat) 
medcost_tbl <- stateval_tbl(multinom3_exdata$costs$medical,
                            dist = "gamma",
                            hesim_data = hesim_dat)  

## -----------------------------------------------------------------------------
n_samples <- 100

## -----------------------------------------------------------------------------
tintervals <- time_intervals(unique(transitions_data[, .(year_cat)])
                             [, time_start := c(0, 2, 6)])
transmod_data <- expand(hesim_dat, times = tintervals)
transmod <- create_CohortDtstmTrans(transfits,
                                    input_data = transmod_data,
                                    trans_mat = tmat,
                                    n = n_samples,
                                    point_estimate = FALSE)

## -----------------------------------------------------------------------------
# Utility
utilitymod <- create_StateVals(utility_tbl, n = n_samples)

# Costs
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples)
costmods <- list(Drug = drugcostmod,
                 Medical = medcostmod)

## -----------------------------------------------------------------------------
econmod <- CohortDtstm$new(trans_model = transmod,
                           utility_model = utilitymod,
                           cost_models = costmods)

## ----stateprobs, fig.width = 7, fig.height = 4--------------------------------
econmod$sim_stateprobs(n_cycles = 20)

# Plot
library("ggplot2")
theme_set(theme_bw())
stateprob_summary <- econmod$stateprobs_[, .(prob_mean = mean(prob),
                                              prob_lower = quantile(prob, .025),
                                              prob_upper = quantile(prob, .975)),
                                          by = c("strategy_id", "state_id", "t")]
stateprob_summary[, strategy_name := factor(strategy_id,
                                          labels = hesim_dat$strategies$strategy_name)]
ggplot(stateprob_summary, aes(x = t, y = prob_mean)) +
  geom_line(aes(col = strategy_name)) +
  geom_ribbon(aes(x = t, ymin = prob_lower, ymax = prob_upper,
                  fill = strategy_name), alpha = .3) +
  facet_wrap(~factor(state_id, labels = c(hesim_dat$states$state_name, "Death"))) +
  xlab("Year") + ylab("Probability") +
  scale_fill_discrete("Strategy") + scale_color_discrete("Strategy")

## -----------------------------------------------------------------------------
econmod$sim_qalys(dr = .03)
econmod$sim_costs(dr = .03)

## -----------------------------------------------------------------------------
ce_sim <- econmod$summarize(by_grp = TRUE)
cea_pw_out <- cea_pw(ce_sim, comparator = 1, dr_qalys = .03, dr_costs = .03)

## ----icerHeterogeneity, fig.width = 7, fig.height = 4-------------------------
icers <- merge(cea_pw_out$summary,
               hesim_dat$patients[, .(grp_id, age, female)],
               by = "grp_id")
icers[, gender := factor(female, 
                         levels = c(0, 1),
                         labels = c("Male", "Female"))]

# Plot of ICER by demographics
ggplot(icers, aes(x = age, y = as.numeric(gsub(",", "", icer)), col = gender)) +
  geom_point() +
  xlab("Age") + ylab("Incremental cost-effectiveness ratio") +
  scale_y_continuous(label = scales::dollar_format()) +
  scale_colour_discrete(name = "")

