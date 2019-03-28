## ----echo = FALSE, message = FALSE, warning = FALSE----------------------
psm <- c("N-state partitioned survival model (PSM)", "[Psm](../reference/Psm.html)")
ictstm <- c("Individual-level continuous time state transition model (iCTSTM)", "[IndivCtstm](../reference/IndivCtstm.html)")
tbl <- rbind(psm, ictstm)
colnames(tbl) <- c("Economic model", "R6 class")
knitr::kable(tbl, row.names = FALSE)

## ---- out.width = "600px", echo = FALSE----------------------------------
knitr::include_graphics("econ-eval-process-hesim.png")

## ----warning = FALSE, message = FALSE------------------------------------
library("hesim")
library("data.table")
strategies <- data.table(strategy_id = c(1, 2))
n_patients <- 1000
patients <- data.table(patient_id = 1:n_patients,
                          age = rnorm(n_patients, mean = 45, sd = 7),
                          female = rbinom(n_patients, size = 1, prob = .51))
states <- data.table(state_id = c(1, 2),
                     state_name = c("Healthy", "Sick")) # Non-death health states
tmat <- rbind(c(NA, 1, 2),
              c(3, NA, 4),
              c(NA, NA, NA))
colnames(tmat) <- rownames(tmat) <- c("Healthy", "Sick", "Dead")
transitions <- create_trans_dt(tmat)
transitions[, trans := factor(transition_id)]
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients, 
                        states = states,
                        transitions = transitions)
print(hesim_dat)

## ----echo = FALSE, message = FALSE, warning = FALSE----------------------
psm <- c("[Psm](../reference/Psm.html)", "Independent survival models", "[params_surv_list](../reference/params_surv_list.html)", "[hesim::flexsurvreg_list](../reference/flexsurvreg_list.html)")
ictstm <- c("[IndivCtstm](../reference/IndivCtstm.html)", "Multi-state model", "[params_surv](../reference/params_surv.html) or [params_surv_list](../reference/params_surv_list.html)", "[flexsurv::flexsurvreg](https://www.rdocumentation.org/packages/flexsurv/versions/1.0.0/topics/flexsurvreg) or [hesim::flexsurvreg_list](../reference/flexsurvreg_list.html)")
tbl <- rbind(psm, ictstm)
colnames(tbl) <- c("Economic model (R6 class)", "Statistical model", "Parameter object", "Model fit object")
knitr::kable(tbl, row.names = FALSE)

## ---- message = FALSE, warning = FALSE-----------------------------------
library("flexsurv")
mstate_data <- data.table(ctstm3_exdata$transitions)
mstate_data[, trans := factor(trans)]
fit_wei <- flexsurv::flexsurvreg(Surv(years, status) ~ trans + 
                                                       factor(strategy_id):trans +
                                                       age:trans + 
                                                       female: trans +
                                                       shape(trans), 
                                 data = mstate_data, 
                                 dist = "weibull")

## ----echo = FALSE, message = FALSE, warning = FALSE----------------------
means <- c("Mean model", "[params_mean](../reference/params_mean.html)", "[stateval_tbl](../reference/stateval_tbl.html)")
lm <- c("Linear model", "[params_lm](../reference/params_lm.html)", "[stats::lm](https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/lm)")
tbl <- rbind(means, lm)
colnames(tbl) <- c("Statistical model", "Parameter object", "Model fit object")
knitr::kable(tbl, row.names = FALSE)

## ------------------------------------------------------------------------
# Utility
utility_tbl <- stateval_tbl(data.table(state_id = states$state_id,
                                       mean = ctstm3_exdata$utility$mean,
                                       se = ctstm3_exdata$utility$se),
                            dist = "beta",
                            hesim_data = hesim_dat)

# Costs
drugcost_tbl <- stateval_tbl(data.table(strategy_id = strategies$strategy_id,
                                       est = ctstm3_exdata$costs$drugs$costs),
                            dist = "fixed",
                            hesim_data = hesim_dat) 
medcost_tbl <- stateval_tbl(data.table(state_id = states$state_id,
                                       mean = ctstm3_exdata$costs$medical$mean,
                                       se = ctstm3_exdata$costs$medical$se),
                            dist = "gamma",
                            hesim_data = hesim_dat)  

## ----echo = FALSE, message = FALSE, warning = FALSE----------------------
psm <- c("[Psm](../reference/Psm.html)", "[PsmCurves](../reference/PsmCurves.html)",
         "[StateVals](../reference/StateVals.html)", "[StateVals](../reference/StateVals.html)")
ictstm <- c("[IndivCtstm](../reference/IndivCtstm.html)", "[IndivCtstmTrans](../reference/IndivCtstmTrans.html)",
         "[StateVals](../reference/StateVals.html)", "[StateVals](../reference/StateVals.html)")
tbl <- rbind(psm, ictstm)
colnames(tbl) <- c("Economic model", "Disease model", "Utility model", "Cost model(s)")
knitr::kable(tbl, row.names = FALSE)

## ------------------------------------------------------------------------
n_samples <- 1000

## ----warning = FALSE, message = FALSE------------------------------------
transmod_data <- expand(hesim_dat, 
                        by = c("strategies", "patients", "transitions"))
head(transmod_data)
attr(transmod_data, "id_vars")

## ------------------------------------------------------------------------
transmod <- create_IndivCtstmTrans(fit_wei, transmod_data,
                                   trans_mat = tmat, n = n_samples)
class(transmod)

## ------------------------------------------------------------------------
# Utility
utilitymod <- create_StateVals(utility_tbl, n = n_samples)

# Costs
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples)
costmods <- list(drugs = drugcostmod,
                 medical = medcostmod)

## ------------------------------------------------------------------------
ictstm <- IndivCtstm$new(trans_model = transmod,
                         utility_model = utilitymod,
                         cost_models = costmods)

## ----echo = FALSE, message = FALSE, warning = FALSE----------------------
psm_methods <- c("[Psm](../reference/Psm.html)", "$sim_survival() and $sim_stateprobs()", "$sim_qalys()", "$sim_costs()")
ictstm_methods <- c("[IndivCtstm](../reference/IndivCtstm.html)", "$sim_disease() and $sim_stateprobs()", "$sim_qalys()", "$sim_costs()")
tbl <- rbind(psm_methods, ictstm_methods)
colnames(tbl) <- c("Economic model (R6 class)", "Disease progression", "QALYs", "Costs")
knitr::kable(tbl, row.names = FALSE)

## ------------------------------------------------------------------------
ictstm$sim_disease()
head(ictstm$disprog_)

## ------------------------------------------------------------------------
ictstm$sim_stateprobs(t = c(0:10))
head(ictstm$stateprobs_)

## ------------------------------------------------------------------------
# QALYs
ictstm$sim_qalys(dr = .03)
head(ictstm$qalys_)

# Costs
ictstm$sim_costs(dr = .03)
head(ictstm$costs_)

## ------------------------------------------------------------------------
ce <- ictstm$summarize()
print(ce)

## ------------------------------------------------------------------------
icea <- icea(ce, dr_qalys = .03, dr_costs = .03)
icea_pw <- icea_pw(ce, dr_qalys = .03, dr_costs = .03, comparator = 1)

## ----ceac_plot, warning = FALSE, message = FALSE-------------------------
library("ggplot2")
ggplot2::ggplot(icea_pw$ceac, aes(x = k, y = prob, col = factor(strategy_id))) +
  geom_line() + xlab("Willingness to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, 200000, 100000), label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Strategy") + 
  theme_minimal()

