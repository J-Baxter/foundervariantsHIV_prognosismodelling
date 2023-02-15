nlform <- bf(cum ~ theta + sex + 1|age,
             theta ~ , omega ~ 1, theta ~ 1, nl = TRUE)

nlprior <- c(prior(normal(5000, 1000), nlpar = "ult"),
             prior(normal(1, 2), nlpar = "omega"),
             prior(normal(45, 10), nlpar = "theta"))

fit_loss1 <- brm(formula = nlform, data = loss, family = gaussian(),
                 prior = nlprior, control = list(adapt_delta = 0.9))


CD4_decline = alpha * SPVL_log10^2 + age + sex
SPVL_log10 = func(variants)
variants = sample according to empirical distirbution & normal distirbution with mean (T)