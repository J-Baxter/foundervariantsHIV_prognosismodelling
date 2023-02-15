


nlform <- bf(CD4_decline ~ alpha * SPVL_log10^2 + sex,
             SPVL_log10 ~ Additive(variants) + E,
             sex ~ 1, 
             alpha ~ 1,
             E ~ 1,
             variants ~ rnorm(theta, T, x),
             x ~ 1,
             nl = TRUE)



nlprior <- c(prior(normal(,), nlpar = "alpha"),
             prior(normal(,), nlpar = "omega"),
             prior(normal(,), nlpar = "theta"))

fit_loss1 <- brm(formula = nlform, data = loss, family = gaussian(),
                 prior = nlprior, control = list(adapt_delta = 0.9))


CD4_decline = alpha * SPVL_log10^2 + age + sex
SPVL_log10 = func(variants)
variants = & normal distirbution with mean (transmitter_vl)

transmitter_vl (fixed)