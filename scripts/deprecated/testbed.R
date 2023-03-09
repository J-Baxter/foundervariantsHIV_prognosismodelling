#H2 model
library(caret)
NPAIRS = 100

rnormt <- function(n, range, mu, s = 1) {
  
  # range is a vector of two values
  
  F.a <- pnorm(min(range), mean = mu, sd = s)
  F.b <- pnorm(max(range), mean = mu, sd = s)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qnorm(u, mean = mu, sd = s)
  
}


test <- cbind.data.frame(gud_ = sample(c('pos','neg'), NPAIRS,  prob= c(0.3,0.7), rep = T),
                         pos_ = sample(c(rep('index', NPAIRS/2), rep('recipient', NPAIRS/2)), NPAIRS, rep = F),
                         age = sample(seq(15, 70, by= 1), NPAIRS, rep = T))  %>% 
  group_by(pos_) %>% 
  mutate(couple = row_number()) %>%
  ungroup() %>%
  group_by(couple) %>%
  group_modify(function(.x, .y) .x %>% mutate(log10spvl_mean = rnormt(1, c(2.6, 7.14), mu = 4.39, s = 0.84))) %>%
  mutate(sex = sample(c('M', 'F'),2, rep = F)) %>%
  mutate(sub_ = sample(c('A', 'C', 'D', 'R'), 1, rep = T)) %>%
  ungroup() %>%
  mutate(age_group = cut(age , breaks = c(-Inf,24,29,39, Inf), labels = c('A_1', 'A_2', 'A_3', 'A_4'))) 


test_dummy <-  caret::dummyVars("~ sex + gud_ + sub_ + pos_ + log10spvl_mean + age_group + couple", data= test) %>% predict(.,newdata = test) %>% 
  data.frame() %>%
  rename(A_1 = age_group.A_1, A_2 = age_group.A_2,A_3 = age_group.A_3,A_4 = age_group.A_4) 

results <- test_dummy %>%
   mutate(log10spvl_est = log10spvl_mean - 0.15*sexM - 0.15*A_1 - 0.32*A_2 + 0.35*A_3 + 0.73*A_4 + 
            0.42*gud_pos + 0.62*gud_neg +
            0.068*sub_A + 1.85*sub_C + 0.39*sub_D + 0.78*sub_R + 
            0.42*pos_index + 0*pos_recipient) %>%
  dplyr::select(log10spvl_est)%>%cbind.data.frame(test) %>%
  pivot_wider(values_from = !couple, names_from = pos_, names_glue =  "{pos_}_{.value}") %>%
  dplyr::select(!ends_with('pos_')  & !ends_with('group')) #names with prefix pos_
  

ggplot(results) + 
  geom_point(aes(x = index_log10spvl_est, y = recipient_log10spvl_est))


           
