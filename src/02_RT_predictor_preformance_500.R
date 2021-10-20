library(dplyr)
library(ggplot2)
library(gridExtra)
library(readxl)
library(stringr)
library(tidyr)
library(vroom)

regression_stats <- function(obs, pred){
  # lm
  data = data.frame(pred = pred, 
                    obs = obs)
  pred.lm = lm(pred ~ obs, data = data)
  
  # metrics
  # correlation coefficients
  pc = cor(obs, pred, method = "pearson")
  sm = cor(obs, pred, method = "spearman")
  
  # mean squared error
  mse = (obs - pred)^2 %>% mean() %>% round(4)
  # root mean squared error
  rmse = sqrt(mse) %>% round(4)
  # mean absolute deviation
  mae = (obs - pred) %>% abs() %>% mean() %>% round(4)
  
  # sumarize
  all.metrics = c(summary(pred.lm)$r.squared, pc, mse, rmse, mae)
  names(all.metrics) = c("Rsquared", "PCC", "MSE", "RMSE", "MAE")
  
  return(all.metrics)
}

# Ground truth
files_test <- list.files("results_train_500/test/achrom/", pattern = ".csv", full.names = T) %>%
  lapply(vroom, show_col_types = FALSE)
names(files_test) <- list.files("results_train_500/test/achrom/", pattern = ".csv", full.names = F) %>%
  str_remove(".csv")

# Achrom
pred_achrom <- list.files("results_train_500/predict/achrom/", pattern = ".csv", full.names = T) %>%
  lapply(vroom, show_col_types = FALSE)
names(pred_achrom) <- list.files("results_train_500/predict/achrom//", pattern = ".csv", full.names = F, recursive = T) %>%
  str_remove_all(pattern = ".csv")
pred_achrom <- pred_achrom %>%
  bind_rows(.id = "file") %>%
  select(-"...1")  %>%
  rename(RT_pred = RT) %>%
  split(~file)
RT_predictors <- c("achrom")

# Deep LC
pred_DeepLC <- list.files("results_train_500/predict/DeepLC/", pattern = ".csv", full.names = T) %>%
  lapply(vroom, show_col_types = FALSE)
names(pred_DeepLC) <- list.files("results_train_500/predict/DeepLC/", pattern = ".csv", full.names = F, recursive = T) %>%
  str_remove_all(pattern = ".csv")
pred_DeepLC <- pred_DeepLC %>%
  bind_rows(.id = "file") %>%
  select(-"...1") %>%
  select(-modifications) %>%
  rename(Mascot_seq = seq,
         RT_pred = predicted_tr)  %>%
  split(~file)
RT_predictors <- c(RT_predictors, "DeepLC")

# AutoRT
pred_AutoRT <- list.files("results_train_500/predict/AutoRT/", pattern = "test.tsv", full.names = T, recursive = T) %>%
  lapply(vroom, show_col_types = FALSE)
names(pred_AutoRT) <- list.files("results_train_500/predict/AutoRT/", pattern = "test.tsv", full.names = F, recursive = T) %>%
  str_remove_all(pattern = "/test.tsv")

pred <- lapply(pred_AutoRT, function(x) predict(lm(y ~ y_pred, data = x))) %>%
  lapply(as_tibble) %>%
  bind_rows()

pred_AutoRT <- pred_AutoRT %>%
  bind_rows(.id = "file") %>%
  rename(Mascot_seq = x)  %>%
  mutate(RT_pred = pred$value) %>%
  select(-y) %>%
  select(-y_pred) %>%
  split(~file)

RT_predictors <- c(RT_predictors, "AutoRT")


# Prosit
pred_Prosit <- list.files("results_train_500/predict/Prosit/", pattern = "prositPred", full.names = T, recursive = T) %>%
  lapply(vroom, show_col_types = FALSE)
names(pred_Prosit) <- list.files("results_train_500/predict/Prosit/", pattern = "prositPred", full.names = F, recursive = T) %>%
  str_remove_all(pattern = "prositPred") %>%
  str_remove_all(pattern = ".csv")

pred <- lapply(pred_Prosit, function(x) predict(lm(RT ~ prositIRT, data = x))) %>%
  lapply(as_tibble) %>%
  bind_rows()

pred_Prosit <- pred_Prosit %>%
  bind_rows(.id = "file") %>%
  rename(RT_pred = prositIRT)  %>%
  select(-RT) %>%
  mutate(RT_pred = pred$value) %>%
  split(~file)
RT_predictors <- c(RT_predictors, "Prosit")

# Evaluate each sampling iteration separately
split_structure <- files_test %>%
  names() %>%
  as_tibble() %>%
  separate(value, sep = "[.]", into = c("date", "SPL", "sample")) %>%
  mutate(sample = as.integer(str_remove(sample, "sample"))) %>%
  select(date, SPL) %>%
  unique()

RT_predictors <- paste0("pred_", RT_predictors)

# ---------------------------------------- (1) Pre-processing ------------------------------------------------
out <- tibble()
for (i in seq_along(split_structure$SPL)) {
  
  SPL <- split_structure$SPL[i]
  keep <- files_test[grep(SPL, names(files_test))]
  
  for (j in seq_along(keep)) {
    
    for (pred in RT_predictors) {
      
      dt = keep[[j]] %>%
        left_join(get(pred)[[names(keep)[j]]]) %>%
        na.omit()
      
      out_tmp <- regression_stats(obs = dt$RT, 
                                  pred = dt$RT_pred) %>%
        t() %>%
        as_tibble() %>%
        mutate(SPL = SPL,
               sample = names(keep)[j],
               predictor = pred)
      out <- rbind(out, out_tmp)
    }
  }
}
out <- out %>%
  mutate(predictor = str_remove(predictor, "pred_"),
         sample = str_split_fixed(sample, ".sample", 2)[,2]) %>%
  pivot_longer(cols = c("Rsquared", "PCC", "MSE", "RMSE", "MAE"), names_to = "metric")


metrics = c("Rsquared", "PCC", "MSE", "RMSE", "MAE")
gg <- list()
for (i in seq_along(metrics)) {
  
  keep_col = metrics[[i]]
  gg[[i]] <- out %>%
    filter(metric == keep_col) %>%
    ggplot(aes(y=value, x=as.factor(SPL))) + 
    geom_boxplot(aes(fill=as.factor(predictor)), stat="boxplot", position="dodge", alpha=0.5, width=0.2) + 
    theme_bw() + 
    theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
    guides(fill=guide_legend(title="predictor")) + 
    xlab("SPL") + 
    ylab(keep_col)
  names(gg)[i] = keep_col
  
}
gg[1]
gg[2]
gg[3]
gg[4]
gg[5]

gg$performance <- arrangeGrob(grobs = gg[c(1,2,5)], ncol = 1)
plot(gg$performance)

### Scatterplots
RT_predictors <- c("files_test", RT_predictors)
RT_predictors

dt <- lapply(RT_predictors, get) %>%
  lapply(bind_rows, .id = "file") 
names(dt) <- RT_predictors
dt <- dt[2:length(dt)] %>%
  bind_rows(.id = "method") %>%
  left_join(dt$files_test) %>%
  mutate(SPL = str_split_fixed(file, pattern = fixed("."), n = 3)[,2]) %>%
  select(-file)

library(ggpubr)
gg$scatter_correlation <- ggscatter(dt, x="RT", y="RT_pred", 
                  fill = "lightgrey", 
                  color = "black", shape = 21, size = 2, alpha = 0.3,
                  add.params = list(color = "orange", fill = "blue"), # Customize reg. line
                  facet.by = c("SPL", "method"),
                  title = paste("Predicted vs Observed retention time"), 
                  conf.int = T, 
                  conf.int.level = 1 - 10^-15,
                  # conf.int.level = 0.95,
                  add = "reg.line",
                  cor.coef = T, 
                  ggtheme = theme_bw())
plot(gg$scatter_correlation)

for (i in 1:5) {
  ggsave(gg[[i]], filename = paste0("results_train_500/plots/",names(gg)[i],".pdf"), width = 8, height = 6, dpi = "retina")
  ggsave(gg[[i]], filename = paste0("results_train_500/plots/",names(gg)[i],".png"), width = 8, height = 6, dpi = "retina")
}

for (i in 6:7) {
  ggsave(gg[[i]], filename = paste0("results_train_500/plots/",names(gg)[i],".pdf"), width = 16, height = 9, dpi = "retina")
  ggsave(gg[[i]], filename = paste0("results_train_500/plots/",names(gg)[i],".png"), width = 16, height = 9, dpi = "retina")
}
