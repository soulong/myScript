library(rio)
options(rio.import.class = "tbl")
library(tidyverse)
library(tidymodels)
library(vip)

# rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd()
setwd("C:\\Users\\haohe\\Desktop")


auc <- import("elife-78810-supp4-v2.xlsx") %>% print()
expr <- import("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv") %>% print()
mut <- import("OmicsSomaticMutationsMatrixHotspot.csv") %>% print()

yap_related_genes <- c("YAP1", "WWTR1", 
                       "TEAD1", "TEAD2", "TEAD3", "TEAD4", 
                       "VGLL1", "VGLL2", "VGLL3", "VGLL4"
                       # "NF2", "FAT1", "GNAQ", "GNA11"
)
expr_related <- expr %>%
  # select(ModelID, matches(str_c("^", yap_related_genes))) %>%
  select(ModelID, where(is.numeric)) %>% 
  rename_with(\(x) str_replace(x, " (.*)", "")) %>%
  reframe(across(everything(), sum), .by = ModelID) %>%
  print()
mut_related <- mut %>%
  filter(IsDefaultEntryForModel=='Yes') %>% 
  select(ModelID, matches(str_c("^", yap_related_genes))) %>%
  rename_with(\(x) str_replace(x, " (.*)", "_mut")) %>%
  reframe(across(everything(), sum), .by = ModelID) %>%
  mutate(across(!ModelID, \(x) magrittr::is_greater_than(x, 0))) %>% 
  print()
# filter(mut, ModelID == 'ACH-001517') %>% glimpse() %>% head()
data <- auc %>%
  select(1:5) %>%
  rename(ModelID = Depmap_ID) %>%
  left_join(expr_related) %>%
  # left_join(mut_related) %>%
  mutate(AUC = as.numeric(AUC)) %>%
  filter(!is.na(AUC)) %>%
  glimpse()

data$AUC %>% hist(breaks=200)


# Split data into training and testing sets
set.seed(123)
data_split <- data %>% 
  # filter(AUC < 0.99) %>% 
  mutate(AUC = -log10(AUC)) %>% 
  initial_split(prop = 0.8, strata = AUC)
print(data_split)
train_data <- training(data_split)
test_data <- testing(data_split)

# Specify recipe (only use numeric predictors for modeling)
rec <- recipe(AUC ~ ., data = train_data) %>%
  update_role(where(is.character), new_role = "meta") %>%
  # step_zv(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_nzv(all_predictors()) %>% 
  print()

# Specify model (random forest) with feature importance and multi-cpu support
mod <- rand_forest() %>%
  set_engine("ranger", importance = "permutation", num.threads = parallel::detectCores()) %>%
  set_mode("regression")

# Create workflow
wf <- workflow() %>%
  add_recipe(rec) %>%
  add_model(mod)

# Fit model
fit <- wf %>% 
  fit(data = train_data)

# Predict on test set
preds <- augment(fit, test_data)
glimpse(preds)

# Evaluate performance
metrics(preds, truth = AUC, estimate = .pred)

preds %>% 
  ggplot(aes(AUC, .pred)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggpubr::stat_cor() +
  coord_equal() +
  # coord_fixed() +
  theme_bw()

# Model feature importance analysis
# vi(fit)
# fit$fit$fit$fit$variable.importance
vip(fit, num_features=50) +
  ggtitle("Feature Importance (Permutation)") +
  theme_bw()

# Get ranked feature importance as a tibble
importance_tbl <- vi(fit) %>%
  print()

