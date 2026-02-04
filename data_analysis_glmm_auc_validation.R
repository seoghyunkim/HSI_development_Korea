

#' # Session information
sessionInfo()

#' # Load library
rm(list=ls(all.names=T)) # clear all data
pacman::p_load(tidyverse, # load package
               glmmTMB,
               jtools,
               ggpubr,
               MuMIn,
               DHARMa,
               car,
               dplyr,
               broom.mixed,
               performance,
               pROC)



#' # Normalization function
#'- This function is to normalize probability value to HSI
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}


#' # Read data
#'# Read data 
#'read species data
HSI_DB <- read.csv("./data_outcome/HSI_DB.csv", header=TRUE) %>%
  select(-X) %>% # remove "X" column
  filter(!site == "YEONGSAN_RIVER_MANGWOL") %>% # remove site which is out of streamnetwork
  filter(!site == "YEONGSAN_RIVER_GUJEONG") %>% # remove site which is out of streamnetwork
  filter(!site == "DONGGANG_ST1-2") %>%
  mutate(species = str_replace_all(species, pattern = "\\s", replacement = "_") ) # replace space to underscore in species name

#' read cluster data
site_group <- read.csv("./data_outcome/site_group_dat_k.csv", header = TRUE) %>%
  select(-X) %>% # 'X' delete
  rename(site_group = Cluster) %>% # 'Cluster' rename 'site_group'
  mutate(
    site_group = case_when(
      site_group == 1 ~ "GRP1", # change from numeric to categorical name
      site_group == 2 ~ "GRP2" # change from numeric to categorical name
    )
  )


#'- merge data
HSI_DB2 <- HSI_DB %>%
  left_join(site_group, by = "site") %>%
  drop_na(site_group) %>%
  drop_na(substrate2) %>%
  select(basin, site_group, site_by_season, site, date, season, 
         sample_number, species, individual_number, 
         depth, velocity, substrate = substrate2, MBSNCD)%>%
  mutate(substrate = case_when(substrate == "bedrock" ~ "boulder",
                               TRUE ~ as.character(substrate))) %>% # fix typo
  drop_na()


##### check the data #####

length(unique(HSI_DB2$site) ) # site number n = 195
length(unique(HSI_DB2$basin) ) # basin number n = 6
list(unique(HSI_DB2$basin))
length(unique(HSI_DB2$species) ) # species number n = 98
sum(HSI_DB2$individual_number) # n = 88534 individuals
set.seed(1234)


##### Zacco_platypus model #####

#'## (1) Zacco_platypus model

#'### data management
Zacco_platypus_dat <- HSI_DB2 %>%
  mutate(Zacco_platypus_PA = ifelse(species == "Zacco_platypus", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Zacco_platypus_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Zacco_platypus_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Zacco_platypus_PA, depth, velocity, substrate, MBSNCD)

Zacco_platypus_Pre <- HSI_DB2 %>%
  mutate(Zacco_platypus_PA = ifelse(species == "Zacco_platypus", 1, 0)) %>%
  aggregate(Zacco_platypus_PA ~ site, "sum") %>%
  mutate(Zacco_platypus_Pre = if_else(Zacco_platypus_PA > 0, "Yes", "No")) %>%
  select(site, Zacco_platypus_Pre)

Zacco_platypus_dat <- Zacco_platypus_dat %>%
  left_join(Zacco_platypus_Pre, by = "site") %>%
  filter(Zacco_platypus_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Zacco_platypus_glmm01 <- glmmTMB(Zacco_platypus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Zacco_platypus_dat, na.action ="na.fail")

summary(Zacco_platypus_glmm01)



# prediction
Zacco_platypus_glmm01_global <- predict(Zacco_platypus_glmm01, newdata = Zacco_platypus_dat, type = "response")

# ROC production
Zacco_platypus_glmm01_global_roc <- roc(Zacco_platypus_dat$Zacco_platypus_PA, Zacco_platypus_glmm01_global)

# AUC
General_AUC_Zacco_platypus_global<- auc(Zacco_platypus_glmm01_global_roc)


#####cv_Zacco_platypus_global#####

Zacco_platypus_glmm01_indices <- sample(1:nrow(Zacco_platypus_dat))
Zacco_platypus_glmm01_fold_size <- floor(nrow(Zacco_platypus_dat) / 10)
Zacco_platypus_glmm01_cv_folds <- split(Zacco_platypus_glmm01_indices, 
                                        rep(1:10, each = Zacco_platypus_glmm01_fold_size, 
                                            length.out = length(Zacco_platypus_glmm01_indices)))

Zacco_platypus_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Zacco_platypus_glmm01_train_indices <- setdiff(1:nrow(Zacco_platypus_dat), Zacco_platypus_glmm01_cv_folds[[i]])
  Zacco_platypus_glmm01_validation_indices <- Zacco_platypus_glmm01_cv_folds[[i]]
  
  Zacco_platypus_glmm01_train_data <- Zacco_platypus_dat[Zacco_platypus_glmm01_train_indices, ]
  Zacco_platypus_glmm01_validation_data <- Zacco_platypus_dat[Zacco_platypus_glmm01_validation_indices, ]
  
  # GLMM training
  Zacco_platypus_glmm01_model <- glmmTMB(Zacco_platypus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                           (1|site) + (1|season), 
                                         family = binomial(link = "logit"), data = Zacco_platypus_glmm01_train_data)
  
  # prediction
  Zacco_platypus_glmm01_pred_probs <- predict(Zacco_platypus_glmm01_model, newdata = Zacco_platypus_glmm01_validation_data, type = "response")
  
  # ROC production
  Zacco_platypus_glmm01_roc_obj <- roc(Zacco_platypus_glmm01_validation_data$Zacco_platypus_PA, Zacco_platypus_glmm01_pred_probs)
  
  # fold -  AUC
  Zacco_platypus_glmm01_auc_values[i] <- auc(Zacco_platypus_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Zacco_platypus_glmm01_auc_values[i], 4)))
}

# mean AUC
Zacco_platypus_glmm01_mean_auc <- mean(Zacco_platypus_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Zacco_platypus_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Zacco_platypus_global <- data.frame(
  Type = c("General_AUC_Zacco_platypus_global", "CV_AUC_Zacco_platypus_global"),
  AUC = c(General_AUC_Zacco_platypus_global, Zacco_platypus_glmm01_mean_auc)
)







#'## (1) Zacco_platypus model
##### Zacco_platypus group1 model #####

#'## select Group1

Zacco_platypus_dat02 <- Zacco_platypus_dat %>%
  filter(site_group == "GRP1")

#'##### glmm
#'- combination of site and season are used as a random effect ``
Zacco_platypus_glmm02 <- glmmTMB(Zacco_platypus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Zacco_platypus_dat02, na.action ="na.fail")

summary(Zacco_platypus_glmm02)



# prediction
Zacco_platypus_glmm02_grp1 <- predict(Zacco_platypus_glmm02, newdata = Zacco_platypus_dat02, type = "response")

# ROC production
Zacco_platypus_glmm02_grp1_roc <- roc(Zacco_platypus_dat02$Zacco_platypus_PA, Zacco_platypus_glmm02_grp1)

# AUC
General_AUC_Zacco_platypus_grp1 <-auc(Zacco_platypus_glmm02_grp1_roc)


#####cv_Zacco_platypus_grp1#####

Zacco_platypus_glmm02_indices <- sample(1:nrow(Zacco_platypus_dat02))
Zacco_platypus_glmm02_fold_size <- floor(nrow(Zacco_platypus_dat02) / 10)
Zacco_platypus_glmm02_cv_folds <- split(Zacco_platypus_glmm02_indices, 
                                        rep(1:10, each = Zacco_platypus_glmm02_fold_size, 
                                            length.out = length(Zacco_platypus_glmm02_indices)))

Zacco_platypus_glmm02_auc_values <- numeric(10)

for (i in 1:10) {
  Zacco_platypus_glmm02_train_indices <- setdiff(1:nrow(Zacco_platypus_dat02), Zacco_platypus_glmm02_cv_folds[[i]])
  Zacco_platypus_glmm02_validation_indices <- Zacco_platypus_glmm02_cv_folds[[i]]
  
  Zacco_platypus_glmm02_train_dat02a <- Zacco_platypus_dat02[Zacco_platypus_glmm02_train_indices, ]
  Zacco_platypus_glmm02_validation_dat02a <- Zacco_platypus_dat02[Zacco_platypus_glmm02_validation_indices, ]
  
  # GLMM training
  Zacco_platypus_glmm02_model <- glmmTMB(Zacco_platypus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                           (1|site) + (1|season), 
                                         family = binomial(link = "logit"), data = Zacco_platypus_glmm02_train_dat02a)
  
  # prediction
  Zacco_platypus_glmm02_pred_probs <- predict(Zacco_platypus_glmm02_model, newdata = Zacco_platypus_glmm02_validation_dat02a, type = "response")
  
  # ROC production
  Zacco_platypus_glmm02_roc_obj <- roc(Zacco_platypus_glmm02_validation_dat02a$Zacco_platypus_PA, Zacco_platypus_glmm02_pred_probs)
  
  # fold -  AUC
  Zacco_platypus_glmm02_auc_values[i] <- auc(Zacco_platypus_glmm02_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Zacco_platypus_glmm02_auc_values[i], 4)))
}

# mean AUC
Zacco_platypus_glmm02_mean_auc <- mean(Zacco_platypus_glmm02_auc_values)
print(paste0("Mean AUC across folds: ", round(Zacco_platypus_glmm02_mean_auc, 4)))



#####AUC-cvAUC_grp1#####
AUC_cvAUC_Zacco_platypus_grp1 <- data.frame(
  Type = c("General_AUC_Zacco_platypus_grp1", "CV_AUC_Zacco_platypus_grp1"),
  AUC = c(General_AUC_Zacco_platypus_grp1, Zacco_platypus_glmm02_mean_auc)
)














##### Zacco_platypus group2 model #####

#'## select Group2
Zacco_platypus_dat03 <- Zacco_platypus_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Zacco_platypus_glmm03 <- glmmTMB(Zacco_platypus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Zacco_platypus_dat03, na.action ="na.fail")

summary(Zacco_platypus_glmm03)

# prediction
Zacco_platypus_glmm03_grp2 <- predict(Zacco_platypus_glmm03, newdata = Zacco_platypus_dat03, type = "response")

# ROC production
Zacco_platypus_glmm03_grp2_roc <- roc(Zacco_platypus_dat03$Zacco_platypus_PA, Zacco_platypus_glmm03_grp2)

# AUC
General_AUC_Zacco_platypus_grp2 <- auc(Zacco_platypus_glmm03_grp2_roc)




#####cv_Zacco_platypus_grp2#####

Zacco_platypus_glmm03_indices <- sample(1:nrow(Zacco_platypus_dat03))
Zacco_platypus_glmm03_fold_size <- floor(nrow(Zacco_platypus_dat03) / 10)
Zacco_platypus_glmm03_cv_folds <- split(Zacco_platypus_glmm03_indices, 
                                        rep(1:10, each = Zacco_platypus_glmm03_fold_size, 
                                            length.out = length(Zacco_platypus_glmm03_indices)))

Zacco_platypus_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Zacco_platypus_glmm03_train_indices <- setdiff(1:nrow(Zacco_platypus_dat03), Zacco_platypus_glmm03_cv_folds[[i]])
  Zacco_platypus_glmm03_validation_indices <- Zacco_platypus_glmm03_cv_folds[[i]]
  
  Zacco_platypus_glmm03_train_dat03 <- Zacco_platypus_dat03[Zacco_platypus_glmm03_train_indices, ]
  Zacco_platypus_glmm03_validation_dat03 <- Zacco_platypus_dat03[Zacco_platypus_glmm03_validation_indices, ]
  
  # GLMM training
  Zacco_platypus_glmm03_model <- glmmTMB(Zacco_platypus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                           (1|site) + (1|season), 
                                         family = binomial(link = "logit"), data = Zacco_platypus_glmm03_train_dat03)
  
  # prediction
  Zacco_platypus_glmm03_pred_probs <- predict(Zacco_platypus_glmm03_model, newdata = Zacco_platypus_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Zacco_platypus_glmm03_roc_obj <- roc(Zacco_platypus_glmm03_validation_dat03$Zacco_platypus_PA, Zacco_platypus_glmm03_pred_probs)
  
  # fold -  AUC
  Zacco_platypus_glmm03_auc_values[i] <- auc(Zacco_platypus_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Zacco_platypus_glmm03_auc_values[i], 4)))
}

# mean AUC
Zacco_platypus_glmm03_mean_auc <- mean(Zacco_platypus_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Zacco_platypus_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Zacco_platypus_grp2 <- data.frame(
  Type = c("General_AUC_Zacco_platypus_grp2", "CV_AUC_Zacco_platypus_grp2"),
  AUC = c(General_AUC_Zacco_platypus_grp2, Zacco_platypus_glmm03_mean_auc)
)


#####Zacco_platypus_AUC_merge#####


Zacco_platypus_AUC <- AUC_cvAUC_Zacco_platypus_global%>%
  rbind(AUC_cvAUC_Zacco_platypus_grp1)%>%
  rbind(AUC_cvAUC_Zacco_platypus_grp2)






##### Zacco_koreanus model #####

#'## (1) Zacco_koreanus model

#'### data management
Zacco_koreanus_dat <- HSI_DB2 %>%
  mutate(Zacco_koreanus_PA = ifelse(species == "Zacco_koreanus", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Zacco_koreanus_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Zacco_koreanus_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Zacco_koreanus_PA, depth, velocity, substrate, MBSNCD)

Zacco_koreanus_Pre <- HSI_DB2 %>%
  mutate(Zacco_koreanus_PA = ifelse(species == "Zacco_koreanus", 1, 0)) %>%
  aggregate(Zacco_koreanus_PA ~ site, "sum") %>%
  mutate(Zacco_koreanus_Pre = if_else(Zacco_koreanus_PA > 0, "Yes", "No")) %>%
  select(site, Zacco_koreanus_Pre)

Zacco_koreanus_dat <- Zacco_koreanus_dat %>%
  left_join(Zacco_koreanus_Pre, by = "site") %>%
  filter(Zacco_koreanus_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Zacco_koreanus_glmm01 <- glmmTMB(Zacco_koreanus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Zacco_koreanus_dat, na.action ="na.fail")

summary(Zacco_koreanus_glmm01)



# prediction
Zacco_koreanus_glmm01_global <- predict(Zacco_koreanus_glmm01, newdata = Zacco_koreanus_dat, type = "response")

# ROC production
Zacco_koreanus_glmm01_global_roc <- roc(Zacco_koreanus_dat$Zacco_koreanus_PA, Zacco_koreanus_glmm01_global)

# AUC
General_AUC_Zacco_koreanus_global<- auc(Zacco_koreanus_glmm01_global_roc)


#####cv_Zacco_koreanus_global#####

Zacco_koreanus_glmm01_indices <- sample(1:nrow(Zacco_koreanus_dat))
Zacco_koreanus_glmm01_fold_size <- floor(nrow(Zacco_koreanus_dat) / 10)
Zacco_koreanus_glmm01_cv_folds <- split(Zacco_koreanus_glmm01_indices, 
                                        rep(1:10, each = Zacco_koreanus_glmm01_fold_size, 
                                            length.out = length(Zacco_koreanus_glmm01_indices)))

Zacco_koreanus_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Zacco_koreanus_glmm01_train_indices <- setdiff(1:nrow(Zacco_koreanus_dat), Zacco_koreanus_glmm01_cv_folds[[i]])
  Zacco_koreanus_glmm01_validation_indices <- Zacco_koreanus_glmm01_cv_folds[[i]]
  
  Zacco_koreanus_glmm01_train_data <- Zacco_koreanus_dat[Zacco_koreanus_glmm01_train_indices, ]
  Zacco_koreanus_glmm01_validation_data <- Zacco_koreanus_dat[Zacco_koreanus_glmm01_validation_indices, ]
  
  # GLMM training
  Zacco_koreanus_glmm01_model <- glmmTMB(Zacco_koreanus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                           (1|site) + (1|season), 
                                         family = binomial(link = "logit"), data = Zacco_koreanus_glmm01_train_data)
  
  # prediction
  Zacco_koreanus_glmm01_pred_probs <- predict(Zacco_koreanus_glmm01_model, newdata = Zacco_koreanus_glmm01_validation_data, type = "response")
  
  # ROC production
  Zacco_koreanus_glmm01_roc_obj <- roc(Zacco_koreanus_glmm01_validation_data$Zacco_koreanus_PA, Zacco_koreanus_glmm01_pred_probs)
  
  # fold -  AUC
  Zacco_koreanus_glmm01_auc_values[i] <- auc(Zacco_koreanus_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Zacco_koreanus_glmm01_auc_values[i], 4)))
}

# mean AUC
Zacco_koreanus_glmm01_mean_auc <- mean(Zacco_koreanus_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Zacco_koreanus_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Zacco_koreanus_global <- data.frame(
  Type = c("General_AUC_Zacco_koreanus_global", "CV_AUC_Zacco_koreanus_global"),
  AUC = c(General_AUC_Zacco_koreanus_global, Zacco_koreanus_glmm01_mean_auc)
)

















##### Zacco_koreanus group2 model #####

#'## select Group2
Zacco_koreanus_dat03 <- Zacco_koreanus_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Zacco_koreanus_glmm03 <- glmmTMB(Zacco_koreanus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Zacco_koreanus_dat03, na.action ="na.fail")

summary(Zacco_koreanus_glmm03)

# prediction
Zacco_koreanus_glmm03_grp2 <- predict(Zacco_koreanus_glmm03, newdata = Zacco_koreanus_dat03, type = "response")

# ROC production
Zacco_koreanus_glmm03_grp2_roc <- roc(Zacco_koreanus_dat03$Zacco_koreanus_PA, Zacco_koreanus_glmm03_grp2)

# AUC
General_AUC_Zacco_koreanus_grp2 <- auc(Zacco_koreanus_glmm03_grp2_roc)




#####cv_Zacco_koreanus_grp2#####

Zacco_koreanus_glmm03_indices <- sample(1:nrow(Zacco_koreanus_dat03))
Zacco_koreanus_glmm03_fold_size <- floor(nrow(Zacco_koreanus_dat03) / 10)
Zacco_koreanus_glmm03_cv_folds <- split(Zacco_koreanus_glmm03_indices, 
                                        rep(1:10, each = Zacco_koreanus_glmm03_fold_size, 
                                            length.out = length(Zacco_koreanus_glmm03_indices)))

Zacco_koreanus_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Zacco_koreanus_glmm03_train_indices <- setdiff(1:nrow(Zacco_koreanus_dat03), Zacco_koreanus_glmm03_cv_folds[[i]])
  Zacco_koreanus_glmm03_validation_indices <- Zacco_koreanus_glmm03_cv_folds[[i]]
  
  Zacco_koreanus_glmm03_train_dat03 <- Zacco_koreanus_dat03[Zacco_koreanus_glmm03_train_indices, ]
  Zacco_koreanus_glmm03_validation_dat03 <- Zacco_koreanus_dat03[Zacco_koreanus_glmm03_validation_indices, ]
  
  # GLMM training
  Zacco_koreanus_glmm03_model <- glmmTMB(Zacco_koreanus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                           (1|site) + (1|season), 
                                         family = binomial(link = "logit"), data = Zacco_koreanus_glmm03_train_dat03)
  
  # prediction
  Zacco_koreanus_glmm03_pred_probs <- predict(Zacco_koreanus_glmm03_model, newdata = Zacco_koreanus_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Zacco_koreanus_glmm03_roc_obj <- roc(Zacco_koreanus_glmm03_validation_dat03$Zacco_koreanus_PA, Zacco_koreanus_glmm03_pred_probs)
  
  # fold -  AUC
  Zacco_koreanus_glmm03_auc_values[i] <- auc(Zacco_koreanus_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Zacco_koreanus_glmm03_auc_values[i], 4)))
}

# mean AUC
Zacco_koreanus_glmm03_mean_auc <- mean(Zacco_koreanus_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Zacco_koreanus_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Zacco_koreanus_grp2 <- data.frame(
  Type = c("General_AUC_Zacco_koreanus_grp2", "CV_AUC_Zacco_koreanus_grp2"),
  AUC = c(General_AUC_Zacco_koreanus_grp2, Zacco_koreanus_glmm03_mean_auc)
)


#####Zacco_koreanus_AUC_merge#####


Zacco_koreanus_AUC <- AUC_cvAUC_Zacco_koreanus_global%>%
  rbind(AUC_cvAUC_Zacco_koreanus_grp2)






##### Pungtungia_herzi model #####

#'## (1) Pungtungia_herzi model

#'### data management
Pungtungia_herzi_dat <- HSI_DB2 %>%
  mutate(Pungtungia_herzi_PA = ifelse(species == "Pungtungia_herzi", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Pungtungia_herzi_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Pungtungia_herzi_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Pungtungia_herzi_PA, depth, velocity, substrate, MBSNCD)

Pungtungia_herzi_Pre <- HSI_DB2 %>%
  mutate(Pungtungia_herzi_PA = ifelse(species == "Pungtungia_herzi", 1, 0)) %>%
  aggregate(Pungtungia_herzi_PA ~ site, "sum") %>%
  mutate(Pungtungia_herzi_Pre = if_else(Pungtungia_herzi_PA > 0, "Yes", "No")) %>%
  select(site, Pungtungia_herzi_Pre)

Pungtungia_herzi_dat <- Pungtungia_herzi_dat %>%
  left_join(Pungtungia_herzi_Pre, by = "site") %>%
  filter(Pungtungia_herzi_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Pungtungia_herzi_glmm01 <- glmmTMB(Pungtungia_herzi_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                     (1|site) + (1|season), family=binomial(link = "logit"), data= Pungtungia_herzi_dat, na.action ="na.fail")

summary(Pungtungia_herzi_glmm01)



# prediction
Pungtungia_herzi_glmm01_global <- predict(Pungtungia_herzi_glmm01, newdata = Pungtungia_herzi_dat, type = "response")

# ROC production
Pungtungia_herzi_glmm01_global_roc <- roc(Pungtungia_herzi_dat$Pungtungia_herzi_PA, Pungtungia_herzi_glmm01_global)

# AUC
General_AUC_Pungtungia_herzi_global<- auc(Pungtungia_herzi_glmm01_global_roc)


#####cv_Pungtungia_herzi_global#####

Pungtungia_herzi_glmm01_indices <- sample(1:nrow(Pungtungia_herzi_dat))
Pungtungia_herzi_glmm01_fold_size <- floor(nrow(Pungtungia_herzi_dat) / 10)
Pungtungia_herzi_glmm01_cv_folds <- split(Pungtungia_herzi_glmm01_indices, 
                                          rep(1:10, each = Pungtungia_herzi_glmm01_fold_size, 
                                              length.out = length(Pungtungia_herzi_glmm01_indices)))

Pungtungia_herzi_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Pungtungia_herzi_glmm01_train_indices <- setdiff(1:nrow(Pungtungia_herzi_dat), Pungtungia_herzi_glmm01_cv_folds[[i]])
  Pungtungia_herzi_glmm01_validation_indices <- Pungtungia_herzi_glmm01_cv_folds[[i]]
  
  Pungtungia_herzi_glmm01_train_data <- Pungtungia_herzi_dat[Pungtungia_herzi_glmm01_train_indices, ]
  Pungtungia_herzi_glmm01_validation_data <- Pungtungia_herzi_dat[Pungtungia_herzi_glmm01_validation_indices, ]
  
  Pungtungia_herzi_glmm01_train_data$substrate <- factor(Pungtungia_herzi_glmm01_train_data$substrate)
  
  Pungtungia_herzi_glmm01_validation_data$substrate <- factor(
    Pungtungia_herzi_glmm01_validation_data$substrate,
    levels = levels(Pungtungia_herzi_glmm01_train_data$substrate)
  )
  
  
  
  # GLMM training
  Pungtungia_herzi_glmm01_model <- glmmTMB(Pungtungia_herzi_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                             (1|site) + (1|season), 
                                           family = binomial(link = "logit"), data = Pungtungia_herzi_glmm01_train_data)
  
  # prediction
  Pungtungia_herzi_glmm01_pred_probs <- predict(Pungtungia_herzi_glmm01_model, newdata = Pungtungia_herzi_glmm01_validation_data, type = "response")
  
  # ROC production
  Pungtungia_herzi_glmm01_roc_obj <- roc(Pungtungia_herzi_glmm01_validation_data$Pungtungia_herzi_PA, Pungtungia_herzi_glmm01_pred_probs)
  
  # fold -  AUC
  Pungtungia_herzi_glmm01_auc_values[i] <- auc(Pungtungia_herzi_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Pungtungia_herzi_glmm01_auc_values[i], 4)))
}

# mean AUC
Pungtungia_herzi_glmm01_mean_auc <- mean(Pungtungia_herzi_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Pungtungia_herzi_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Pungtungia_herzi_global <- data.frame(
  Type = c("General_AUC_Pungtungia_herzi_global", "CV_AUC_Pungtungia_herzi_global"),
  AUC = c(General_AUC_Pungtungia_herzi_global, Pungtungia_herzi_glmm01_mean_auc)
)

















##### Pungtungia_herzi group2 model #####

#'## select Group2
Pungtungia_herzi_dat03 <- Pungtungia_herzi_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Pungtungia_herzi_glmm03 <- glmmTMB(Pungtungia_herzi_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                     (1|site) + (1|season), family=binomial(link = "logit"), data= Pungtungia_herzi_dat03, na.action ="na.fail")

summary(Pungtungia_herzi_glmm03)

# prediction
Pungtungia_herzi_glmm03_grp2 <- predict(Pungtungia_herzi_glmm03, newdata = Pungtungia_herzi_dat03, type = "response")

# ROC production
Pungtungia_herzi_glmm03_grp2_roc <- roc(Pungtungia_herzi_dat03$Pungtungia_herzi_PA, Pungtungia_herzi_glmm03_grp2)

# AUC
General_AUC_Pungtungia_herzi_grp2 <- auc(Pungtungia_herzi_glmm03_grp2_roc)




#####cv_Pungtungia_herzi_grp2#####

Pungtungia_herzi_glmm03_indices <- sample(1:nrow(Pungtungia_herzi_dat03))
Pungtungia_herzi_glmm03_fold_size <- floor(nrow(Pungtungia_herzi_dat03) / 10)
Pungtungia_herzi_glmm03_cv_folds <- split(Pungtungia_herzi_glmm03_indices, 
                                          rep(1:10, each = Pungtungia_herzi_glmm03_fold_size, 
                                              length.out = length(Pungtungia_herzi_glmm03_indices)))

Pungtungia_herzi_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Pungtungia_herzi_glmm03_train_indices <- setdiff(1:nrow(Pungtungia_herzi_dat03), Pungtungia_herzi_glmm03_cv_folds[[i]])
  Pungtungia_herzi_glmm03_validation_indices <- Pungtungia_herzi_glmm03_cv_folds[[i]]
  
  Pungtungia_herzi_glmm03_train_dat03 <- Pungtungia_herzi_dat03[Pungtungia_herzi_glmm03_train_indices, ]
  Pungtungia_herzi_glmm03_validation_dat03 <- Pungtungia_herzi_dat03[Pungtungia_herzi_glmm03_validation_indices, ]
 
  
  
  # GLMM training
  Pungtungia_herzi_glmm03_model <- glmmTMB(Pungtungia_herzi_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                             (1|site) + (1|season), 
                                           family = binomial(link = "logit"), data = Pungtungia_herzi_glmm03_train_dat03)
  
  # prediction
  Pungtungia_herzi_glmm03_pred_probs <- predict(Pungtungia_herzi_glmm03_model, newdata = Pungtungia_herzi_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Pungtungia_herzi_glmm03_roc_obj <- roc(Pungtungia_herzi_glmm03_validation_dat03$Pungtungia_herzi_PA, Pungtungia_herzi_glmm03_pred_probs)
  
  # fold -  AUC
  Pungtungia_herzi_glmm03_auc_values[i] <- auc(Pungtungia_herzi_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Pungtungia_herzi_glmm03_auc_values[i], 4)))
}

# mean AUC
Pungtungia_herzi_glmm03_mean_auc <- mean(Pungtungia_herzi_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Pungtungia_herzi_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Pungtungia_herzi_grp2 <- data.frame(
  Type = c("General_AUC_Pungtungia_herzi_grp2", "CV_AUC_Pungtungia_herzi_grp2"),
  AUC = c(General_AUC_Pungtungia_herzi_grp2, Pungtungia_herzi_glmm03_mean_auc)
)


#####Pungtungia_herzi_AUC_merge#####


Pungtungia_herzi_AUC <- AUC_cvAUC_Pungtungia_herzi_global%>%
  rbind(AUC_cvAUC_Pungtungia_herzi_grp2)









##### Pseudopungtungia_nigra model #####

#'## (1) Pseudopungtungia_nigra model

#'### data management
Pseudopungtungia_nigra_dat <- HSI_DB2 %>%
  mutate(Pseudopungtungia_nigra_PA = ifelse(species == "Pseudopungtungia_nigra", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Pseudopungtungia_nigra_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Pseudopungtungia_nigra_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Pseudopungtungia_nigra_PA, depth, velocity, substrate, MBSNCD)

Pseudopungtungia_nigra_Pre <- HSI_DB2 %>%
  mutate(Pseudopungtungia_nigra_PA = ifelse(species == "Pseudopungtungia_nigra", 1, 0)) %>%
  aggregate(Pseudopungtungia_nigra_PA ~ site, "sum") %>%
  mutate(Pseudopungtungia_nigra_Pre = if_else(Pseudopungtungia_nigra_PA > 0, "Yes", "No")) %>%
  select(site, Pseudopungtungia_nigra_Pre)

Pseudopungtungia_nigra_dat <- Pseudopungtungia_nigra_dat %>%
  left_join(Pseudopungtungia_nigra_Pre, by = "site") %>%
  filter(Pseudopungtungia_nigra_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Pseudopungtungia_nigra_glmm01 <- glmmTMB(Pseudopungtungia_nigra_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                           (1|site) + (1|season), family=binomial(link = "logit"), data= Pseudopungtungia_nigra_dat, na.action ="na.fail")

summary(Pseudopungtungia_nigra_glmm01)



# prediction
Pseudopungtungia_nigra_glmm01_global <- predict(Pseudopungtungia_nigra_glmm01, newdata = Pseudopungtungia_nigra_dat, type = "response")

# ROC production
Pseudopungtungia_nigra_glmm01_global_roc <- roc(Pseudopungtungia_nigra_dat$Pseudopungtungia_nigra_PA, Pseudopungtungia_nigra_glmm01_global)

# AUC
General_AUC_Pseudopungtungia_nigra_global<- auc(Pseudopungtungia_nigra_glmm01_global_roc)


#####cv_Pseudopungtungia_nigra_global#####

Pseudopungtungia_nigra_glmm01_indices <- sample(1:nrow(Pseudopungtungia_nigra_dat))
Pseudopungtungia_nigra_glmm01_fold_size <- floor(nrow(Pseudopungtungia_nigra_dat) / 10)
Pseudopungtungia_nigra_glmm01_cv_folds <- split(Pseudopungtungia_nigra_glmm01_indices, 
                                                rep(1:10, each = Pseudopungtungia_nigra_glmm01_fold_size, 
                                                    length.out = length(Pseudopungtungia_nigra_glmm01_indices)))

Pseudopungtungia_nigra_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Pseudopungtungia_nigra_glmm01_train_indices <- setdiff(1:nrow(Pseudopungtungia_nigra_dat), Pseudopungtungia_nigra_glmm01_cv_folds[[i]])
  Pseudopungtungia_nigra_glmm01_validation_indices <- Pseudopungtungia_nigra_glmm01_cv_folds[[i]]
  
  Pseudopungtungia_nigra_glmm01_train_data <- Pseudopungtungia_nigra_dat[Pseudopungtungia_nigra_glmm01_train_indices, ]
  Pseudopungtungia_nigra_glmm01_validation_data <- Pseudopungtungia_nigra_dat[Pseudopungtungia_nigra_glmm01_validation_indices, ]
  
  Pseudopungtungia_nigra_glmm01_train_data$substrate <- factor(Pseudopungtungia_nigra_glmm01_train_data$substrate)
  
  Pseudopungtungia_nigra_glmm01_validation_data$substrate <- factor(
    Pseudopungtungia_nigra_glmm01_validation_data$substrate,
    levels = levels(Pseudopungtungia_nigra_glmm01_train_data$substrate)
  )
  
  # GLMM training
  Pseudopungtungia_nigra_glmm01_model <- glmmTMB(Pseudopungtungia_nigra_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                                   (1|site) + (1|season), 
                                                 family = binomial(link = "logit"), data = Pseudopungtungia_nigra_glmm01_train_data)
  
  # prediction
  Pseudopungtungia_nigra_glmm01_pred_probs <- predict(Pseudopungtungia_nigra_glmm01_model, newdata = Pseudopungtungia_nigra_glmm01_validation_data, type = "response")
  
  # ROC production
  Pseudopungtungia_nigra_glmm01_roc_obj <- roc(Pseudopungtungia_nigra_glmm01_validation_data$Pseudopungtungia_nigra_PA, Pseudopungtungia_nigra_glmm01_pred_probs)
  
  # fold -  AUC
  Pseudopungtungia_nigra_glmm01_auc_values[i] <- auc(Pseudopungtungia_nigra_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Pseudopungtungia_nigra_glmm01_auc_values[i], 4)))
}

# mean AUC
Pseudopungtungia_nigra_glmm01_mean_auc <- mean(Pseudopungtungia_nigra_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Pseudopungtungia_nigra_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Pseudopungtungia_nigra_global <- data.frame(
  Type = c("General_AUC_Pseudopungtungia_nigra_global", "CV_AUC_Pseudopungtungia_nigra_global"),
  AUC = c(General_AUC_Pseudopungtungia_nigra_global, Pseudopungtungia_nigra_glmm01_mean_auc)
)

  
  
##### Pseudopungtungia_nigra group2 model #####
  
#'## select Group2
Pseudopungtungia_nigra_dat03 <- Pseudopungtungia_nigra_dat %>%
    filter(site_group == "GRP2")
  
  
#'##### glmm
#'- combination of site and season are used as a random effect 
Pseudopungtungia_nigra_glmm03 <- glmmTMB(Pseudopungtungia_nigra_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                             (1|site) + (1|season), family=binomial(link = "logit"), data= Pseudopungtungia_nigra_dat03, na.action ="na.fail")
  
summary(Pseudopungtungia_nigra_glmm03)
  
# prediction
Pseudopungtungia_nigra_glmm03_grp2 <- predict(Pseudopungtungia_nigra_glmm03, newdata = Pseudopungtungia_nigra_dat03, type = "response")
  
# ROC production
Pseudopungtungia_nigra_glmm03_grp2_roc <- roc(Pseudopungtungia_nigra_dat03$Pseudopungtungia_nigra_PA, Pseudopungtungia_nigra_glmm03_grp2)
  
# AUC
General_AUC_Pseudopungtungia_nigra_grp2 <- auc(Pseudopungtungia_nigra_glmm03_grp2_roc)
  
  
  
  
#####cv_Pseudopungtungia_nigra_grp2#####
  
Pseudopungtungia_nigra_glmm03_indices <- sample(1:nrow(Pseudopungtungia_nigra_dat03))
Pseudopungtungia_nigra_glmm03_fold_size <- floor(nrow(Pseudopungtungia_nigra_dat03) / 10)
Pseudopungtungia_nigra_glmm03_cv_folds <- split(Pseudopungtungia_nigra_glmm03_indices, 
                                                  rep(1:10, each = Pseudopungtungia_nigra_glmm03_fold_size, 
                                                      length.out = length(Pseudopungtungia_nigra_glmm03_indices)))
  
Pseudopungtungia_nigra_glmm03_auc_values <- numeric(10)
  
for (i in 1:10) {
Pseudopungtungia_nigra_glmm03_train_indices <- setdiff(1:nrow(Pseudopungtungia_nigra_dat03), Pseudopungtungia_nigra_glmm03_cv_folds[[i]])
Pseudopungtungia_nigra_glmm03_validation_indices <- Pseudopungtungia_nigra_glmm03_cv_folds[[i]]
    
Pseudopungtungia_nigra_glmm03_train_dat03 <- Pseudopungtungia_nigra_dat03[Pseudopungtungia_nigra_glmm03_train_indices, ]
Pseudopungtungia_nigra_glmm03_validation_dat03 <- Pseudopungtungia_nigra_dat03[Pseudopungtungia_nigra_glmm03_validation_indices, ]
    
    
Pseudopungtungia_nigra_glmm03_train_dat03$substrate <- factor(Pseudopungtungia_nigra_glmm03_train_dat03$substrate)
    
Pseudopungtungia_nigra_glmm03_validation_dat03$substrate <- factor(
      Pseudopungtungia_nigra_glmm03_validation_dat03$substrate,
      levels = levels(Pseudopungtungia_nigra_glmm03_train_dat03$substrate)
    )
    
    
# GLMM training
Pseudopungtungia_nigra_glmm03_model <- glmmTMB(Pseudopungtungia_nigra_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                     (1|site) + (1|season), 
                                                   family = binomial(link = "logit"), data = Pseudopungtungia_nigra_glmm03_train_dat03)
    
# prediction
Pseudopungtungia_nigra_glmm03_pred_probs <- predict(Pseudopungtungia_nigra_glmm03_model, newdata = Pseudopungtungia_nigra_glmm03_validation_dat03, type = "response")
    
# ROC production
Pseudopungtungia_nigra_glmm03_roc_obj <- roc(Pseudopungtungia_nigra_glmm03_validation_dat03$Pseudopungtungia_nigra_PA, Pseudopungtungia_nigra_glmm03_pred_probs)
    
# fold -  AUC
Pseudopungtungia_nigra_glmm03_auc_values[i] <- auc(Pseudopungtungia_nigra_glmm03_roc_obj)
    
print(paste0("Fold ", i, " AUC: ", round(Pseudopungtungia_nigra_glmm03_auc_values[i], 4)))
  }
  
# mean AUC
  Pseudopungtungia_nigra_glmm03_mean_auc <- mean(Pseudopungtungia_nigra_glmm03_auc_values)
  print(paste0("Mean AUC across folds: ", round(Pseudopungtungia_nigra_glmm03_mean_auc, 4)))
  
  
  
#####AUC-cvAUC_grp2#####
  
  
  
AUC_cvAUC_Pseudopungtungia_nigra_grp2 <- data.frame(
    Type = c("General_AUC_Pseudopungtungia_nigra_grp2", "CV_AUC_Pseudopungtungia_nigra_grp2"),
    AUC = c(General_AUC_Pseudopungtungia_nigra_grp2, Pseudopungtungia_nigra_glmm03_mean_auc)
  )
  
  
  
  
#####Pseudopungtungia_nigra_AUC_merge#####
  
  
Pseudopungtungia_nigra_AUC <- AUC_cvAUC_Pseudopungtungia_nigra_global%>%
    rbind(AUC_cvAUC_Pseudopungtungia_nigra_grp2)
  

  
  
  
  
  
  



  
  
  
  
  
  
  
  
  
##### Coreoleuciscus_splendidus model #####

#'## (1) Coreoleuciscus_splendidus model

#'### data management
Coreoleuciscus_splendidus_dat <- HSI_DB2 %>%
  mutate(Coreoleuciscus_splendidus_PA = ifelse(species == "Coreoleuciscus_splendidus", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Coreoleuciscus_splendidus_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Coreoleuciscus_splendidus_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Coreoleuciscus_splendidus_PA, depth, velocity, substrate, MBSNCD)

Coreoleuciscus_splendidus_Pre <- HSI_DB2 %>%
  mutate(Coreoleuciscus_splendidus_PA = ifelse(species == "Coreoleuciscus_splendidus", 1, 0)) %>%
  aggregate(Coreoleuciscus_splendidus_PA ~ site, "sum") %>%
  mutate(Coreoleuciscus_splendidus_Pre = if_else(Coreoleuciscus_splendidus_PA > 0, "Yes", "No")) %>%
  select(site, Coreoleuciscus_splendidus_Pre)

Coreoleuciscus_splendidus_dat <- Coreoleuciscus_splendidus_dat %>%
  left_join(Coreoleuciscus_splendidus_Pre, by = "site") %>%
  filter(Coreoleuciscus_splendidus_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Coreoleuciscus_splendidus_glmm01 <- glmmTMB(Coreoleuciscus_splendidus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                              (1|site) + (1|season), family=binomial(link = "logit"), data= Coreoleuciscus_splendidus_dat, na.action ="na.fail")

summary(Coreoleuciscus_splendidus_glmm01)



# prediction
Coreoleuciscus_splendidus_glmm01_global <- predict(Coreoleuciscus_splendidus_glmm01, newdata = Coreoleuciscus_splendidus_dat, type = "response")

# ROC production
Coreoleuciscus_splendidus_glmm01_global_roc <- roc(Coreoleuciscus_splendidus_dat$Coreoleuciscus_splendidus_PA, Coreoleuciscus_splendidus_glmm01_global)

# AUC
General_AUC_Coreoleuciscus_splendidus_global<- auc(Coreoleuciscus_splendidus_glmm01_global_roc)


#####cv_Coreoleuciscus_splendidus_global#####

Coreoleuciscus_splendidus_glmm01_indices <- sample(1:nrow(Coreoleuciscus_splendidus_dat))
Coreoleuciscus_splendidus_glmm01_fold_size <- floor(nrow(Coreoleuciscus_splendidus_dat) / 10)
Coreoleuciscus_splendidus_glmm01_cv_folds <- split(Coreoleuciscus_splendidus_glmm01_indices, 
                                                   rep(1:10, each = Coreoleuciscus_splendidus_glmm01_fold_size, 
                                                       length.out = length(Coreoleuciscus_splendidus_glmm01_indices)))

Coreoleuciscus_splendidus_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Coreoleuciscus_splendidus_glmm01_train_indices <- setdiff(1:nrow(Coreoleuciscus_splendidus_dat), Coreoleuciscus_splendidus_glmm01_cv_folds[[i]])
  Coreoleuciscus_splendidus_glmm01_validation_indices <- Coreoleuciscus_splendidus_glmm01_cv_folds[[i]]
  
  Coreoleuciscus_splendidus_glmm01_train_data <- Coreoleuciscus_splendidus_dat[Coreoleuciscus_splendidus_glmm01_train_indices, ]
  Coreoleuciscus_splendidus_glmm01_validation_data <- Coreoleuciscus_splendidus_dat[Coreoleuciscus_splendidus_glmm01_validation_indices, ]
  
  Coreoleuciscus_splendidus_glmm01_train_data$substrate <- factor(Coreoleuciscus_splendidus_glmm01_train_data$substrate)
  
  Coreoleuciscus_splendidus_glmm01_validation_data$substrate <- factor(
    Coreoleuciscus_splendidus_glmm01_validation_data$substrate,
    levels = levels(Coreoleuciscus_splendidus_glmm01_train_data$substrate)
  )
  
  # GLMM training
  Coreoleuciscus_splendidus_glmm01_model <- glmmTMB(Coreoleuciscus_splendidus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                      (1|site) + (1|season), 
                                                    family = binomial(link = "logit"), data = Coreoleuciscus_splendidus_glmm01_train_data)
  
  # prediction
  Coreoleuciscus_splendidus_glmm01_pred_probs <- predict(Coreoleuciscus_splendidus_glmm01_model, newdata = Coreoleuciscus_splendidus_glmm01_validation_data, type = "response")
  
  # ROC production
  Coreoleuciscus_splendidus_glmm01_roc_obj <- roc(Coreoleuciscus_splendidus_glmm01_validation_data$Coreoleuciscus_splendidus_PA, Coreoleuciscus_splendidus_glmm01_pred_probs)
  
  # fold -  AUC
  Coreoleuciscus_splendidus_glmm01_auc_values[i] <- auc(Coreoleuciscus_splendidus_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Coreoleuciscus_splendidus_glmm01_auc_values[i], 4)))
}

# mean AUC
Coreoleuciscus_splendidus_glmm01_mean_auc <- mean(Coreoleuciscus_splendidus_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Coreoleuciscus_splendidus_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Coreoleuciscus_splendidus_global <- data.frame(
  Type = c("General_AUC_Coreoleuciscus_splendidus_global", "CV_AUC_Coreoleuciscus_splendidus_global"),
  AUC = c(General_AUC_Coreoleuciscus_splendidus_global, Coreoleuciscus_splendidus_glmm01_mean_auc)
)



##### Coreoleuciscus_splendidus group2 model #####

#'## select Group2
Coreoleuciscus_splendidus_dat03 <- Coreoleuciscus_splendidus_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Coreoleuciscus_splendidus_glmm03 <- glmmTMB(Coreoleuciscus_splendidus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                              (1|site) + (1|season), family=binomial(link = "logit"), data= Coreoleuciscus_splendidus_dat03, na.action ="na.fail")

summary(Coreoleuciscus_splendidus_glmm03)

# prediction
Coreoleuciscus_splendidus_glmm03_grp2 <- predict(Coreoleuciscus_splendidus_glmm03, newdata = Coreoleuciscus_splendidus_dat03, type = "response")

# ROC production
Coreoleuciscus_splendidus_glmm03_grp2_roc <- roc(Coreoleuciscus_splendidus_dat03$Coreoleuciscus_splendidus_PA, Coreoleuciscus_splendidus_glmm03_grp2)

# AUC
General_AUC_Coreoleuciscus_splendidus_grp2 <- auc(Coreoleuciscus_splendidus_glmm03_grp2_roc)




#####cv_Coreoleuciscus_splendidus_grp2#####

Coreoleuciscus_splendidus_glmm03_indices <- sample(1:nrow(Coreoleuciscus_splendidus_dat03))
Coreoleuciscus_splendidus_glmm03_fold_size <- floor(nrow(Coreoleuciscus_splendidus_dat03) / 10)
Coreoleuciscus_splendidus_glmm03_cv_folds <- split(Coreoleuciscus_splendidus_glmm03_indices, 
                                                   rep(1:10, each = Coreoleuciscus_splendidus_glmm03_fold_size, 
                                                       length.out = length(Coreoleuciscus_splendidus_glmm03_indices)))

Coreoleuciscus_splendidus_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Coreoleuciscus_splendidus_glmm03_train_indices <- setdiff(1:nrow(Coreoleuciscus_splendidus_dat03), Coreoleuciscus_splendidus_glmm03_cv_folds[[i]])
  Coreoleuciscus_splendidus_glmm03_validation_indices <- Coreoleuciscus_splendidus_glmm03_cv_folds[[i]]
  
  Coreoleuciscus_splendidus_glmm03_train_dat03 <- Coreoleuciscus_splendidus_dat03[Coreoleuciscus_splendidus_glmm03_train_indices, ]
  Coreoleuciscus_splendidus_glmm03_validation_dat03 <- Coreoleuciscus_splendidus_dat03[Coreoleuciscus_splendidus_glmm03_validation_indices, ]
  
  
  Coreoleuciscus_splendidus_glmm03_train_dat03$substrate <- factor(Coreoleuciscus_splendidus_glmm03_train_dat03$substrate)
  
  Coreoleuciscus_splendidus_glmm03_validation_dat03$substrate <- factor(
    Coreoleuciscus_splendidus_glmm03_validation_dat03$substrate,
    levels = levels(Coreoleuciscus_splendidus_glmm03_train_dat03$substrate)
  )
  
  
  # GLMM training
  Coreoleuciscus_splendidus_glmm03_model <- glmmTMB(Coreoleuciscus_splendidus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                      (1|site) + (1|season), 
                                                    family = binomial(link = "logit"), data = Coreoleuciscus_splendidus_glmm03_train_dat03)
  
  # prediction
  Coreoleuciscus_splendidus_glmm03_pred_probs <- predict(Coreoleuciscus_splendidus_glmm03_model, newdata = Coreoleuciscus_splendidus_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Coreoleuciscus_splendidus_glmm03_roc_obj <- roc(Coreoleuciscus_splendidus_glmm03_validation_dat03$Coreoleuciscus_splendidus_PA, Coreoleuciscus_splendidus_glmm03_pred_probs)
  
  # fold -  AUC
  Coreoleuciscus_splendidus_glmm03_auc_values[i] <- auc(Coreoleuciscus_splendidus_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Coreoleuciscus_splendidus_glmm03_auc_values[i], 4)))
}

# mean AUC
Coreoleuciscus_splendidus_glmm03_mean_auc <- mean(Coreoleuciscus_splendidus_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Coreoleuciscus_splendidus_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Coreoleuciscus_splendidus_grp2 <- data.frame(
  Type = c("General_AUC_Coreoleuciscus_splendidus_grp2", "CV_AUC_Coreoleuciscus_splendidus_grp2"),
  AUC = c(General_AUC_Coreoleuciscus_splendidus_grp2, Coreoleuciscus_splendidus_glmm03_mean_auc)
)


  
  
  
#####Coreoleuciscus_splendidus_AUC_merge#####


Coreoleuciscus_splendidus_AUC <- AUC_cvAUC_Coreoleuciscus_splendidus_global%>%
  rbind(AUC_cvAUC_Coreoleuciscus_splendidus_grp2)
  
  

  
  
  ##### Opsariichthys_uncirostris_amurensis model #####
  
  #'## (1) Opsariichthys_uncirostris_amurensis model
  
  #'### data management
  Opsariichthys_uncirostris_amurensis_dat <- HSI_DB2 %>%
    mutate(Opsariichthys_uncirostris_amurensis_PA = ifelse(species == "Opsariichthys_uncirostris_amurensis", 1, 0)) %>% # create column based on P/A of target species
    distinct(site_by_season, sample_number, Opsariichthys_uncirostris_amurensis_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
    arrange(basin, site_by_season, sample_number, desc(Opsariichthys_uncirostris_amurensis_PA)) %>% # arrange order
    group_by(site_by_season, sample_number) %>% 
    slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
    select(basin, site_group, site_by_season, site, date, season, sample_number, Opsariichthys_uncirostris_amurensis_PA, depth, velocity, substrate, MBSNCD)
  
  Opsariichthys_uncirostris_amurensis_Pre <- HSI_DB2 %>%
    mutate(Opsariichthys_uncirostris_amurensis_PA = ifelse(species == "Opsariichthys_uncirostris_amurensis", 1, 0)) %>%
    aggregate(Opsariichthys_uncirostris_amurensis_PA ~ site, "sum") %>%
    mutate(Opsariichthys_uncirostris_amurensis_Pre = if_else(Opsariichthys_uncirostris_amurensis_PA > 0, "Yes", "No")) %>%
    select(site, Opsariichthys_uncirostris_amurensis_Pre)
  
  Opsariichthys_uncirostris_amurensis_dat <- Opsariichthys_uncirostris_amurensis_dat %>%
    left_join(Opsariichthys_uncirostris_amurensis_Pre, by = "site") %>%
    filter(Opsariichthys_uncirostris_amurensis_Pre == "Yes")
  
  
  #'### Global Model
  #'##### glmm
  #'- combination of site and season are used as a random effect 
  Opsariichthys_uncirostris_amurensis_glmm01 <- glmmTMB(Opsariichthys_uncirostris_amurensis_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                                          (1|site) + (1|season), family=binomial(link = "logit"), data= Opsariichthys_uncirostris_amurensis_dat, na.action ="na.fail")
  
  summary(Opsariichthys_uncirostris_amurensis_glmm01)
  
  
  
  # prediction
  Opsariichthys_uncirostris_amurensis_glmm01_global <- predict(Opsariichthys_uncirostris_amurensis_glmm01, newdata = Opsariichthys_uncirostris_amurensis_dat, type = "response")
  
  # ROC production
  Opsariichthys_uncirostris_amurensis_glmm01_global_roc <- roc(Opsariichthys_uncirostris_amurensis_dat$Opsariichthys_uncirostris_amurensis_PA, Opsariichthys_uncirostris_amurensis_glmm01_global)
  
  # AUC
  General_AUC_Opsariichthys_uncirostris_amurensis_global<- auc(Opsariichthys_uncirostris_amurensis_glmm01_global_roc)
  
  
  #####cv_Opsariichthys_uncirostris_amurensis_global#####
  
  Opsariichthys_uncirostris_amurensis_glmm01_indices <- sample(1:nrow(Opsariichthys_uncirostris_amurensis_dat))
  Opsariichthys_uncirostris_amurensis_glmm01_fold_size <- floor(nrow(Opsariichthys_uncirostris_amurensis_dat) / 10)
  Opsariichthys_uncirostris_amurensis_glmm01_cv_folds <- split(Opsariichthys_uncirostris_amurensis_glmm01_indices, 
                                                               rep(1:10, each = Opsariichthys_uncirostris_amurensis_glmm01_fold_size, 
                                                                   length.out = length(Opsariichthys_uncirostris_amurensis_glmm01_indices)))
  
  Opsariichthys_uncirostris_amurensis_glmm01_auc_values <- numeric(10)
  
  for (i in 1:10) {
    Opsariichthys_uncirostris_amurensis_glmm01_train_indices <- setdiff(1:nrow(Opsariichthys_uncirostris_amurensis_dat), Opsariichthys_uncirostris_amurensis_glmm01_cv_folds[[i]])
    Opsariichthys_uncirostris_amurensis_glmm01_validation_indices <- Opsariichthys_uncirostris_amurensis_glmm01_cv_folds[[i]]
    
    Opsariichthys_uncirostris_amurensis_glmm01_train_data <- Opsariichthys_uncirostris_amurensis_dat[Opsariichthys_uncirostris_amurensis_glmm01_train_indices, ]
    Opsariichthys_uncirostris_amurensis_glmm01_validation_data <- Opsariichthys_uncirostris_amurensis_dat[Opsariichthys_uncirostris_amurensis_glmm01_validation_indices, ]
    
    # GLMM training
    Opsariichthys_uncirostris_amurensis_glmm01_model <- glmmTMB(Opsariichthys_uncirostris_amurensis_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                                  (1|site) + (1|season), 
                                                                family = binomial(link = "logit"), data = Opsariichthys_uncirostris_amurensis_glmm01_train_data)
    
    # prediction
    Opsariichthys_uncirostris_amurensis_glmm01_pred_probs <- predict(Opsariichthys_uncirostris_amurensis_glmm01_model, newdata = Opsariichthys_uncirostris_amurensis_glmm01_validation_data, type = "response")
    
    # ROC production
    Opsariichthys_uncirostris_amurensis_glmm01_roc_obj <- roc(Opsariichthys_uncirostris_amurensis_glmm01_validation_data$Opsariichthys_uncirostris_amurensis_PA, Opsariichthys_uncirostris_amurensis_glmm01_pred_probs)
    
    # fold -  AUC
    Opsariichthys_uncirostris_amurensis_glmm01_auc_values[i] <- auc(Opsariichthys_uncirostris_amurensis_glmm01_roc_obj)
    
    print(paste0("Fold ", i, " AUC: ", round(Opsariichthys_uncirostris_amurensis_glmm01_auc_values[i], 4)))
  }
  
  # mean AUC
  Opsariichthys_uncirostris_amurensis_glmm01_mean_auc <- mean(Opsariichthys_uncirostris_amurensis_glmm01_auc_values)
  print(paste0("Mean AUC across folds: ", round(Opsariichthys_uncirostris_amurensis_glmm01_mean_auc, 4)))
  
  
  
  #####AUC-cvAUC_global#####
  
  
  
  AUC_cvAUC_Opsariichthys_uncirostris_amurensis_global <- data.frame(
    Type = c("General_AUC_Opsariichthys_uncirostris_amurensis_global", "CV_AUC_Opsariichthys_uncirostris_amurensis_global"),
    AUC = c(General_AUC_Opsariichthys_uncirostris_amurensis_global, Opsariichthys_uncirostris_amurensis_glmm01_mean_auc)
  )
  
  
  
  
  
  
  
  #'## (1) Opsariichthys_uncirostris_amurensis model
  ##### Opsariichthys_uncirostris_amurensis group1 model #####
  
  #'## select Group1
  
  Opsariichthys_uncirostris_amurensis_dat02 <- Opsariichthys_uncirostris_amurensis_dat %>%
    filter(site_group == "GRP1")
  
  #'##### glmm
  #'- combination of site and season are used as a random effect ``
  Opsariichthys_uncirostris_amurensis_glmm02 <- glmmTMB(Opsariichthys_uncirostris_amurensis_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                                          (1|site) + (1|season), family=binomial(link = "logit"), data= Opsariichthys_uncirostris_amurensis_dat02, na.action ="na.fail")
  
  summary(Opsariichthys_uncirostris_amurensis_glmm02)
  
  
  
  # prediction
  Opsariichthys_uncirostris_amurensis_glmm02_grp1 <- predict(Opsariichthys_uncirostris_amurensis_glmm02, newdata = Opsariichthys_uncirostris_amurensis_dat02, type = "response")
  
  # ROC production
  Opsariichthys_uncirostris_amurensis_glmm02_grp1_roc <- roc(Opsariichthys_uncirostris_amurensis_dat02$Opsariichthys_uncirostris_amurensis_PA, Opsariichthys_uncirostris_amurensis_glmm02_grp1)
  
  # AUC
  General_AUC_Opsariichthys_uncirostris_amurensis_grp1 <-auc(Opsariichthys_uncirostris_amurensis_glmm02_grp1_roc)
  
  
  #####cv_Opsariichthys_uncirostris_amurensis_grp1#####
  
  Opsariichthys_uncirostris_amurensis_glmm02_indices <- sample(1:nrow(Opsariichthys_uncirostris_amurensis_dat02))
  Opsariichthys_uncirostris_amurensis_glmm02_fold_size <- floor(nrow(Opsariichthys_uncirostris_amurensis_dat02) / 10)
  Opsariichthys_uncirostris_amurensis_glmm02_cv_folds <- split(Opsariichthys_uncirostris_amurensis_glmm02_indices, 
                                                               rep(1:10, each = Opsariichthys_uncirostris_amurensis_glmm02_fold_size, 
                                                                   length.out = length(Opsariichthys_uncirostris_amurensis_glmm02_indices)))
  
  Opsariichthys_uncirostris_amurensis_glmm02_auc_values <- numeric(10)
  
  for (i in 1:10) {
    Opsariichthys_uncirostris_amurensis_glmm02_train_indices <- setdiff(1:nrow(Opsariichthys_uncirostris_amurensis_dat02), Opsariichthys_uncirostris_amurensis_glmm02_cv_folds[[i]])
    Opsariichthys_uncirostris_amurensis_glmm02_validation_indices <- Opsariichthys_uncirostris_amurensis_glmm02_cv_folds[[i]]
    
    Opsariichthys_uncirostris_amurensis_glmm02_train_dat02a <- Opsariichthys_uncirostris_amurensis_dat02[Opsariichthys_uncirostris_amurensis_glmm02_train_indices, ]
    Opsariichthys_uncirostris_amurensis_glmm02_validation_dat02a <- Opsariichthys_uncirostris_amurensis_dat02[Opsariichthys_uncirostris_amurensis_glmm02_validation_indices, ]
    
    # GLMM training
    Opsariichthys_uncirostris_amurensis_glmm02_model <- glmmTMB(Opsariichthys_uncirostris_amurensis_PA ~ poly(depth, 1) + poly(velocity,1) + substrate + 
                                                                  (1|site) + (1|season), 
                                                                family = binomial(link = "logit"), data = Opsariichthys_uncirostris_amurensis_glmm02_train_dat02a)
    
    # prediction
    Opsariichthys_uncirostris_amurensis_glmm02_pred_probs <- predict(Opsariichthys_uncirostris_amurensis_glmm02_model, newdata = Opsariichthys_uncirostris_amurensis_glmm02_validation_dat02a, type = "response")
    
    # ROC production
    Opsariichthys_uncirostris_amurensis_glmm02_roc_obj <- roc(Opsariichthys_uncirostris_amurensis_glmm02_validation_dat02a$Opsariichthys_uncirostris_amurensis_PA, Opsariichthys_uncirostris_amurensis_glmm02_pred_probs)
    
    # fold -  AUC
    Opsariichthys_uncirostris_amurensis_glmm02_auc_values[i] <- auc(Opsariichthys_uncirostris_amurensis_glmm02_roc_obj)
    
    print(paste0("Fold ", i, " AUC: ", round(Opsariichthys_uncirostris_amurensis_glmm02_auc_values[i], 4)))
  }
  
  # mean AUC
  Opsariichthys_uncirostris_amurensis_glmm02_mean_auc <- mean(Opsariichthys_uncirostris_amurensis_glmm02_auc_values)
  print(paste0("Mean AUC across folds: ", round(Opsariichthys_uncirostris_amurensis_glmm02_mean_auc, 4)))
  
  
  
  #####AUC-cvAUC_grp1#####
  AUC_cvAUC_Opsariichthys_uncirostris_amurensis_grp1 <- data.frame(
    Type = c("General_AUC_Opsariichthys_uncirostris_amurensis_grp1", "CV_AUC_Opsariichthys_uncirostris_amurensis_grp1"),
    AUC = c(General_AUC_Opsariichthys_uncirostris_amurensis_grp1, Opsariichthys_uncirostris_amurensis_glmm02_mean_auc)
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ##### Opsariichthys_uncirostris_amurensis group2 model #####
  
  #'## select Group2
  Opsariichthys_uncirostris_amurensis_dat03 <- Opsariichthys_uncirostris_amurensis_dat %>%
    filter(site_group == "GRP2")
  
  
  #'##### glmm
  #'- combination of site and season are used as a random effect 
  Opsariichthys_uncirostris_amurensis_glmm03 <- glmmTMB(Opsariichthys_uncirostris_amurensis_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                          (1|site) + (1|season), family=binomial(link = "logit"), data= Opsariichthys_uncirostris_amurensis_dat03, na.action ="na.fail")
  
  summary(Opsariichthys_uncirostris_amurensis_glmm03)
  
  # prediction
  Opsariichthys_uncirostris_amurensis_glmm03_grp2 <- predict(Opsariichthys_uncirostris_amurensis_glmm03, newdata = Opsariichthys_uncirostris_amurensis_dat03, type = "response")
  
  # ROC production
  Opsariichthys_uncirostris_amurensis_glmm03_grp2_roc <- roc(Opsariichthys_uncirostris_amurensis_dat03$Opsariichthys_uncirostris_amurensis_PA, Opsariichthys_uncirostris_amurensis_glmm03_grp2)
  
  # AUC
  General_AUC_Opsariichthys_uncirostris_amurensis_grp2 <- auc(Opsariichthys_uncirostris_amurensis_glmm03_grp2_roc)
  
  
  
  
  #####cv_Opsariichthys_uncirostris_amurensis_grp2#####
  
  Opsariichthys_uncirostris_amurensis_glmm03_indices <- sample(1:nrow(Opsariichthys_uncirostris_amurensis_dat03))
  Opsariichthys_uncirostris_amurensis_glmm03_fold_size <- floor(nrow(Opsariichthys_uncirostris_amurensis_dat03) / 10)
  Opsariichthys_uncirostris_amurensis_glmm03_cv_folds <- split(Opsariichthys_uncirostris_amurensis_glmm03_indices, 
                                                               rep(1:10, each = Opsariichthys_uncirostris_amurensis_glmm03_fold_size, 
                                                                   length.out = length(Opsariichthys_uncirostris_amurensis_glmm03_indices)))
  
  Opsariichthys_uncirostris_amurensis_glmm03_auc_values <- numeric(10)
  
  for (i in 1:10) {
    Opsariichthys_uncirostris_amurensis_glmm03_train_indices <- setdiff(1:nrow(Opsariichthys_uncirostris_amurensis_dat03), Opsariichthys_uncirostris_amurensis_glmm03_cv_folds[[i]])
    Opsariichthys_uncirostris_amurensis_glmm03_validation_indices <- Opsariichthys_uncirostris_amurensis_glmm03_cv_folds[[i]]
    
    Opsariichthys_uncirostris_amurensis_glmm03_train_dat03 <- Opsariichthys_uncirostris_amurensis_dat03[Opsariichthys_uncirostris_amurensis_glmm03_train_indices, ]
    Opsariichthys_uncirostris_amurensis_glmm03_validation_dat03 <- Opsariichthys_uncirostris_amurensis_dat03[Opsariichthys_uncirostris_amurensis_glmm03_validation_indices, ]
    
    # GLMM training
    Opsariichthys_uncirostris_amurensis_glmm03_model <- glmmTMB(Opsariichthys_uncirostris_amurensis_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                                  (1|site) + (1|season), 
                                                                family = binomial(link = "logit"), data = Opsariichthys_uncirostris_amurensis_glmm03_train_dat03)
    
    # prediction
    Opsariichthys_uncirostris_amurensis_glmm03_pred_probs <- predict(Opsariichthys_uncirostris_amurensis_glmm03_model, newdata = Opsariichthys_uncirostris_amurensis_glmm03_validation_dat03, type = "response")
    
    # ROC production
    Opsariichthys_uncirostris_amurensis_glmm03_roc_obj <- roc(Opsariichthys_uncirostris_amurensis_glmm03_validation_dat03$Opsariichthys_uncirostris_amurensis_PA, Opsariichthys_uncirostris_amurensis_glmm03_pred_probs)
    
    # fold -  AUC
    Opsariichthys_uncirostris_amurensis_glmm03_auc_values[i] <- auc(Opsariichthys_uncirostris_amurensis_glmm03_roc_obj)
    
    print(paste0("Fold ", i, " AUC: ", round(Opsariichthys_uncirostris_amurensis_glmm03_auc_values[i], 4)))
  }
  
  # mean AUC
  Opsariichthys_uncirostris_amurensis_glmm03_mean_auc <- mean(Opsariichthys_uncirostris_amurensis_glmm03_auc_values)
  print(paste0("Mean AUC across folds: ", round(Opsariichthys_uncirostris_amurensis_glmm03_mean_auc, 4)))
  
  
  
  #####AUC-cvAUC_grp2#####
  
  
  
  AUC_cvAUC_Opsariichthys_uncirostris_amurensis_grp2 <- data.frame(
    Type = c("General_AUC_Opsariichthys_uncirostris_amurensis_grp2", "CV_AUC_Opsariichthys_uncirostris_amurensis_grp2"),
    AUC = c(General_AUC_Opsariichthys_uncirostris_amurensis_grp2, Opsariichthys_uncirostris_amurensis_glmm03_mean_auc)
  )
  
  
  #####Opsariichthys_uncirostris_amurensis_AUC_merge#####
  
  
  Opsariichthys_uncirostris_amurensis_AUC <- AUC_cvAUC_Opsariichthys_uncirostris_amurensis_global%>%
    rbind(AUC_cvAUC_Opsariichthys_uncirostris_amurensis_grp1)%>%
    rbind(AUC_cvAUC_Opsariichthys_uncirostris_amurensis_grp2)
  
  

##### Pseudogobio_esocinus model #####

#'## (1) Pseudogobio_esocinus model

#'### data management
Pseudogobio_esocinus_dat <- HSI_DB2 %>%
  mutate(Pseudogobio_esocinus_PA = ifelse(species == "Pseudogobio_esocinus", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Pseudogobio_esocinus_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Pseudogobio_esocinus_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Pseudogobio_esocinus_PA, depth, velocity, substrate, MBSNCD)

Pseudogobio_esocinus_Pre <- HSI_DB2 %>%
  mutate(Pseudogobio_esocinus_PA = ifelse(species == "Pseudogobio_esocinus", 1, 0)) %>%
  aggregate(Pseudogobio_esocinus_PA ~ site, "sum") %>%
  mutate(Pseudogobio_esocinus_Pre = if_else(Pseudogobio_esocinus_PA > 0, "Yes", "No")) %>%
  select(site, Pseudogobio_esocinus_Pre)

Pseudogobio_esocinus_dat <- Pseudogobio_esocinus_dat %>%
  left_join(Pseudogobio_esocinus_Pre, by = "site") %>%
  filter(Pseudogobio_esocinus_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Pseudogobio_esocinus_glmm01 <- glmmTMB(Pseudogobio_esocinus_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Pseudogobio_esocinus_dat, na.action ="na.fail")

summary(Pseudogobio_esocinus_glmm01)



# prediction
Pseudogobio_esocinus_glmm01_global <- predict(Pseudogobio_esocinus_glmm01, newdata = Pseudogobio_esocinus_dat, type = "response")

# ROC production
Pseudogobio_esocinus_glmm01_global_roc <- roc(Pseudogobio_esocinus_dat$Pseudogobio_esocinus_PA, Pseudogobio_esocinus_glmm01_global)

# AUC
General_AUC_Pseudogobio_esocinus_global<- auc(Pseudogobio_esocinus_glmm01_global_roc)


#####cv_Pseudogobio_esocinus_global#####

Pseudogobio_esocinus_glmm01_indices <- sample(1:nrow(Pseudogobio_esocinus_dat))
Pseudogobio_esocinus_glmm01_fold_size <- floor(nrow(Pseudogobio_esocinus_dat) / 10)
Pseudogobio_esocinus_glmm01_cv_folds <- split(Pseudogobio_esocinus_glmm01_indices, 
                                        rep(1:10, each = Pseudogobio_esocinus_glmm01_fold_size, 
                                            length.out = length(Pseudogobio_esocinus_glmm01_indices)))

Pseudogobio_esocinus_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Pseudogobio_esocinus_glmm01_train_indices <- setdiff(1:nrow(Pseudogobio_esocinus_dat), Pseudogobio_esocinus_glmm01_cv_folds[[i]])
  Pseudogobio_esocinus_glmm01_validation_indices <- Pseudogobio_esocinus_glmm01_cv_folds[[i]]
  
  Pseudogobio_esocinus_glmm01_train_data <- Pseudogobio_esocinus_dat[Pseudogobio_esocinus_glmm01_train_indices, ]
  Pseudogobio_esocinus_glmm01_validation_data <- Pseudogobio_esocinus_dat[Pseudogobio_esocinus_glmm01_validation_indices, ]
  
  # GLMM training
  Pseudogobio_esocinus_glmm01_model <- glmmTMB(Pseudogobio_esocinus_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                           (1|site) + (1|season), 
                                         family = binomial(link = "logit"), data = Pseudogobio_esocinus_glmm01_train_data)
  
  # prediction
  Pseudogobio_esocinus_glmm01_pred_probs <- predict(Pseudogobio_esocinus_glmm01_model, newdata = Pseudogobio_esocinus_glmm01_validation_data, type = "response")
  
  # ROC production
  Pseudogobio_esocinus_glmm01_roc_obj <- roc(Pseudogobio_esocinus_glmm01_validation_data$Pseudogobio_esocinus_PA, Pseudogobio_esocinus_glmm01_pred_probs)
  
  # fold -  AUC
  Pseudogobio_esocinus_glmm01_auc_values[i] <- auc(Pseudogobio_esocinus_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Pseudogobio_esocinus_glmm01_auc_values[i], 4)))
}

# mean AUC
Pseudogobio_esocinus_glmm01_mean_auc <- mean(Pseudogobio_esocinus_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Pseudogobio_esocinus_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Pseudogobio_esocinus_global <- data.frame(
  Type = c("General_AUC_Pseudogobio_esocinus_global", "CV_AUC_Pseudogobio_esocinus_global"),
  AUC = c(General_AUC_Pseudogobio_esocinus_global, Pseudogobio_esocinus_glmm01_mean_auc)
)







#'## (1) Pseudogobio_esocinus model
##### Pseudogobio_esocinus group1 model #####

#'## select Group1

Pseudogobio_esocinus_dat02 <- Pseudogobio_esocinus_dat %>%
  filter(site_group == "GRP1")

#'##### glmm
#'- combination of site and season are used as a random effect ``
Pseudogobio_esocinus_glmm02 <- glmmTMB(Pseudogobio_esocinus_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Pseudogobio_esocinus_dat02, na.action ="na.fail")

summary(Pseudogobio_esocinus_glmm02)



# prediction
Pseudogobio_esocinus_glmm02_grp1 <- predict(Pseudogobio_esocinus_glmm02, newdata = Pseudogobio_esocinus_dat02, type = "response")

# ROC production
Pseudogobio_esocinus_glmm02_grp1_roc <- roc(Pseudogobio_esocinus_dat02$Pseudogobio_esocinus_PA, Pseudogobio_esocinus_glmm02_grp1)

# AUC
General_AUC_Pseudogobio_esocinus_grp1 <-auc(Pseudogobio_esocinus_glmm02_grp1_roc)


#####cv_Pseudogobio_esocinus_grp1#####

Pseudogobio_esocinus_glmm02_indices <- sample(1:nrow(Pseudogobio_esocinus_dat02))
Pseudogobio_esocinus_glmm02_fold_size <- floor(nrow(Pseudogobio_esocinus_dat02) / 10)
Pseudogobio_esocinus_glmm02_cv_folds <- split(Pseudogobio_esocinus_glmm02_indices, 
                                        rep(1:10, each = Pseudogobio_esocinus_glmm02_fold_size, 
                                            length.out = length(Pseudogobio_esocinus_glmm02_indices)))

Pseudogobio_esocinus_glmm02_auc_values <- numeric(10)

for (i in 1:10) {
  Pseudogobio_esocinus_glmm02_train_indices <- setdiff(1:nrow(Pseudogobio_esocinus_dat02), Pseudogobio_esocinus_glmm02_cv_folds[[i]])
  Pseudogobio_esocinus_glmm02_validation_indices <- Pseudogobio_esocinus_glmm02_cv_folds[[i]]
  
  Pseudogobio_esocinus_glmm02_train_dat02a <- Pseudogobio_esocinus_dat02[Pseudogobio_esocinus_glmm02_train_indices, ]
  Pseudogobio_esocinus_glmm02_validation_dat02a <- Pseudogobio_esocinus_dat02[Pseudogobio_esocinus_glmm02_validation_indices, ]
  
  # GLMM training
  Pseudogobio_esocinus_glmm02_model <- glmmTMB(Pseudogobio_esocinus_PA ~ poly(depth, 2) + poly(velocity,1) + substrate + 
                                           (1|site) + (1|season), 
                                         family = binomial(link = "logit"), data = Pseudogobio_esocinus_glmm02_train_dat02a)
  
  # prediction
  Pseudogobio_esocinus_glmm02_pred_probs <- predict(Pseudogobio_esocinus_glmm02_model, newdata = Pseudogobio_esocinus_glmm02_validation_dat02a, type = "response")
  
  # ROC production
  Pseudogobio_esocinus_glmm02_roc_obj <- roc(Pseudogobio_esocinus_glmm02_validation_dat02a$Pseudogobio_esocinus_PA, Pseudogobio_esocinus_glmm02_pred_probs)
  
  # fold -  AUC
  Pseudogobio_esocinus_glmm02_auc_values[i] <- auc(Pseudogobio_esocinus_glmm02_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Pseudogobio_esocinus_glmm02_auc_values[i], 4)))
}

# mean AUC
Pseudogobio_esocinus_glmm02_mean_auc <- mean(Pseudogobio_esocinus_glmm02_auc_values)
print(paste0("Mean AUC across folds: ", round(Pseudogobio_esocinus_glmm02_mean_auc, 4)))



#####AUC-cvAUC_grp1#####
AUC_cvAUC_Pseudogobio_esocinus_grp1 <- data.frame(
  Type = c("General_AUC_Pseudogobio_esocinus_grp1", "CV_AUC_Pseudogobio_esocinus_grp1"),
  AUC = c(General_AUC_Pseudogobio_esocinus_grp1, Pseudogobio_esocinus_glmm02_mean_auc)
)














##### Pseudogobio_esocinus group2 model #####

#'## select Group2
Pseudogobio_esocinus_dat03 <- Pseudogobio_esocinus_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Pseudogobio_esocinus_glmm03 <- glmmTMB(Pseudogobio_esocinus_PA ~ poly(depth, 1) + poly(velocity, 1) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Pseudogobio_esocinus_dat03, na.action ="na.fail")

summary(Pseudogobio_esocinus_glmm03)

# prediction
Pseudogobio_esocinus_glmm03_grp2 <- predict(Pseudogobio_esocinus_glmm03, newdata = Pseudogobio_esocinus_dat03, type = "response")

# ROC production
Pseudogobio_esocinus_glmm03_grp2_roc <- roc(Pseudogobio_esocinus_dat03$Pseudogobio_esocinus_PA, Pseudogobio_esocinus_glmm03_grp2)

# AUC
General_AUC_Pseudogobio_esocinus_grp2 <- auc(Pseudogobio_esocinus_glmm03_grp2_roc)




#####cv_Pseudogobio_esocinus_grp2#####

Pseudogobio_esocinus_glmm03_indices <- sample(1:nrow(Pseudogobio_esocinus_dat03))
Pseudogobio_esocinus_glmm03_fold_size <- floor(nrow(Pseudogobio_esocinus_dat03) / 10)
Pseudogobio_esocinus_glmm03_cv_folds <- split(Pseudogobio_esocinus_glmm03_indices, 
                                        rep(1:10, each = Pseudogobio_esocinus_glmm03_fold_size, 
                                            length.out = length(Pseudogobio_esocinus_glmm03_indices)))

Pseudogobio_esocinus_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Pseudogobio_esocinus_glmm03_train_indices <- setdiff(1:nrow(Pseudogobio_esocinus_dat03), Pseudogobio_esocinus_glmm03_cv_folds[[i]])
  Pseudogobio_esocinus_glmm03_validation_indices <- Pseudogobio_esocinus_glmm03_cv_folds[[i]]
  
  Pseudogobio_esocinus_glmm03_train_dat03 <- Pseudogobio_esocinus_dat03[Pseudogobio_esocinus_glmm03_train_indices, ]
  Pseudogobio_esocinus_glmm03_validation_dat03 <- Pseudogobio_esocinus_dat03[Pseudogobio_esocinus_glmm03_validation_indices, ]
  
  # GLMM training
  Pseudogobio_esocinus_glmm03_model <- glmmTMB(Pseudogobio_esocinus_PA ~ poly(depth, 1) + poly(velocity, 1) + substrate + 
                                           (1|site) + (1|season), 
                                         family = binomial(link = "logit"), data = Pseudogobio_esocinus_glmm03_train_dat03)
  
  # prediction
  Pseudogobio_esocinus_glmm03_pred_probs <- predict(Pseudogobio_esocinus_glmm03_model, newdata = Pseudogobio_esocinus_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Pseudogobio_esocinus_glmm03_roc_obj <- roc(Pseudogobio_esocinus_glmm03_validation_dat03$Pseudogobio_esocinus_PA, Pseudogobio_esocinus_glmm03_pred_probs)
  
  # fold -  AUC
  Pseudogobio_esocinus_glmm03_auc_values[i] <- auc(Pseudogobio_esocinus_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Pseudogobio_esocinus_glmm03_auc_values[i], 4)))
}

# mean AUC
Pseudogobio_esocinus_glmm03_mean_auc <- mean(Pseudogobio_esocinus_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Pseudogobio_esocinus_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Pseudogobio_esocinus_grp2 <- data.frame(
  Type = c("General_AUC_Pseudogobio_esocinus_grp2", "CV_AUC_Pseudogobio_esocinus_grp2"),
  AUC = c(General_AUC_Pseudogobio_esocinus_grp2, Pseudogobio_esocinus_glmm03_mean_auc)
)


#####Pseudogobio_esocinus_AUC_merge#####


Pseudogobio_esocinus_AUC <- AUC_cvAUC_Pseudogobio_esocinus_global%>%
  rbind(AUC_cvAUC_Pseudogobio_esocinus_grp1)%>%
  rbind(AUC_cvAUC_Pseudogobio_esocinus_grp2)







##### Microphysogobio_yaluensis model #####

#'## (1) Microphysogobio_yaluensis model

#'### data management
Microphysogobio_yaluensis_dat <- HSI_DB2 %>%
  mutate(Microphysogobio_yaluensis_PA = ifelse(species == "Microphysogobio_yaluensis", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Microphysogobio_yaluensis_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Microphysogobio_yaluensis_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Microphysogobio_yaluensis_PA, depth, velocity, substrate, MBSNCD)

Microphysogobio_yaluensis_Pre <- HSI_DB2 %>%
  mutate(Microphysogobio_yaluensis_PA = ifelse(species == "Microphysogobio_yaluensis", 1, 0)) %>%
  aggregate(Microphysogobio_yaluensis_PA ~ site, "sum") %>%
  mutate(Microphysogobio_yaluensis_Pre = if_else(Microphysogobio_yaluensis_PA > 0, "Yes", "No")) %>%
  select(site, Microphysogobio_yaluensis_Pre)

Microphysogobio_yaluensis_dat <- Microphysogobio_yaluensis_dat %>%
  left_join(Microphysogobio_yaluensis_Pre, by = "site") %>%
  filter(Microphysogobio_yaluensis_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Microphysogobio_yaluensis_glmm01 <- glmmTMB(Microphysogobio_yaluensis_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                              (1|site) + (1|season), family=binomial(link = "logit"), data= Microphysogobio_yaluensis_dat, na.action ="na.fail")

summary(Microphysogobio_yaluensis_glmm01)



# prediction
Microphysogobio_yaluensis_glmm01_global <- predict(Microphysogobio_yaluensis_glmm01, newdata = Microphysogobio_yaluensis_dat, type = "response")

# ROC production
Microphysogobio_yaluensis_glmm01_global_roc <- roc(Microphysogobio_yaluensis_dat$Microphysogobio_yaluensis_PA, Microphysogobio_yaluensis_glmm01_global)

# AUC
General_AUC_Microphysogobio_yaluensis_global<- auc(Microphysogobio_yaluensis_glmm01_global_roc)


#####cv_Microphysogobio_yaluensis_global#####

Microphysogobio_yaluensis_glmm01_indices <- sample(1:nrow(Microphysogobio_yaluensis_dat))
Microphysogobio_yaluensis_glmm01_fold_size <- floor(nrow(Microphysogobio_yaluensis_dat) / 10)
Microphysogobio_yaluensis_glmm01_cv_folds <- split(Microphysogobio_yaluensis_glmm01_indices, 
                                                   rep(1:10, each = Microphysogobio_yaluensis_glmm01_fold_size, 
                                                       length.out = length(Microphysogobio_yaluensis_glmm01_indices)))

Microphysogobio_yaluensis_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Microphysogobio_yaluensis_glmm01_train_indices <- setdiff(1:nrow(Microphysogobio_yaluensis_dat), Microphysogobio_yaluensis_glmm01_cv_folds[[i]])
  Microphysogobio_yaluensis_glmm01_validation_indices <- Microphysogobio_yaluensis_glmm01_cv_folds[[i]]
  
  Microphysogobio_yaluensis_glmm01_train_data <- Microphysogobio_yaluensis_dat[Microphysogobio_yaluensis_glmm01_train_indices, ]
  Microphysogobio_yaluensis_glmm01_validation_data <- Microphysogobio_yaluensis_dat[Microphysogobio_yaluensis_glmm01_validation_indices, ]
  
  Microphysogobio_yaluensis_glmm01_train_data$substrate <- factor(Microphysogobio_yaluensis_glmm01_train_data$substrate)
  
  Microphysogobio_yaluensis_glmm01_validation_data$substrate <- factor(
    Microphysogobio_yaluensis_glmm01_validation_data$substrate,
    levels = levels(Microphysogobio_yaluensis_glmm01_train_data$substrate)
  )
  
  # GLMM training
  Microphysogobio_yaluensis_glmm01_model <- glmmTMB(Microphysogobio_yaluensis_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                                      (1|site) + (1|season), 
                                                    family = binomial(link = "logit"), data = Microphysogobio_yaluensis_glmm01_train_data)
  
  # prediction
  Microphysogobio_yaluensis_glmm01_pred_probs <- predict(Microphysogobio_yaluensis_glmm01_model, newdata = Microphysogobio_yaluensis_glmm01_validation_data, type = "response")
  
  # ROC production
  Microphysogobio_yaluensis_glmm01_roc_obj <- roc(Microphysogobio_yaluensis_glmm01_validation_data$Microphysogobio_yaluensis_PA, Microphysogobio_yaluensis_glmm01_pred_probs)
  
  # fold -  AUC
  Microphysogobio_yaluensis_glmm01_auc_values[i] <- auc(Microphysogobio_yaluensis_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Microphysogobio_yaluensis_glmm01_auc_values[i], 4)))
}

# mean AUC
Microphysogobio_yaluensis_glmm01_mean_auc <- mean(Microphysogobio_yaluensis_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Microphysogobio_yaluensis_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Microphysogobio_yaluensis_global <- data.frame(
  Type = c("General_AUC_Microphysogobio_yaluensis_global", "CV_AUC_Microphysogobio_yaluensis_global"),
  AUC = c(General_AUC_Microphysogobio_yaluensis_global, Microphysogobio_yaluensis_glmm01_mean_auc)
)

















##### Microphysogobio_yaluensis group2 model #####

#'## select Group2
Microphysogobio_yaluensis_dat03 <- Microphysogobio_yaluensis_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Microphysogobio_yaluensis_glmm03 <- glmmTMB(Microphysogobio_yaluensis_PA ~ poly(depth, 1) + poly(velocity,2) + substrate + 
                                              (1|site) + (1|season) , family=binomial(link = "logit"), data= Microphysogobio_yaluensis_dat03, na.action ="na.fail")

summary(Microphysogobio_yaluensis_glmm03)

# prediction
Microphysogobio_yaluensis_glmm03_grp2 <- predict(Microphysogobio_yaluensis_glmm03, newdata = Microphysogobio_yaluensis_dat03, type = "response")

# ROC production
Microphysogobio_yaluensis_glmm03_grp2_roc <- roc(Microphysogobio_yaluensis_dat03$Microphysogobio_yaluensis_PA, Microphysogobio_yaluensis_glmm03_grp2)

# AUC
General_AUC_Microphysogobio_yaluensis_grp2 <- auc(Microphysogobio_yaluensis_glmm03_grp2_roc)




#####cv_Microphysogobio_yaluensis_grp2#####

Microphysogobio_yaluensis_glmm03_indices <- sample(1:nrow(Microphysogobio_yaluensis_dat03))
Microphysogobio_yaluensis_glmm03_fold_size <- floor(nrow(Microphysogobio_yaluensis_dat03) / 10)
Microphysogobio_yaluensis_glmm03_cv_folds <- split(Microphysogobio_yaluensis_glmm03_indices, 
                                                   rep(1:10, each = Microphysogobio_yaluensis_glmm03_fold_size, 
                                                       length.out = length(Microphysogobio_yaluensis_glmm03_indices)))

Microphysogobio_yaluensis_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Microphysogobio_yaluensis_glmm03_train_indices <- setdiff(1:nrow(Microphysogobio_yaluensis_dat03), Microphysogobio_yaluensis_glmm03_cv_folds[[i]])
  Microphysogobio_yaluensis_glmm03_validation_indices <- Microphysogobio_yaluensis_glmm03_cv_folds[[i]]
  
  Microphysogobio_yaluensis_glmm03_train_dat03 <- Microphysogobio_yaluensis_dat03[Microphysogobio_yaluensis_glmm03_train_indices, ]
  Microphysogobio_yaluensis_glmm03_validation_dat03 <- Microphysogobio_yaluensis_dat03[Microphysogobio_yaluensis_glmm03_validation_indices, ]
  
  # GLMM training
  Microphysogobio_yaluensis_glmm03_model <- glmmTMB(Microphysogobio_yaluensis_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                                      (1|site) + (1|season), 
                                                    family = binomial(link = "logit"), data = Microphysogobio_yaluensis_glmm03_train_dat03)
  
  # prediction
  Microphysogobio_yaluensis_glmm03_pred_probs <- predict(Microphysogobio_yaluensis_glmm03_model, newdata = Microphysogobio_yaluensis_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Microphysogobio_yaluensis_glmm03_roc_obj <- roc(Microphysogobio_yaluensis_glmm03_validation_dat03$Microphysogobio_yaluensis_PA, Microphysogobio_yaluensis_glmm03_pred_probs)
  
  # fold -  AUC
  Microphysogobio_yaluensis_glmm03_auc_values[i] <- auc(Microphysogobio_yaluensis_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Microphysogobio_yaluensis_glmm03_auc_values[i], 4)))
}

# mean AUC
Microphysogobio_yaluensis_glmm03_mean_auc <- mean(Microphysogobio_yaluensis_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Microphysogobio_yaluensis_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Microphysogobio_yaluensis_grp2 <- data.frame(
  Type = c("General_AUC_Microphysogobio_yaluensis_grp2", "CV_AUC_Microphysogobio_yaluensis_grp2"),
  AUC = c(General_AUC_Microphysogobio_yaluensis_grp2, Microphysogobio_yaluensis_glmm03_mean_auc)
)


#####Microphysogobio_yaluensis_AUC_merge#####


Microphysogobio_yaluensis_AUC <- AUC_cvAUC_Microphysogobio_yaluensis_global%>%
  rbind(AUC_cvAUC_Microphysogobio_yaluensis_grp2)









##### Rhinogobius_brunneus model #####

#'## (1) Rhinogobius_brunneus model

#'### data management
Rhinogobius_brunneus_dat <- HSI_DB2 %>%
  mutate(Rhinogobius_brunneus_PA = ifelse(species == "Rhinogobius_brunneus", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Rhinogobius_brunneus_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Rhinogobius_brunneus_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Rhinogobius_brunneus_PA, depth, velocity, substrate, MBSNCD)

Rhinogobius_brunneus_Pre <- HSI_DB2 %>%
  mutate(Rhinogobius_brunneus_PA = ifelse(species == "Rhinogobius_brunneus", 1, 0)) %>%
  aggregate(Rhinogobius_brunneus_PA ~ site, "sum") %>%
  mutate(Rhinogobius_brunneus_Pre = if_else(Rhinogobius_brunneus_PA > 0, "Yes", "No")) %>%
  select(site, Rhinogobius_brunneus_Pre)

Rhinogobius_brunneus_dat <- Rhinogobius_brunneus_dat %>%
  left_join(Rhinogobius_brunneus_Pre, by = "site") %>%
  filter(Rhinogobius_brunneus_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Rhinogobius_brunneus_glmm01 <- glmmTMB(Rhinogobius_brunneus_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                         (1|site) + (1|season), family=binomial(link = "logit"), data= Rhinogobius_brunneus_dat, na.action ="na.fail")

summary(Rhinogobius_brunneus_glmm01)



# prediction
Rhinogobius_brunneus_glmm01_global <- predict(Rhinogobius_brunneus_glmm01, newdata = Rhinogobius_brunneus_dat, type = "response")

# ROC production
Rhinogobius_brunneus_glmm01_global_roc <- roc(Rhinogobius_brunneus_dat$Rhinogobius_brunneus_PA, Rhinogobius_brunneus_glmm01_global)

# AUC
General_AUC_Rhinogobius_brunneus_global<- auc(Rhinogobius_brunneus_glmm01_global_roc)


#####cv_Rhinogobius_brunneus_global#####

Rhinogobius_brunneus_glmm01_indices <- sample(1:nrow(Rhinogobius_brunneus_dat))
Rhinogobius_brunneus_glmm01_fold_size <- floor(nrow(Rhinogobius_brunneus_dat) / 10)
Rhinogobius_brunneus_glmm01_cv_folds <- split(Rhinogobius_brunneus_glmm01_indices, 
                                              rep(1:10, each = Rhinogobius_brunneus_glmm01_fold_size, 
                                                  length.out = length(Rhinogobius_brunneus_glmm01_indices)))

Rhinogobius_brunneus_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Rhinogobius_brunneus_glmm01_train_indices <- setdiff(1:nrow(Rhinogobius_brunneus_dat), Rhinogobius_brunneus_glmm01_cv_folds[[i]])
  Rhinogobius_brunneus_glmm01_validation_indices <- Rhinogobius_brunneus_glmm01_cv_folds[[i]]
  
  Rhinogobius_brunneus_glmm01_train_data <- Rhinogobius_brunneus_dat[Rhinogobius_brunneus_glmm01_train_indices, ]
  Rhinogobius_brunneus_glmm01_validation_data <- Rhinogobius_brunneus_dat[Rhinogobius_brunneus_glmm01_validation_indices, ]
  
  # GLMM training
  Rhinogobius_brunneus_glmm01_model <- glmmTMB(Rhinogobius_brunneus_PA ~ poly(depth, 1) + poly(velocity, 1) + substrate + 
                                                 (1|site) + (1|season), 
                                               family = binomial(link = "logit"), data = Rhinogobius_brunneus_glmm01_train_data)
  
  # prediction
  Rhinogobius_brunneus_glmm01_pred_probs <- predict(Rhinogobius_brunneus_glmm01_model, newdata = Rhinogobius_brunneus_glmm01_validation_data, type = "response")
  
  # ROC production
  Rhinogobius_brunneus_glmm01_roc_obj <- roc(Rhinogobius_brunneus_glmm01_validation_data$Rhinogobius_brunneus_PA, Rhinogobius_brunneus_glmm01_pred_probs)
  
  # fold -  AUC
  Rhinogobius_brunneus_glmm01_auc_values[i] <- auc(Rhinogobius_brunneus_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Rhinogobius_brunneus_glmm01_auc_values[i], 4)))
}

# mean AUC
Rhinogobius_brunneus_glmm01_mean_auc <- mean(Rhinogobius_brunneus_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Rhinogobius_brunneus_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Rhinogobius_brunneus_global <- data.frame(
  Type = c("General_AUC_Rhinogobius_brunneus_global", "CV_AUC_Rhinogobius_brunneus_global"),
  AUC = c(General_AUC_Rhinogobius_brunneus_global, Rhinogobius_brunneus_glmm01_mean_auc)
)







#'## (1) Rhinogobius_brunneus model
##### Rhinogobius_brunneus group1 model #####

#'## select Group1

Rhinogobius_brunneus_dat02 <- Rhinogobius_brunneus_dat %>%
  filter(site_group == "GRP1")

#'##### glmm
#'- combination of site and season are used as a random effect ``
Rhinogobius_brunneus_glmm02 <- glmmTMB(Rhinogobius_brunneus_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                         (1|site) + (1|season), family=binomial(link = "logit"), data= Rhinogobius_brunneus_dat02, na.action ="na.fail")

summary(Rhinogobius_brunneus_glmm02)



# prediction
Rhinogobius_brunneus_glmm02_grp1 <- predict(Rhinogobius_brunneus_glmm02, newdata = Rhinogobius_brunneus_dat02, type = "response")

# ROC production
Rhinogobius_brunneus_glmm02_grp1_roc <- roc(Rhinogobius_brunneus_dat02$Rhinogobius_brunneus_PA, Rhinogobius_brunneus_glmm02_grp1)

# AUC
General_AUC_Rhinogobius_brunneus_grp1 <-auc(Rhinogobius_brunneus_glmm02_grp1_roc)


#####cv_Rhinogobius_brunneus_grp1#####

Rhinogobius_brunneus_glmm02_indices <- sample(1:nrow(Rhinogobius_brunneus_dat02))
Rhinogobius_brunneus_glmm02_fold_size <- floor(nrow(Rhinogobius_brunneus_dat02) / 10)
Rhinogobius_brunneus_glmm02_cv_folds <- split(Rhinogobius_brunneus_glmm02_indices, 
                                              rep(1:10, each = Rhinogobius_brunneus_glmm02_fold_size, 
                                                  length.out = length(Rhinogobius_brunneus_glmm02_indices)))

Rhinogobius_brunneus_glmm02_auc_values <- numeric(10)

for (i in 1:10) {
  Rhinogobius_brunneus_glmm02_train_indices <- setdiff(1:nrow(Rhinogobius_brunneus_dat02), Rhinogobius_brunneus_glmm02_cv_folds[[i]])
  Rhinogobius_brunneus_glmm02_validation_indices <- Rhinogobius_brunneus_glmm02_cv_folds[[i]]
  
  Rhinogobius_brunneus_glmm02_train_dat02a <- Rhinogobius_brunneus_dat02[Rhinogobius_brunneus_glmm02_train_indices, ]
  Rhinogobius_brunneus_glmm02_validation_dat02a <- Rhinogobius_brunneus_dat02[Rhinogobius_brunneus_glmm02_validation_indices, ]
  
  # GLMM training
  Rhinogobius_brunneus_glmm02_model <- glmmTMB(Rhinogobius_brunneus_PA ~ poly(depth, 1) + poly(velocity,2) + substrate + 
                                                 (1|site) + (1|season), 
                                               family = binomial(link = "logit"), data = Rhinogobius_brunneus_glmm02_train_dat02a)
  
  # prediction
  Rhinogobius_brunneus_glmm02_pred_probs <- predict(Rhinogobius_brunneus_glmm02_model, newdata = Rhinogobius_brunneus_glmm02_validation_dat02a, type = "response")
  
  # ROC production
  Rhinogobius_brunneus_glmm02_roc_obj <- roc(Rhinogobius_brunneus_glmm02_validation_dat02a$Rhinogobius_brunneus_PA, Rhinogobius_brunneus_glmm02_pred_probs)
  
  # fold -  AUC
  Rhinogobius_brunneus_glmm02_auc_values[i] <- auc(Rhinogobius_brunneus_glmm02_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Rhinogobius_brunneus_glmm02_auc_values[i], 4)))
}

# mean AUC
Rhinogobius_brunneus_glmm02_mean_auc <- mean(Rhinogobius_brunneus_glmm02_auc_values)
print(paste0("Mean AUC across folds: ", round(Rhinogobius_brunneus_glmm02_mean_auc, 4)))



#####AUC-cvAUC_grp1#####
AUC_cvAUC_Rhinogobius_brunneus_grp1 <- data.frame(
  Type = c("General_AUC_Rhinogobius_brunneus_grp1", "CV_AUC_Rhinogobius_brunneus_grp1"),
  AUC = c(General_AUC_Rhinogobius_brunneus_grp1, Rhinogobius_brunneus_glmm02_mean_auc)
)














##### Rhinogobius_brunneus group2 model #####

#'## select Group2
Rhinogobius_brunneus_dat03 <- Rhinogobius_brunneus_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Rhinogobius_brunneus_glmm03 <- glmmTMB(Rhinogobius_brunneus_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                         (1|site) + (1|season), family=binomial(link = "logit"), data= Rhinogobius_brunneus_dat03, na.action ="na.fail")

summary(Rhinogobius_brunneus_glmm03)

# prediction
Rhinogobius_brunneus_glmm03_grp2 <- predict(Rhinogobius_brunneus_glmm03, newdata = Rhinogobius_brunneus_dat03, type = "response")

# ROC production
Rhinogobius_brunneus_glmm03_grp2_roc <- roc(Rhinogobius_brunneus_dat03$Rhinogobius_brunneus_PA, Rhinogobius_brunneus_glmm03_grp2)

# AUC
General_AUC_Rhinogobius_brunneus_grp2 <- auc(Rhinogobius_brunneus_glmm03_grp2_roc)




#####cv_Rhinogobius_brunneus_grp2#####

Rhinogobius_brunneus_glmm03_indices <- sample(1:nrow(Rhinogobius_brunneus_dat03))
Rhinogobius_brunneus_glmm03_fold_size <- floor(nrow(Rhinogobius_brunneus_dat03) / 10)
Rhinogobius_brunneus_glmm03_cv_folds <- split(Rhinogobius_brunneus_glmm03_indices, 
                                              rep(1:10, each = Rhinogobius_brunneus_glmm03_fold_size, 
                                                  length.out = length(Rhinogobius_brunneus_glmm03_indices)))

Rhinogobius_brunneus_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Rhinogobius_brunneus_glmm03_train_indices <- setdiff(1:nrow(Rhinogobius_brunneus_dat03), Rhinogobius_brunneus_glmm03_cv_folds[[i]])
  Rhinogobius_brunneus_glmm03_validation_indices <- Rhinogobius_brunneus_glmm03_cv_folds[[i]]
  
  Rhinogobius_brunneus_glmm03_train_dat03 <- Rhinogobius_brunneus_dat03[Rhinogobius_brunneus_glmm03_train_indices, ]
  Rhinogobius_brunneus_glmm03_validation_dat03 <- Rhinogobius_brunneus_dat03[Rhinogobius_brunneus_glmm03_validation_indices, ]
  
  # GLMM training
  Rhinogobius_brunneus_glmm03_model <- glmmTMB(Rhinogobius_brunneus_PA ~ poly(depth, 1) + poly(velocity, 1) + substrate + 
                                                 (1|site) + (1|season), 
                                               family = binomial(link = "logit"), data = Rhinogobius_brunneus_glmm03_train_dat03)
  
  # prediction
  Rhinogobius_brunneus_glmm03_pred_probs <- predict(Rhinogobius_brunneus_glmm03_model, newdata = Rhinogobius_brunneus_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Rhinogobius_brunneus_glmm03_roc_obj <- roc(Rhinogobius_brunneus_glmm03_validation_dat03$Rhinogobius_brunneus_PA, Rhinogobius_brunneus_glmm03_pred_probs)
  
  # fold -  AUC
  Rhinogobius_brunneus_glmm03_auc_values[i] <- auc(Rhinogobius_brunneus_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Rhinogobius_brunneus_glmm03_auc_values[i], 4)))
}

# mean AUC
Rhinogobius_brunneus_glmm03_mean_auc <- mean(Rhinogobius_brunneus_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Rhinogobius_brunneus_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Rhinogobius_brunneus_grp2 <- data.frame(
  Type = c("General_AUC_Rhinogobius_brunneus_grp2", "CV_AUC_Rhinogobius_brunneus_grp2"),
  AUC = c(General_AUC_Rhinogobius_brunneus_grp2, Rhinogobius_brunneus_glmm03_mean_auc)
)


#####Rhinogobius_brunneus_AUC_merge#####


Rhinogobius_brunneus_AUC <- AUC_cvAUC_Rhinogobius_brunneus_global%>%
  rbind(AUC_cvAUC_Rhinogobius_brunneus_grp1)%>%
  rbind(AUC_cvAUC_Rhinogobius_brunneus_grp2)






##### Hemibarbus_labeo model #####

#'## (1) Hemibarbus_labeo model

#'### data management
Hemibarbus_labeo_dat <- HSI_DB2 %>%
  mutate(Hemibarbus_labeo_PA = ifelse(species == "Hemibarbus_labeo", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Hemibarbus_labeo_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Hemibarbus_labeo_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Hemibarbus_labeo_PA, depth, velocity, substrate, MBSNCD)

Hemibarbus_labeo_Pre <- HSI_DB2 %>%
  mutate(Hemibarbus_labeo_PA = ifelse(species == "Hemibarbus_labeo", 1, 0)) %>%
  aggregate(Hemibarbus_labeo_PA ~ site, "sum") %>%
  mutate(Hemibarbus_labeo_Pre = if_else(Hemibarbus_labeo_PA > 0, "Yes", "No")) %>%
  select(site, Hemibarbus_labeo_Pre)

Hemibarbus_labeo_dat <- Hemibarbus_labeo_dat %>%
  left_join(Hemibarbus_labeo_Pre, by = "site") %>%
  filter(Hemibarbus_labeo_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Hemibarbus_labeo_glmm01 <- glmmTMB(Hemibarbus_labeo_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                     (1|site) + (1|season), family=binomial(link = "logit"), data= Hemibarbus_labeo_dat, na.action ="na.fail")

summary(Hemibarbus_labeo_glmm01)



# prediction
Hemibarbus_labeo_glmm01_global <- predict(Hemibarbus_labeo_glmm01, newdata = Hemibarbus_labeo_dat, type = "response")

# ROC production
Hemibarbus_labeo_glmm01_global_roc <- roc(Hemibarbus_labeo_dat$Hemibarbus_labeo_PA, Hemibarbus_labeo_glmm01_global)

# AUC
General_AUC_Hemibarbus_labeo_global<- auc(Hemibarbus_labeo_glmm01_global_roc)


#####cv_Hemibarbus_labeo_global#####

Hemibarbus_labeo_glmm01_indices <- sample(1:nrow(Hemibarbus_labeo_dat))
Hemibarbus_labeo_glmm01_fold_size <- floor(nrow(Hemibarbus_labeo_dat) / 10)
Hemibarbus_labeo_glmm01_cv_folds <- split(Hemibarbus_labeo_glmm01_indices, 
                                          rep(1:10, each = Hemibarbus_labeo_glmm01_fold_size, 
                                              length.out = length(Hemibarbus_labeo_glmm01_indices)))

Hemibarbus_labeo_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Hemibarbus_labeo_glmm01_train_indices <- setdiff(1:nrow(Hemibarbus_labeo_dat), Hemibarbus_labeo_glmm01_cv_folds[[i]])
  Hemibarbus_labeo_glmm01_validation_indices <- Hemibarbus_labeo_glmm01_cv_folds[[i]]
  
  Hemibarbus_labeo_glmm01_train_data <- Hemibarbus_labeo_dat[Hemibarbus_labeo_glmm01_train_indices, ]
  Hemibarbus_labeo_glmm01_validation_data <- Hemibarbus_labeo_dat[Hemibarbus_labeo_glmm01_validation_indices, ]
  
  # GLMM training
  Hemibarbus_labeo_glmm01_model <- glmmTMB(Hemibarbus_labeo_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                             (1|site) + (1|season), 
                                           family = binomial(link = "logit"), data = Hemibarbus_labeo_glmm01_train_data)
  
  # prediction
  Hemibarbus_labeo_glmm01_pred_probs <- predict(Hemibarbus_labeo_glmm01_model, newdata = Hemibarbus_labeo_glmm01_validation_data, type = "response")
  
  # ROC production
  Hemibarbus_labeo_glmm01_roc_obj <- roc(Hemibarbus_labeo_glmm01_validation_data$Hemibarbus_labeo_PA, Hemibarbus_labeo_glmm01_pred_probs)
  
  # fold -  AUC
  Hemibarbus_labeo_glmm01_auc_values[i] <- auc(Hemibarbus_labeo_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Hemibarbus_labeo_glmm01_auc_values[i], 4)))
}

# mean AUC
Hemibarbus_labeo_glmm01_mean_auc <- mean(Hemibarbus_labeo_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Hemibarbus_labeo_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Hemibarbus_labeo_global <- data.frame(
  Type = c("General_AUC_Hemibarbus_labeo_global", "CV_AUC_Hemibarbus_labeo_global"),
  AUC = c(General_AUC_Hemibarbus_labeo_global, Hemibarbus_labeo_glmm01_mean_auc)
)







#'## (1) Hemibarbus_labeo model
##### Hemibarbus_labeo group1 model #####

#'## select Group1

Hemibarbus_labeo_dat02 <- Hemibarbus_labeo_dat %>%
  filter(site_group == "GRP1")

#'##### glmm
#'- combination of site and season are used as a random effect ``
Hemibarbus_labeo_glmm02 <- glmmTMB(Hemibarbus_labeo_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                     (1|site) + (1|season), family=binomial(link = "logit"), data= Hemibarbus_labeo_dat02, na.action ="na.fail")

summary(Hemibarbus_labeo_glmm02)



# prediction
Hemibarbus_labeo_glmm02_grp1 <- predict(Hemibarbus_labeo_glmm02, newdata = Hemibarbus_labeo_dat02, type = "response")

# ROC production
Hemibarbus_labeo_glmm02_grp1_roc <- roc(Hemibarbus_labeo_dat02$Hemibarbus_labeo_PA, Hemibarbus_labeo_glmm02_grp1)

# AUC
General_AUC_Hemibarbus_labeo_grp1 <-auc(Hemibarbus_labeo_glmm02_grp1_roc)


#####cv_Hemibarbus_labeo_grp1#####

Hemibarbus_labeo_glmm02_indices <- sample(1:nrow(Hemibarbus_labeo_dat02))
Hemibarbus_labeo_glmm02_fold_size <- floor(nrow(Hemibarbus_labeo_dat02) / 10)
Hemibarbus_labeo_glmm02_cv_folds <- split(Hemibarbus_labeo_glmm02_indices, 
                                          rep(1:10, each = Hemibarbus_labeo_glmm02_fold_size, 
                                              length.out = length(Hemibarbus_labeo_glmm02_indices)))

Hemibarbus_labeo_glmm02_auc_values <- numeric(10)

for (i in 1:10) {
  Hemibarbus_labeo_glmm02_train_indices <- setdiff(1:nrow(Hemibarbus_labeo_dat02), Hemibarbus_labeo_glmm02_cv_folds[[i]])
  Hemibarbus_labeo_glmm02_validation_indices <- Hemibarbus_labeo_glmm02_cv_folds[[i]]
  
  Hemibarbus_labeo_glmm02_train_dat02a <- Hemibarbus_labeo_dat02[Hemibarbus_labeo_glmm02_train_indices, ]
  Hemibarbus_labeo_glmm02_validation_dat02a <- Hemibarbus_labeo_dat02[Hemibarbus_labeo_glmm02_validation_indices, ]
  
  # GLMM training
  Hemibarbus_labeo_glmm02_model <- glmmTMB(Hemibarbus_labeo_PA ~ poly(depth, 2) + poly(velocity,1) + substrate + 
                                             (1|site) + (1|season), 
                                           family = binomial(link = "logit"), data = Hemibarbus_labeo_glmm02_train_dat02a)
  
  # prediction
  Hemibarbus_labeo_glmm02_pred_probs <- predict(Hemibarbus_labeo_glmm02_model, newdata = Hemibarbus_labeo_glmm02_validation_dat02a, type = "response")
  
  # ROC production
  Hemibarbus_labeo_glmm02_roc_obj <- roc(Hemibarbus_labeo_glmm02_validation_dat02a$Hemibarbus_labeo_PA, Hemibarbus_labeo_glmm02_pred_probs)
  
  # fold -  AUC
  Hemibarbus_labeo_glmm02_auc_values[i] <- auc(Hemibarbus_labeo_glmm02_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Hemibarbus_labeo_glmm02_auc_values[i], 4)))
}

# mean AUC
Hemibarbus_labeo_glmm02_mean_auc <- mean(Hemibarbus_labeo_glmm02_auc_values)
print(paste0("Mean AUC across folds: ", round(Hemibarbus_labeo_glmm02_mean_auc, 4)))



#####AUC-cvAUC_grp1#####
AUC_cvAUC_Hemibarbus_labeo_grp1 <- data.frame(
  Type = c("General_AUC_Hemibarbus_labeo_grp1", "CV_AUC_Hemibarbus_labeo_grp1"),
  AUC = c(General_AUC_Hemibarbus_labeo_grp1, Hemibarbus_labeo_glmm02_mean_auc)
)














##### Hemibarbus_labeo group2 model #####

#'## select Group2
Hemibarbus_labeo_dat03 <- Hemibarbus_labeo_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Hemibarbus_labeo_glmm03 <- glmmTMB(Hemibarbus_labeo_PA ~ poly(depth, 1) + poly(velocity, 1) + substrate + 
                                     (1|site) + (1|season), family=binomial(link = "logit"), data= Hemibarbus_labeo_dat03, na.action ="na.fail")

summary(Hemibarbus_labeo_glmm03)

# prediction
Hemibarbus_labeo_glmm03_grp2 <- predict(Hemibarbus_labeo_glmm03, newdata = Hemibarbus_labeo_dat03, type = "response")

# ROC production
Hemibarbus_labeo_glmm03_grp2_roc <- roc(Hemibarbus_labeo_dat03$Hemibarbus_labeo_PA, Hemibarbus_labeo_glmm03_grp2)

# AUC
General_AUC_Hemibarbus_labeo_grp2 <- auc(Hemibarbus_labeo_glmm03_grp2_roc)




#####cv_Hemibarbus_labeo_grp2#####

Hemibarbus_labeo_glmm03_indices <- sample(1:nrow(Hemibarbus_labeo_dat03))
Hemibarbus_labeo_glmm03_fold_size <- floor(nrow(Hemibarbus_labeo_dat03) / 10)
Hemibarbus_labeo_glmm03_cv_folds <- split(Hemibarbus_labeo_glmm03_indices, 
                                          rep(1:10, each = Hemibarbus_labeo_glmm03_fold_size, 
                                              length.out = length(Hemibarbus_labeo_glmm03_indices)))

Hemibarbus_labeo_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Hemibarbus_labeo_glmm03_train_indices <- setdiff(1:nrow(Hemibarbus_labeo_dat03), Hemibarbus_labeo_glmm03_cv_folds[[i]])
  Hemibarbus_labeo_glmm03_validation_indices <- Hemibarbus_labeo_glmm03_cv_folds[[i]]
  
  Hemibarbus_labeo_glmm03_train_dat03 <- Hemibarbus_labeo_dat03[Hemibarbus_labeo_glmm03_train_indices, ]
  Hemibarbus_labeo_glmm03_validation_dat03 <- Hemibarbus_labeo_dat03[Hemibarbus_labeo_glmm03_validation_indices, ]
  
  # GLMM training
  Hemibarbus_labeo_glmm03_model <- glmmTMB(Hemibarbus_labeo_PA ~ poly(depth, 1) + poly(velocity, 1) + substrate + 
                                             (1|site) + (1|season), 
                                           family = binomial(link = "logit"), data = Hemibarbus_labeo_glmm03_train_dat03)
  
  # prediction
  Hemibarbus_labeo_glmm03_pred_probs <- predict(Hemibarbus_labeo_glmm03_model, newdata = Hemibarbus_labeo_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Hemibarbus_labeo_glmm03_roc_obj <- roc(Hemibarbus_labeo_glmm03_validation_dat03$Hemibarbus_labeo_PA, Hemibarbus_labeo_glmm03_pred_probs)
  
  # fold -  AUC
  Hemibarbus_labeo_glmm03_auc_values[i] <- auc(Hemibarbus_labeo_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Hemibarbus_labeo_glmm03_auc_values[i], 4)))
}

# mean AUC
Hemibarbus_labeo_glmm03_mean_auc <- mean(Hemibarbus_labeo_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Hemibarbus_labeo_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Hemibarbus_labeo_grp2 <- data.frame(
  Type = c("General_AUC_Hemibarbus_labeo_grp2", "CV_AUC_Hemibarbus_labeo_grp2"),
  AUC = c(General_AUC_Hemibarbus_labeo_grp2, Hemibarbus_labeo_glmm03_mean_auc)
)


#####Hemibarbus_labeo_AUC_merge#####


Hemibarbus_labeo_AUC <- AUC_cvAUC_Hemibarbus_labeo_global%>%
  rbind(AUC_cvAUC_Hemibarbus_labeo_grp1)%>%
  rbind(AUC_cvAUC_Hemibarbus_labeo_grp2)


















##### Coreoperca_herzi model #####

#'## (1) Coreoperca_herzi model

#'### data management
Coreoperca_herzi_dat <- HSI_DB2 %>%
  mutate(Coreoperca_herzi_PA = ifelse(species == "Coreoperca_herzi", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Coreoperca_herzi_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Coreoperca_herzi_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Coreoperca_herzi_PA, depth, velocity, substrate, MBSNCD)

Coreoperca_herzi_Pre <- HSI_DB2 %>%
  mutate(Coreoperca_herzi_PA = ifelse(species == "Coreoperca_herzi", 1, 0)) %>%
  aggregate(Coreoperca_herzi_PA ~ site, "sum") %>%
  mutate(Coreoperca_herzi_Pre = if_else(Coreoperca_herzi_PA > 0, "Yes", "No")) %>%
  select(site, Coreoperca_herzi_Pre)

Coreoperca_herzi_dat <- Coreoperca_herzi_dat %>%
  left_join(Coreoperca_herzi_Pre, by = "site") %>%
  filter(Coreoperca_herzi_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Coreoperca_herzi_glmm01 <- glmmTMB(Coreoperca_herzi_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                     (1|site) + (1|season), family=binomial(link = "logit"), data= Coreoperca_herzi_dat, na.action ="na.fail")

summary(Coreoperca_herzi_glmm01)



# prediction
Coreoperca_herzi_glmm01_global <- predict(Coreoperca_herzi_glmm01, newdata = Coreoperca_herzi_dat, type = "response")

# ROC production
Coreoperca_herzi_glmm01_global_roc <- roc(Coreoperca_herzi_dat$Coreoperca_herzi_PA, Coreoperca_herzi_glmm01_global)

# AUC
General_AUC_Coreoperca_herzi_global<- auc(Coreoperca_herzi_glmm01_global_roc)


#####cv_Coreoperca_herzi_global#####

Coreoperca_herzi_glmm01_indices <- sample(1:nrow(Coreoperca_herzi_dat))
Coreoperca_herzi_glmm01_fold_size <- floor(nrow(Coreoperca_herzi_dat) / 10)
Coreoperca_herzi_glmm01_cv_folds <- split(Coreoperca_herzi_glmm01_indices, 
                                          rep(1:10, each = Coreoperca_herzi_glmm01_fold_size, 
                                              length.out = length(Coreoperca_herzi_glmm01_indices)))

Coreoperca_herzi_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Coreoperca_herzi_glmm01_train_indices <- setdiff(1:nrow(Coreoperca_herzi_dat), Coreoperca_herzi_glmm01_cv_folds[[i]])
  Coreoperca_herzi_glmm01_validation_indices <- Coreoperca_herzi_glmm01_cv_folds[[i]]
  
  Coreoperca_herzi_glmm01_train_data <- Coreoperca_herzi_dat[Coreoperca_herzi_glmm01_train_indices, ]
  Coreoperca_herzi_glmm01_validation_data <- Coreoperca_herzi_dat[Coreoperca_herzi_glmm01_validation_indices, ]
  
  Coreoperca_herzi_glmm01_train_data$substrate <- factor(Coreoperca_herzi_glmm01_train_data$substrate)
  
  Coreoperca_herzi_glmm01_validation_data$substrate <- factor(
    Coreoperca_herzi_glmm01_validation_data$substrate,
    levels = levels(Coreoperca_herzi_glmm01_train_data$substrate)
  )
  
  # GLMM training
  Coreoperca_herzi_glmm01_model <- glmmTMB(Coreoperca_herzi_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                             (1|site) + (1|season), 
                                           family = binomial(link = "logit"), data = Coreoperca_herzi_glmm01_train_data)
  
  # prediction
  Coreoperca_herzi_glmm01_pred_probs <- predict(Coreoperca_herzi_glmm01_model, newdata = Coreoperca_herzi_glmm01_validation_data, type = "response")
  
  # ROC production
  Coreoperca_herzi_glmm01_roc_obj <- roc(Coreoperca_herzi_glmm01_validation_data$Coreoperca_herzi_PA, Coreoperca_herzi_glmm01_pred_probs)
  
  # fold -  AUC
  Coreoperca_herzi_glmm01_auc_values[i] <- auc(Coreoperca_herzi_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Coreoperca_herzi_glmm01_auc_values[i], 4)))
}

# mean AUC
Coreoperca_herzi_glmm01_mean_auc <- mean(Coreoperca_herzi_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Coreoperca_herzi_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Coreoperca_herzi_global <- data.frame(
  Type = c("General_AUC_Coreoperca_herzi_global", "CV_AUC_Coreoperca_herzi_global"),
  AUC = c(General_AUC_Coreoperca_herzi_global, Coreoperca_herzi_glmm01_mean_auc)
)

















##### Coreoperca_herzi group2 model #####

#'## select Group2
Coreoperca_herzi_dat03 <- Coreoperca_herzi_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Coreoperca_herzi_glmm03 <- glmmTMB(Coreoperca_herzi_PA ~ poly(depth, 2) + poly(velocity,1) + substrate + 
                                     (1|site) + (1|season) , family=binomial(link = "logit"), data= Coreoperca_herzi_dat03, na.action ="na.fail")

summary(Coreoperca_herzi_glmm03)

# prediction
Coreoperca_herzi_glmm03_grp2 <- predict(Coreoperca_herzi_glmm03, newdata = Coreoperca_herzi_dat03, type = "response")

# ROC production
Coreoperca_herzi_glmm03_grp2_roc <- roc(Coreoperca_herzi_dat03$Coreoperca_herzi_PA, Coreoperca_herzi_glmm03_grp2)

# AUC
General_AUC_Coreoperca_herzi_grp2 <- auc(Coreoperca_herzi_glmm03_grp2_roc)




#####cv_Coreoperca_herzi_grp2#####

Coreoperca_herzi_glmm03_indices <- sample(1:nrow(Coreoperca_herzi_dat03))
Coreoperca_herzi_glmm03_fold_size <- floor(nrow(Coreoperca_herzi_dat03) / 10)
Coreoperca_herzi_glmm03_cv_folds <- split(Coreoperca_herzi_glmm03_indices, 
                                          rep(1:10, each = Coreoperca_herzi_glmm03_fold_size, 
                                              length.out = length(Coreoperca_herzi_glmm03_indices)))

Coreoperca_herzi_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Coreoperca_herzi_glmm03_train_indices <- setdiff(1:nrow(Coreoperca_herzi_dat03), Coreoperca_herzi_glmm03_cv_folds[[i]])
  Coreoperca_herzi_glmm03_validation_indices <- Coreoperca_herzi_glmm03_cv_folds[[i]]
  
  Coreoperca_herzi_glmm03_train_dat03 <- Coreoperca_herzi_dat03[Coreoperca_herzi_glmm03_train_indices, ]
  Coreoperca_herzi_glmm03_validation_dat03 <- Coreoperca_herzi_dat03[Coreoperca_herzi_glmm03_validation_indices, ]
  
  # GLMM training
  Coreoperca_herzi_glmm03_model <- glmmTMB(Coreoperca_herzi_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                             (1|site) + (1|season), 
                                           family = binomial(link = "logit"), data = Coreoperca_herzi_glmm03_train_dat03)
  
  # prediction
  Coreoperca_herzi_glmm03_pred_probs <- predict(Coreoperca_herzi_glmm03_model, newdata = Coreoperca_herzi_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Coreoperca_herzi_glmm03_roc_obj <- roc(Coreoperca_herzi_glmm03_validation_dat03$Coreoperca_herzi_PA, Coreoperca_herzi_glmm03_pred_probs)
  
  # fold -  AUC
  Coreoperca_herzi_glmm03_auc_values[i] <- auc(Coreoperca_herzi_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Coreoperca_herzi_glmm03_auc_values[i], 4)))
}

# mean AUC
Coreoperca_herzi_glmm03_mean_auc <- mean(Coreoperca_herzi_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Coreoperca_herzi_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Coreoperca_herzi_grp2 <- data.frame(
  Type = c("General_AUC_Coreoperca_herzi_grp2", "CV_AUC_Coreoperca_herzi_grp2"),
  AUC = c(General_AUC_Coreoperca_herzi_grp2, Coreoperca_herzi_glmm03_mean_auc)
)


#####Coreoperca_herzi_AUC_merge#####


Coreoperca_herzi_AUC <- AUC_cvAUC_Coreoperca_herzi_global%>%
  rbind(AUC_cvAUC_Coreoperca_herzi_grp2)














##### Squalidus_chankaensis_tsuchigae model #####

#'## (1) Squalidus_chankaensis_tsuchigae model

#'### data management
Squalidus_chankaensis_tsuchigae_dat <- HSI_DB2 %>%
  mutate(Squalidus_chankaensis_tsuchigae_PA = ifelse(species == "Squalidus_chankaensis_tsuchigae", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Squalidus_chankaensis_tsuchigae_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Squalidus_chankaensis_tsuchigae_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Squalidus_chankaensis_tsuchigae_PA, depth, velocity, substrate, MBSNCD)

Squalidus_chankaensis_tsuchigae_Pre <- HSI_DB2 %>%
  mutate(Squalidus_chankaensis_tsuchigae_PA = ifelse(species == "Squalidus_chankaensis_tsuchigae", 1, 0)) %>%
  aggregate(Squalidus_chankaensis_tsuchigae_PA ~ site, "sum") %>%
  mutate(Squalidus_chankaensis_tsuchigae_Pre = if_else(Squalidus_chankaensis_tsuchigae_PA > 0, "Yes", "No")) %>%
  select(site, Squalidus_chankaensis_tsuchigae_Pre)

Squalidus_chankaensis_tsuchigae_dat <- Squalidus_chankaensis_tsuchigae_dat %>%
  left_join(Squalidus_chankaensis_tsuchigae_Pre, by = "site") %>%
  filter(Squalidus_chankaensis_tsuchigae_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Squalidus_chankaensis_tsuchigae_glmm01 <- glmmTMB(Squalidus_chankaensis_tsuchigae_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                    (1|site) + (1|season), family=binomial(link = "logit"), data= Squalidus_chankaensis_tsuchigae_dat, na.action ="na.fail")

summary(Squalidus_chankaensis_tsuchigae_glmm01)



# prediction
Squalidus_chankaensis_tsuchigae_glmm01_global <- predict(Squalidus_chankaensis_tsuchigae_glmm01, newdata = Squalidus_chankaensis_tsuchigae_dat, type = "response")

# ROC production
Squalidus_chankaensis_tsuchigae_glmm01_global_roc <- roc(Squalidus_chankaensis_tsuchigae_dat$Squalidus_chankaensis_tsuchigae_PA, Squalidus_chankaensis_tsuchigae_glmm01_global)

# AUC
General_AUC_Squalidus_chankaensis_tsuchigae_global<- auc(Squalidus_chankaensis_tsuchigae_glmm01_global_roc)


#####cv_Squalidus_chankaensis_tsuchigae_global#####

Squalidus_chankaensis_tsuchigae_glmm01_indices <- sample(1:nrow(Squalidus_chankaensis_tsuchigae_dat))
Squalidus_chankaensis_tsuchigae_glmm01_fold_size <- floor(nrow(Squalidus_chankaensis_tsuchigae_dat) / 10)
Squalidus_chankaensis_tsuchigae_glmm01_cv_folds <- split(Squalidus_chankaensis_tsuchigae_glmm01_indices, 
                                                         rep(1:10, each = Squalidus_chankaensis_tsuchigae_glmm01_fold_size, 
                                                             length.out = length(Squalidus_chankaensis_tsuchigae_glmm01_indices)))

Squalidus_chankaensis_tsuchigae_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Squalidus_chankaensis_tsuchigae_glmm01_train_indices <- setdiff(1:nrow(Squalidus_chankaensis_tsuchigae_dat), Squalidus_chankaensis_tsuchigae_glmm01_cv_folds[[i]])
  Squalidus_chankaensis_tsuchigae_glmm01_validation_indices <- Squalidus_chankaensis_tsuchigae_glmm01_cv_folds[[i]]
  
  Squalidus_chankaensis_tsuchigae_glmm01_train_data <- Squalidus_chankaensis_tsuchigae_dat[Squalidus_chankaensis_tsuchigae_glmm01_train_indices, ]
  Squalidus_chankaensis_tsuchigae_glmm01_validation_data <- Squalidus_chankaensis_tsuchigae_dat[Squalidus_chankaensis_tsuchigae_glmm01_validation_indices, ]
  
  # GLMM training
  Squalidus_chankaensis_tsuchigae_glmm01_model <- glmmTMB(Squalidus_chankaensis_tsuchigae_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                            (1|site) + (1|season), 
                                                          family = binomial(link = "logit"), data = Squalidus_chankaensis_tsuchigae_glmm01_train_data)
  
  # prediction
  Squalidus_chankaensis_tsuchigae_glmm01_pred_probs <- predict(Squalidus_chankaensis_tsuchigae_glmm01_model, newdata = Squalidus_chankaensis_tsuchigae_glmm01_validation_data, type = "response")
  
  # ROC production
  Squalidus_chankaensis_tsuchigae_glmm01_roc_obj <- roc(Squalidus_chankaensis_tsuchigae_glmm01_validation_data$Squalidus_chankaensis_tsuchigae_PA, Squalidus_chankaensis_tsuchigae_glmm01_pred_probs)
  
  # fold -  AUC
  Squalidus_chankaensis_tsuchigae_glmm01_auc_values[i] <- auc(Squalidus_chankaensis_tsuchigae_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Squalidus_chankaensis_tsuchigae_glmm01_auc_values[i], 4)))
}

# mean AUC
Squalidus_chankaensis_tsuchigae_glmm01_mean_auc <- mean(Squalidus_chankaensis_tsuchigae_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Squalidus_chankaensis_tsuchigae_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Squalidus_chankaensis_tsuchigae_global <- data.frame(
  Type = c("General_AUC_Squalidus_chankaensis_tsuchigae_global", "CV_AUC_Squalidus_chankaensis_tsuchigae_global"),
  AUC = c(General_AUC_Squalidus_chankaensis_tsuchigae_global, Squalidus_chankaensis_tsuchigae_glmm01_mean_auc)
)







#'## (1) Squalidus_chankaensis_tsuchigae model
##### Squalidus_chankaensis_tsuchigae group1 model #####

#'## select Group1

Squalidus_chankaensis_tsuchigae_dat02 <- Squalidus_chankaensis_tsuchigae_dat %>%
  filter(site_group == "GRP1")

#'##### glmm
#'- combination of site and season are used as a random effect ``
Squalidus_chankaensis_tsuchigae_glmm02 <- glmmTMB(Squalidus_chankaensis_tsuchigae_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                    (1|site) + (1|season), family=binomial(link = "logit"), data= Squalidus_chankaensis_tsuchigae_dat02, na.action ="na.fail")

summary(Squalidus_chankaensis_tsuchigae_glmm02)



# prediction
Squalidus_chankaensis_tsuchigae_glmm02_grp1 <- predict(Squalidus_chankaensis_tsuchigae_glmm02, newdata = Squalidus_chankaensis_tsuchigae_dat02, type = "response")

# ROC production
Squalidus_chankaensis_tsuchigae_glmm02_grp1_roc <- roc(Squalidus_chankaensis_tsuchigae_dat02$Squalidus_chankaensis_tsuchigae_PA, Squalidus_chankaensis_tsuchigae_glmm02_grp1)

# AUC
General_AUC_Squalidus_chankaensis_tsuchigae_grp1 <-auc(Squalidus_chankaensis_tsuchigae_glmm02_grp1_roc)


#####cv_Squalidus_chankaensis_tsuchigae_grp1#####

Squalidus_chankaensis_tsuchigae_glmm02_indices <- sample(1:nrow(Squalidus_chankaensis_tsuchigae_dat02))
Squalidus_chankaensis_tsuchigae_glmm02_fold_size <- floor(nrow(Squalidus_chankaensis_tsuchigae_dat02) / 10)
Squalidus_chankaensis_tsuchigae_glmm02_cv_folds <- split(Squalidus_chankaensis_tsuchigae_glmm02_indices, 
                                                         rep(1:10, each = Squalidus_chankaensis_tsuchigae_glmm02_fold_size, 
                                                             length.out = length(Squalidus_chankaensis_tsuchigae_glmm02_indices)))

Squalidus_chankaensis_tsuchigae_glmm02_auc_values <- numeric(10)

for (i in 1:10) {
  Squalidus_chankaensis_tsuchigae_glmm02_train_indices <- setdiff(1:nrow(Squalidus_chankaensis_tsuchigae_dat02), Squalidus_chankaensis_tsuchigae_glmm02_cv_folds[[i]])
  Squalidus_chankaensis_tsuchigae_glmm02_validation_indices <- Squalidus_chankaensis_tsuchigae_glmm02_cv_folds[[i]]
  
  Squalidus_chankaensis_tsuchigae_glmm02_train_dat02a <- Squalidus_chankaensis_tsuchigae_dat02[Squalidus_chankaensis_tsuchigae_glmm02_train_indices, ]
  Squalidus_chankaensis_tsuchigae_glmm02_validation_dat02a <- Squalidus_chankaensis_tsuchigae_dat02[Squalidus_chankaensis_tsuchigae_glmm02_validation_indices, ]
  
  # GLMM training
  Squalidus_chankaensis_tsuchigae_glmm02_model <- glmmTMB(Squalidus_chankaensis_tsuchigae_PA ~ poly(depth, 2) + poly(velocity,2) + substrate + 
                                                            (1|site) + (1|season), 
                                                          family = binomial(link = "logit"), data = Squalidus_chankaensis_tsuchigae_glmm02_train_dat02a)
  
  # prediction
  Squalidus_chankaensis_tsuchigae_glmm02_pred_probs <- predict(Squalidus_chankaensis_tsuchigae_glmm02_model, newdata = Squalidus_chankaensis_tsuchigae_glmm02_validation_dat02a, type = "response")
  
  # ROC production
  Squalidus_chankaensis_tsuchigae_glmm02_roc_obj <- roc(Squalidus_chankaensis_tsuchigae_glmm02_validation_dat02a$Squalidus_chankaensis_tsuchigae_PA, Squalidus_chankaensis_tsuchigae_glmm02_pred_probs)
  
  # fold -  AUC
  Squalidus_chankaensis_tsuchigae_glmm02_auc_values[i] <- auc(Squalidus_chankaensis_tsuchigae_glmm02_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Squalidus_chankaensis_tsuchigae_glmm02_auc_values[i], 4)))
}

# mean AUC
Squalidus_chankaensis_tsuchigae_glmm02_mean_auc <- mean(Squalidus_chankaensis_tsuchigae_glmm02_auc_values)
print(paste0("Mean AUC across folds: ", round(Squalidus_chankaensis_tsuchigae_glmm02_mean_auc, 4)))



#####AUC-cvAUC_grp1#####
AUC_cvAUC_Squalidus_chankaensis_tsuchigae_grp1 <- data.frame(
  Type = c("General_AUC_Squalidus_chankaensis_tsuchigae_grp1", "CV_AUC_Squalidus_chankaensis_tsuchigae_grp1"),
  AUC = c(General_AUC_Squalidus_chankaensis_tsuchigae_grp1, Squalidus_chankaensis_tsuchigae_glmm02_mean_auc)
)














##### Squalidus_chankaensis_tsuchigae group2 model #####

#'## select Group2
Squalidus_chankaensis_tsuchigae_dat03 <- Squalidus_chankaensis_tsuchigae_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Squalidus_chankaensis_tsuchigae_glmm03 <- glmmTMB(Squalidus_chankaensis_tsuchigae_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                                    (1|site) + (1|season), family=binomial(link = "logit"), data= Squalidus_chankaensis_tsuchigae_dat03, na.action ="na.fail")

summary(Squalidus_chankaensis_tsuchigae_glmm03)

# prediction
Squalidus_chankaensis_tsuchigae_glmm03_grp2 <- predict(Squalidus_chankaensis_tsuchigae_glmm03, newdata = Squalidus_chankaensis_tsuchigae_dat03, type = "response")

# ROC production
Squalidus_chankaensis_tsuchigae_glmm03_grp2_roc <- roc(Squalidus_chankaensis_tsuchigae_dat03$Squalidus_chankaensis_tsuchigae_PA, Squalidus_chankaensis_tsuchigae_glmm03_grp2)

# AUC
General_AUC_Squalidus_chankaensis_tsuchigae_grp2 <- auc(Squalidus_chankaensis_tsuchigae_glmm03_grp2_roc)




#####cv_Squalidus_chankaensis_tsuchigae_grp2#####

Squalidus_chankaensis_tsuchigae_glmm03_indices <- sample(1:nrow(Squalidus_chankaensis_tsuchigae_dat03))
Squalidus_chankaensis_tsuchigae_glmm03_fold_size <- floor(nrow(Squalidus_chankaensis_tsuchigae_dat03) / 10)
Squalidus_chankaensis_tsuchigae_glmm03_cv_folds <- split(Squalidus_chankaensis_tsuchigae_glmm03_indices, 
                                                         rep(1:10, each = Squalidus_chankaensis_tsuchigae_glmm03_fold_size, 
                                                             length.out = length(Squalidus_chankaensis_tsuchigae_glmm03_indices)))

Squalidus_chankaensis_tsuchigae_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Squalidus_chankaensis_tsuchigae_glmm03_train_indices <- setdiff(1:nrow(Squalidus_chankaensis_tsuchigae_dat03), Squalidus_chankaensis_tsuchigae_glmm03_cv_folds[[i]])
  Squalidus_chankaensis_tsuchigae_glmm03_validation_indices <- Squalidus_chankaensis_tsuchigae_glmm03_cv_folds[[i]]
  
  Squalidus_chankaensis_tsuchigae_glmm03_train_dat03 <- Squalidus_chankaensis_tsuchigae_dat03[Squalidus_chankaensis_tsuchigae_glmm03_train_indices, ]
  Squalidus_chankaensis_tsuchigae_glmm03_validation_dat03 <- Squalidus_chankaensis_tsuchigae_dat03[Squalidus_chankaensis_tsuchigae_glmm03_validation_indices, ]
  
  # GLMM training
  Squalidus_chankaensis_tsuchigae_glmm03_model <- glmmTMB(Squalidus_chankaensis_tsuchigae_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                                            (1|site) + (1|season), 
                                                          family = binomial(link = "logit"), data = Squalidus_chankaensis_tsuchigae_glmm03_train_dat03)
  
  # prediction
  Squalidus_chankaensis_tsuchigae_glmm03_pred_probs <- predict(Squalidus_chankaensis_tsuchigae_glmm03_model, newdata = Squalidus_chankaensis_tsuchigae_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Squalidus_chankaensis_tsuchigae_glmm03_roc_obj <- roc(Squalidus_chankaensis_tsuchigae_glmm03_validation_dat03$Squalidus_chankaensis_tsuchigae_PA, Squalidus_chankaensis_tsuchigae_glmm03_pred_probs)
  
  # fold -  AUC
  Squalidus_chankaensis_tsuchigae_glmm03_auc_values[i] <- auc(Squalidus_chankaensis_tsuchigae_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Squalidus_chankaensis_tsuchigae_glmm03_auc_values[i], 4)))
}

# mean AUC
Squalidus_chankaensis_tsuchigae_glmm03_mean_auc <- mean(Squalidus_chankaensis_tsuchigae_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Squalidus_chankaensis_tsuchigae_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Squalidus_chankaensis_tsuchigae_grp2 <- data.frame(
  Type = c("General_AUC_Squalidus_chankaensis_tsuchigae_grp2", "CV_AUC_Squalidus_chankaensis_tsuchigae_grp2"),
  AUC = c(General_AUC_Squalidus_chankaensis_tsuchigae_grp2, Squalidus_chankaensis_tsuchigae_glmm03_mean_auc)
)


#####Squalidus_chankaensis_tsuchigae_AUC_merge#####


Squalidus_chankaensis_tsuchigae_AUC <- AUC_cvAUC_Squalidus_chankaensis_tsuchigae_global%>%
  rbind(AUC_cvAUC_Squalidus_chankaensis_tsuchigae_grp1)%>%
  rbind(AUC_cvAUC_Squalidus_chankaensis_tsuchigae_grp2)













##### Microphysogobio_longidorsalis model #####

#'## (1) Microphysogobio_longidorsalis model

#'### data management
Microphysogobio_longidorsalis_dat <- HSI_DB2 %>%
  mutate(Microphysogobio_longidorsalis_PA = ifelse(species == "Microphysogobio_longidorsalis", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Microphysogobio_longidorsalis_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Microphysogobio_longidorsalis_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Microphysogobio_longidorsalis_PA, depth, velocity, substrate, MBSNCD)

Microphysogobio_longidorsalis_Pre <- HSI_DB2 %>%
  mutate(Microphysogobio_longidorsalis_PA = ifelse(species == "Microphysogobio_longidorsalis", 1, 0)) %>%
  aggregate(Microphysogobio_longidorsalis_PA ~ site, "sum") %>%
  mutate(Microphysogobio_longidorsalis_Pre = if_else(Microphysogobio_longidorsalis_PA > 0, "Yes", "No")) %>%
  select(site, Microphysogobio_longidorsalis_Pre)

Microphysogobio_longidorsalis_dat <- Microphysogobio_longidorsalis_dat %>%
  left_join(Microphysogobio_longidorsalis_Pre, by = "site") %>%
  filter(Microphysogobio_longidorsalis_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Microphysogobio_longidorsalis_glmm01 <- glmmTMB(Microphysogobio_longidorsalis_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                  (1|site) + (1|season), family=binomial(link = "logit"), data= Microphysogobio_longidorsalis_dat, na.action ="na.fail")

summary(Microphysogobio_longidorsalis_glmm01)



# prediction
Microphysogobio_longidorsalis_glmm01_global <- predict(Microphysogobio_longidorsalis_glmm01, newdata = Microphysogobio_longidorsalis_dat, type = "response")

# ROC production
Microphysogobio_longidorsalis_glmm01_global_roc <- roc(Microphysogobio_longidorsalis_dat$Microphysogobio_longidorsalis_PA, Microphysogobio_longidorsalis_glmm01_global)

# AUC
General_AUC_Microphysogobio_longidorsalis_global<- auc(Microphysogobio_longidorsalis_glmm01_global_roc)


#####cv_Microphysogobio_longidorsalis_global#####

Microphysogobio_longidorsalis_glmm01_indices <- sample(1:nrow(Microphysogobio_longidorsalis_dat))
Microphysogobio_longidorsalis_glmm01_fold_size <- floor(nrow(Microphysogobio_longidorsalis_dat) / 10)
Microphysogobio_longidorsalis_glmm01_cv_folds <- split(Microphysogobio_longidorsalis_glmm01_indices, 
                                                       rep(1:10, each = Microphysogobio_longidorsalis_glmm01_fold_size, 
                                                           length.out = length(Microphysogobio_longidorsalis_glmm01_indices)))

Microphysogobio_longidorsalis_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Microphysogobio_longidorsalis_glmm01_train_indices <- setdiff(1:nrow(Microphysogobio_longidorsalis_dat), Microphysogobio_longidorsalis_glmm01_cv_folds[[i]])
  Microphysogobio_longidorsalis_glmm01_validation_indices <- Microphysogobio_longidorsalis_glmm01_cv_folds[[i]]
  
  Microphysogobio_longidorsalis_glmm01_train_data <- Microphysogobio_longidorsalis_dat[Microphysogobio_longidorsalis_glmm01_train_indices, ]
  Microphysogobio_longidorsalis_glmm01_validation_data <- Microphysogobio_longidorsalis_dat[Microphysogobio_longidorsalis_glmm01_validation_indices, ]
  
  Microphysogobio_longidorsalis_glmm01_train_data$substrate <- factor(Microphysogobio_longidorsalis_glmm01_train_data$substrate)
  
  Microphysogobio_longidorsalis_glmm01_validation_data$substrate <- factor(
    Microphysogobio_longidorsalis_glmm01_validation_data$substrate,
    levels = levels(Microphysogobio_longidorsalis_glmm01_train_data$substrate)
  )
  
  # GLMM training
  Microphysogobio_longidorsalis_glmm01_model <- glmmTMB(Microphysogobio_longidorsalis_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                          (1|site) + (1|season), 
                                                        family = binomial(link = "logit"), data = Microphysogobio_longidorsalis_glmm01_train_data)
  
  # prediction
  Microphysogobio_longidorsalis_glmm01_pred_probs <- predict(Microphysogobio_longidorsalis_glmm01_model, newdata = Microphysogobio_longidorsalis_glmm01_validation_data, type = "response")
  
  # ROC production
  Microphysogobio_longidorsalis_glmm01_roc_obj <- roc(Microphysogobio_longidorsalis_glmm01_validation_data$Microphysogobio_longidorsalis_PA, Microphysogobio_longidorsalis_glmm01_pred_probs)
  
  # fold -  AUC
  Microphysogobio_longidorsalis_glmm01_auc_values[i] <- auc(Microphysogobio_longidorsalis_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Microphysogobio_longidorsalis_glmm01_auc_values[i], 4)))
}

# mean AUC
Microphysogobio_longidorsalis_glmm01_mean_auc <- mean(Microphysogobio_longidorsalis_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Microphysogobio_longidorsalis_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Microphysogobio_longidorsalis_global <- data.frame(
  Type = c("General_AUC_Microphysogobio_longidorsalis_global", "CV_AUC_Microphysogobio_longidorsalis_global"),
  AUC = c(General_AUC_Microphysogobio_longidorsalis_global, Microphysogobio_longidorsalis_glmm01_mean_auc)
)

















##### Microphysogobio_longidorsalis group2 model #####

#'## select Group2
Microphysogobio_longidorsalis_dat03 <- Microphysogobio_longidorsalis_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Microphysogobio_longidorsalis_glmm03 <- glmmTMB(Microphysogobio_longidorsalis_PA ~ poly(depth, 2) + poly(velocity,2) + substrate + 
                                                  (1|site) + (1|season) , family=binomial(link = "logit"), data= Microphysogobio_longidorsalis_dat03, na.action ="na.fail")

summary(Microphysogobio_longidorsalis_glmm03)

# prediction
Microphysogobio_longidorsalis_glmm03_grp2 <- predict(Microphysogobio_longidorsalis_glmm03, newdata = Microphysogobio_longidorsalis_dat03, type = "response")

# ROC production
Microphysogobio_longidorsalis_glmm03_grp2_roc <- roc(Microphysogobio_longidorsalis_dat03$Microphysogobio_longidorsalis_PA, Microphysogobio_longidorsalis_glmm03_grp2)

# AUC
General_AUC_Microphysogobio_longidorsalis_grp2 <- auc(Microphysogobio_longidorsalis_glmm03_grp2_roc)




#####cv_Microphysogobio_longidorsalis_grp2#####

Microphysogobio_longidorsalis_glmm03_indices <- sample(1:nrow(Microphysogobio_longidorsalis_dat03))
Microphysogobio_longidorsalis_glmm03_fold_size <- floor(nrow(Microphysogobio_longidorsalis_dat03) / 10)
Microphysogobio_longidorsalis_glmm03_cv_folds <- split(Microphysogobio_longidorsalis_glmm03_indices, 
                                                       rep(1:10, each = Microphysogobio_longidorsalis_glmm03_fold_size, 
                                                           length.out = length(Microphysogobio_longidorsalis_glmm03_indices)))

Microphysogobio_longidorsalis_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Microphysogobio_longidorsalis_glmm03_train_indices <- setdiff(1:nrow(Microphysogobio_longidorsalis_dat03), Microphysogobio_longidorsalis_glmm03_cv_folds[[i]])
  Microphysogobio_longidorsalis_glmm03_validation_indices <- Microphysogobio_longidorsalis_glmm03_cv_folds[[i]]
  
  Microphysogobio_longidorsalis_glmm03_train_dat03 <- Microphysogobio_longidorsalis_dat03[Microphysogobio_longidorsalis_glmm03_train_indices, ]
  Microphysogobio_longidorsalis_glmm03_validation_dat03 <- Microphysogobio_longidorsalis_dat03[Microphysogobio_longidorsalis_glmm03_validation_indices, ]
  
  # GLMM training
  Microphysogobio_longidorsalis_glmm03_model <- glmmTMB(Microphysogobio_longidorsalis_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                          (1|site) + (1|season), 
                                                        family = binomial(link = "logit"), data = Microphysogobio_longidorsalis_glmm03_train_dat03)
  
  # prediction
  Microphysogobio_longidorsalis_glmm03_pred_probs <- predict(Microphysogobio_longidorsalis_glmm03_model, newdata = Microphysogobio_longidorsalis_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Microphysogobio_longidorsalis_glmm03_roc_obj <- roc(Microphysogobio_longidorsalis_glmm03_validation_dat03$Microphysogobio_longidorsalis_PA, Microphysogobio_longidorsalis_glmm03_pred_probs)
  
  # fold -  AUC
  Microphysogobio_longidorsalis_glmm03_auc_values[i] <- auc(Microphysogobio_longidorsalis_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Microphysogobio_longidorsalis_glmm03_auc_values[i], 4)))
}

# mean AUC
Microphysogobio_longidorsalis_glmm03_mean_auc <- mean(Microphysogobio_longidorsalis_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Microphysogobio_longidorsalis_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Microphysogobio_longidorsalis_grp2 <- data.frame(
  Type = c("General_AUC_Microphysogobio_longidorsalis_grp2", "CV_AUC_Microphysogobio_longidorsalis_grp2"),
  AUC = c(General_AUC_Microphysogobio_longidorsalis_grp2, Microphysogobio_longidorsalis_glmm03_mean_auc)
)


#####Microphysogobio_longidorsalis_AUC_merge#####


Microphysogobio_longidorsalis_AUC <- AUC_cvAUC_Microphysogobio_longidorsalis_global%>%
  rbind(AUC_cvAUC_Microphysogobio_longidorsalis_grp2)














##### Hemibarbus_longirostris model #####

#'## (1) Hemibarbus_longirostris model

#'### data management
Hemibarbus_longirostris_dat <- HSI_DB2 %>%
  mutate(Hemibarbus_longirostris_PA = ifelse(species == "Hemibarbus_longirostris", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Hemibarbus_longirostris_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Hemibarbus_longirostris_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Hemibarbus_longirostris_PA, depth, velocity, substrate, MBSNCD)

Hemibarbus_longirostris_Pre <- HSI_DB2 %>%
  mutate(Hemibarbus_longirostris_PA = ifelse(species == "Hemibarbus_longirostris", 1, 0)) %>%
  aggregate(Hemibarbus_longirostris_PA ~ site, "sum") %>%
  mutate(Hemibarbus_longirostris_Pre = if_else(Hemibarbus_longirostris_PA > 0, "Yes", "No")) %>%
  select(site, Hemibarbus_longirostris_Pre)

Hemibarbus_longirostris_dat <- Hemibarbus_longirostris_dat %>%
  left_join(Hemibarbus_longirostris_Pre, by = "site") %>%
  filter(Hemibarbus_longirostris_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Hemibarbus_longirostris_glmm01 <- glmmTMB(Hemibarbus_longirostris_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                            (1|site) + (1|season), family=binomial(link = "logit"), data= Hemibarbus_longirostris_dat, na.action ="na.fail")

summary(Hemibarbus_longirostris_glmm01)



# prediction
Hemibarbus_longirostris_glmm01_global <- predict(Hemibarbus_longirostris_glmm01, newdata = Hemibarbus_longirostris_dat, type = "response")

# ROC production
Hemibarbus_longirostris_glmm01_global_roc <- roc(Hemibarbus_longirostris_dat$Hemibarbus_longirostris_PA, Hemibarbus_longirostris_glmm01_global)

# AUC
General_AUC_Hemibarbus_longirostris_global<- auc(Hemibarbus_longirostris_glmm01_global_roc)


#####cv_Hemibarbus_longirostris_global#####

Hemibarbus_longirostris_glmm01_indices <- sample(1:nrow(Hemibarbus_longirostris_dat))
Hemibarbus_longirostris_glmm01_fold_size <- floor(nrow(Hemibarbus_longirostris_dat) / 10)
Hemibarbus_longirostris_glmm01_cv_folds <- split(Hemibarbus_longirostris_glmm01_indices, 
                                                 rep(1:10, each = Hemibarbus_longirostris_glmm01_fold_size, 
                                                     length.out = length(Hemibarbus_longirostris_glmm01_indices)))

Hemibarbus_longirostris_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Hemibarbus_longirostris_glmm01_train_indices <- setdiff(1:nrow(Hemibarbus_longirostris_dat), Hemibarbus_longirostris_glmm01_cv_folds[[i]])
  Hemibarbus_longirostris_glmm01_validation_indices <- Hemibarbus_longirostris_glmm01_cv_folds[[i]]
  
  Hemibarbus_longirostris_glmm01_train_data <- Hemibarbus_longirostris_dat[Hemibarbus_longirostris_glmm01_train_indices, ]
  Hemibarbus_longirostris_glmm01_validation_data <- Hemibarbus_longirostris_dat[Hemibarbus_longirostris_glmm01_validation_indices, ]
  
  # GLMM training
  Hemibarbus_longirostris_glmm01_model <- glmmTMB(Hemibarbus_longirostris_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                    (1|site) + (1|season), 
                                                  family = binomial(link = "logit"), data = Hemibarbus_longirostris_glmm01_train_data)
  
  # prediction
  Hemibarbus_longirostris_glmm01_pred_probs <- predict(Hemibarbus_longirostris_glmm01_model, newdata = Hemibarbus_longirostris_glmm01_validation_data, type = "response")
  
  # ROC production
  Hemibarbus_longirostris_glmm01_roc_obj <- roc(Hemibarbus_longirostris_glmm01_validation_data$Hemibarbus_longirostris_PA, Hemibarbus_longirostris_glmm01_pred_probs)
  
  # fold -  AUC
  Hemibarbus_longirostris_glmm01_auc_values[i] <- auc(Hemibarbus_longirostris_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Hemibarbus_longirostris_glmm01_auc_values[i], 4)))
}

# mean AUC
Hemibarbus_longirostris_glmm01_mean_auc <- mean(Hemibarbus_longirostris_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Hemibarbus_longirostris_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Hemibarbus_longirostris_global <- data.frame(
  Type = c("General_AUC_Hemibarbus_longirostris_global", "CV_AUC_Hemibarbus_longirostris_global"),
  AUC = c(General_AUC_Hemibarbus_longirostris_global, Hemibarbus_longirostris_glmm01_mean_auc)
)







#'## (1) Hemibarbus_longirostris model
##### Hemibarbus_longirostris group1 model #####

#'## select Group1

Hemibarbus_longirostris_dat02 <- Hemibarbus_longirostris_dat %>%
  filter(site_group == "GRP1")

#'##### glmm
#'- combination of site and season are used as a random effect ``
Hemibarbus_longirostris_glmm02 <- glmmTMB(Hemibarbus_longirostris_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                            (1|site) + (1|season), family=binomial(link = "logit"), data= Hemibarbus_longirostris_dat02, na.action ="na.fail")

summary(Hemibarbus_longirostris_glmm02)



# prediction
Hemibarbus_longirostris_glmm02_grp1 <- predict(Hemibarbus_longirostris_glmm02, newdata = Hemibarbus_longirostris_dat02, type = "response")

# ROC production
Hemibarbus_longirostris_glmm02_grp1_roc <- roc(Hemibarbus_longirostris_dat02$Hemibarbus_longirostris_PA, Hemibarbus_longirostris_glmm02_grp1)

# AUC
General_AUC_Hemibarbus_longirostris_grp1 <-auc(Hemibarbus_longirostris_glmm02_grp1_roc)


#####cv_Hemibarbus_longirostris_grp1#####

Hemibarbus_longirostris_glmm02_indices <- sample(1:nrow(Hemibarbus_longirostris_dat02))
Hemibarbus_longirostris_glmm02_fold_size <- floor(nrow(Hemibarbus_longirostris_dat02) / 10)
Hemibarbus_longirostris_glmm02_cv_folds <- split(Hemibarbus_longirostris_glmm02_indices, 
                                                 rep(1:10, each = Hemibarbus_longirostris_glmm02_fold_size, 
                                                     length.out = length(Hemibarbus_longirostris_glmm02_indices)))

Hemibarbus_longirostris_glmm02_auc_values <- numeric(10)

for (i in 1:10) {
  Hemibarbus_longirostris_glmm02_train_indices <- setdiff(1:nrow(Hemibarbus_longirostris_dat02), Hemibarbus_longirostris_glmm02_cv_folds[[i]])
  Hemibarbus_longirostris_glmm02_validation_indices <- Hemibarbus_longirostris_glmm02_cv_folds[[i]]
  
  Hemibarbus_longirostris_glmm02_train_dat02a <- Hemibarbus_longirostris_dat02[Hemibarbus_longirostris_glmm02_train_indices, ]
  Hemibarbus_longirostris_glmm02_validation_dat02a <- Hemibarbus_longirostris_dat02[Hemibarbus_longirostris_glmm02_validation_indices, ]
  
  # GLMM training
  Hemibarbus_longirostris_glmm02_model <- glmmTMB(Hemibarbus_longirostris_PA ~ poly(depth, 1) + poly(velocity,2) + substrate + 
                                                    (1|site) + (1|season), 
                                                  family = binomial(link = "logit"), data = Hemibarbus_longirostris_glmm02_train_dat02a)
  
  # prediction
  Hemibarbus_longirostris_glmm02_pred_probs <- predict(Hemibarbus_longirostris_glmm02_model, newdata = Hemibarbus_longirostris_glmm02_validation_dat02a, type = "response")
  
  # ROC production
  Hemibarbus_longirostris_glmm02_roc_obj <- roc(Hemibarbus_longirostris_glmm02_validation_dat02a$Hemibarbus_longirostris_PA, Hemibarbus_longirostris_glmm02_pred_probs)
  
  # fold -  AUC
  Hemibarbus_longirostris_glmm02_auc_values[i] <- auc(Hemibarbus_longirostris_glmm02_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Hemibarbus_longirostris_glmm02_auc_values[i], 4)))
}

# mean AUC
Hemibarbus_longirostris_glmm02_mean_auc <- mean(Hemibarbus_longirostris_glmm02_auc_values)
print(paste0("Mean AUC across folds: ", round(Hemibarbus_longirostris_glmm02_mean_auc, 4)))



#####AUC-cvAUC_grp1#####
AUC_cvAUC_Hemibarbus_longirostris_grp1 <- data.frame(
  Type = c("General_AUC_Hemibarbus_longirostris_grp1", "CV_AUC_Hemibarbus_longirostris_grp1"),
  AUC = c(General_AUC_Hemibarbus_longirostris_grp1, Hemibarbus_longirostris_glmm02_mean_auc)
)














##### Hemibarbus_longirostris group2 model #####

#'## select Group2
Hemibarbus_longirostris_dat03 <- Hemibarbus_longirostris_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Hemibarbus_longirostris_glmm03 <- glmmTMB(Hemibarbus_longirostris_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                            (1|site) + (1|season), family=binomial(link = "logit"), data= Hemibarbus_longirostris_dat03, na.action ="na.fail")

summary(Hemibarbus_longirostris_glmm03)

# prediction
Hemibarbus_longirostris_glmm03_grp2 <- predict(Hemibarbus_longirostris_glmm03, newdata = Hemibarbus_longirostris_dat03, type = "response")

# ROC production
Hemibarbus_longirostris_glmm03_grp2_roc <- roc(Hemibarbus_longirostris_dat03$Hemibarbus_longirostris_PA, Hemibarbus_longirostris_glmm03_grp2)

# AUC
General_AUC_Hemibarbus_longirostris_grp2 <- auc(Hemibarbus_longirostris_glmm03_grp2_roc)




#####cv_Hemibarbus_longirostris_grp2#####

Hemibarbus_longirostris_glmm03_indices <- sample(1:nrow(Hemibarbus_longirostris_dat03))
Hemibarbus_longirostris_glmm03_fold_size <- floor(nrow(Hemibarbus_longirostris_dat03) / 10)
Hemibarbus_longirostris_glmm03_cv_folds <- split(Hemibarbus_longirostris_glmm03_indices, 
                                                 rep(1:10, each = Hemibarbus_longirostris_glmm03_fold_size, 
                                                     length.out = length(Hemibarbus_longirostris_glmm03_indices)))

Hemibarbus_longirostris_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Hemibarbus_longirostris_glmm03_train_indices <- setdiff(1:nrow(Hemibarbus_longirostris_dat03), Hemibarbus_longirostris_glmm03_cv_folds[[i]])
  Hemibarbus_longirostris_glmm03_validation_indices <- Hemibarbus_longirostris_glmm03_cv_folds[[i]]
  
  Hemibarbus_longirostris_glmm03_train_dat03 <- Hemibarbus_longirostris_dat03[Hemibarbus_longirostris_glmm03_train_indices, ]
  Hemibarbus_longirostris_glmm03_validation_dat03 <- Hemibarbus_longirostris_dat03[Hemibarbus_longirostris_glmm03_validation_indices, ]
  
  # GLMM training
  Hemibarbus_longirostris_glmm03_model <- glmmTMB(Hemibarbus_longirostris_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                    (1|site) + (1|season), 
                                                  family = binomial(link = "logit"), data = Hemibarbus_longirostris_glmm03_train_dat03)
  
  # prediction
  Hemibarbus_longirostris_glmm03_pred_probs <- predict(Hemibarbus_longirostris_glmm03_model, newdata = Hemibarbus_longirostris_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Hemibarbus_longirostris_glmm03_roc_obj <- roc(Hemibarbus_longirostris_glmm03_validation_dat03$Hemibarbus_longirostris_PA, Hemibarbus_longirostris_glmm03_pred_probs)
  
  # fold -  AUC
  Hemibarbus_longirostris_glmm03_auc_values[i] <- auc(Hemibarbus_longirostris_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Hemibarbus_longirostris_glmm03_auc_values[i], 4)))
}

# mean AUC
Hemibarbus_longirostris_glmm03_mean_auc <- mean(Hemibarbus_longirostris_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Hemibarbus_longirostris_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Hemibarbus_longirostris_grp2 <- data.frame(
  Type = c("General_AUC_Hemibarbus_longirostris_grp2", "CV_AUC_Hemibarbus_longirostris_grp2"),
  AUC = c(General_AUC_Hemibarbus_longirostris_grp2, Hemibarbus_longirostris_glmm03_mean_auc)
)


#####Hemibarbus_longirostris_AUC_merge#####


Hemibarbus_longirostris_AUC <- AUC_cvAUC_Hemibarbus_longirostris_global%>%
  rbind(AUC_cvAUC_Hemibarbus_longirostris_grp1)%>%
  rbind(AUC_cvAUC_Hemibarbus_longirostris_grp2)





















##### Sarcocheilichthys_variegatus_wakiyae model #####

#'## (1) Sarcocheilichthys_variegatus_wakiyae model

#'### data management
Sarcocheilichthys_variegatus_wakiyae_dat <- HSI_DB2 %>%
  mutate(Sarcocheilichthys_variegatus_wakiyae_PA = ifelse(species == "Sarcocheilichthys_variegatus_wakiyae", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Sarcocheilichthys_variegatus_wakiyae_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Sarcocheilichthys_variegatus_wakiyae_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Sarcocheilichthys_variegatus_wakiyae_PA, depth, velocity, substrate, MBSNCD)

Sarcocheilichthys_variegatus_wakiyae_Pre <- HSI_DB2 %>%
  mutate(Sarcocheilichthys_variegatus_wakiyae_PA = ifelse(species == "Sarcocheilichthys_variegatus_wakiyae", 1, 0)) %>%
  aggregate(Sarcocheilichthys_variegatus_wakiyae_PA ~ site, "sum") %>%
  mutate(Sarcocheilichthys_variegatus_wakiyae_Pre = if_else(Sarcocheilichthys_variegatus_wakiyae_PA > 0, "Yes", "No")) %>%
  select(site, Sarcocheilichthys_variegatus_wakiyae_Pre)

Sarcocheilichthys_variegatus_wakiyae_dat <- Sarcocheilichthys_variegatus_wakiyae_dat %>%
  left_join(Sarcocheilichthys_variegatus_wakiyae_Pre, by = "site") %>%
  filter(Sarcocheilichthys_variegatus_wakiyae_Pre == "Yes")


#'### Global Model
#'##### glmm
#'- combination of site and season are used as a random effect 
Sarcocheilichthys_variegatus_wakiyae_glmm01 <- glmmTMB(Sarcocheilichthys_variegatus_wakiyae_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                         (1|site) + (1|season), family=binomial(link = "logit"), data= Sarcocheilichthys_variegatus_wakiyae_dat, na.action ="na.fail")

summary(Sarcocheilichthys_variegatus_wakiyae_glmm01)



# prediction
Sarcocheilichthys_variegatus_wakiyae_glmm01_global <- predict(Sarcocheilichthys_variegatus_wakiyae_glmm01, newdata = Sarcocheilichthys_variegatus_wakiyae_dat, type = "response")

# ROC production
Sarcocheilichthys_variegatus_wakiyae_glmm01_global_roc <- roc(Sarcocheilichthys_variegatus_wakiyae_dat$Sarcocheilichthys_variegatus_wakiyae_PA, Sarcocheilichthys_variegatus_wakiyae_glmm01_global)

# AUC
General_AUC_Sarcocheilichthys_variegatus_wakiyae_global<- auc(Sarcocheilichthys_variegatus_wakiyae_glmm01_global_roc)


#####cv_Sarcocheilichthys_variegatus_wakiyae_global#####

Sarcocheilichthys_variegatus_wakiyae_glmm01_indices <- sample(1:nrow(Sarcocheilichthys_variegatus_wakiyae_dat))
Sarcocheilichthys_variegatus_wakiyae_glmm01_fold_size <- floor(nrow(Sarcocheilichthys_variegatus_wakiyae_dat) / 10)
Sarcocheilichthys_variegatus_wakiyae_glmm01_cv_folds <- split(Sarcocheilichthys_variegatus_wakiyae_glmm01_indices, 
                                                              rep(1:10, each = Sarcocheilichthys_variegatus_wakiyae_glmm01_fold_size, 
                                                                  length.out = length(Sarcocheilichthys_variegatus_wakiyae_glmm01_indices)))

Sarcocheilichthys_variegatus_wakiyae_glmm01_auc_values <- numeric(10)

for (i in 1:10) {
  Sarcocheilichthys_variegatus_wakiyae_glmm01_train_indices <- setdiff(1:nrow(Sarcocheilichthys_variegatus_wakiyae_dat), Sarcocheilichthys_variegatus_wakiyae_glmm01_cv_folds[[i]])
  Sarcocheilichthys_variegatus_wakiyae_glmm01_validation_indices <- Sarcocheilichthys_variegatus_wakiyae_glmm01_cv_folds[[i]]
  
  Sarcocheilichthys_variegatus_wakiyae_glmm01_train_data <- Sarcocheilichthys_variegatus_wakiyae_dat[Sarcocheilichthys_variegatus_wakiyae_glmm01_train_indices, ]
  Sarcocheilichthys_variegatus_wakiyae_glmm01_validation_data <- Sarcocheilichthys_variegatus_wakiyae_dat[Sarcocheilichthys_variegatus_wakiyae_glmm01_validation_indices, ]
  
  Sarcocheilichthys_variegatus_wakiyae_glmm01_train_data$substrate <- factor(Sarcocheilichthys_variegatus_wakiyae_glmm01_train_data$substrate)
  
  Sarcocheilichthys_variegatus_wakiyae_glmm01_validation_data$substrate <- factor(
    Sarcocheilichthys_variegatus_wakiyae_glmm01_validation_data$substrate,
    levels = levels(Sarcocheilichthys_variegatus_wakiyae_glmm01_train_data$substrate)
  )
  
  # GLMM training
  Sarcocheilichthys_variegatus_wakiyae_glmm01_model <- glmmTMB(Sarcocheilichthys_variegatus_wakiyae_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                                 (1|site) + (1|season), 
                                                               family = binomial(link = "logit"), data = Sarcocheilichthys_variegatus_wakiyae_glmm01_train_data)
  
  # prediction
  Sarcocheilichthys_variegatus_wakiyae_glmm01_pred_probs <- predict(Sarcocheilichthys_variegatus_wakiyae_glmm01_model, newdata = Sarcocheilichthys_variegatus_wakiyae_glmm01_validation_data, type = "response")
  
  # ROC production
  Sarcocheilichthys_variegatus_wakiyae_glmm01_roc_obj <- roc(Sarcocheilichthys_variegatus_wakiyae_glmm01_validation_data$Sarcocheilichthys_variegatus_wakiyae_PA, Sarcocheilichthys_variegatus_wakiyae_glmm01_pred_probs)
  
  # fold -  AUC
  Sarcocheilichthys_variegatus_wakiyae_glmm01_auc_values[i] <- auc(Sarcocheilichthys_variegatus_wakiyae_glmm01_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Sarcocheilichthys_variegatus_wakiyae_glmm01_auc_values[i], 4)))
}

# mean AUC
Sarcocheilichthys_variegatus_wakiyae_glmm01_mean_auc <- mean(Sarcocheilichthys_variegatus_wakiyae_glmm01_auc_values)
print(paste0("Mean AUC across folds: ", round(Sarcocheilichthys_variegatus_wakiyae_glmm01_mean_auc, 4)))



#####AUC-cvAUC_global#####



AUC_cvAUC_Sarcocheilichthys_variegatus_wakiyae_global <- data.frame(
  Type = c("General_AUC_Sarcocheilichthys_variegatus_wakiyae_global", "CV_AUC_Sarcocheilichthys_variegatus_wakiyae_global"),
  AUC = c(General_AUC_Sarcocheilichthys_variegatus_wakiyae_global, Sarcocheilichthys_variegatus_wakiyae_glmm01_mean_auc)
)

















##### Sarcocheilichthys_variegatus_wakiyae group2 model #####

#'## select Group2
Sarcocheilichthys_variegatus_wakiyae_dat03 <- Sarcocheilichthys_variegatus_wakiyae_dat %>%
  filter(site_group == "GRP2")


#'##### glmm
#'- combination of site and season are used as a random effect 
Sarcocheilichthys_variegatus_wakiyae_glmm03 <- glmmTMB(Sarcocheilichthys_variegatus_wakiyae_PA ~ poly(depth, 2) + poly(velocity,2) + substrate + 
                                                         (1|site) + (1|season) , family=binomial(link = "logit"), data= Sarcocheilichthys_variegatus_wakiyae_dat03, na.action ="na.fail")

summary(Sarcocheilichthys_variegatus_wakiyae_glmm03)

# prediction
Sarcocheilichthys_variegatus_wakiyae_glmm03_grp2 <- predict(Sarcocheilichthys_variegatus_wakiyae_glmm03, newdata = Sarcocheilichthys_variegatus_wakiyae_dat03, type = "response")

# ROC production
Sarcocheilichthys_variegatus_wakiyae_glmm03_grp2_roc <- roc(Sarcocheilichthys_variegatus_wakiyae_dat03$Sarcocheilichthys_variegatus_wakiyae_PA, Sarcocheilichthys_variegatus_wakiyae_glmm03_grp2)

# AUC
General_AUC_Sarcocheilichthys_variegatus_wakiyae_grp2 <- auc(Sarcocheilichthys_variegatus_wakiyae_glmm03_grp2_roc)




#####cv_Sarcocheilichthys_variegatus_wakiyae_grp2#####

Sarcocheilichthys_variegatus_wakiyae_glmm03_indices <- sample(1:nrow(Sarcocheilichthys_variegatus_wakiyae_dat03))
Sarcocheilichthys_variegatus_wakiyae_glmm03_fold_size <- floor(nrow(Sarcocheilichthys_variegatus_wakiyae_dat03) / 10)
Sarcocheilichthys_variegatus_wakiyae_glmm03_cv_folds <- split(Sarcocheilichthys_variegatus_wakiyae_glmm03_indices, 
                                                              rep(1:10, each = Sarcocheilichthys_variegatus_wakiyae_glmm03_fold_size, 
                                                                  length.out = length(Sarcocheilichthys_variegatus_wakiyae_glmm03_indices)))

Sarcocheilichthys_variegatus_wakiyae_glmm03_auc_values <- numeric(10)

for (i in 1:10) {
  Sarcocheilichthys_variegatus_wakiyae_glmm03_train_indices <- setdiff(1:nrow(Sarcocheilichthys_variegatus_wakiyae_dat03), Sarcocheilichthys_variegatus_wakiyae_glmm03_cv_folds[[i]])
  Sarcocheilichthys_variegatus_wakiyae_glmm03_validation_indices <- Sarcocheilichthys_variegatus_wakiyae_glmm03_cv_folds[[i]]
  
  Sarcocheilichthys_variegatus_wakiyae_glmm03_train_dat03 <- Sarcocheilichthys_variegatus_wakiyae_dat03[Sarcocheilichthys_variegatus_wakiyae_glmm03_train_indices, ]
  Sarcocheilichthys_variegatus_wakiyae_glmm03_validation_dat03 <- Sarcocheilichthys_variegatus_wakiyae_dat03[Sarcocheilichthys_variegatus_wakiyae_glmm03_validation_indices, ]
  
  # GLMM training
  Sarcocheilichthys_variegatus_wakiyae_glmm03_model <- glmmTMB(Sarcocheilichthys_variegatus_wakiyae_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                                 (1|site) + (1|season), 
                                                               family = binomial(link = "logit"), data = Sarcocheilichthys_variegatus_wakiyae_glmm03_train_dat03)
  
  # prediction
  Sarcocheilichthys_variegatus_wakiyae_glmm03_pred_probs <- predict(Sarcocheilichthys_variegatus_wakiyae_glmm03_model, newdata = Sarcocheilichthys_variegatus_wakiyae_glmm03_validation_dat03, type = "response")
  
  # ROC production
  Sarcocheilichthys_variegatus_wakiyae_glmm03_roc_obj <- roc(Sarcocheilichthys_variegatus_wakiyae_glmm03_validation_dat03$Sarcocheilichthys_variegatus_wakiyae_PA, Sarcocheilichthys_variegatus_wakiyae_glmm03_pred_probs)
  
  # fold -  AUC
  Sarcocheilichthys_variegatus_wakiyae_glmm03_auc_values[i] <- auc(Sarcocheilichthys_variegatus_wakiyae_glmm03_roc_obj)
  
  print(paste0("Fold ", i, " AUC: ", round(Sarcocheilichthys_variegatus_wakiyae_glmm03_auc_values[i], 4)))
}

# mean AUC
Sarcocheilichthys_variegatus_wakiyae_glmm03_mean_auc <- mean(Sarcocheilichthys_variegatus_wakiyae_glmm03_auc_values)
print(paste0("Mean AUC across folds: ", round(Sarcocheilichthys_variegatus_wakiyae_glmm03_mean_auc, 4)))



#####AUC-cvAUC_grp2#####



AUC_cvAUC_Sarcocheilichthys_variegatus_wakiyae_grp2 <- data.frame(
  Type = c("General_AUC_Sarcocheilichthys_variegatus_wakiyae_grp2", "CV_AUC_Sarcocheilichthys_variegatus_wakiyae_grp2"),
  AUC = c(General_AUC_Sarcocheilichthys_variegatus_wakiyae_grp2, Sarcocheilichthys_variegatus_wakiyae_glmm03_mean_auc)
)


#####Sarcocheilichthys_variegatus_wakiyae_AUC_merge#####


Sarcocheilichthys_variegatus_wakiyae_AUC <- AUC_cvAUC_Sarcocheilichthys_variegatus_wakiyae_global%>%
  rbind(AUC_cvAUC_Sarcocheilichthys_variegatus_wakiyae_grp2)
















##### all_species AUC-CV Auc Value #####


Result_CV_General_AUC <- rbind(Zacco_platypus_AUC, Zacco_koreanus_AUC, Pungtungia_herzi_AUC, Pseudopungtungia_nigra_AUC, Coreoleuciscus_splendidus_AUC, 
Opsariichthys_uncirostris_amurensis_AUC, Pseudogobio_esocinus_AUC, Microphysogobio_yaluensis_AUC, Rhinogobius_brunneus_AUC, 
Hemibarbus_labeo_AUC, Coreoperca_herzi_AUC, Squalidus_chankaensis_tsuchigae_AUC, Microphysogobio_longidorsalis_AUC, 
Hemibarbus_longirostris_AUC, Sarcocheilichthys_variegatus_wakiyae_AUC)

Result_CV_General_AUC

write.csv(Result_CV_General_AUC, "./data_outcome/AUC_result.csv")


