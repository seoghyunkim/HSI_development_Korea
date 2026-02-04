

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
               performance)



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




length(unique(HSI_DB2$site) ) # site number n = 195
length(unique(HSI_DB2$basin) ) # basin number n = 6
list(unique(HSI_DB2$basin))
length(unique(HSI_DB2$species) ) # species number n = 98
sum(HSI_DB2$individual_number) # n = 88534 individuals




##### Zacco_platypus model #####

#'## (2) Zacco_platypus model

#'### data management
Zacco_platypus_dat <- HSI_DB2 %>%
  mutate(Zacco_platypus_PA = ifelse(species == "Zacco_platypus", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Zacco_platypus_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Zacco_platypus_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Zacco_platypus_PA, depth, velocity, substrate, MBSNCD) 

head(Zacco_platypus_dat)


#'## remove site where species did not occur
Zacco_platypus_Pre <- HSI_DB2 %>%
  mutate(Zacco_platypus_PA = ifelse(species == "Zacco_platypus", 1, 0)) %>%
  aggregate(Zacco_platypus_PA ~ site, "sum") %>%
  mutate(Zacco_platypus_Pre = if_else(Zacco_platypus_PA > 0, "Yes", "No")) %>%
  select(site, Zacco_platypus_Pre)

Zacco_platypus_dat <- Zacco_platypus_dat %>%
  left_join(Zacco_platypus_Pre, by = "site") %>%
  filter(Zacco_platypus_Pre == "Yes") %>%
  as_tibble()





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

#'##### glm
Zacco_platypus_glm01 <- glm(Zacco_platypus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                            family=binomial(link = "logit"), data= Zacco_platypus_dat, na.action ="na.fail")

summary(Zacco_platypus_glm01)

vif(Zacco_platypus_glm01)



#'##### glmm
#'- combination of site and season are used as a random effect 
Zacco_platypus_glmm01 <- glmmTMB(Zacco_platypus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Zacco_platypus_dat, na.action ="na.fail")

summary(Zacco_platypus_glmm01)


#'- Model validation
r.squaredGLMM(Zacco_platypus_glmm01)
cor.test(predict(Zacco_platypus_glmm01, type = "response", allow.new.levels=TRUE), Zacco_platypus_dat$Zacco_platypus_PA)





#'- Zacco_platypus depth
Zacco_platypus_glmm01_dep_fit <- effect_plot(Zacco_platypus_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Zacco_platypus_glmm01_dep_fit_PA <- as.data.frame(Zacco_platypus_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Zacco_platypus_PA"]])
Zacco_platypus_glmm01_dep_fit_dep <- as.data.frame(Zacco_platypus_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Zacco_platypus_dep_res <- cbind(Zacco_platypus_glmm01_dep_fit_PA, Zacco_platypus_glmm01_dep_fit_dep)
colnames(Zacco_platypus_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Zacco_platypus_dep_res <- Zacco_platypus_dep_res %>%
  mutate(HSI_dep = scale_values(prob))
#'- plot
Zacco_platypus_depth_fig <- ggplot(Zacco_platypus_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Zacco_platypus depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_platypus_depth_fig

#'- Zacco_platypus velocity
Zacco_platypus_glmm01_vel_fit <- effect_plot(Zacco_platypus_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Zacco_platypus_glmm01_vel_fit_PA <- as.data.frame(Zacco_platypus_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Zacco_platypus_PA"]])
Zacco_platypus_glmm01_vel_fit_vel <- as.data.frame(Zacco_platypus_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Zacco_platypus_vel_res <- cbind(Zacco_platypus_glmm01_vel_fit_PA, Zacco_platypus_glmm01_vel_fit_vel)
colnames(Zacco_platypus_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Zacco_platypus_vel_res <- Zacco_platypus_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Zacco_platypus_velocity_fig <- ggplot(Zacco_platypus_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Zacco_platypus Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_platypus_velocity_fig

#'- Zacco_platypus substrate
Zacco_platypus_glmm01_sub_fit <- effect_plot(Zacco_platypus_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Zacco_platypus_glmm01_sub_fit_PA <- as.data.frame(Zacco_platypus_glmm01_sub_fit[["data"]][["Zacco_platypus_PA"]])
Zacco_platypus_glmm01_sub_fit_sub <- as.data.frame(Zacco_platypus_glmm01_sub_fit[["data"]][["substrate"]])
Zacco_platypus_sub_res <- cbind(Zacco_platypus_glmm01_sub_fit_PA, Zacco_platypus_glmm01_sub_fit_sub)
colnames(Zacco_platypus_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Zacco_platypus_sub_res <- Zacco_platypus_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Zacco_platypus_sub_res$sub <- factor(Zacco_platypus_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                            'pebble', 'cobble', 'boulder'))

#'- plot
Zacco_platypus_substrate_fig <- ggplot(Zacco_platypus_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Zacco_platypus Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_platypus_substrate_fig



#'## (1) Zacco_platypus model
##### Zacco_platypus group1 model #####

#'## select Group1

Zacco_platypus_dat02 <- Zacco_platypus_dat %>%
  filter(site_group == "GRP1")

#'##### glm
Zacco_platypus_glm02 <- glm(Zacco_platypus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                            family=binomial(link = "logit"), data= Zacco_platypus_dat02, na.action ="na.fail")

summary(Zacco_platypus_glm02)

vif(Zacco_platypus_glm02)


#'##### glmm
#'- combination of site and season are used as a random effect ``
Zacco_platypus_glmm02 <- glmmTMB(Zacco_platypus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Zacco_platypus_dat02, na.action ="na.fail")

summary(Zacco_platypus_glmm02)



#'- Model validation
r.squaredGLMM(Zacco_platypus_glmm02)
cor.test(predict(Zacco_platypus_glmm02, type = "response", allow.new.levels=TRUE), Zacco_platypus_dat02$Zacco_platypus_PA)

#'- Zacco_platypus depth
Zacco_platypus_glmm02_dep_fit <- effect_plot(Zacco_platypus_glmm02, pred = depth, interval = TRUE, plot.points = FALSE)
Zacco_platypus_glmm02_dep_fit_PA <- as.data.frame(Zacco_platypus_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["Zacco_platypus_PA"]])
Zacco_platypus_glmm02_dep_fit_depth <- as.data.frame(Zacco_platypus_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Zacco_platypus_dep_res_grp1 <- cbind(Zacco_platypus_glmm02_dep_fit_PA, Zacco_platypus_glmm02_dep_fit_depth)
colnames(Zacco_platypus_dep_res_grp1) <- c("prob", "dep")

#'- normalize 0 to 1
Zacco_platypus_dep_res_grp1 <- Zacco_platypus_dep_res_grp1 %>%
  mutate(HSI_dep_grp1 = scale_values(prob))

#'- plot
Zacco_platypus_depth_fig1 <- ggplot(Zacco_platypus_dep_res_grp1, aes(x = dep, y = HSI_dep_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Zacco_platypus depth HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_platypus_depth_fig1

#'- Zacco_platypus velocity
Zacco_platypus_glmm02_vel_fit <- effect_plot(Zacco_platypus_glmm02, pred = velocity, interval = TRUE, plot.points = FALSE)
Zacco_platypus_glmm02_vel_fit_PA <- as.data.frame(Zacco_platypus_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["Zacco_platypus_PA"]])
Zacco_platypus_glmm02_vel_fit_vel <- as.data.frame(Zacco_platypus_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Zacco_platypus_vel_res_grp1 <- cbind(Zacco_platypus_glmm02_vel_fit_PA, Zacco_platypus_glmm02_vel_fit_vel)
colnames(Zacco_platypus_vel_res_grp1) <- c("prob", "vel")

#'- normalize 0 to 1
Zacco_platypus_vel_res_grp1 <- Zacco_platypus_vel_res_grp1 %>%
  mutate(HSI_vel_grp1 = scale_values(prob))

#'- plot
Zacco_platypus_velocity_fig1 <- ggplot(Zacco_platypus_vel_res_grp1, aes(x = vel, y = HSI_vel_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Zacco_platypus Velocity HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_platypus_velocity_fig1

#'- Zacco_platypus substrate
Zacco_platypus_glmm02_sub_fit <- effect_plot(Zacco_platypus_glmm02, pred = substrate, interval = FALSE, plot.points = FALSE)
Zacco_platypus_glmm02_sub_fit_PA <- as.data.frame(Zacco_platypus_glmm02_sub_fit[["data"]][["Zacco_platypus_PA"]])
Zacco_platypus_glmm02_sub_fit_sub <- as.data.frame(Zacco_platypus_glmm02_sub_fit[["data"]][["substrate"]])
Zacco_platypus_sub_res_grp1 <- cbind(Zacco_platypus_glmm02_sub_fit_PA, Zacco_platypus_glmm02_sub_fit_sub)
colnames(Zacco_platypus_sub_res_grp1) <- c("prob", "sub")

#'- normalize 0 to 1
Zacco_platypus_sub_res_grp1 <- Zacco_platypus_sub_res_grp1 %>%
  mutate(HSI_sub_grp1 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Zacco_platypus_sub_res_grp1$sub <- factor(Zacco_platypus_sub_res_grp1$sub, levels = c('silt', 'sand', 'gravel',
                                                                                      'pebble', 'cobble', 'boulder'))

#'- plot
Zacco_platypus_substrate_fig1 <- ggplot(Zacco_platypus_sub_res_grp1, aes(x = sub, y = HSI_sub_grp1)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Zacco_platypus Substrate HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_platypus_substrate_fig1



##### Zacco_platypus group2 model #####

#'## select Group2
Zacco_platypus_dat03 <- Zacco_platypus_dat %>%
  filter(site_group == "GRP2")

#'##### glm
Zacco_platypus_glm03 <- glm(Zacco_platypus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                            family=binomial(link = "logit"), data= Zacco_platypus_dat03, na.action ="na.fail")

summary(Zacco_platypus_glm03)

vif(Zacco_platypus_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect 
Zacco_platypus_glmm03 <- glmmTMB(Zacco_platypus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Zacco_platypus_dat03, na.action ="na.fail")

summary(Zacco_platypus_glmm03)




#'- Model validation
r.squaredGLMM(Zacco_platypus_glmm03)
cor.test(predict(Zacco_platypus_glmm03, type = "response", allow.new.levels=TRUE), Zacco_platypus_dat03$Zacco_platypus_PA)

#'- Zacco_platypus depth
Zacco_platypus_glmm03_dep_fit <- effect_plot(Zacco_platypus_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Zacco_platypus_glmm03_dep_fit_PA <- as.data.frame(Zacco_platypus_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Zacco_platypus_PA"]])
Zacco_platypus_glmm03_dep_fit_depth <- as.data.frame(Zacco_platypus_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Zacco_platypus_dep_res_grp2 <- cbind(Zacco_platypus_glmm03_dep_fit_PA, Zacco_platypus_glmm03_dep_fit_depth)
colnames(Zacco_platypus_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Zacco_platypus_dep_res_grp2 <- Zacco_platypus_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Zacco_platypus_depth_fig2 <- ggplot(Zacco_platypus_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Zacco_platypus depth HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_platypus_depth_fig2

#'- Zacco_platypus velocity
Zacco_platypus_glmm03_vel_fit <- effect_plot(Zacco_platypus_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Zacco_platypus_glmm03_vel_fit_PA <- as.data.frame(Zacco_platypus_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Zacco_platypus_PA"]])
Zacco_platypus_glmm03_vel_fit_vel <- as.data.frame(Zacco_platypus_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Zacco_platypus_vel_res_grp2 <- cbind(Zacco_platypus_glmm03_vel_fit_PA, Zacco_platypus_glmm03_vel_fit_vel)
colnames(Zacco_platypus_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Zacco_platypus_vel_res_grp2 <- Zacco_platypus_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Zacco_platypus_velocity_fig2 <- ggplot(Zacco_platypus_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Zacco_platypus Velocity HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_platypus_velocity_fig2

#'- Zacco_platypus substrate
Zacco_platypus_glmm03_sub_fit <- effect_plot(Zacco_platypus_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Zacco_platypus_glmm03_sub_fit_PA <- as.data.frame(Zacco_platypus_glmm03_sub_fit[["data"]][["Zacco_platypus_PA"]])
Zacco_platypus_glmm03_sub_fit_sub <- as.data.frame(Zacco_platypus_glmm03_sub_fit[["data"]][["substrate"]])
Zacco_platypus_sub_res_grp2 <- cbind(Zacco_platypus_glmm03_sub_fit_PA, Zacco_platypus_glmm03_sub_fit_sub)
colnames(Zacco_platypus_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Zacco_platypus_sub_res_grp2 <- Zacco_platypus_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Zacco_platypus_sub_res_grp2$sub <- factor(Zacco_platypus_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                      'pebble', 'cobble', 'boulder'))

#'- plot
Zacco_platypus_substrate_fig2 <- ggplot(Zacco_platypus_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Zacco platypus Substrate HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_platypus_substrate_fig2







#'## Combined figures

#'- merge mulitple effect plot
Zacco_platypus_dep_res_1 <- Zacco_platypus_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)

Zacco_platypus_dep_res_grp1_1 <- Zacco_platypus_dep_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp1)

Zacco_platypus_dep_res_grp2_1 <- Zacco_platypus_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)


Zacco_platypus_dep_res_comb <- rbind(Zacco_platypus_dep_res_1, Zacco_platypus_dep_res_grp1_1,
                                     Zacco_platypus_dep_res_grp2_1)


#'- depth plot
Zacco_platypus_depth_fig4 <- ggplot(Zacco_platypus_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)", limits = c(0,2.0) ,breaks = seq(0, 3, 0.5)) +     
  labs(title= "Zacco platypus Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Zacco_platypus_depth_fig4 <-Zacco_platypus_depth_fig4+scale_color_manual(values=c("black","blue3", "red"))
Zacco_platypus_depth_fig4




#'- merge mulitple effect plot
Zacco_platypus_vel_res_1 <- Zacco_platypus_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)

Zacco_platypus_vel_res_grp1_1 <- Zacco_platypus_vel_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp1)

Zacco_platypus_vel_res_grp2_1 <- Zacco_platypus_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)




Zacco_platypus_vel_res_comb <- rbind(Zacco_platypus_vel_res_1, Zacco_platypus_vel_res_grp1_1,
                                     Zacco_platypus_vel_res_grp2_1)


#'- velocity plot
Zacco_platypus_velocity_fig4 <- ggplot(Zacco_platypus_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)", limits = c(0, 3.1), breaks = seq(0, 3, 0.5)) +     
  labs(title= "Zacco platypus Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Zacco_platypus_velocity_fig4 <- Zacco_platypus_velocity_fig4 + scale_color_manual(values=c("black","blue3", "red"))
Zacco_platypus_velocity_fig4



#'- merge mulitple effect plot
Zacco_platypus_sub_res_1 <- Zacco_platypus_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)

Zacco_platypus_sub_res_grp1_1 <- Zacco_platypus_sub_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp1)

Zacco_platypus_sub_res_grp2_1 <- Zacco_platypus_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)





Zacco_platypus_sub_res_comb <- rbind(Zacco_platypus_sub_res_1, Zacco_platypus_sub_res_grp1_1,
                                     Zacco_platypus_sub_res_grp2_1)

#'- substrate plot
Zacco_platypus_substrate_fig4 <- ggplot(Zacco_platypus_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Zacco platypus Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Zacco_platypus_substrate_fig4 <- Zacco_platypus_substrate_fig4 + scale_fill_manual(values=c("black","blue3", "red"))
Zacco_platypus_substrate_fig4


#'## combine figures

ggarrange(Zacco_platypus_depth_fig4, Zacco_platypus_velocity_fig4, Zacco_platypus_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 






#------------------------------------------------------------------------------------------------------------------------------------------#





##### Zacco_koreanus model #####

#'## (1) Zacco_koreanus model

#'### data management
Zacco_koreanus_dat <- HSI_DB2 %>%
  mutate(Zacco_koreanus_PA = ifelse(species == "Zacco_koreanus", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Zacco_koreanus_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Zacco_koreanus_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Zacco_koreanus_PA, depth, velocity, substrate, MBSNCD) 

head(Zacco_koreanus_dat)


#'## remove site where species did not occur
Zacco_koreanus_Pre <- HSI_DB2 %>%
  mutate(Zacco_koreanus_PA = ifelse(species == "Zacco_koreanus", 1, 0)) %>%
  aggregate(Zacco_koreanus_PA ~ site, "sum") %>%
  mutate(Zacco_koreanus_Pre = if_else(Zacco_koreanus_PA > 0, "Yes", "No")) %>%
  select(site, Zacco_koreanus_Pre)

Zacco_koreanus_dat <- Zacco_koreanus_dat %>%
  left_join(Zacco_koreanus_Pre, by = "site") %>%
  filter(Zacco_koreanus_Pre == "Yes") %>%
  as_tibble()




#'### Global Model

#'##### glm
Zacco_koreanus_glm01 <- glm(Zacco_koreanus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                            family=binomial(link = "logit"), data= Zacco_koreanus_dat, na.action ="na.fail")

summary(Zacco_koreanus_glm01)


vif(Zacco_koreanus_glm01)

#'##### glmm
#'- combination of site and season are used as a random effect 
Zacco_koreanus_glmm01 <- glmmTMB(Zacco_koreanus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Zacco_koreanus_dat, na.action ="na.fail")

summary(Zacco_koreanus_glmm01)






#'- Model validation
r.squaredGLMM(Zacco_koreanus_glmm01)
cor.test(predict(Zacco_koreanus_glmm01, type = "response", allow.new.levels=TRUE), Zacco_koreanus_dat$Zacco_koreanus_PA)

#'- Zacco_koreanus depth
Zacco_koreanus_glmm01_dep_fit <- effect_plot(Zacco_koreanus_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Zacco_koreanus_glmm01_dep_fit_PA <- as.data.frame(Zacco_koreanus_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Zacco_koreanus_PA"]])
Zacco_koreanus_glmm01_dep_fit_dep <- as.data.frame(Zacco_koreanus_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Zacco_koreanus_dep_res <- cbind(Zacco_koreanus_glmm01_dep_fit_PA, Zacco_koreanus_glmm01_dep_fit_dep)
colnames(Zacco_koreanus_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Zacco_koreanus_dep_res <- Zacco_koreanus_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Zacco_koreanus_depth_fig <- ggplot(Zacco_koreanus_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Zacco_koreanus depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_koreanus_depth_fig

#'- Zacco_koreanus velocity
Zacco_koreanus_glmm01_vel_fit <- effect_plot(Zacco_koreanus_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Zacco_koreanus_glmm01_vel_fit_PA <- as.data.frame(Zacco_koreanus_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Zacco_koreanus_PA"]])
Zacco_koreanus_glmm01_vel_fit_vel <- as.data.frame(Zacco_koreanus_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Zacco_koreanus_vel_res <- cbind(Zacco_koreanus_glmm01_vel_fit_PA, Zacco_koreanus_glmm01_vel_fit_vel)
colnames(Zacco_koreanus_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Zacco_koreanus_vel_res <- Zacco_koreanus_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Zacco_koreanus_velocity_fig <- ggplot(Zacco_koreanus_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Zacco_koreanus Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_koreanus_velocity_fig

#'- Zacco_koreanus substrate
Zacco_koreanus_glmm01_sub_fit <- effect_plot(Zacco_koreanus_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Zacco_koreanus_glmm01_sub_fit_PA <- as.data.frame(Zacco_koreanus_glmm01_sub_fit[["data"]][["Zacco_koreanus_PA"]])
Zacco_koreanus_glmm01_sub_fit_sub <- as.data.frame(Zacco_koreanus_glmm01_sub_fit[["data"]][["substrate"]])
Zacco_koreanus_sub_res <- cbind(Zacco_koreanus_glmm01_sub_fit_PA, Zacco_koreanus_glmm01_sub_fit_sub)
colnames(Zacco_koreanus_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Zacco_koreanus_sub_res <- Zacco_koreanus_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Zacco_koreanus_sub_res$sub <- factor(Zacco_koreanus_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                            'pebble', 'cobble', 'boulder'))

#'- plot
Zacco_koreanus_substrate_fig <- ggplot(Zacco_koreanus_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Zacco_koreanus Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_koreanus_substrate_fig




##### Zacco_koreanus group2 model #####

#'## select Group2
Zacco_koreanus_dat03 <- Zacco_koreanus_dat %>%
  filter(site_group == "GRP2")

#'##### glm
Zacco_koreanus_glm03 <- glm(Zacco_koreanus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                            family=binomial(link = "logit"), data= Zacco_koreanus_dat03, na.action ="na.fail")

summary(Zacco_koreanus_glm03)
vif(Zacco_koreanus_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect 
Zacco_koreanus_glmm03 <- glmmTMB(Zacco_koreanus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                   (1|site) + (1|season), family=binomial(link = "logit"), data= Zacco_koreanus_dat03, na.action ="na.fail")

summary(Zacco_koreanus_glmm03)






#'- Model validation
r.squaredGLMM(Zacco_koreanus_glmm03)
cor.test(predict(Zacco_koreanus_glmm03, type = "response", allow.new.levels=TRUE), Zacco_koreanus_dat03$Zacco_koreanus_PA)

#'- Zacco_koreanus depth
Zacco_koreanus_glmm03_dep_fit <- effect_plot(Zacco_koreanus_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Zacco_koreanus_glmm03_dep_fit_PA <- as.data.frame(Zacco_koreanus_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Zacco_koreanus_PA"]])
Zacco_koreanus_glmm03_dep_fit_depth <- as.data.frame(Zacco_koreanus_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Zacco_koreanus_dep_res_grp2 <- cbind(Zacco_koreanus_glmm03_dep_fit_PA, Zacco_koreanus_glmm03_dep_fit_depth)
colnames(Zacco_koreanus_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Zacco_koreanus_dep_res_grp2 <- Zacco_koreanus_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Zacco_koreanus_depth_fig2 <- ggplot(Zacco_koreanus_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2.5, 0.5)) +     
  labs(title= "Zacco_koreanus depth HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_koreanus_depth_fig2

#'- Zacco_koreanus velocity
Zacco_koreanus_glmm03_vel_fit <- effect_plot(Zacco_koreanus_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Zacco_koreanus_glmm03_vel_fit_PA <- as.data.frame(Zacco_koreanus_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Zacco_koreanus_PA"]])
Zacco_koreanus_glmm03_vel_fit_vel <- as.data.frame(Zacco_koreanus_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Zacco_koreanus_vel_res_grp2 <- cbind(Zacco_koreanus_glmm03_vel_fit_PA, Zacco_koreanus_glmm03_vel_fit_vel)
colnames(Zacco_koreanus_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Zacco_koreanus_vel_res_grp2 <- Zacco_koreanus_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Zacco_koreanus_velocity_fig2 <- ggplot(Zacco_koreanus_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Zacco_koreanus Velocity HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_koreanus_velocity_fig2

#'- Zacco_koreanus substrate
Zacco_koreanus_glmm03_sub_fit <- effect_plot(Zacco_koreanus_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Zacco_koreanus_glmm03_sub_fit_PA <- as.data.frame(Zacco_koreanus_glmm03_sub_fit[["data"]][["Zacco_koreanus_PA"]])
Zacco_koreanus_glmm03_sub_fit_sub <- as.data.frame(Zacco_koreanus_glmm03_sub_fit[["data"]][["substrate"]])
Zacco_koreanus_sub_res_grp2 <- cbind(Zacco_koreanus_glmm03_sub_fit_PA, Zacco_koreanus_glmm03_sub_fit_sub)
colnames(Zacco_koreanus_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Zacco_koreanus_sub_res_grp2 <- Zacco_koreanus_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Zacco_koreanus_sub_res_grp2$sub <- factor(Zacco_koreanus_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                      'pebble', 'cobble', 'boulder'))

#'- plot
Zacco_koreanus_substrate_fig2 <- ggplot(Zacco_koreanus_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Zacco_koreanus Substrate HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Zacco_koreanus_substrate_fig2




#'## Combined figures

#'- merge mulitple effect plot
Zacco_koreanus_dep_res_1 <- Zacco_koreanus_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)


Zacco_koreanus_dep_res_grp2_1 <- Zacco_koreanus_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)


Zacco_koreanus_dep_res_comb <- rbind(Zacco_koreanus_dep_res_1, 
                                     Zacco_koreanus_dep_res_grp2_1)


#'- depth plot
Zacco_koreanus_depth_fig4 <- ggplot(Zacco_koreanus_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)", limits = c(0,2.0), breaks = seq(0, 2, 0.5)) +     
  labs(title= "Zacco koreanus Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Zacco_koreanus_depth_fig4 <-Zacco_koreanus_depth_fig4+scale_color_manual(values=c("black", "red"))
Zacco_koreanus_depth_fig4




#'- merge mulitple effect plot
Zacco_koreanus_vel_res_1 <- Zacco_koreanus_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)



Zacco_koreanus_vel_res_grp2_1 <- Zacco_koreanus_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)



Zacco_koreanus_vel_res_comb <- rbind(Zacco_koreanus_vel_res_1,
                                     Zacco_koreanus_vel_res_grp2_1)


#'- velocity plot
Zacco_koreanus_velocity_fig4 <- ggplot(Zacco_koreanus_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)", limits = c(0, 3.1), breaks = seq(0, 3, 0.5)) +     
  labs(title= "Zacco koreanus Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Zacco_koreanus_velocity_fig4 <- Zacco_koreanus_velocity_fig4 + scale_color_manual(values=c("black", "red"))
Zacco_koreanus_velocity_fig4



#'- merge mulitple effect plot
Zacco_koreanus_sub_res_1 <- Zacco_koreanus_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)



Zacco_koreanus_sub_res_grp2_1 <- Zacco_koreanus_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)



Zacco_koreanus_sub_res_comb <- rbind(Zacco_koreanus_sub_res_1, 
                                     Zacco_koreanus_sub_res_grp2_1)

#'- substrate plot
Zacco_koreanus_substrate_fig4 <- ggplot(Zacco_koreanus_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Zacco koreanus Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Zacco_koreanus_substrate_fig4 <- Zacco_koreanus_substrate_fig4 + scale_fill_manual(values=c("black", "red"))
Zacco_koreanus_substrate_fig4


#'## combine figures

ggarrange(Zacco_koreanus_depth_fig4, Zacco_koreanus_velocity_fig4, Zacco_koreanus_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 




#------------------------------------------------------------------------------------------------------------------------------------------#




##### Pungtungia_herzi model #####



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

#'##### glm
Pungtungia_herzi_glm01 <- glm(Pungtungia_herzi_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                              family=binomial(link = "logit"), data= Pungtungia_herzi_dat, na.action ="na.fail")

summary(Pungtungia_herzi_glm01)

vif(Pungtungia_herzi_glm01)
#'##### glmm
#'- combination of site and season are used as a random effect 
Pungtungia_herzi_glmm01 <- glmmTMB(Pungtungia_herzi_PA ~ poly(depth, 2) + poly(velocity,1) + substrate + 
                                     (1|site) + (1|season), family=binomial(link = "logit"), data= Pungtungia_herzi_dat, na.action ="na.fail")

summary(Pungtungia_herzi_glmm01)


#'- Model validation
r.squaredGLMM(Pungtungia_herzi_glmm01)
cor.test(predict(Pungtungia_herzi_glmm01, type = "response", allow.new.levels=TRUE), Pungtungia_herzi_dat$Pungtungia_herzi_PA)

#'- Pungtungia_herzi depth
Pungtungia_herzi_glmm01_dep_fit <- effect_plot(Pungtungia_herzi_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Pungtungia_herzi_glmm01_dep_fit_PA <- as.data.frame(Pungtungia_herzi_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Pungtungia_herzi_PA"]])
Pungtungia_herzi_glmm01_dep_fit_dep <- as.data.frame(Pungtungia_herzi_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Pungtungia_herzi_dep_res <- cbind(Pungtungia_herzi_glmm01_dep_fit_PA, Pungtungia_herzi_glmm01_dep_fit_dep)
colnames(Pungtungia_herzi_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Pungtungia_herzi_dep_res <- Pungtungia_herzi_dep_res %>%
  mutate(HSI_dep = scale_values(prob))



#'- plot
Pungtungia_herzi_depth_fig <- ggplot(Pungtungia_herzi_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pungtungia_herzi depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pungtungia_herzi_depth_fig

#'- Pungtungia_herzi velocity
Pungtungia_herzi_glmm01_vel_fit <- effect_plot(Pungtungia_herzi_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Pungtungia_herzi_glmm01_vel_fit_PA <- as.data.frame(Pungtungia_herzi_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Pungtungia_herzi_PA"]])
Pungtungia_herzi_glmm01_vel_fit_vel <- as.data.frame(Pungtungia_herzi_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Pungtungia_herzi_vel_res <- cbind(Pungtungia_herzi_glmm01_vel_fit_PA, Pungtungia_herzi_glmm01_vel_fit_vel)
colnames(Pungtungia_herzi_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Pungtungia_herzi_vel_res <- Pungtungia_herzi_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Pungtungia_herzi_velocity_fig <- ggplot(Pungtungia_herzi_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pungtungia_herzi Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pungtungia_herzi_velocity_fig

#'- Pungtungia_herzi substrate
Pungtungia_herzi_glmm01_sub_fit <- effect_plot(Pungtungia_herzi_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Pungtungia_herzi_glmm01_sub_fit_PA <- as.data.frame(Pungtungia_herzi_glmm01_sub_fit[["data"]][["Pungtungia_herzi_PA"]])
Pungtungia_herzi_glmm01_sub_fit_sub <- as.data.frame(Pungtungia_herzi_glmm01_sub_fit[["data"]][["substrate"]])
Pungtungia_herzi_sub_res <- cbind(Pungtungia_herzi_glmm01_sub_fit_PA, Pungtungia_herzi_glmm01_sub_fit_sub)
colnames(Pungtungia_herzi_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Pungtungia_herzi_sub_res <- Pungtungia_herzi_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Pungtungia_herzi_sub_res$sub <- factor(Pungtungia_herzi_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                'pebble', 'cobble', 'boulder'))

#'- plot
Pungtungia_herzi_substrate_fig <- ggplot(Pungtungia_herzi_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Pungtungia_herzi Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pungtungia_herzi_substrate_fig




##### Pungtungia_herzi group2 model #####

#'## select Group2
Pungtungia_herzi_dat03 <- Pungtungia_herzi_dat %>%
  filter(site_group == "GRP2")

#'##### glm
Pungtungia_herzi_glm03 <- glm(Pungtungia_herzi_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                              family=binomial(link = "logit"), data= Pungtungia_herzi_dat03, na.action ="na.fail")

summary(Pungtungia_herzi_glm03)

vif(Pungtungia_herzi_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect 
Pungtungia_herzi_glmm03 <- glmmTMB(Pungtungia_herzi_PA ~ poly(depth, 2) + poly(velocity,2 ) + substrate + 
                                     (1|site) + (1|season), family=binomial(link = "logit"), data= Pungtungia_herzi_dat03, na.action ="na.fail")

summary(Pungtungia_herzi_glmm03)






#'- Model validation
r.squaredGLMM(Pungtungia_herzi_glmm03)
cor.test(predict(Pungtungia_herzi_glmm03, type = "response", allow.new.levels=TRUE), Pungtungia_herzi_dat03$Pungtungia_herzi_PA)

#'- Pungtungia_herzi depth
Pungtungia_herzi_glmm03_dep_fit <- effect_plot(Pungtungia_herzi_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Pungtungia_herzi_glmm03_dep_fit_PA <- as.data.frame(Pungtungia_herzi_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Pungtungia_herzi_PA"]])
Pungtungia_herzi_glmm03_dep_fit_depth <- as.data.frame(Pungtungia_herzi_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Pungtungia_herzi_dep_res_grp2 <- cbind(Pungtungia_herzi_glmm03_dep_fit_PA, Pungtungia_herzi_glmm03_dep_fit_depth)
colnames(Pungtungia_herzi_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Pungtungia_herzi_dep_res_grp2 <- Pungtungia_herzi_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Pungtungia_herzi_depth_fig2 <- ggplot(Pungtungia_herzi_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pungtungia_herzi depth HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pungtungia_herzi_depth_fig2

#'- Pungtungia_herzi velocity
Pungtungia_herzi_glmm03_vel_fit <- effect_plot(Pungtungia_herzi_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Pungtungia_herzi_glmm03_vel_fit_PA <- as.data.frame(Pungtungia_herzi_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Pungtungia_herzi_PA"]])
Pungtungia_herzi_glmm03_vel_fit_vel <- as.data.frame(Pungtungia_herzi_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Pungtungia_herzi_vel_res_grp2 <- cbind(Pungtungia_herzi_glmm03_vel_fit_PA, Pungtungia_herzi_glmm03_vel_fit_vel)
colnames(Pungtungia_herzi_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Pungtungia_herzi_vel_res_grp2 <- Pungtungia_herzi_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Pungtungia_herzi_velocity_fig2 <- ggplot(Pungtungia_herzi_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pungtungia_herzi Velocity HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pungtungia_herzi_velocity_fig2

#'- Pungtungia_herzi substrate
Pungtungia_herzi_glmm03_sub_fit <- effect_plot(Pungtungia_herzi_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Pungtungia_herzi_glmm03_sub_fit_PA <- as.data.frame(Pungtungia_herzi_glmm03_sub_fit[["data"]][["Pungtungia_herzi_PA"]])
Pungtungia_herzi_glmm03_sub_fit_sub <- as.data.frame(Pungtungia_herzi_glmm03_sub_fit[["data"]][["substrate"]])
Pungtungia_herzi_sub_res_grp2 <- cbind(Pungtungia_herzi_glmm03_sub_fit_PA, Pungtungia_herzi_glmm03_sub_fit_sub)
colnames(Pungtungia_herzi_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Pungtungia_herzi_sub_res_grp2 <- Pungtungia_herzi_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Pungtungia_herzi_sub_res_grp2$sub <- factor(Pungtungia_herzi_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                          'pebble', 'cobble', 'boulder'))

#'- plot
Pungtungia_herzi_substrate_fig2 <- ggplot(Pungtungia_herzi_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Pungtungia_herzi Substrate HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pungtungia_herzi_substrate_fig2







#'## Combined figures

#'- merge mulitple effect plot
Pungtungia_herzi_dep_res_1 <- Pungtungia_herzi_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)



Pungtungia_herzi_dep_res_grp2_1 <- Pungtungia_herzi_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)


Pungtungia_herzi_dep_res_comb <- rbind(Pungtungia_herzi_dep_res_1, 
                                       Pungtungia_herzi_dep_res_grp2_1)


#'- depth plot
Pungtungia_herzi_depth_fig4 <- ggplot(Pungtungia_herzi_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) +
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)", limits = c(0,2.0) ,breaks = seq(0, 2, 0.5)) +
  labs(title= "Pungtungia herzi Depth") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "top",
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))
Pungtungia_herzi_depth_fig4 <-Pungtungia_herzi_depth_fig4+scale_color_manual(values=c("black", "red"))
Pungtungia_herzi_depth_fig4




#'- merge mulitple effect plot
Pungtungia_herzi_vel_res_1 <- Pungtungia_herzi_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)



Pungtungia_herzi_vel_res_grp2_1 <- Pungtungia_herzi_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)




Pungtungia_herzi_vel_res_comb <- rbind(Pungtungia_herzi_vel_res_1, 
                                       Pungtungia_herzi_vel_res_grp2_1)


#'- velocity plot
Pungtungia_herzi_velocity_fig4 <- ggplot(Pungtungia_herzi_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",limits = c(0, 3.1) , breaks = seq(0, 3, 0.5)) +     
  labs(title= "Pungtungia herzi Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Pungtungia_herzi_velocity_fig4 <- Pungtungia_herzi_velocity_fig4 + scale_color_manual(values=c("black", "red"))
Pungtungia_herzi_velocity_fig4



#'- merge mulitple effect plot
Pungtungia_herzi_sub_res_1 <- Pungtungia_herzi_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)


Pungtungia_herzi_sub_res_grp2_1 <- Pungtungia_herzi_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)




Pungtungia_herzi_sub_res_comb <- rbind(Pungtungia_herzi_sub_res_1, 
                                       Pungtungia_herzi_sub_res_grp2_1)

#'- substrate plot
Pungtungia_herzi_substrate_fig4 <- ggplot(Pungtungia_herzi_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Pungtungia herzi Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Pungtungia_herzi_substrate_fig4 <- Pungtungia_herzi_substrate_fig4 + scale_fill_manual(values=c("black", "red"))
Pungtungia_herzi_substrate_fig4


#'## combine figures

ggarrange(Pungtungia_herzi_depth_fig4, Pungtungia_herzi_velocity_fig4, Pungtungia_herzi_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 





#------------------------------------------------------------------------------------------------------------------------------------------#





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
  filter(Pseudopungtungia_nigra_Pre == "Yes")%>%
  as_tibble()


#'### Global Model

#'##### glm
Pseudopungtungia_nigra_glm01 <- glm(Pseudopungtungia_nigra_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                    family=binomial(link = "logit"), data= Pseudopungtungia_nigra_dat, na.action ="na.fail")

summary(Pseudopungtungia_nigra_glm01)
vif(Pseudopungtungia_nigra_glm01)

#'##### glmm
#'- combination of site and season are used as a random effect 
Pseudopungtungia_nigra_glmm01 <- glmmTMB(Pseudopungtungia_nigra_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                           (1|site) + (1|season), family=binomial(link = "logit"), data= Pseudopungtungia_nigra_dat, na.action ="na.fail")

summary(Pseudopungtungia_nigra_glmm01)











#'- Model validation
r.squaredGLMM(Pseudopungtungia_nigra_glmm01)
cor.test(predict(Pseudopungtungia_nigra_glmm01, type = "response", allow.new.levels=TRUE), Pseudopungtungia_nigra_dat$Pseudopungtungia_nigra_PA)

#'- Pseudopungtungia_nigra depth
Pseudopungtungia_nigra_glmm01_dep_fit <- effect_plot(Pseudopungtungia_nigra_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Pseudopungtungia_nigra_glmm01_dep_fit_PA <- as.data.frame(Pseudopungtungia_nigra_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Pseudopungtungia_nigra_PA"]])
Pseudopungtungia_nigra_glmm01_dep_fit_dep <- as.data.frame(Pseudopungtungia_nigra_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Pseudopungtungia_nigra_dep_res <- cbind(Pseudopungtungia_nigra_glmm01_dep_fit_PA, Pseudopungtungia_nigra_glmm01_dep_fit_dep)
colnames(Pseudopungtungia_nigra_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Pseudopungtungia_nigra_dep_res <- Pseudopungtungia_nigra_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Pseudopungtungia_nigra_depth_fig <- ggplot(Pseudopungtungia_nigra_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudopungtungia_nigra depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudopungtungia_nigra_depth_fig

#'- Pseudopungtungia_nigra velocity
Pseudopungtungia_nigra_glmm01_vel_fit <- effect_plot(Pseudopungtungia_nigra_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Pseudopungtungia_nigra_glmm01_vel_fit_PA <- as.data.frame(Pseudopungtungia_nigra_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Pseudopungtungia_nigra_PA"]])
Pseudopungtungia_nigra_glmm01_vel_fit_vel <- as.data.frame(Pseudopungtungia_nigra_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Pseudopungtungia_nigra_vel_res <- cbind(Pseudopungtungia_nigra_glmm01_vel_fit_PA, Pseudopungtungia_nigra_glmm01_vel_fit_vel)
colnames(Pseudopungtungia_nigra_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Pseudopungtungia_nigra_vel_res <- Pseudopungtungia_nigra_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Pseudopungtungia_nigra_velocity_fig <- ggplot(Pseudopungtungia_nigra_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudopungtungia_nigra Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudopungtungia_nigra_velocity_fig

#'- Pseudopungtungia_nigra substrate
Pseudopungtungia_nigra_glmm01_sub_fit <- effect_plot(Pseudopungtungia_nigra_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Pseudopungtungia_nigra_glmm01_sub_fit_PA <- as.data.frame(Pseudopungtungia_nigra_glmm01_sub_fit[["data"]][["Pseudopungtungia_nigra_PA"]])
Pseudopungtungia_nigra_glmm01_sub_fit_sub <- as.data.frame(Pseudopungtungia_nigra_glmm01_sub_fit[["data"]][["substrate"]])
Pseudopungtungia_nigra_sub_res <- cbind(Pseudopungtungia_nigra_glmm01_sub_fit_PA, Pseudopungtungia_nigra_glmm01_sub_fit_sub)
colnames(Pseudopungtungia_nigra_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Pseudopungtungia_nigra_sub_res <- Pseudopungtungia_nigra_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Pseudopungtungia_nigra_sub_res$sub <- factor(Pseudopungtungia_nigra_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                            'pebble', 'cobble', 'boulder'))

#'- plot
Pseudopungtungia_nigra_substrate_fig <- ggplot(Pseudopungtungia_nigra_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Pseudopungtungia_nigra Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudopungtungia_nigra_substrate_fig






##### Pseudopungtungia_nigra group2 model #####

#'## select Group2
Pseudopungtungia_nigra_dat03 <- Pseudopungtungia_nigra_dat %>%
  filter(site_group == "GRP2")

#'##### glm
Pseudopungtungia_nigra_glm03 <- glm(Pseudopungtungia_nigra_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                    family=binomial(link = "logit"), data= Pseudopungtungia_nigra_dat03, na.action ="na.fail")

summary(Pseudopungtungia_nigra_glm03)

vif(Pseudopungtungia_nigra_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect 
Pseudopungtungia_nigra_glmm03 <- glmmTMB(Pseudopungtungia_nigra_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                           (1|site) + (1|season), family=binomial(link = "logit"), data= Pseudopungtungia_nigra_dat03, na.action ="na.fail")

summary(Pseudopungtungia_nigra_glmm03)









#'- Model validation
r.squaredGLMM(Pseudopungtungia_nigra_glmm03)
cor.test(predict(Pseudopungtungia_nigra_glmm03, type = "response", allow.new.levels=TRUE), Pseudopungtungia_nigra_dat03$Pseudopungtungia_nigra_PA)

#'- Pseudopungtungia_nigra depth
Pseudopungtungia_nigra_glmm03_dep_fit <- effect_plot(Pseudopungtungia_nigra_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Pseudopungtungia_nigra_glmm03_dep_fit_PA <- as.data.frame(Pseudopungtungia_nigra_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Pseudopungtungia_nigra_PA"]])
Pseudopungtungia_nigra_glmm03_dep_fit_depth <- as.data.frame(Pseudopungtungia_nigra_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Pseudopungtungia_nigra_dep_res_grp2 <- cbind(Pseudopungtungia_nigra_glmm03_dep_fit_PA, Pseudopungtungia_nigra_glmm03_dep_fit_depth)
colnames(Pseudopungtungia_nigra_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Pseudopungtungia_nigra_dep_res_grp2 <- Pseudopungtungia_nigra_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Pseudopungtungia_nigra_depth_fig2 <- ggplot(Pseudopungtungia_nigra_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudopungtungia_nigra depth HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudopungtungia_nigra_depth_fig2

#'- Pseudopungtungia_nigra velocity
Pseudopungtungia_nigra_glmm03_vel_fit <- effect_plot(Pseudopungtungia_nigra_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Pseudopungtungia_nigra_glmm03_vel_fit_PA <- as.data.frame(Pseudopungtungia_nigra_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Pseudopungtungia_nigra_PA"]])
Pseudopungtungia_nigra_glmm03_vel_fit_vel <- as.data.frame(Pseudopungtungia_nigra_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Pseudopungtungia_nigra_vel_res_grp2 <- cbind(Pseudopungtungia_nigra_glmm03_vel_fit_PA, Pseudopungtungia_nigra_glmm03_vel_fit_vel)
colnames(Pseudopungtungia_nigra_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Pseudopungtungia_nigra_vel_res_grp2 <- Pseudopungtungia_nigra_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Pseudopungtungia_nigra_velocity_fig2 <- ggplot(Pseudopungtungia_nigra_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudopungtungia_nigra Velocity HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudopungtungia_nigra_velocity_fig2

#'- Pseudopungtungia_nigra substrate
Pseudopungtungia_nigra_glmm03_sub_fit <- effect_plot(Pseudopungtungia_nigra_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Pseudopungtungia_nigra_glmm03_sub_fit_PA <- as.data.frame(Pseudopungtungia_nigra_glmm03_sub_fit[["data"]][["Pseudopungtungia_nigra_PA"]])
Pseudopungtungia_nigra_glmm03_sub_fit_sub <- as.data.frame(Pseudopungtungia_nigra_glmm03_sub_fit[["data"]][["substrate"]])
Pseudopungtungia_nigra_sub_res_grp2 <- cbind(Pseudopungtungia_nigra_glmm03_sub_fit_PA, Pseudopungtungia_nigra_glmm03_sub_fit_sub)
colnames(Pseudopungtungia_nigra_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Pseudopungtungia_nigra_sub_res_grp2 <- Pseudopungtungia_nigra_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Pseudopungtungia_nigra_sub_res_grp2$sub <- factor(Pseudopungtungia_nigra_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                      'pebble', 'cobble', 'boulder'))

#'- plot
Pseudopungtungia_nigra_substrate_fig2 <- ggplot(Pseudopungtungia_nigra_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Pseudopungtungia_nigra Substrate HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudopungtungia_nigra_substrate_fig2



#'## Combined figures

#'- merge mulitple effect plot
Pseudopungtungia_nigra_dep_res_1 <- Pseudopungtungia_nigra_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)


Pseudopungtungia_nigra_dep_res_grp2_1 <- Pseudopungtungia_nigra_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)


Pseudopungtungia_nigra_dep_res_comb <- rbind(Pseudopungtungia_nigra_dep_res_1,
                                             Pseudopungtungia_nigra_dep_res_grp2_1)


#'- depth plot
Pseudopungtungia_nigra_depth_fig4 <- ggplot(Pseudopungtungia_nigra_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudopungtungia_nigra Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Pseudopungtungia_nigra_depth_fig4 <-Pseudopungtungia_nigra_depth_fig4+scale_color_manual(values=c("black", "red"))
Pseudopungtungia_nigra_depth_fig4




#'- merge mulitple effect plot
Pseudopungtungia_nigra_vel_res_1 <- Pseudopungtungia_nigra_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)


Pseudopungtungia_nigra_vel_res_grp2_1 <- Pseudopungtungia_nigra_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)



Pseudopungtungia_nigra_vel_res_comb <- rbind(Pseudopungtungia_nigra_vel_res_1, 
                                             Pseudopungtungia_nigra_vel_res_grp2_1)


#'- velocity plot
Pseudopungtungia_nigra_velocity_fig4 <- ggplot(Pseudopungtungia_nigra_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudopungtungia_nigra Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Pseudopungtungia_nigra_velocity_fig4 <- Pseudopungtungia_nigra_velocity_fig4 + scale_color_manual(values=c("black","red"))
Pseudopungtungia_nigra_velocity_fig4



#'- merge mulitple effect plot
Pseudopungtungia_nigra_sub_res_1 <- Pseudopungtungia_nigra_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)


Pseudopungtungia_nigra_sub_res_grp2_1 <- Pseudopungtungia_nigra_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)



Pseudopungtungia_nigra_sub_res_comb <- rbind(Pseudopungtungia_nigra_sub_res_1,
                                             Pseudopungtungia_nigra_sub_res_grp2_1 )

#'- substrate plot
Pseudopungtungia_nigra_substrate_fig4 <- ggplot(Pseudopungtungia_nigra_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Pseudopungtungia_nigra Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Pseudopungtungia_nigra_substrate_fig4 <- Pseudopungtungia_nigra_substrate_fig4 + scale_fill_manual(values=c("black","red"))
Pseudopungtungia_nigra_substrate_fig4


#'## combine figures

ggarrange(Pseudopungtungia_nigra_depth_fig4, Pseudopungtungia_nigra_velocity_fig4, Pseudopungtungia_nigra_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 







#------------------------------------------------------------------------------------------------------------------------------------------#













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
  filter(Coreoleuciscus_splendidus_Pre == "Yes")%>%
  ungroup()



#'## (2) Coreoleuciscus_splendidus model

#'### data management
Coreoleuciscus_splendidus_dat <- HSI_DB2 %>%
  mutate(Coreoleuciscus_splendidus_PA = ifelse(species == "Coreoleuciscus_splendidus", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Coreoleuciscus_splendidus_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Coreoleuciscus_splendidus_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  filter(!basin == "Nakdong_River") %>%
  filter(!basin == "Seomjin_River") %>%
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Coreoleuciscus_splendidus_PA, depth, velocity, substrate, MBSNCD) 


#'## remove site where species did not occur
Coreoleuciscus_splendidus_Pre <- HSI_DB2 %>%
  mutate(Coreoleuciscus_splendidus_PA = ifelse(species == "Coreoleuciscus_splendidus", 1, 0)) %>%
  aggregate(Coreoleuciscus_splendidus_PA ~ site, "sum") %>%
  mutate(Coreoleuciscus_splendidus_Pre = if_else(Coreoleuciscus_splendidus_PA > 0, "Yes", "No")) %>%
  select(site, Coreoleuciscus_splendidus_Pre)

Coreoleuciscus_splendidus_dat <- Coreoleuciscus_splendidus_dat %>%
  left_join(Coreoleuciscus_splendidus_Pre, by = "site") %>%
  filter(Coreoleuciscus_splendidus_Pre == "Yes") %>%
  as_tibble()






Coreoleuciscus_splendidus_zero_dat <- Coreoleuciscus_splendidus_dat%>%
  select(basin, site_group, site_by_season, site, date, season, MBSNCD, Coreoleuciscus_splendidus_Pre)%>%
  distinct()





#'### Global Model

#'##### glm
Coreoleuciscus_splendidus_glm01 <- glm(Coreoleuciscus_splendidus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                       family=binomial(link = "logit"), data= Coreoleuciscus_splendidus_dat, na.action ="na.fail")

summary(Coreoleuciscus_splendidus_glm01)

vif(Coreoleuciscus_splendidus_glm01)

#'##### glmm
#'- combination of site and season are used as a random effect 
Coreoleuciscus_splendidus_glmm01 <- glmmTMB(Coreoleuciscus_splendidus_PA ~ poly(depth, 2)  + poly(velocity, 2) + substrate + 
                                              (1|site) + (1|season), family=binomial(link = "logit"), data= Coreoleuciscus_splendidus_dat, na.action ="na.fail")

summary(Coreoleuciscus_splendidus_glmm01)



#'- Model validation
r.squaredGLMM(Coreoleuciscus_splendidus_glmm01)
cor.test(predict(Coreoleuciscus_splendidus_glmm01, type = "response", allow.new.levels=TRUE), Coreoleuciscus_splendidus_dat$Coreoleuciscus_splendidus_PA)

#'- Coreoleuciscus_splendidus depth
Coreoleuciscus_splendidus_glmm01_dep_fit <- effect_plot(Coreoleuciscus_splendidus_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Coreoleuciscus_splendidus_glmm01_dep_fit_PA <- as.data.frame(Coreoleuciscus_splendidus_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Coreoleuciscus_splendidus_PA"]])
Coreoleuciscus_splendidus_glmm01_dep_fit_dep <- as.data.frame(Coreoleuciscus_splendidus_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Coreoleuciscus_splendidus_dep_res <- cbind(Coreoleuciscus_splendidus_glmm01_dep_fit_PA, Coreoleuciscus_splendidus_glmm01_dep_fit_dep)
colnames(Coreoleuciscus_splendidus_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Coreoleuciscus_splendidus_dep_res <- Coreoleuciscus_splendidus_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Coreoleuciscus_splendidus_depth_fig <- ggplot(Coreoleuciscus_splendidus_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Coreoleuciscus_splendidus depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Coreoleuciscus_splendidus_depth_fig

#'- Coreoleuciscus_splendidus velocity
Coreoleuciscus_splendidus_glmm01_vel_fit <- effect_plot(Coreoleuciscus_splendidus_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Coreoleuciscus_splendidus_glmm01_vel_fit_PA <- as.data.frame(Coreoleuciscus_splendidus_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Coreoleuciscus_splendidus_PA"]])
Coreoleuciscus_splendidus_glmm01_vel_fit_vel <- as.data.frame(Coreoleuciscus_splendidus_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Coreoleuciscus_splendidus_vel_res <- cbind(Coreoleuciscus_splendidus_glmm01_vel_fit_PA, Coreoleuciscus_splendidus_glmm01_vel_fit_vel)
colnames(Coreoleuciscus_splendidus_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Coreoleuciscus_splendidus_vel_res <- Coreoleuciscus_splendidus_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Coreoleuciscus_splendidus_velocity_fig <- ggplot(Coreoleuciscus_splendidus_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Coreoleuciscus_splendidus Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Coreoleuciscus_splendidus_velocity_fig

#'- Coreoleuciscus_splendidus substrate
Coreoleuciscus_splendidus_glmm01_sub_fit <- effect_plot(Coreoleuciscus_splendidus_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Coreoleuciscus_splendidus_glmm01_sub_fit_PA <- as.data.frame(Coreoleuciscus_splendidus_glmm01_sub_fit[["data"]][["Coreoleuciscus_splendidus_PA"]])
Coreoleuciscus_splendidus_glmm01_sub_fit_sub <- as.data.frame(Coreoleuciscus_splendidus_glmm01_sub_fit[["data"]][["substrate"]])
Coreoleuciscus_splendidus_sub_res <- cbind(Coreoleuciscus_splendidus_glmm01_sub_fit_PA, Coreoleuciscus_splendidus_glmm01_sub_fit_sub)
colnames(Coreoleuciscus_splendidus_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Coreoleuciscus_splendidus_sub_res <- Coreoleuciscus_splendidus_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Coreoleuciscus_splendidus_sub_res$sub <- factor(Coreoleuciscus_splendidus_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                  'pebble', 'cobble', 'boulder'))

#'- plot
Coreoleuciscus_splendidus_substrate_fig <- ggplot(Coreoleuciscus_splendidus_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Coreoleuciscus_splendidus Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Coreoleuciscus_splendidus_substrate_fig






##### Coreoleuciscus_splendidus group2 model #####

#'## select Group2
Coreoleuciscus_splendidus_dat03 <- Coreoleuciscus_splendidus_dat %>%
  filter(site_group == "GRP2")

#'##### glm
Coreoleuciscus_splendidus_glm03 <- glm(Coreoleuciscus_splendidus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                       family=binomial(link = "logit"), data= Coreoleuciscus_splendidus_dat03, na.action ="na.fail")

summary(Coreoleuciscus_splendidus_glm03)

vif(Coreoleuciscus_splendidus_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect 
Coreoleuciscus_splendidus_glmm03 <- glmmTMB(Coreoleuciscus_splendidus_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + (1|site) + (1|season), family=binomial(link = "logit"), data= Coreoleuciscus_splendidus_dat03, na.action ="na.fail")

summary(Coreoleuciscus_splendidus_glmm03)




#'- Model validation
r.squaredGLMM(Coreoleuciscus_splendidus_glmm03)
cor.test(predict(Coreoleuciscus_splendidus_glmm03, type = "response", allow.new.levels=TRUE), Coreoleuciscus_splendidus_dat03$Coreoleuciscus_splendidus_PA)

#'- Coreoleuciscus_splendidus depth
Coreoleuciscus_splendidus_glmm03_dep_fit <- effect_plot(Coreoleuciscus_splendidus_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Coreoleuciscus_splendidus_glmm03_dep_fit_PA <- as.data.frame(Coreoleuciscus_splendidus_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Coreoleuciscus_splendidus_PA"]])
Coreoleuciscus_splendidus_glmm03_dep_fit_depth <- as.data.frame(Coreoleuciscus_splendidus_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Coreoleuciscus_splendidus_dep_res_grp2 <- cbind(Coreoleuciscus_splendidus_glmm03_dep_fit_PA, Coreoleuciscus_splendidus_glmm03_dep_fit_depth)
colnames(Coreoleuciscus_splendidus_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Coreoleuciscus_splendidus_dep_res_grp2 <- Coreoleuciscus_splendidus_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Coreoleuciscus_splendidus_depth_fig2 <- ggplot(Coreoleuciscus_splendidus_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Coreoleuciscus_splendidus depth HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Coreoleuciscus_splendidus_depth_fig2

#'- Coreoleuciscus_splendidus velocity
Coreoleuciscus_splendidus_glmm03_vel_fit <- effect_plot(Coreoleuciscus_splendidus_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Coreoleuciscus_splendidus_glmm03_vel_fit_PA <- as.data.frame(Coreoleuciscus_splendidus_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Coreoleuciscus_splendidus_PA"]])
Coreoleuciscus_splendidus_glmm03_vel_fit_vel <- as.data.frame(Coreoleuciscus_splendidus_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Coreoleuciscus_splendidus_vel_res_grp2 <- cbind(Coreoleuciscus_splendidus_glmm03_vel_fit_PA, Coreoleuciscus_splendidus_glmm03_vel_fit_vel)
colnames(Coreoleuciscus_splendidus_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Coreoleuciscus_splendidus_vel_res_grp2 <- Coreoleuciscus_splendidus_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Coreoleuciscus_splendidus_velocity_fig2 <- ggplot(Coreoleuciscus_splendidus_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Coreoleuciscus_splendidus Velocity HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Coreoleuciscus_splendidus_velocity_fig2

#'- Coreoleuciscus_splendidus substrate
Coreoleuciscus_splendidus_glmm03_sub_fit <- effect_plot(Coreoleuciscus_splendidus_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Coreoleuciscus_splendidus_glmm03_sub_fit_PA <- as.data.frame(Coreoleuciscus_splendidus_glmm03_sub_fit[["data"]][["Coreoleuciscus_splendidus_PA"]])
Coreoleuciscus_splendidus_glmm03_sub_fit_sub <- as.data.frame(Coreoleuciscus_splendidus_glmm03_sub_fit[["data"]][["substrate"]])
Coreoleuciscus_splendidus_sub_res_grp2 <- cbind(Coreoleuciscus_splendidus_glmm03_sub_fit_PA, Coreoleuciscus_splendidus_glmm03_sub_fit_sub)
colnames(Coreoleuciscus_splendidus_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Coreoleuciscus_splendidus_sub_res_grp2 <- Coreoleuciscus_splendidus_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Coreoleuciscus_splendidus_sub_res_grp2$sub <- factor(Coreoleuciscus_splendidus_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                            'pebble', 'cobble', 'boulder'))

#'- plot
Coreoleuciscus_splendidus_substrate_fig2 <- ggplot(Coreoleuciscus_splendidus_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Coreoleuciscus_splendidus Substrate HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Coreoleuciscus_splendidus_substrate_fig2







#'## Combined figures

#'- merge mulitple effect plot
Coreoleuciscus_splendidus_dep_res_1 <- Coreoleuciscus_splendidus_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)

Coreoleuciscus_splendidus_dep_res_grp2_1 <- Coreoleuciscus_splendidus_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)

Coreoleuciscus_splendidus_dep_res_comb <- rbind(Coreoleuciscus_splendidus_dep_res_1, 
                                                Coreoleuciscus_splendidus_dep_res_grp2_1)

#'- depth plot
Coreoleuciscus_splendidus_depth_fig4 <- ggplot(Coreoleuciscus_splendidus_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",limits = c(0,2.0) , breaks = seq(0, 2.5, 0.5)) +     
  labs(title= "Coreoleuciscus splendidus Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Coreoleuciscus_splendidus_depth_fig4 <-Coreoleuciscus_splendidus_depth_fig4+scale_color_manual(values=c("black","red"))
Coreoleuciscus_splendidus_depth_fig4




#'- merge mulitple effect plot
Coreoleuciscus_splendidus_vel_res_1 <- Coreoleuciscus_splendidus_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)


Coreoleuciscus_splendidus_vel_res_grp2_1 <- Coreoleuciscus_splendidus_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)


Coreoleuciscus_splendidus_vel_res_comb <- rbind(Coreoleuciscus_splendidus_vel_res_1,
                                                Coreoleuciscus_splendidus_vel_res_grp2_1)


#'- velocity plot
Coreoleuciscus_splendidus_velocity_fig4 <- ggplot(Coreoleuciscus_splendidus_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)", limits = c(0,3.1), breaks = seq(0, 3, 0.5)) +     
  labs(title= "Coreoleuciscus splendidus Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Coreoleuciscus_splendidus_velocity_fig4 <- Coreoleuciscus_splendidus_velocity_fig4 + scale_color_manual(values=c("black", "red"))
Coreoleuciscus_splendidus_velocity_fig4



#'- merge mulitple effect plot
Coreoleuciscus_splendidus_sub_res_1 <- Coreoleuciscus_splendidus_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)


Coreoleuciscus_splendidus_sub_res_grp2_1 <- Coreoleuciscus_splendidus_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)


Coreoleuciscus_splendidus_sub_res_comb <- rbind(Coreoleuciscus_splendidus_sub_res_1, 
                                                Coreoleuciscus_splendidus_sub_res_grp2_1)

#'- substrate plot
Coreoleuciscus_splendidus_substrate_fig4 <- ggplot(Coreoleuciscus_splendidus_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Coreoleuciscus splendidus Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Coreoleuciscus_splendidus_substrate_fig4 <- Coreoleuciscus_splendidus_substrate_fig4 + scale_fill_manual(values=c("black", "red"))
Coreoleuciscus_splendidus_substrate_fig4


#'## combine figures

ggarrange(Coreoleuciscus_splendidus_depth_fig4, Coreoleuciscus_splendidus_velocity_fig4, Coreoleuciscus_splendidus_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 








##### Opsariichthys_uncirostris_amurensis model #####

#'## (2) Opsariichthys_uncirostris_amurensis model

#'### data management
Opsariichthys_uncirostris_amurensis_dat <- HSI_DB2 %>%
  mutate(Opsariichthys_uncirostris_amurensis_PA = ifelse(species == "Opsariichthys_uncirostris_amurensis", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Opsariichthys_uncirostris_amurensis_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Opsariichthys_uncirostris_amurensis_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Opsariichthys_uncirostris_amurensis_PA, depth, velocity, substrate, MBSNCD) 

head(Opsariichthys_uncirostris_amurensis_dat)


#'## remove site where species did not occur
Opsariichthys_uncirostris_amurensis_Pre <- HSI_DB2 %>%
  mutate(Opsariichthys_uncirostris_amurensis_PA = ifelse(species == "Opsariichthys_uncirostris_amurensis", 1, 0)) %>%
  aggregate(Opsariichthys_uncirostris_amurensis_PA ~ site, "sum") %>%
  mutate(Opsariichthys_uncirostris_amurensis_Pre = if_else(Opsariichthys_uncirostris_amurensis_PA > 0, "Yes", "No")) %>%
  select(site, Opsariichthys_uncirostris_amurensis_Pre)

Opsariichthys_uncirostris_amurensis_dat <- Opsariichthys_uncirostris_amurensis_dat %>%
  left_join(Opsariichthys_uncirostris_amurensis_Pre, by = "site") %>%
  filter(Opsariichthys_uncirostris_amurensis_Pre == "Yes") %>%
  as_tibble()


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

#'##### glm
Opsariichthys_uncirostris_amurensis_glm01 <- glm(Opsariichthys_uncirostris_amurensis_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                                 family=binomial, data= Opsariichthys_uncirostris_amurensis_dat, na.action ="na.fail")

summary(Opsariichthys_uncirostris_amurensis_glm01)

vif(Opsariichthys_uncirostris_amurensis_glm01)



#'##### glmm
#'- combination of site and season are used as a random effect 
Opsariichthys_uncirostris_amurensis_glmm01 <- glmmTMB(Opsariichthys_uncirostris_amurensis_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                                        (1|site) + (1|season), family=binomial, data= Opsariichthys_uncirostris_amurensis_dat, na.action ="na.fail")

summary(Opsariichthys_uncirostris_amurensis_glmm01)




#'- Model validation
r.squaredGLMM(Opsariichthys_uncirostris_amurensis_glmm01)
cor.test(predict(Opsariichthys_uncirostris_amurensis_glmm01, type = "response", allow.new.levels=TRUE), Opsariichthys_uncirostris_amurensis_dat$Opsariichthys_uncirostris_amurensis_PA)

#'- Opsariichthys_uncirostris_amurensis depth
Opsariichthys_uncirostris_amurensis_glmm01_dep_fit <- effect_plot(Opsariichthys_uncirostris_amurensis_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Opsariichthys_uncirostris_amurensis_glmm01_dep_fit_PA <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Opsariichthys_uncirostris_amurensis_PA"]])
Opsariichthys_uncirostris_amurensis_glmm01_dep_fit_dep <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Opsariichthys_uncirostris_amurensis_dep_res <- cbind(Opsariichthys_uncirostris_amurensis_glmm01_dep_fit_PA, Opsariichthys_uncirostris_amurensis_glmm01_dep_fit_dep)
colnames(Opsariichthys_uncirostris_amurensis_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Opsariichthys_uncirostris_amurensis_dep_res <- Opsariichthys_uncirostris_amurensis_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Opsariichthys_uncirostris_amurensis_depth_fig <- ggplot(Opsariichthys_uncirostris_amurensis_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Opsariichthys_uncirostris_amurensis depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Opsariichthys_uncirostris_amurensis_depth_fig

#'- Opsariichthys_uncirostris_amurensis velocity
Opsariichthys_uncirostris_amurensis_glmm01_vel_fit <- effect_plot(Opsariichthys_uncirostris_amurensis_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Opsariichthys_uncirostris_amurensis_glmm01_vel_fit_PA <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Opsariichthys_uncirostris_amurensis_PA"]])
Opsariichthys_uncirostris_amurensis_glmm01_vel_fit_vel <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Opsariichthys_uncirostris_amurensis_vel_res <- cbind(Opsariichthys_uncirostris_amurensis_glmm01_vel_fit_PA, Opsariichthys_uncirostris_amurensis_glmm01_vel_fit_vel)
colnames(Opsariichthys_uncirostris_amurensis_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Opsariichthys_uncirostris_amurensis_vel_res <- Opsariichthys_uncirostris_amurensis_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Opsariichthys_uncirostris_amurensis_velocity_fig <- ggplot(Opsariichthys_uncirostris_amurensis_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Opsariichthys_uncirostris_amurensis Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Opsariichthys_uncirostris_amurensis_velocity_fig

#'- Opsariichthys_uncirostris_amurensis substrate
Opsariichthys_uncirostris_amurensis_glmm01_sub_fit <- effect_plot(Opsariichthys_uncirostris_amurensis_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Opsariichthys_uncirostris_amurensis_glmm01_sub_fit_PA <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm01_sub_fit[["data"]][["Opsariichthys_uncirostris_amurensis_PA"]])
Opsariichthys_uncirostris_amurensis_glmm01_sub_fit_sub <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm01_sub_fit[["data"]][["substrate"]])
Opsariichthys_uncirostris_amurensis_sub_res <- cbind(Opsariichthys_uncirostris_amurensis_glmm01_sub_fit_PA, Opsariichthys_uncirostris_amurensis_glmm01_sub_fit_sub)
colnames(Opsariichthys_uncirostris_amurensis_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Opsariichthys_uncirostris_amurensis_sub_res <- Opsariichthys_uncirostris_amurensis_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Opsariichthys_uncirostris_amurensis_sub_res$sub <- factor(Opsariichthys_uncirostris_amurensis_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                                      'pebble', 'cobble', 'boulder'))

#'- plot
Opsariichthys_uncirostris_amurensis_substrate_fig <- ggplot(Opsariichthys_uncirostris_amurensis_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Opsariichthys_uncirostris_amurensis Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Opsariichthys_uncirostris_amurensis_substrate_fig


#'## (1) Opsariichthys_uncirostris_amurensis model
##### Opsariichthys_uncirostris_amurensis group1 model #####

#'## select Group1

Opsariichthys_uncirostris_amurensis_dat02 <- Opsariichthys_uncirostris_amurensis_dat %>%
  filter(site_group == "GRP1")

#'##### glm
Opsariichthys_uncirostris_amurensis_glm02 <- glm(Opsariichthys_uncirostris_amurensis_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                                 family=binomial(link = "logit"), data= Opsariichthys_uncirostris_amurensis_dat02, na.action ="na.fail")

summary(Opsariichthys_uncirostris_amurensis_glm02)

vif(Opsariichthys_uncirostris_amurensis_glm02)
#'##### glmm
#'- combination of site and season are used as a random effect ``
Opsariichthys_uncirostris_amurensis_glmm02 <- glmmTMB(Opsariichthys_uncirostris_amurensis_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                                        (1|site) + (1|season), family=binomial(link = "logit"), data= Opsariichthys_uncirostris_amurensis_dat02, na.action ="na.fail")

summary(Opsariichthys_uncirostris_amurensis_glmm02)



#'- Model validation
r.squaredGLMM(Opsariichthys_uncirostris_amurensis_glmm02)
cor.test(predict(Opsariichthys_uncirostris_amurensis_glmm02, type = "response", allow.new.levels=TRUE), Opsariichthys_uncirostris_amurensis_dat02$Opsariichthys_uncirostris_amurensis_PA)

#'- Opsariichthys_uncirostris_amurensis depth
Opsariichthys_uncirostris_amurensis_glmm02_dep_fit <- effect_plot(Opsariichthys_uncirostris_amurensis_glmm02, pred = depth, interval = TRUE, plot.points = FALSE)
Opsariichthys_uncirostris_amurensis_glmm02_dep_fit_PA <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["Opsariichthys_uncirostris_amurensis_PA"]])
Opsariichthys_uncirostris_amurensis_glmm02_dep_fit_depth <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Opsariichthys_uncirostris_amurensis_dep_res_grp1 <- cbind(Opsariichthys_uncirostris_amurensis_glmm02_dep_fit_PA, Opsariichthys_uncirostris_amurensis_glmm02_dep_fit_depth)
colnames(Opsariichthys_uncirostris_amurensis_dep_res_grp1) <- c("prob", "dep")

#'- normalize 0 to 1
Opsariichthys_uncirostris_amurensis_dep_res_grp1 <- Opsariichthys_uncirostris_amurensis_dep_res_grp1 %>%
  mutate(HSI_dep_grp1 = scale_values(prob))

#'- plot
Opsariichthys_uncirostris_amurensis_depth_fig1 <- ggplot(Opsariichthys_uncirostris_amurensis_dep_res_grp1, aes(x = dep, y = HSI_dep_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Opsariichthys_uncirostris_amurensis depth HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Opsariichthys_uncirostris_amurensis_depth_fig1

#'- Opsariichthys_uncirostris_amurensis velocity
Opsariichthys_uncirostris_amurensis_glmm02_vel_fit <- effect_plot(Opsariichthys_uncirostris_amurensis_glmm02, pred = velocity, interval = TRUE, plot.points = FALSE)
Opsariichthys_uncirostris_amurensis_glmm02_vel_fit_PA <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["Opsariichthys_uncirostris_amurensis_PA"]])
Opsariichthys_uncirostris_amurensis_glmm02_vel_fit_vel <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Opsariichthys_uncirostris_amurensis_vel_res_grp1 <- cbind(Opsariichthys_uncirostris_amurensis_glmm02_vel_fit_PA, Opsariichthys_uncirostris_amurensis_glmm02_vel_fit_vel)
colnames(Opsariichthys_uncirostris_amurensis_vel_res_grp1) <- c("prob", "vel")

#'- normalize 0 to 1
Opsariichthys_uncirostris_amurensis_vel_res_grp1 <- Opsariichthys_uncirostris_amurensis_vel_res_grp1 %>%
  mutate(HSI_vel_grp1 = scale_values(prob))

#'- plot
Opsariichthys_uncirostris_amurensis_velocity_fig1 <- ggplot(Opsariichthys_uncirostris_amurensis_vel_res_grp1, aes(x = vel, y = HSI_vel_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Opsariichthys_uncirostris_amurensis Velocity HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Opsariichthys_uncirostris_amurensis_velocity_fig1

#'- Opsariichthys_uncirostris_amurensis substrate
Opsariichthys_uncirostris_amurensis_glmm02_sub_fit <- effect_plot(Opsariichthys_uncirostris_amurensis_glmm02, pred = substrate, interval = FALSE, plot.points = FALSE)
Opsariichthys_uncirostris_amurensis_glmm02_sub_fit_PA <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm02_sub_fit[["data"]][["Opsariichthys_uncirostris_amurensis_PA"]])
Opsariichthys_uncirostris_amurensis_glmm02_sub_fit_sub <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm02_sub_fit[["data"]][["substrate"]])
Opsariichthys_uncirostris_amurensis_sub_res_grp1 <- cbind(Opsariichthys_uncirostris_amurensis_glmm02_sub_fit_PA, Opsariichthys_uncirostris_amurensis_glmm02_sub_fit_sub)
colnames(Opsariichthys_uncirostris_amurensis_sub_res_grp1) <- c("prob", "sub")

#'- normalize 0 to 1
Opsariichthys_uncirostris_amurensis_sub_res_grp1 <- Opsariichthys_uncirostris_amurensis_sub_res_grp1 %>%
  mutate(HSI_sub_grp1 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Opsariichthys_uncirostris_amurensis_sub_res_grp1$sub <- factor(Opsariichthys_uncirostris_amurensis_sub_res_grp1$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                                                'pebble', 'cobble', 'boulder'))

#'- plot
Opsariichthys_uncirostris_amurensis_substrate_fig1 <- ggplot(Opsariichthys_uncirostris_amurensis_sub_res_grp1, aes(x = sub, y = HSI_sub_grp1)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Opsariichthys_uncirostris_amurensis Substrate HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Opsariichthys_uncirostris_amurensis_substrate_fig1



##### Opsariichthys_uncirostris_amurensis group2 model #####

#'## select Group2
Opsariichthys_uncirostris_amurensis_dat03 <- Opsariichthys_uncirostris_amurensis_dat %>%
  filter(site_group == "GRP2")

#'##### glm
Opsariichthys_uncirostris_amurensis_glm03 <- glm(Opsariichthys_uncirostris_amurensis_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                                 family=binomial(link = "logit"), data= Opsariichthys_uncirostris_amurensis_dat03, na.action ="na.fail")

summary(Opsariichthys_uncirostris_amurensis_glm03)

vif(Opsariichthys_uncirostris_amurensis_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect 
Opsariichthys_uncirostris_amurensis_glmm03 <- glmmTMB(Opsariichthys_uncirostris_amurensis_PA ~ poly(depth, 1) + poly(velocity, 1) + substrate + 
                                                        (1|site) + (1|season), family=binomial(link = "logit"), data= Opsariichthys_uncirostris_amurensis_dat03, na.action ="na.fail")

summary(Opsariichthys_uncirostris_amurensis_glmm03)

#'- Model validation
r.squaredGLMM(Opsariichthys_uncirostris_amurensis_glmm03)
cor.test(predict(Opsariichthys_uncirostris_amurensis_glmm03, type = "response", allow.new.levels=TRUE), Opsariichthys_uncirostris_amurensis_dat03$Opsariichthys_uncirostris_amurensis_PA)

#'- Opsariichthys_uncirostris_amurensis depth
Opsariichthys_uncirostris_amurensis_glmm03_dep_fit <- effect_plot(Opsariichthys_uncirostris_amurensis_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Opsariichthys_uncirostris_amurensis_glmm03_dep_fit_PA <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Opsariichthys_uncirostris_amurensis_PA"]])
Opsariichthys_uncirostris_amurensis_glmm03_dep_fit_depth <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Opsariichthys_uncirostris_amurensis_dep_res_grp2 <- cbind(Opsariichthys_uncirostris_amurensis_glmm03_dep_fit_PA, Opsariichthys_uncirostris_amurensis_glmm03_dep_fit_depth)
colnames(Opsariichthys_uncirostris_amurensis_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Opsariichthys_uncirostris_amurensis_dep_res_grp2 <- Opsariichthys_uncirostris_amurensis_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Opsariichthys_uncirostris_amurensis_depth_fig2 <- ggplot(Opsariichthys_uncirostris_amurensis_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Opsariichthys_uncirostris_amurensis depth HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Opsariichthys_uncirostris_amurensis_depth_fig2

#'- Opsariichthys_uncirostris_amurensis velocity
Opsariichthys_uncirostris_amurensis_glmm03_vel_fit <- effect_plot(Opsariichthys_uncirostris_amurensis_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Opsariichthys_uncirostris_amurensis_glmm03_vel_fit_PA <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Opsariichthys_uncirostris_amurensis_PA"]])
Opsariichthys_uncirostris_amurensis_glmm03_vel_fit_vel <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Opsariichthys_uncirostris_amurensis_vel_res_grp2 <- cbind(Opsariichthys_uncirostris_amurensis_glmm03_vel_fit_PA, Opsariichthys_uncirostris_amurensis_glmm03_vel_fit_vel)
colnames(Opsariichthys_uncirostris_amurensis_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Opsariichthys_uncirostris_amurensis_vel_res_grp2 <- Opsariichthys_uncirostris_amurensis_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Opsariichthys_uncirostris_amurensis_velocity_fig2 <- ggplot(Opsariichthys_uncirostris_amurensis_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Opsariichthys_uncirostris_amurensis Velocity HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Opsariichthys_uncirostris_amurensis_velocity_fig2

#'- Opsariichthys_uncirostris_amurensis substrate
Opsariichthys_uncirostris_amurensis_glmm03_sub_fit <- effect_plot(Opsariichthys_uncirostris_amurensis_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Opsariichthys_uncirostris_amurensis_glmm03_sub_fit_PA <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm03_sub_fit[["data"]][["Opsariichthys_uncirostris_amurensis_PA"]])
Opsariichthys_uncirostris_amurensis_glmm03_sub_fit_sub <- as.data.frame(Opsariichthys_uncirostris_amurensis_glmm03_sub_fit[["data"]][["substrate"]])
Opsariichthys_uncirostris_amurensis_sub_res_grp2 <- cbind(Opsariichthys_uncirostris_amurensis_glmm03_sub_fit_PA, Opsariichthys_uncirostris_amurensis_glmm03_sub_fit_sub)
colnames(Opsariichthys_uncirostris_amurensis_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Opsariichthys_uncirostris_amurensis_sub_res_grp2 <- Opsariichthys_uncirostris_amurensis_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Opsariichthys_uncirostris_amurensis_sub_res_grp2$sub <- factor(Opsariichthys_uncirostris_amurensis_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                                                'pebble', 'cobble', 'boulder'))

#'- plot
Opsariichthys_uncirostris_amurensis_substrate_fig2 <- ggplot(Opsariichthys_uncirostris_amurensis_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Opsariichthys uncirostris amurensis Substrate HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Opsariichthys_uncirostris_amurensis_substrate_fig2








#'## Combined figures

#'- merge mulitple effect plot
Opsariichthys_uncirostris_amurensis_dep_res_1 <- Opsariichthys_uncirostris_amurensis_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)

Opsariichthys_uncirostris_amurensis_dep_res_grp1_1 <- Opsariichthys_uncirostris_amurensis_dep_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp1)

Opsariichthys_uncirostris_amurensis_dep_res_grp2_1 <- Opsariichthys_uncirostris_amurensis_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)



Opsariichthys_uncirostris_amurensis_dep_res_comb <- rbind(Opsariichthys_uncirostris_amurensis_dep_res_1, 
                                                          Opsariichthys_uncirostris_amurensis_dep_res_grp1_1, Opsariichthys_uncirostris_amurensis_dep_res_grp2_1)


#'- depth plot
Opsariichthys_uncirostris_amurensis_depth_fig4 <- ggplot(Opsariichthys_uncirostris_amurensis_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2.5, 0.5)) +     
  labs(title= "Opsariichthys uncirostris amurensis Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Opsariichthys_uncirostris_amurensis_depth_fig4 <-Opsariichthys_uncirostris_amurensis_depth_fig4+scale_color_manual(values=c("black", "blue", "red"))
Opsariichthys_uncirostris_amurensis_depth_fig4












#'- merge mulitple effect plot
Opsariichthys_uncirostris_amurensis_vel_res_1 <- Opsariichthys_uncirostris_amurensis_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)

Opsariichthys_uncirostris_amurensis_vel_res_grp1_1 <- Opsariichthys_uncirostris_amurensis_vel_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp1)

Opsariichthys_uncirostris_amurensis_vel_res_grp2_1 <- Opsariichthys_uncirostris_amurensis_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)




Opsariichthys_uncirostris_amurensis_vel_res_comb <- rbind(Opsariichthys_uncirostris_amurensis_vel_res_1, 
                                                          Opsariichthys_uncirostris_amurensis_vel_res_grp1_1, Opsariichthys_uncirostris_amurensis_vel_res_grp2_1)


#'- velocity plot
Opsariichthys_uncirostris_amurensis_velocity_fig4 <- ggplot(Opsariichthys_uncirostris_amurensis_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2.5, 0.5)) +     
  labs(title= "Opsariichthys uncirostris amurensis Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Opsariichthys_uncirostris_amurensis_velocity_fig4 <- Opsariichthys_uncirostris_amurensis_velocity_fig4 + scale_color_manual(values=c("black", "blue", "red"))
Opsariichthys_uncirostris_amurensis_velocity_fig4



#'- merge mulitple effect plot
Opsariichthys_uncirostris_amurensis_sub_res_1 <- Opsariichthys_uncirostris_amurensis_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)

Opsariichthys_uncirostris_amurensis_sub_res_grp1_1 <- Opsariichthys_uncirostris_amurensis_sub_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp1)

Opsariichthys_uncirostris_amurensis_sub_res_grp2_1 <- Opsariichthys_uncirostris_amurensis_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)




Opsariichthys_uncirostris_amurensis_sub_res_comb <- rbind(Opsariichthys_uncirostris_amurensis_sub_res_1, 
                                                          Opsariichthys_uncirostris_amurensis_sub_res_grp1_1, Opsariichthys_uncirostris_amurensis_sub_res_grp2_1)

#'- substrate plot
Opsariichthys_uncirostris_amurensis_substrate_fig4 <- ggplot(Opsariichthys_uncirostris_amurensis_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Opsariichthys uncirostris_amurensis Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Opsariichthys_uncirostris_amurensis_substrate_fig4 <- Opsariichthys_uncirostris_amurensis_substrate_fig4 + scale_fill_manual(values=c("black", "blue", "red"))
Opsariichthys_uncirostris_amurensis_substrate_fig4


#'## combine figures

ggarrange(Opsariichthys_uncirostris_amurensis_depth_fig4, Opsariichthys_uncirostris_amurensis_velocity_fig4, Opsariichthys_uncirostris_amurensis_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 




##### Pseudogobio_esocinus model #####

#'## (2) Pseudogobio_esocinus model

#'### data management
Pseudogobio_esocinus_dat <- HSI_DB2 %>%
  mutate(Pseudogobio_esocinus_PA = ifelse(species == "Pseudogobio_esocinus", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Pseudogobio_esocinus_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Pseudogobio_esocinus_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Pseudogobio_esocinus_PA, depth, velocity, substrate, MBSNCD) 

head(Pseudogobio_esocinus_dat)


#'## remove site where species did not occur
Pseudogobio_esocinus_Pre <- HSI_DB2 %>%
  mutate(Pseudogobio_esocinus_PA = ifelse(species == "Pseudogobio_esocinus", 1, 0)) %>%
  aggregate(Pseudogobio_esocinus_PA ~ site, "sum") %>%
  mutate(Pseudogobio_esocinus_Pre = if_else(Pseudogobio_esocinus_PA > 0, "Yes", "No")) %>%
  select(site, Pseudogobio_esocinus_Pre)

Pseudogobio_esocinus_dat <- Pseudogobio_esocinus_dat %>%
  left_join(Pseudogobio_esocinus_Pre, by = "site") %>%
  filter(Pseudogobio_esocinus_Pre == "Yes") %>%
  as_tibble()



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

#'##### glm
Pseudogobio_esocinus_glm01 <- glm(Pseudogobio_esocinus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                  family=binomial, data= Pseudogobio_esocinus_dat, na.action ="na.fail")

summary(Pseudogobio_esocinus_glm01)

vif(Pseudogobio_esocinus_glm01)



#'##### glmm
#'- combination of site and season are used as a random effect 
Pseudogobio_esocinus_glmm01 <- glmmTMB(Pseudogobio_esocinus_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                         (1|site) + (1|season), family=binomial, data= Pseudogobio_esocinus_dat, na.action ="na.fail")

summary(Pseudogobio_esocinus_glmm01)




#'- Model validation
r.squaredGLMM(Pseudogobio_esocinus_glmm01)
cor.test(predict(Pseudogobio_esocinus_glmm01, type = "response", allow.new.levels=TRUE), Pseudogobio_esocinus_dat$Pseudogobio_esocinus_PA)

#'- Pseudogobio_esocinus depth
Pseudogobio_esocinus_glmm01_dep_fit <- effect_plot(Pseudogobio_esocinus_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Pseudogobio_esocinus_glmm01_dep_fit_PA <- as.data.frame(Pseudogobio_esocinus_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Pseudogobio_esocinus_PA"]])
Pseudogobio_esocinus_glmm01_dep_fit_dep <- as.data.frame(Pseudogobio_esocinus_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Pseudogobio_esocinus_dep_res <- cbind(Pseudogobio_esocinus_glmm01_dep_fit_PA, Pseudogobio_esocinus_glmm01_dep_fit_dep)
colnames(Pseudogobio_esocinus_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Pseudogobio_esocinus_dep_res <- Pseudogobio_esocinus_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Pseudogobio_esocinus_depth_fig <- ggplot(Pseudogobio_esocinus_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudogobio_esocinus depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudogobio_esocinus_depth_fig

#'- Pseudogobio_esocinus velocity
Pseudogobio_esocinus_glmm01_vel_fit <- effect_plot(Pseudogobio_esocinus_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Pseudogobio_esocinus_glmm01_vel_fit_PA <- as.data.frame(Pseudogobio_esocinus_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Pseudogobio_esocinus_PA"]])
Pseudogobio_esocinus_glmm01_vel_fit_vel <- as.data.frame(Pseudogobio_esocinus_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Pseudogobio_esocinus_vel_res <- cbind(Pseudogobio_esocinus_glmm01_vel_fit_PA, Pseudogobio_esocinus_glmm01_vel_fit_vel)
colnames(Pseudogobio_esocinus_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Pseudogobio_esocinus_vel_res <- Pseudogobio_esocinus_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Pseudogobio_esocinus_velocity_fig <- ggplot(Pseudogobio_esocinus_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudogobio_esocinus Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudogobio_esocinus_velocity_fig

#'- Pseudogobio_esocinus substrate
Pseudogobio_esocinus_glmm01_sub_fit <- effect_plot(Pseudogobio_esocinus_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Pseudogobio_esocinus_glmm01_sub_fit_PA <- as.data.frame(Pseudogobio_esocinus_glmm01_sub_fit[["data"]][["Pseudogobio_esocinus_PA"]])
Pseudogobio_esocinus_glmm01_sub_fit_sub <- as.data.frame(Pseudogobio_esocinus_glmm01_sub_fit[["data"]][["substrate"]])
Pseudogobio_esocinus_sub_res <- cbind(Pseudogobio_esocinus_glmm01_sub_fit_PA, Pseudogobio_esocinus_glmm01_sub_fit_sub)
colnames(Pseudogobio_esocinus_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Pseudogobio_esocinus_sub_res <- Pseudogobio_esocinus_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Pseudogobio_esocinus_sub_res$sub <- factor(Pseudogobio_esocinus_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                        'pebble', 'cobble', 'boulder'))

#'- plot
Pseudogobio_esocinus_substrate_fig <- ggplot(Pseudogobio_esocinus_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Pseudogobio_esocinus Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudogobio_esocinus_substrate_fig


#'## (1) Pseudogobio_esocinus model
##### Pseudogobio_esocinus group1 model #####

#'## select Group1

Pseudogobio_esocinus_dat02 <- Pseudogobio_esocinus_dat %>%
  filter(site_group == "GRP1")

#'##### glm
Pseudogobio_esocinus_glm02 <- glm(Pseudogobio_esocinus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                  family=binomial(link = "logit"), data= Pseudogobio_esocinus_dat02, na.action ="na.fail")

summary(Pseudogobio_esocinus_glm02)

vif(Pseudogobio_esocinus_glm02)
#'##### glmm
#'- combination of site and season are used as a random effect ``
Pseudogobio_esocinus_glmm02 <- glmmTMB(Pseudogobio_esocinus_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                         (1|site) + (1|season), family=binomial(link = "logit"), data= Pseudogobio_esocinus_dat02, na.action ="na.fail")

summary(Pseudogobio_esocinus_glmm02)



#'- Model validation
r.squaredGLMM(Pseudogobio_esocinus_glmm02)
cor.test(predict(Pseudogobio_esocinus_glmm02, type = "response", allow.new.levels=TRUE), Pseudogobio_esocinus_dat02$Pseudogobio_esocinus_PA)

#'- Pseudogobio_esocinus depth
Pseudogobio_esocinus_glmm02_dep_fit <- effect_plot(Pseudogobio_esocinus_glmm02, pred = depth, interval = TRUE, plot.points = FALSE)
Pseudogobio_esocinus_glmm02_dep_fit_PA <- as.data.frame(Pseudogobio_esocinus_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["Pseudogobio_esocinus_PA"]])
Pseudogobio_esocinus_glmm02_dep_fit_depth <- as.data.frame(Pseudogobio_esocinus_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Pseudogobio_esocinus_dep_res_grp1 <- cbind(Pseudogobio_esocinus_glmm02_dep_fit_PA, Pseudogobio_esocinus_glmm02_dep_fit_depth)
colnames(Pseudogobio_esocinus_dep_res_grp1) <- c("prob", "dep")

#'- normalize 0 to 1
Pseudogobio_esocinus_dep_res_grp1 <- Pseudogobio_esocinus_dep_res_grp1 %>%
  mutate(HSI_dep_grp1 = scale_values(prob))

#'- plot
Pseudogobio_esocinus_depth_fig1 <- ggplot(Pseudogobio_esocinus_dep_res_grp1, aes(x = dep, y = HSI_dep_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudogobio_esocinus depth HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudogobio_esocinus_depth_fig1

#'- Pseudogobio_esocinus velocity
Pseudogobio_esocinus_glmm02_vel_fit <- effect_plot(Pseudogobio_esocinus_glmm02, pred = velocity, interval = TRUE, plot.points = FALSE)
Pseudogobio_esocinus_glmm02_vel_fit_PA <- as.data.frame(Pseudogobio_esocinus_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["Pseudogobio_esocinus_PA"]])
Pseudogobio_esocinus_glmm02_vel_fit_vel <- as.data.frame(Pseudogobio_esocinus_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Pseudogobio_esocinus_vel_res_grp1 <- cbind(Pseudogobio_esocinus_glmm02_vel_fit_PA, Pseudogobio_esocinus_glmm02_vel_fit_vel)
colnames(Pseudogobio_esocinus_vel_res_grp1) <- c("prob", "vel")

#'- normalize 0 to 1
Pseudogobio_esocinus_vel_res_grp1 <- Pseudogobio_esocinus_vel_res_grp1 %>%
  mutate(HSI_vel_grp1 = scale_values(prob))

#'- plot
Pseudogobio_esocinus_velocity_fig1 <- ggplot(Pseudogobio_esocinus_vel_res_grp1, aes(x = vel, y = HSI_vel_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudogobio_esocinus Velocity HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudogobio_esocinus_velocity_fig1

#'- Pseudogobio_esocinus substrate
Pseudogobio_esocinus_glmm02_sub_fit <- effect_plot(Pseudogobio_esocinus_glmm02, pred = substrate, interval = FALSE, plot.points = FALSE)
Pseudogobio_esocinus_glmm02_sub_fit_PA <- as.data.frame(Pseudogobio_esocinus_glmm02_sub_fit[["data"]][["Pseudogobio_esocinus_PA"]])
Pseudogobio_esocinus_glmm02_sub_fit_sub <- as.data.frame(Pseudogobio_esocinus_glmm02_sub_fit[["data"]][["substrate"]])
Pseudogobio_esocinus_sub_res_grp1 <- cbind(Pseudogobio_esocinus_glmm02_sub_fit_PA, Pseudogobio_esocinus_glmm02_sub_fit_sub)
colnames(Pseudogobio_esocinus_sub_res_grp1) <- c("prob", "sub")

#'- normalize 0 to 1
Pseudogobio_esocinus_sub_res_grp1 <- Pseudogobio_esocinus_sub_res_grp1 %>%
  mutate(HSI_sub_grp1 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Pseudogobio_esocinus_sub_res_grp1$sub <- factor(Pseudogobio_esocinus_sub_res_grp1$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                  'pebble', 'cobble', 'boulder'))

#'- plot
Pseudogobio_esocinus_substrate_fig1 <- ggplot(Pseudogobio_esocinus_sub_res_grp1, aes(x = sub, y = HSI_sub_grp1)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Pseudogobio_esocinus Substrate HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudogobio_esocinus_substrate_fig1



##### Pseudogobio_esocinus group2 model #####

#'## select Group2
Pseudogobio_esocinus_dat03 <- Pseudogobio_esocinus_dat %>%
  filter(site_group == "GRP2")

#'##### glm
Pseudogobio_esocinus_glm03 <- glm(Pseudogobio_esocinus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                  family=binomial(link = "logit"), data= Pseudogobio_esocinus_dat03, na.action ="na.fail")

summary(Pseudogobio_esocinus_glm03)

vif(Pseudogobio_esocinus_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect 
Pseudogobio_esocinus_glmm03 <- glmmTMB(Pseudogobio_esocinus_PA ~ poly(depth, 1) + poly(velocity, 1) + substrate + 
                                         (1|site) + (1|season), family=binomial(link = "logit"), data= Pseudogobio_esocinus_dat03, na.action ="na.fail")

summary(Pseudogobio_esocinus_glmm03)

#'- Model validation
r.squaredGLMM(Pseudogobio_esocinus_glmm03)
cor.test(predict(Pseudogobio_esocinus_glmm03, type = "response", allow.new.levels=TRUE), Pseudogobio_esocinus_dat03$Pseudogobio_esocinus_PA)

#'- Pseudogobio_esocinus depth
Pseudogobio_esocinus_glmm03_dep_fit <- effect_plot(Pseudogobio_esocinus_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Pseudogobio_esocinus_glmm03_dep_fit_PA <- as.data.frame(Pseudogobio_esocinus_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Pseudogobio_esocinus_PA"]])
Pseudogobio_esocinus_glmm03_dep_fit_depth <- as.data.frame(Pseudogobio_esocinus_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Pseudogobio_esocinus_dep_res_grp2 <- cbind(Pseudogobio_esocinus_glmm03_dep_fit_PA, Pseudogobio_esocinus_glmm03_dep_fit_depth)
colnames(Pseudogobio_esocinus_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Pseudogobio_esocinus_dep_res_grp2 <- Pseudogobio_esocinus_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Pseudogobio_esocinus_depth_fig2 <- ggplot(Pseudogobio_esocinus_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudogobio_esocinus depth HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudogobio_esocinus_depth_fig2

#'- Pseudogobio_esocinus velocity
Pseudogobio_esocinus_glmm03_vel_fit <- effect_plot(Pseudogobio_esocinus_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Pseudogobio_esocinus_glmm03_vel_fit_PA <- as.data.frame(Pseudogobio_esocinus_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Pseudogobio_esocinus_PA"]])
Pseudogobio_esocinus_glmm03_vel_fit_vel <- as.data.frame(Pseudogobio_esocinus_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Pseudogobio_esocinus_vel_res_grp2 <- cbind(Pseudogobio_esocinus_glmm03_vel_fit_PA, Pseudogobio_esocinus_glmm03_vel_fit_vel)
colnames(Pseudogobio_esocinus_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Pseudogobio_esocinus_vel_res_grp2 <- Pseudogobio_esocinus_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Pseudogobio_esocinus_velocity_fig2 <- ggplot(Pseudogobio_esocinus_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudogobio_esocinus Velocity HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudogobio_esocinus_velocity_fig2

#'- Pseudogobio_esocinus substrate
Pseudogobio_esocinus_glmm03_sub_fit <- effect_plot(Pseudogobio_esocinus_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Pseudogobio_esocinus_glmm03_sub_fit_PA <- as.data.frame(Pseudogobio_esocinus_glmm03_sub_fit[["data"]][["Pseudogobio_esocinus_PA"]])
Pseudogobio_esocinus_glmm03_sub_fit_sub <- as.data.frame(Pseudogobio_esocinus_glmm03_sub_fit[["data"]][["substrate"]])
Pseudogobio_esocinus_sub_res_grp2 <- cbind(Pseudogobio_esocinus_glmm03_sub_fit_PA, Pseudogobio_esocinus_glmm03_sub_fit_sub)
colnames(Pseudogobio_esocinus_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Pseudogobio_esocinus_sub_res_grp2 <- Pseudogobio_esocinus_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Pseudogobio_esocinus_sub_res_grp2$sub <- factor(Pseudogobio_esocinus_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                  'pebble', 'cobble', 'boulder'))

#'- plot
Pseudogobio_esocinus_substrate_fig2 <- ggplot(Pseudogobio_esocinus_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Pseudogobio_esocinus Substrate HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Pseudogobio_esocinus_substrate_fig2



#'## Combined figures

#'- merge mulitple effect plot
Pseudogobio_esocinus_dep_res_1 <- Pseudogobio_esocinus_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)

Pseudogobio_esocinus_dep_res_grp1_1 <- Pseudogobio_esocinus_dep_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp1)

Pseudogobio_esocinus_dep_res_grp2_1 <- Pseudogobio_esocinus_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)


Pseudogobio_esocinus_dep_res_comb <- rbind(Pseudogobio_esocinus_dep_res_1, Pseudogobio_esocinus_dep_res_grp1_1,
                                           Pseudogobio_esocinus_dep_res_grp2_1)


#'- depth plot
Pseudogobio_esocinus_depth_fig4 <- ggplot(Pseudogobio_esocinus_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudogobio_esocinus Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Pseudogobio_esocinus_depth_fig4 <-Pseudogobio_esocinus_depth_fig4+scale_color_manual(values=c("black","blue3", "red", "green3"))
Pseudogobio_esocinus_depth_fig4




#'- merge mulitple effect plot
Pseudogobio_esocinus_vel_res_1 <- Pseudogobio_esocinus_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)

Pseudogobio_esocinus_vel_res_grp1_1 <- Pseudogobio_esocinus_vel_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp1)

Pseudogobio_esocinus_vel_res_grp2_1 <- Pseudogobio_esocinus_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)


Pseudogobio_esocinus_vel_res_comb <- rbind(Pseudogobio_esocinus_vel_res_1, Pseudogobio_esocinus_vel_res_grp1_1,
                                           Pseudogobio_esocinus_vel_res_grp2_1)


#'- velocity plot
Pseudogobio_esocinus_velocity_fig4 <- ggplot(Pseudogobio_esocinus_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Pseudogobio_esocinus Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Pseudogobio_esocinus_velocity_fig4 <- Pseudogobio_esocinus_velocity_fig4 + scale_color_manual(values=c("black","blue3", "red", "green3", "orange3"))
Pseudogobio_esocinus_velocity_fig4



#'- merge mulitple effect plot
Pseudogobio_esocinus_sub_res_1 <- Pseudogobio_esocinus_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)

Pseudogobio_esocinus_sub_res_grp1_1 <- Pseudogobio_esocinus_sub_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp1)

Pseudogobio_esocinus_sub_res_grp2_1 <- Pseudogobio_esocinus_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)

Pseudogobio_esocinus_sub_res_comb <- rbind(Pseudogobio_esocinus_sub_res_1, Pseudogobio_esocinus_sub_res_grp1_1,
                                           Pseudogobio_esocinus_sub_res_grp2_1)

#'- substrate plot
Pseudogobio_esocinus_substrate_fig4 <- ggplot(Pseudogobio_esocinus_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Pseudogobio esocinus Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Pseudogobio_esocinus_substrate_fig4 <- Pseudogobio_esocinus_substrate_fig4 + scale_fill_manual(values=c("black","blue3", "red", "green3"))
Pseudogobio_esocinus_substrate_fig4


#'## combine figures

ggarrange(Pseudogobio_esocinus_depth_fig4, Pseudogobio_esocinus_velocity_fig4, Pseudogobio_esocinus_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 







##### Microphysogobio_yaluensis model #####

#'## (2) Microphysogobio_yaluensis model

#'### data management
Microphysogobio_yaluensis_dat <- HSI_DB2 %>%
  mutate(Microphysogobio_yaluensis_PA = ifelse(species == "Microphysogobio_yaluensis", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Microphysogobio_yaluensis_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Microphysogobio_yaluensis_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Microphysogobio_yaluensis_PA, depth, velocity, substrate, MBSNCD) 

head(Microphysogobio_yaluensis_dat)


#'## remove site where species did not occur
Microphysogobio_yaluensis_Pre <- HSI_DB2 %>%
  mutate(Microphysogobio_yaluensis_PA = ifelse(species == "Microphysogobio_yaluensis", 1, 0)) %>%
  aggregate(Microphysogobio_yaluensis_PA ~ site, "sum") %>%
  mutate(Microphysogobio_yaluensis_Pre = if_else(Microphysogobio_yaluensis_PA > 0, "Yes", "No")) %>%
  select(site, Microphysogobio_yaluensis_Pre)

Microphysogobio_yaluensis_dat <- Microphysogobio_yaluensis_dat %>%
  left_join(Microphysogobio_yaluensis_Pre, by = "site") %>%
  filter(Microphysogobio_yaluensis_Pre == "Yes") %>%
  as_tibble()


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

#'##### glm
Microphysogobio_yaluensis_glm01 <- glm(Microphysogobio_yaluensis_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                       family=binomial, data= Microphysogobio_yaluensis_dat, na.action ="na.fail")

summary(Microphysogobio_yaluensis_glm01)

vif(Microphysogobio_yaluensis_glm01)



#'##### glmm
#'- combination of site and season are used as a random effect 
Microphysogobio_yaluensis_glmm01 <- glmmTMB(Microphysogobio_yaluensis_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                              (1|site) + (1|season), family=binomial, data= Microphysogobio_yaluensis_dat, na.action ="na.fail")

summary(Microphysogobio_yaluensis_glmm01)




#'- Model validation
r.squaredGLMM(Microphysogobio_yaluensis_glmm01)
cor.test(predict(Microphysogobio_yaluensis_glmm01, type = "response", allow.new.levels=TRUE), Microphysogobio_yaluensis_dat$Microphysogobio_yaluensis_PA)

#'- Microphysogobio_yaluensis depth
Microphysogobio_yaluensis_glmm01_dep_fit <- effect_plot(Microphysogobio_yaluensis_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Microphysogobio_yaluensis_glmm01_dep_fit_PA <- as.data.frame(Microphysogobio_yaluensis_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Microphysogobio_yaluensis_PA"]])
Microphysogobio_yaluensis_glmm01_dep_fit_dep <- as.data.frame(Microphysogobio_yaluensis_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Microphysogobio_yaluensis_dep_res <- cbind(Microphysogobio_yaluensis_glmm01_dep_fit_PA, Microphysogobio_yaluensis_glmm01_dep_fit_dep)
colnames(Microphysogobio_yaluensis_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Microphysogobio_yaluensis_dep_res <- Microphysogobio_yaluensis_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Microphysogobio_yaluensis_depth_fig <- ggplot(Microphysogobio_yaluensis_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Microphysogobio_yaluensis depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Microphysogobio_yaluensis_depth_fig

#'- Microphysogobio_yaluensis velocity
Microphysogobio_yaluensis_glmm01_vel_fit <- effect_plot(Microphysogobio_yaluensis_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Microphysogobio_yaluensis_glmm01_vel_fit_PA <- as.data.frame(Microphysogobio_yaluensis_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Microphysogobio_yaluensis_PA"]])
Microphysogobio_yaluensis_glmm01_vel_fit_vel <- as.data.frame(Microphysogobio_yaluensis_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Microphysogobio_yaluensis_vel_res <- cbind(Microphysogobio_yaluensis_glmm01_vel_fit_PA, Microphysogobio_yaluensis_glmm01_vel_fit_vel)
colnames(Microphysogobio_yaluensis_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Microphysogobio_yaluensis_vel_res <- Microphysogobio_yaluensis_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Microphysogobio_yaluensis_velocity_fig <- ggplot(Microphysogobio_yaluensis_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Microphysogobio_yaluensis Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Microphysogobio_yaluensis_velocity_fig

#'- Microphysogobio_yaluensis substrate
Microphysogobio_yaluensis_glmm01_sub_fit <- effect_plot(Microphysogobio_yaluensis_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Microphysogobio_yaluensis_glmm01_sub_fit_PA <- as.data.frame(Microphysogobio_yaluensis_glmm01_sub_fit[["data"]][["Microphysogobio_yaluensis_PA"]])
Microphysogobio_yaluensis_glmm01_sub_fit_sub <- as.data.frame(Microphysogobio_yaluensis_glmm01_sub_fit[["data"]][["substrate"]])
Microphysogobio_yaluensis_sub_res <- cbind(Microphysogobio_yaluensis_glmm01_sub_fit_PA, Microphysogobio_yaluensis_glmm01_sub_fit_sub)
colnames(Microphysogobio_yaluensis_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Microphysogobio_yaluensis_sub_res <- Microphysogobio_yaluensis_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Microphysogobio_yaluensis_sub_res$sub <- factor(Microphysogobio_yaluensis_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                  'pebble', 'cobble', 'boulder'))

#'- plot
Microphysogobio_yaluensis_substrate_fig <- ggplot(Microphysogobio_yaluensis_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Microphysogobio_yaluensis Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Microphysogobio_yaluensis_substrate_fig



##### Microphysogobio_yaluensis group2 model #####

#'## select Group2
Microphysogobio_yaluensis_dat03 <- Microphysogobio_yaluensis_dat %>%
  filter(site_group == "GRP2")

#'##### glm
Microphysogobio_yaluensis_glm03 <- glm(Microphysogobio_yaluensis_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                       family=binomial(link = "logit"), data= Microphysogobio_yaluensis_dat03, na.action ="na.fail")

summary(Microphysogobio_yaluensis_glm03)

vif(Microphysogobio_yaluensis_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect 
Microphysogobio_yaluensis_glmm03 <- glmmTMB(Microphysogobio_yaluensis_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                              (1|site) + (1|season), family=binomial(link = "logit"), data= Microphysogobio_yaluensis_dat03, na.action ="na.fail")

summary(Microphysogobio_yaluensis_glmm03)

#'- Model validation
r.squaredGLMM(Microphysogobio_yaluensis_glmm03)
cor.test(predict(Microphysogobio_yaluensis_glmm03, type = "response", allow.new.levels=TRUE), Microphysogobio_yaluensis_dat03$Microphysogobio_yaluensis_PA)

#'- Microphysogobio_yaluensis depth
Microphysogobio_yaluensis_glmm03_dep_fit <- effect_plot(Microphysogobio_yaluensis_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Microphysogobio_yaluensis_glmm03_dep_fit_PA <- as.data.frame(Microphysogobio_yaluensis_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Microphysogobio_yaluensis_PA"]])
Microphysogobio_yaluensis_glmm03_dep_fit_depth <- as.data.frame(Microphysogobio_yaluensis_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Microphysogobio_yaluensis_dep_res_grp2 <- cbind(Microphysogobio_yaluensis_glmm03_dep_fit_PA, Microphysogobio_yaluensis_glmm03_dep_fit_depth)
colnames(Microphysogobio_yaluensis_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Microphysogobio_yaluensis_dep_res_grp2 <- Microphysogobio_yaluensis_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Microphysogobio_yaluensis_depth_fig2 <- ggplot(Microphysogobio_yaluensis_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Microphysogobio_yaluensis depth HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Microphysogobio_yaluensis_depth_fig2

#'- Microphysogobio_yaluensis velocity
Microphysogobio_yaluensis_glmm03_vel_fit <- effect_plot(Microphysogobio_yaluensis_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Microphysogobio_yaluensis_glmm03_vel_fit_PA <- as.data.frame(Microphysogobio_yaluensis_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Microphysogobio_yaluensis_PA"]])
Microphysogobio_yaluensis_glmm03_vel_fit_vel <- as.data.frame(Microphysogobio_yaluensis_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Microphysogobio_yaluensis_vel_res_grp2 <- cbind(Microphysogobio_yaluensis_glmm03_vel_fit_PA, Microphysogobio_yaluensis_glmm03_vel_fit_vel)
colnames(Microphysogobio_yaluensis_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Microphysogobio_yaluensis_vel_res_grp2 <- Microphysogobio_yaluensis_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Microphysogobio_yaluensis_velocity_fig2 <- ggplot(Microphysogobio_yaluensis_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Microphysogobio_yaluensis Velocity HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Microphysogobio_yaluensis_velocity_fig2

#'- Microphysogobio_yaluensis substrate
Microphysogobio_yaluensis_glmm03_sub_fit <- effect_plot(Microphysogobio_yaluensis_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Microphysogobio_yaluensis_glmm03_sub_fit_PA <- as.data.frame(Microphysogobio_yaluensis_glmm03_sub_fit[["data"]][["Microphysogobio_yaluensis_PA"]])
Microphysogobio_yaluensis_glmm03_sub_fit_sub <- as.data.frame(Microphysogobio_yaluensis_glmm03_sub_fit[["data"]][["substrate"]])
Microphysogobio_yaluensis_sub_res_grp2 <- cbind(Microphysogobio_yaluensis_glmm03_sub_fit_PA, Microphysogobio_yaluensis_glmm03_sub_fit_sub)
colnames(Microphysogobio_yaluensis_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Microphysogobio_yaluensis_sub_res_grp2 <- Microphysogobio_yaluensis_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Microphysogobio_yaluensis_sub_res_grp2$sub <- factor(Microphysogobio_yaluensis_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                            'pebble', 'cobble', 'boulder'))

#'- plot
Microphysogobio_yaluensis_substrate_fig2 <- ggplot(Microphysogobio_yaluensis_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Microphysogobio_yaluensis Substrate HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Microphysogobio_yaluensis_substrate_fig2






#'## Combined figures

#'- merge mulitple effect plot
Microphysogobio_yaluensis_dep_res_1 <- Microphysogobio_yaluensis_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)


Microphysogobio_yaluensis_dep_res_grp2_1 <- Microphysogobio_yaluensis_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)



Microphysogobio_yaluensis_dep_res_comb <- rbind(Microphysogobio_yaluensis_dep_res_1,
                                                Microphysogobio_yaluensis_dep_res_grp2_1)


#'- depth plot
Microphysogobio_yaluensis_depth_fig4 <- ggplot(Microphysogobio_yaluensis_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Microphysogobio_yaluensis Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Microphysogobio_yaluensis_depth_fig4 <-Microphysogobio_yaluensis_depth_fig4+scale_color_manual(values=c("black", "red"))
Microphysogobio_yaluensis_depth_fig4




#'- merge mulitple effect plot
Microphysogobio_yaluensis_vel_res_1 <- Microphysogobio_yaluensis_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)

Microphysogobio_yaluensis_vel_res_grp2_1 <- Microphysogobio_yaluensis_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)

Microphysogobio_yaluensis_vel_res_comb <- rbind(Microphysogobio_yaluensis_vel_res_1, 
                                                Microphysogobio_yaluensis_vel_res_grp2_1)


#'- velocity plot
Microphysogobio_yaluensis_velocity_fig4 <- ggplot(Microphysogobio_yaluensis_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Microphysogobio_yaluensis Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Microphysogobio_yaluensis_velocity_fig4 <- Microphysogobio_yaluensis_velocity_fig4 + scale_color_manual(values=c("black", "red"))
Microphysogobio_yaluensis_velocity_fig4



#'- merge mulitple effect plot
Microphysogobio_yaluensis_sub_res_1 <- Microphysogobio_yaluensis_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)

Microphysogobio_yaluensis_sub_res_grp2_1 <- Microphysogobio_yaluensis_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)


Microphysogobio_yaluensis_sub_res_comb <- rbind(Microphysogobio_yaluensis_sub_res_1, 
                                                Microphysogobio_yaluensis_sub_res_grp2_1)

#'- substrate plot
Microphysogobio_yaluensis_substrate_fig4 <- ggplot(Microphysogobio_yaluensis_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Microphysogobio_yaluensis Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Microphysogobio_yaluensis_substrate_fig4 <- Microphysogobio_yaluensis_substrate_fig4 + scale_fill_manual(values=c("black", "red"))
Microphysogobio_yaluensis_substrate_fig4


#'## combine figures

ggarrange(Microphysogobio_yaluensis_depth_fig4, Microphysogobio_yaluensis_velocity_fig4, Microphysogobio_yaluensis_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 






##### Rhinogobius_brunneus model #####

#'## (2) Rhinogobius_brunneus model

#'### data management
Rhinogobius_brunneus_dat <- HSI_DB2 %>%
  mutate(Rhinogobius_brunneus_PA = ifelse(species == "Rhinogobius_brunneus", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Rhinogobius_brunneus_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Rhinogobius_brunneus_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Rhinogobius_brunneus_PA, depth, velocity, substrate, MBSNCD) 

head(Rhinogobius_brunneus_dat)


#'## remove site where species did not occur
Rhinogobius_brunneus_Pre <- HSI_DB2 %>%
  mutate(Rhinogobius_brunneus_PA = ifelse(species == "Rhinogobius_brunneus", 1, 0)) %>%
  aggregate(Rhinogobius_brunneus_PA ~ site, "sum") %>%
  mutate(Rhinogobius_brunneus_Pre = if_else(Rhinogobius_brunneus_PA > 0, "Yes", "No")) %>%
  select(site, Rhinogobius_brunneus_Pre)

Rhinogobius_brunneus_dat <- Rhinogobius_brunneus_dat %>%
  left_join(Rhinogobius_brunneus_Pre, by = "site") %>%
  filter(Rhinogobius_brunneus_Pre == "Yes") %>%
  as_tibble()




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

#'##### glm
Rhinogobius_brunneus_glm01 <- glm(Rhinogobius_brunneus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                  family=binomial, data= Rhinogobius_brunneus_dat, na.action ="na.fail")

summary(Rhinogobius_brunneus_glm01)

vif(Rhinogobius_brunneus_glm01)



#'##### glmm
#'- combination of site and season are used as a random effect 
Rhinogobius_brunneus_glmm01 <- glmmTMB(Rhinogobius_brunneus_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                         (1|site) + (1|season), family=binomial, data= Rhinogobius_brunneus_dat, na.action ="na.fail")

summary(Rhinogobius_brunneus_glmm01)




#'- Model validation
r.squaredGLMM(Rhinogobius_brunneus_glmm01)
cor.test(predict(Rhinogobius_brunneus_glmm01, type = "response", allow.new.levels=TRUE), Rhinogobius_brunneus_dat$Rhinogobius_brunneus_PA)

#'- Rhinogobius_brunneus depth
Rhinogobius_brunneus_glmm01_dep_fit <- effect_plot(Rhinogobius_brunneus_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Rhinogobius_brunneus_glmm01_dep_fit_PA <- as.data.frame(Rhinogobius_brunneus_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Rhinogobius_brunneus_PA"]])
Rhinogobius_brunneus_glmm01_dep_fit_dep <- as.data.frame(Rhinogobius_brunneus_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Rhinogobius_brunneus_dep_res <- cbind(Rhinogobius_brunneus_glmm01_dep_fit_PA, Rhinogobius_brunneus_glmm01_dep_fit_dep)
colnames(Rhinogobius_brunneus_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Rhinogobius_brunneus_dep_res <- Rhinogobius_brunneus_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Rhinogobius_brunneus_depth_fig <- ggplot(Rhinogobius_brunneus_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Rhinogobius_brunneus depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Rhinogobius_brunneus_depth_fig

#'- Rhinogobius_brunneus velocity
Rhinogobius_brunneus_glmm01_vel_fit <- effect_plot(Rhinogobius_brunneus_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Rhinogobius_brunneus_glmm01_vel_fit_PA <- as.data.frame(Rhinogobius_brunneus_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Rhinogobius_brunneus_PA"]])
Rhinogobius_brunneus_glmm01_vel_fit_vel <- as.data.frame(Rhinogobius_brunneus_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Rhinogobius_brunneus_vel_res <- cbind(Rhinogobius_brunneus_glmm01_vel_fit_PA, Rhinogobius_brunneus_glmm01_vel_fit_vel)
colnames(Rhinogobius_brunneus_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Rhinogobius_brunneus_vel_res <- Rhinogobius_brunneus_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Rhinogobius_brunneus_velocity_fig <- ggplot(Rhinogobius_brunneus_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Rhinogobius_brunneus Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Rhinogobius_brunneus_velocity_fig

#'- Rhinogobius_brunneus substrate
Rhinogobius_brunneus_glmm01_sub_fit <- effect_plot(Rhinogobius_brunneus_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Rhinogobius_brunneus_glmm01_sub_fit_PA <- as.data.frame(Rhinogobius_brunneus_glmm01_sub_fit[["data"]][["Rhinogobius_brunneus_PA"]])
Rhinogobius_brunneus_glmm01_sub_fit_sub <- as.data.frame(Rhinogobius_brunneus_glmm01_sub_fit[["data"]][["substrate"]])
Rhinogobius_brunneus_sub_res <- cbind(Rhinogobius_brunneus_glmm01_sub_fit_PA, Rhinogobius_brunneus_glmm01_sub_fit_sub)
colnames(Rhinogobius_brunneus_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Rhinogobius_brunneus_sub_res <- Rhinogobius_brunneus_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Rhinogobius_brunneus_sub_res$sub <- factor(Rhinogobius_brunneus_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                        'pebble', 'cobble', 'boulder'))

#'- plot
Rhinogobius_brunneus_substrate_fig <- ggplot(Rhinogobius_brunneus_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Rhinogobius_brunneus Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Rhinogobius_brunneus_substrate_fig


#'## (1) Rhinogobius_brunneus model
##### Rhinogobius_brunneus group1 model #####

#'## select Group1

Rhinogobius_brunneus_dat02 <- Rhinogobius_brunneus_dat %>%
  filter(site_group == "GRP1")

#'##### glm
Rhinogobius_brunneus_glm02 <- glm(Rhinogobius_brunneus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                  family=binomial(link = "logit"), data= Rhinogobius_brunneus_dat02, na.action ="na.fail")

summary(Rhinogobius_brunneus_glm02)

vif(Rhinogobius_brunneus_glm02)
#'##### glmm
#'- combination of site and season are used as a random effect ``
Rhinogobius_brunneus_glmm02 <- glmmTMB(Rhinogobius_brunneus_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                         (1|site) + (1|season), family=binomial(link = "logit"), data= Rhinogobius_brunneus_dat02, na.action ="na.fail")

summary(Rhinogobius_brunneus_glmm02)



#'- Model validation
r.squaredGLMM(Rhinogobius_brunneus_glmm02)
cor.test(predict(Rhinogobius_brunneus_glmm02, type = "response", allow.new.levels=TRUE), Rhinogobius_brunneus_dat02$Rhinogobius_brunneus_PA)

#'- Rhinogobius_brunneus depth
Rhinogobius_brunneus_glmm02_dep_fit <- effect_plot(Rhinogobius_brunneus_glmm02, pred = depth, interval = TRUE, plot.points = FALSE)
Rhinogobius_brunneus_glmm02_dep_fit_PA <- as.data.frame(Rhinogobius_brunneus_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["Rhinogobius_brunneus_PA"]])
Rhinogobius_brunneus_glmm02_dep_fit_depth <- as.data.frame(Rhinogobius_brunneus_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Rhinogobius_brunneus_dep_res_grp1 <- cbind(Rhinogobius_brunneus_glmm02_dep_fit_PA, Rhinogobius_brunneus_glmm02_dep_fit_depth)
colnames(Rhinogobius_brunneus_dep_res_grp1) <- c("prob", "dep")

#'- normalize 0 to 1
Rhinogobius_brunneus_dep_res_grp1 <- Rhinogobius_brunneus_dep_res_grp1 %>%
  mutate(HSI_dep_grp1 = scale_values(prob))

#'- plot
Rhinogobius_brunneus_depth_fig1 <- ggplot(Rhinogobius_brunneus_dep_res_grp1, aes(x = dep, y = HSI_dep_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Rhinogobius_brunneus depth HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Rhinogobius_brunneus_depth_fig1

#'- Rhinogobius_brunneus velocity
Rhinogobius_brunneus_glmm02_vel_fit <- effect_plot(Rhinogobius_brunneus_glmm02, pred = velocity, interval = TRUE, plot.points = FALSE)
Rhinogobius_brunneus_glmm02_vel_fit_PA <- as.data.frame(Rhinogobius_brunneus_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["Rhinogobius_brunneus_PA"]])
Rhinogobius_brunneus_glmm02_vel_fit_vel <- as.data.frame(Rhinogobius_brunneus_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Rhinogobius_brunneus_vel_res_grp1 <- cbind(Rhinogobius_brunneus_glmm02_vel_fit_PA, Rhinogobius_brunneus_glmm02_vel_fit_vel)
colnames(Rhinogobius_brunneus_vel_res_grp1) <- c("prob", "vel")

#'- normalize 0 to 1
Rhinogobius_brunneus_vel_res_grp1 <- Rhinogobius_brunneus_vel_res_grp1 %>%
  mutate(HSI_vel_grp1 = scale_values(prob))

#'- plot
Rhinogobius_brunneus_velocity_fig1 <- ggplot(Rhinogobius_brunneus_vel_res_grp1, aes(x = vel, y = HSI_vel_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Rhinogobius_brunneus Velocity HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Rhinogobius_brunneus_velocity_fig1

#'- Rhinogobius_brunneus substrate
Rhinogobius_brunneus_glmm02_sub_fit <- effect_plot(Rhinogobius_brunneus_glmm02, pred = substrate, interval = FALSE, plot.points = FALSE)
Rhinogobius_brunneus_glmm02_sub_fit_PA <- as.data.frame(Rhinogobius_brunneus_glmm02_sub_fit[["data"]][["Rhinogobius_brunneus_PA"]])
Rhinogobius_brunneus_glmm02_sub_fit_sub <- as.data.frame(Rhinogobius_brunneus_glmm02_sub_fit[["data"]][["substrate"]])
Rhinogobius_brunneus_sub_res_grp1 <- cbind(Rhinogobius_brunneus_glmm02_sub_fit_PA, Rhinogobius_brunneus_glmm02_sub_fit_sub)
colnames(Rhinogobius_brunneus_sub_res_grp1) <- c("prob", "sub")

#'- normalize 0 to 1
Rhinogobius_brunneus_sub_res_grp1 <- Rhinogobius_brunneus_sub_res_grp1 %>%
  mutate(HSI_sub_grp1 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Rhinogobius_brunneus_sub_res_grp1$sub <- factor(Rhinogobius_brunneus_sub_res_grp1$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                  'pebble', 'cobble', 'boulder'))

#'- plot
Rhinogobius_brunneus_substrate_fig1 <- ggplot(Rhinogobius_brunneus_sub_res_grp1, aes(x = sub, y = HSI_sub_grp1)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Rhinogobius_brunneus Substrate HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Rhinogobius_brunneus_substrate_fig1



##### Rhinogobius_brunneus group2 model #####

#'## select Group2
Rhinogobius_brunneus_dat03 <- Rhinogobius_brunneus_dat %>%
  filter(site_group == "GRP2")

#'##### glm
Rhinogobius_brunneus_glm03 <- glm(Rhinogobius_brunneus_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                  family=binomial(link = "logit"), data= Rhinogobius_brunneus_dat03, na.action ="na.fail")

summary(Rhinogobius_brunneus_glm03)

vif(Rhinogobius_brunneus_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect 
Rhinogobius_brunneus_glmm03 <- glmmTMB(Rhinogobius_brunneus_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                         (1|site) + (1|season), family=binomial(link = "logit"), data= Rhinogobius_brunneus_dat03, na.action ="na.fail")

summary(Rhinogobius_brunneus_glmm03)

#'- Model validation
r.squaredGLMM(Rhinogobius_brunneus_glmm03)
cor.test(predict(Rhinogobius_brunneus_glmm03, type = "response", allow.new.levels=TRUE), Rhinogobius_brunneus_dat03$Rhinogobius_brunneus_PA)

#'- Rhinogobius_brunneus depth
Rhinogobius_brunneus_glmm03_dep_fit <- effect_plot(Rhinogobius_brunneus_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Rhinogobius_brunneus_glmm03_dep_fit_PA <- as.data.frame(Rhinogobius_brunneus_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Rhinogobius_brunneus_PA"]])
Rhinogobius_brunneus_glmm03_dep_fit_depth <- as.data.frame(Rhinogobius_brunneus_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Rhinogobius_brunneus_dep_res_grp2 <- cbind(Rhinogobius_brunneus_glmm03_dep_fit_PA, Rhinogobius_brunneus_glmm03_dep_fit_depth)
colnames(Rhinogobius_brunneus_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Rhinogobius_brunneus_dep_res_grp2 <- Rhinogobius_brunneus_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Rhinogobius_brunneus_depth_fig2 <- ggplot(Rhinogobius_brunneus_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Rhinogobius_brunneus depth HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Rhinogobius_brunneus_depth_fig2

#'- Rhinogobius_brunneus velocity
Rhinogobius_brunneus_glmm03_vel_fit <- effect_plot(Rhinogobius_brunneus_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Rhinogobius_brunneus_glmm03_vel_fit_PA <- as.data.frame(Rhinogobius_brunneus_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Rhinogobius_brunneus_PA"]])
Rhinogobius_brunneus_glmm03_vel_fit_vel <- as.data.frame(Rhinogobius_brunneus_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Rhinogobius_brunneus_vel_res_grp2 <- cbind(Rhinogobius_brunneus_glmm03_vel_fit_PA, Rhinogobius_brunneus_glmm03_vel_fit_vel)
colnames(Rhinogobius_brunneus_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Rhinogobius_brunneus_vel_res_grp2 <- Rhinogobius_brunneus_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Rhinogobius_brunneus_velocity_fig2 <- ggplot(Rhinogobius_brunneus_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Rhinogobius_brunneus Velocity HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Rhinogobius_brunneus_velocity_fig2

#'- Rhinogobius_brunneus substrate
Rhinogobius_brunneus_glmm03_sub_fit <- effect_plot(Rhinogobius_brunneus_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Rhinogobius_brunneus_glmm03_sub_fit_PA <- as.data.frame(Rhinogobius_brunneus_glmm03_sub_fit[["data"]][["Rhinogobius_brunneus_PA"]])
Rhinogobius_brunneus_glmm03_sub_fit_sub <- as.data.frame(Rhinogobius_brunneus_glmm03_sub_fit[["data"]][["substrate"]])
Rhinogobius_brunneus_sub_res_grp2 <- cbind(Rhinogobius_brunneus_glmm03_sub_fit_PA, Rhinogobius_brunneus_glmm03_sub_fit_sub)
colnames(Rhinogobius_brunneus_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Rhinogobius_brunneus_sub_res_grp2 <- Rhinogobius_brunneus_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Rhinogobius_brunneus_sub_res_grp2$sub <- factor(Rhinogobius_brunneus_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                  'pebble', 'cobble', 'boulder'))

#'- plot
Rhinogobius_brunneus_substrate_fig2 <- ggplot(Rhinogobius_brunneus_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Rhinogobius_brunneus Substrate HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Rhinogobius_brunneus_substrate_fig2






#'## Combined figures

#'- merge mulitple effect plot
Rhinogobius_brunneus_dep_res_1 <- Rhinogobius_brunneus_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)

Rhinogobius_brunneus_dep_res_grp1_1 <- Rhinogobius_brunneus_dep_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp1)

Rhinogobius_brunneus_dep_res_grp2_1 <- Rhinogobius_brunneus_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)


Rhinogobius_brunneus_dep_res_comb <- rbind(Rhinogobius_brunneus_dep_res_1, Rhinogobius_brunneus_dep_res_grp1_1,
                                           Rhinogobius_brunneus_dep_res_grp2_1)


#'- depth plot
Rhinogobius_brunneus_depth_fig4 <- ggplot(Rhinogobius_brunneus_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 3, 0.5)) +     
  labs(title= "Rhinogobius brunneus Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Rhinogobius_brunneus_depth_fig4 <-Rhinogobius_brunneus_depth_fig4+scale_color_manual(values=c("black","blue3", "red", "green3"))
Rhinogobius_brunneus_depth_fig4




#'- merge mulitple effect plot
Rhinogobius_brunneus_vel_res_1 <- Rhinogobius_brunneus_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)

Rhinogobius_brunneus_vel_res_grp1_1 <- Rhinogobius_brunneus_vel_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp1)

Rhinogobius_brunneus_vel_res_grp2_1 <- Rhinogobius_brunneus_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)



Rhinogobius_brunneus_vel_res_comb <- rbind(Rhinogobius_brunneus_vel_res_1, Rhinogobius_brunneus_vel_res_grp1_1,
                                           Rhinogobius_brunneus_vel_res_grp2_1)


#'- velocity plot
Rhinogobius_brunneus_velocity_fig4 <- ggplot(Rhinogobius_brunneus_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2.5, 0.5)) +     
  labs(title= "Rhinogobius brunneus Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Rhinogobius_brunneus_velocity_fig4 <- Rhinogobius_brunneus_velocity_fig4 + scale_color_manual(values=c("black","blue3", "red", "green3", "orange3"))
Rhinogobius_brunneus_velocity_fig4



#'- merge mulitple effect plot
Rhinogobius_brunneus_sub_res_1 <- Rhinogobius_brunneus_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)

Rhinogobius_brunneus_sub_res_grp1_1 <- Rhinogobius_brunneus_sub_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp1)

Rhinogobius_brunneus_sub_res_grp2_1 <- Rhinogobius_brunneus_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)



Rhinogobius_brunneus_sub_res_comb <- rbind(Rhinogobius_brunneus_sub_res_1, Rhinogobius_brunneus_sub_res_grp1_1,
                                           Rhinogobius_brunneus_sub_res_grp2_1)

#'- substrate plot
Rhinogobius_brunneus_substrate_fig4 <- ggplot(Rhinogobius_brunneus_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Rhinogobius brunneus Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Rhinogobius_brunneus_substrate_fig4 <- Rhinogobius_brunneus_substrate_fig4 + scale_fill_manual(values=c("black","blue3", "red", "green3"))
Rhinogobius_brunneus_substrate_fig4


#'## combine figures

ggarrange(Rhinogobius_brunneus_depth_fig4, Rhinogobius_brunneus_velocity_fig4, Rhinogobius_brunneus_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 






##### Hemibarbus_labeo model #####

#'## (2) Hemibarbus_labeo model

#'### data management
Hemibarbus_labeo_dat <- HSI_DB2 %>%
  mutate(Hemibarbus_labeo_PA = ifelse(species == "Hemibarbus_labeo", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Hemibarbus_labeo_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Hemibarbus_labeo_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Hemibarbus_labeo_PA, depth, velocity, substrate, MBSNCD) 

head(Hemibarbus_labeo_dat)


#'## remove site where species did not occur
Hemibarbus_labeo_Pre <- HSI_DB2 %>%
  mutate(Hemibarbus_labeo_PA = ifelse(species == "Hemibarbus_labeo", 1, 0)) %>%
  aggregate(Hemibarbus_labeo_PA ~ site, "sum") %>%
  mutate(Hemibarbus_labeo_Pre = if_else(Hemibarbus_labeo_PA > 0, "Yes", "No")) %>%
  select(site, Hemibarbus_labeo_Pre)

Hemibarbus_labeo_dat <- Hemibarbus_labeo_dat %>%
  left_join(Hemibarbus_labeo_Pre, by = "site") %>%
  filter(Hemibarbus_labeo_Pre == "Yes") %>%
  as_tibble()





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

#'##### glm
Hemibarbus_labeo_glm01 <- glm(Hemibarbus_labeo_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                              family=binomial, data= Hemibarbus_labeo_dat, na.action ="na.fail")

summary(Hemibarbus_labeo_glm01)

vif(Hemibarbus_labeo_glm01)



#'##### glmm
#'- combination of site and season are used as a random effect 
Hemibarbus_labeo_glmm01 <- glmmTMB(Hemibarbus_labeo_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                     (1|site) + (1|season), family=binomial, data= Hemibarbus_labeo_dat, na.action ="na.fail")

summary(Hemibarbus_labeo_glmm01)




#'- Model validation
r.squaredGLMM(Hemibarbus_labeo_glmm01)
cor.test(predict(Hemibarbus_labeo_glmm01, type = "response", allow.new.levels=TRUE), Hemibarbus_labeo_dat$Hemibarbus_labeo_PA)

#'- Hemibarbus_labeo depth
Hemibarbus_labeo_glmm01_dep_fit <- effect_plot(Hemibarbus_labeo_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Hemibarbus_labeo_glmm01_dep_fit_PA <- as.data.frame(Hemibarbus_labeo_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_labeo_PA"]])
Hemibarbus_labeo_glmm01_dep_fit_dep <- as.data.frame(Hemibarbus_labeo_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Hemibarbus_labeo_dep_res <- cbind(Hemibarbus_labeo_glmm01_dep_fit_PA, Hemibarbus_labeo_glmm01_dep_fit_dep)
colnames(Hemibarbus_labeo_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Hemibarbus_labeo_dep_res <- Hemibarbus_labeo_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Hemibarbus_labeo_depth_fig <- ggplot(Hemibarbus_labeo_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_labeo depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_labeo_depth_fig

#'- Hemibarbus_labeo velocity
Hemibarbus_labeo_glmm01_vel_fit <- effect_plot(Hemibarbus_labeo_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Hemibarbus_labeo_glmm01_vel_fit_PA <- as.data.frame(Hemibarbus_labeo_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_labeo_PA"]])
Hemibarbus_labeo_glmm01_vel_fit_vel <- as.data.frame(Hemibarbus_labeo_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Hemibarbus_labeo_vel_res <- cbind(Hemibarbus_labeo_glmm01_vel_fit_PA, Hemibarbus_labeo_glmm01_vel_fit_vel)
colnames(Hemibarbus_labeo_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Hemibarbus_labeo_vel_res <- Hemibarbus_labeo_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Hemibarbus_labeo_velocity_fig <- ggplot(Hemibarbus_labeo_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_labeo Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_labeo_velocity_fig

#'- Hemibarbus_labeo substrate
Hemibarbus_labeo_glmm01_sub_fit <- effect_plot(Hemibarbus_labeo_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Hemibarbus_labeo_glmm01_sub_fit_PA <- as.data.frame(Hemibarbus_labeo_glmm01_sub_fit[["data"]][["Hemibarbus_labeo_PA"]])
Hemibarbus_labeo_glmm01_sub_fit_sub <- as.data.frame(Hemibarbus_labeo_glmm01_sub_fit[["data"]][["substrate"]])
Hemibarbus_labeo_sub_res <- cbind(Hemibarbus_labeo_glmm01_sub_fit_PA, Hemibarbus_labeo_glmm01_sub_fit_sub)
colnames(Hemibarbus_labeo_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Hemibarbus_labeo_sub_res <- Hemibarbus_labeo_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Hemibarbus_labeo_sub_res$sub <- factor(Hemibarbus_labeo_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                'pebble', 'cobble', 'boulder'))

#'- plot
Hemibarbus_labeo_substrate_fig <- ggplot(Hemibarbus_labeo_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Hemibarbus_labeo Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_labeo_substrate_fig


#'## (1) Hemibarbus_labeo model
##### Hemibarbus_labeo group1 model #####

#'## select Group1

Hemibarbus_labeo_dat02 <- Hemibarbus_labeo_dat %>%
  filter(site_group == "GRP1")

#'##### glm
Hemibarbus_labeo_glm02 <- glm(Hemibarbus_labeo_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                              family=binomial(link = "logit"), data= Hemibarbus_labeo_dat02, na.action ="na.fail")

summary(Hemibarbus_labeo_glm02)

vif(Hemibarbus_labeo_glm02)
#'##### glmm
#'- combination of site and season are used as a random effect ``
Hemibarbus_labeo_glmm02 <- glmmTMB(Hemibarbus_labeo_PA ~ poly(depth, 2) + poly(velocity, 1) + substrate + 
                                     (1|site) + (1|season), family=binomial(link = "logit"), data= Hemibarbus_labeo_dat02, na.action ="na.fail")

summary(Hemibarbus_labeo_glmm02)



#'- Model validation
r.squaredGLMM(Hemibarbus_labeo_glmm02)
cor.test(predict(Hemibarbus_labeo_glmm02, type = "response", allow.new.levels=TRUE), Hemibarbus_labeo_dat02$Hemibarbus_labeo_PA)

#'- Hemibarbus_labeo depth
Hemibarbus_labeo_glmm02_dep_fit <- effect_plot(Hemibarbus_labeo_glmm02, pred = depth, interval = TRUE, plot.points = FALSE)
Hemibarbus_labeo_glmm02_dep_fit_PA <- as.data.frame(Hemibarbus_labeo_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_labeo_PA"]])
Hemibarbus_labeo_glmm02_dep_fit_depth <- as.data.frame(Hemibarbus_labeo_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Hemibarbus_labeo_dep_res_grp1 <- cbind(Hemibarbus_labeo_glmm02_dep_fit_PA, Hemibarbus_labeo_glmm02_dep_fit_depth)
colnames(Hemibarbus_labeo_dep_res_grp1) <- c("prob", "dep")

#'- normalize 0 to 1
Hemibarbus_labeo_dep_res_grp1 <- Hemibarbus_labeo_dep_res_grp1 %>%
  mutate(HSI_dep_grp1 = scale_values(prob))

#'- plot
Hemibarbus_labeo_depth_fig1 <- ggplot(Hemibarbus_labeo_dep_res_grp1, aes(x = dep, y = HSI_dep_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_labeo depth HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_labeo_depth_fig1

#'- Hemibarbus_labeo velocity
Hemibarbus_labeo_glmm02_vel_fit <- effect_plot(Hemibarbus_labeo_glmm02, pred = velocity, interval = TRUE, plot.points = FALSE)
Hemibarbus_labeo_glmm02_vel_fit_PA <- as.data.frame(Hemibarbus_labeo_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_labeo_PA"]])
Hemibarbus_labeo_glmm02_vel_fit_vel <- as.data.frame(Hemibarbus_labeo_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Hemibarbus_labeo_vel_res_grp1 <- cbind(Hemibarbus_labeo_glmm02_vel_fit_PA, Hemibarbus_labeo_glmm02_vel_fit_vel)
colnames(Hemibarbus_labeo_vel_res_grp1) <- c("prob", "vel")

#'- normalize 0 to 1
Hemibarbus_labeo_vel_res_grp1 <- Hemibarbus_labeo_vel_res_grp1 %>%
  mutate(HSI_vel_grp1 = scale_values(prob))

#'- plot
Hemibarbus_labeo_velocity_fig1 <- ggplot(Hemibarbus_labeo_vel_res_grp1, aes(x = vel, y = HSI_vel_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_labeo Velocity HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_labeo_velocity_fig1

#'- Hemibarbus_labeo substrate
Hemibarbus_labeo_glmm02_sub_fit <- effect_plot(Hemibarbus_labeo_glmm02, pred = substrate, interval = FALSE, plot.points = FALSE)
Hemibarbus_labeo_glmm02_sub_fit_PA <- as.data.frame(Hemibarbus_labeo_glmm02_sub_fit[["data"]][["Hemibarbus_labeo_PA"]])
Hemibarbus_labeo_glmm02_sub_fit_sub <- as.data.frame(Hemibarbus_labeo_glmm02_sub_fit[["data"]][["substrate"]])
Hemibarbus_labeo_sub_res_grp1 <- cbind(Hemibarbus_labeo_glmm02_sub_fit_PA, Hemibarbus_labeo_glmm02_sub_fit_sub)
colnames(Hemibarbus_labeo_sub_res_grp1) <- c("prob", "sub")

#'- normalize 0 to 1
Hemibarbus_labeo_sub_res_grp1 <- Hemibarbus_labeo_sub_res_grp1 %>%
  mutate(HSI_sub_grp1 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Hemibarbus_labeo_sub_res_grp1$sub <- factor(Hemibarbus_labeo_sub_res_grp1$sub, levels = c('silt', 'sand', 'gravel',
                                                                                          'pebble', 'cobble', 'boulder'))

#'- plot
Hemibarbus_labeo_substrate_fig1 <- ggplot(Hemibarbus_labeo_sub_res_grp1, aes(x = sub, y = HSI_sub_grp1)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Hemibarbus_labeo Substrate HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_labeo_substrate_fig1



##### Hemibarbus_labeo group2 model #####

#'## select Group2
Hemibarbus_labeo_dat03 <- Hemibarbus_labeo_dat %>%
  filter(site_group == "GRP2")

#'##### glm
Hemibarbus_labeo_glm03 <- glm(Hemibarbus_labeo_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                              family=binomial(link = "logit"), data= Hemibarbus_labeo_dat03, na.action ="na.fail")

summary(Hemibarbus_labeo_glm03)

vif(Hemibarbus_labeo_glm03)


#'##### glmm
#'- combination of site and season are used as a random effect 
Hemibarbus_labeo_glmm03 <- glmmTMB(Hemibarbus_labeo_PA ~ poly(depth, 1) + poly(velocity, 1) + substrate + 
                                     (1|site) + (1|season), family=binomial(link = "logit"), data= Hemibarbus_labeo_dat03, na.action ="na.fail")

summary(Hemibarbus_labeo_glmm03)

#'- Model validation
r.squaredGLMM(Hemibarbus_labeo_glmm03)
cor.test(predict(Hemibarbus_labeo_glmm03, type = "response", allow.new.levels=TRUE), Hemibarbus_labeo_dat03$Hemibarbus_labeo_PA)

#'- Hemibarbus_labeo depth
Hemibarbus_labeo_glmm03_dep_fit <- effect_plot(Hemibarbus_labeo_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Hemibarbus_labeo_glmm03_dep_fit_PA <- as.data.frame(Hemibarbus_labeo_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_labeo_PA"]])
Hemibarbus_labeo_glmm03_dep_fit_depth <- as.data.frame(Hemibarbus_labeo_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Hemibarbus_labeo_dep_res_grp2 <- cbind(Hemibarbus_labeo_glmm03_dep_fit_PA, Hemibarbus_labeo_glmm03_dep_fit_depth)
colnames(Hemibarbus_labeo_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Hemibarbus_labeo_dep_res_grp2 <- Hemibarbus_labeo_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Hemibarbus_labeo_depth_fig2 <- ggplot(Hemibarbus_labeo_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_labeo depth HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_labeo_depth_fig2

#'- Hemibarbus_labeo velocity
Hemibarbus_labeo_glmm03_vel_fit <- effect_plot(Hemibarbus_labeo_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Hemibarbus_labeo_glmm03_vel_fit_PA <- as.data.frame(Hemibarbus_labeo_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_labeo_PA"]])
Hemibarbus_labeo_glmm03_vel_fit_vel <- as.data.frame(Hemibarbus_labeo_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Hemibarbus_labeo_vel_res_grp2 <- cbind(Hemibarbus_labeo_glmm03_vel_fit_PA, Hemibarbus_labeo_glmm03_vel_fit_vel)
colnames(Hemibarbus_labeo_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Hemibarbus_labeo_vel_res_grp2 <- Hemibarbus_labeo_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Hemibarbus_labeo_velocity_fig2 <- ggplot(Hemibarbus_labeo_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_labeo Velocity HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_labeo_velocity_fig2

#'- Hemibarbus_labeo substrate
Hemibarbus_labeo_glmm03_sub_fit <- effect_plot(Hemibarbus_labeo_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Hemibarbus_labeo_glmm03_sub_fit_PA <- as.data.frame(Hemibarbus_labeo_glmm03_sub_fit[["data"]][["Hemibarbus_labeo_PA"]])
Hemibarbus_labeo_glmm03_sub_fit_sub <- as.data.frame(Hemibarbus_labeo_glmm03_sub_fit[["data"]][["substrate"]])
Hemibarbus_labeo_sub_res_grp2 <- cbind(Hemibarbus_labeo_glmm03_sub_fit_PA, Hemibarbus_labeo_glmm03_sub_fit_sub)
colnames(Hemibarbus_labeo_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Hemibarbus_labeo_sub_res_grp2 <- Hemibarbus_labeo_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Hemibarbus_labeo_sub_res_grp2$sub <- factor(Hemibarbus_labeo_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                          'pebble', 'cobble', 'boulder'))

#'- plot
Hemibarbus_labeo_substrate_fig2 <- ggplot(Hemibarbus_labeo_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Hemibarbus_labeo Substrate HSI Group2") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_labeo_substrate_fig2





#'## Combined figures

#'- merge mulitple effect plot
Hemibarbus_labeo_dep_res_1 <- Hemibarbus_labeo_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)

Hemibarbus_labeo_dep_res_grp1_1 <- Hemibarbus_labeo_dep_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp1)

Hemibarbus_labeo_dep_res_grp2_1 <- Hemibarbus_labeo_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)



Hemibarbus_labeo_dep_res_comb <- rbind(Hemibarbus_labeo_dep_res_1, Hemibarbus_labeo_dep_res_grp1_1,
                                       Hemibarbus_labeo_dep_res_grp2_1)


#'- depth plot
Hemibarbus_labeo_depth_fig4 <- ggplot(Hemibarbus_labeo_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_labeo Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Hemibarbus_labeo_depth_fig4 <-Hemibarbus_labeo_depth_fig4+scale_color_manual(values=c("black","blue3", "red", "green3"))
Hemibarbus_labeo_depth_fig4




#'- merge mulitple effect plot
Hemibarbus_labeo_vel_res_1 <- Hemibarbus_labeo_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)

Hemibarbus_labeo_vel_res_grp1_1 <- Hemibarbus_labeo_vel_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp1)

Hemibarbus_labeo_vel_res_grp2_1 <- Hemibarbus_labeo_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)



Hemibarbus_labeo_vel_res_comb <- rbind(Hemibarbus_labeo_vel_res_1, Hemibarbus_labeo_vel_res_grp1_1,
                                       Hemibarbus_labeo_vel_res_grp2_1)


#'- velocity plot
Hemibarbus_labeo_velocity_fig4 <- ggplot(Hemibarbus_labeo_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_labeo Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Hemibarbus_labeo_velocity_fig4 <- Hemibarbus_labeo_velocity_fig4 + scale_color_manual(values=c("black","blue3", "red", "green3", "orange3"))
Hemibarbus_labeo_velocity_fig4



#'- merge mulitple effect plot
Hemibarbus_labeo_sub_res_1 <- Hemibarbus_labeo_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)

Hemibarbus_labeo_sub_res_grp1_1 <- Hemibarbus_labeo_sub_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp1)

Hemibarbus_labeo_sub_res_grp2_1 <- Hemibarbus_labeo_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)




Hemibarbus_labeo_sub_res_comb <- rbind(Hemibarbus_labeo_sub_res_1, Hemibarbus_labeo_sub_res_grp1_1,
                                       Hemibarbus_labeo_sub_res_grp2_1)

#'- substrate plot
Hemibarbus_labeo_substrate_fig4 <- ggplot(Hemibarbus_labeo_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Hemibarbus labeo Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Hemibarbus_labeo_substrate_fig4 <- Hemibarbus_labeo_substrate_fig4 + scale_fill_manual(values=c("black","blue3", "red", "green3"))
Hemibarbus_labeo_substrate_fig4


#'## combine figures

ggarrange(Hemibarbus_labeo_depth_fig4, Hemibarbus_labeo_velocity_fig4, Hemibarbus_labeo_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 








##### Squalidus_chankaensis_tsuchigae model #####

#'## (2) Squalidus_chankaensis_tsuchigae model

#'### data management
Squalidus_chankaensis_tsuchigae_dat <- HSI_DB2 %>%
  mutate(Squalidus_chankaensis_tsuchigae_PA = ifelse(species == "Squalidus_chankaensis_tsuchigae", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Squalidus_chankaensis_tsuchigae_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Squalidus_chankaensis_tsuchigae_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Squalidus_chankaensis_tsuchigae_PA, depth, velocity, substrate, MBSNCD) 

head(Squalidus_chankaensis_tsuchigae_dat)


#'## remove site where species did not occur
Squalidus_chankaensis_tsuchigae_Pre <- HSI_DB2 %>%
  mutate(Squalidus_chankaensis_tsuchigae_PA = ifelse(species == "Squalidus_chankaensis_tsuchigae", 1, 0)) %>%
  aggregate(Squalidus_chankaensis_tsuchigae_PA ~ site, "sum") %>%
  mutate(Squalidus_chankaensis_tsuchigae_Pre = if_else(Squalidus_chankaensis_tsuchigae_PA > 0, "Yes", "No")) %>%
  select(site, Squalidus_chankaensis_tsuchigae_Pre)

Squalidus_chankaensis_tsuchigae_dat <- Squalidus_chankaensis_tsuchigae_dat %>%
  left_join(Squalidus_chankaensis_tsuchigae_Pre, by = "site") %>%
  filter(Squalidus_chankaensis_tsuchigae_Pre == "Yes") %>%
  as_tibble()



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

#'##### glm
Squalidus_chankaensis_tsuchigae_glm01 <- glm(Squalidus_chankaensis_tsuchigae_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                             family=binomial(link = "logit"), data= Squalidus_chankaensis_tsuchigae_dat, na.action ="na.fail")

summary(Squalidus_chankaensis_tsuchigae_glm01)

vif(Squalidus_chankaensis_tsuchigae_glm01)



#'##### glmm
#'- combination of site and season are used as a random effect 
Squalidus_chankaensis_tsuchigae_glmm01 <- glmmTMB(Squalidus_chankaensis_tsuchigae_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                    (1|MBSNCD/site) + (1|season), family=binomial(link = "logit"), data= Squalidus_chankaensis_tsuchigae_dat, na.action ="na.fail")

summary(Squalidus_chankaensis_tsuchigae_glmm01)










#'- Model validation
r.squaredGLMM(Squalidus_chankaensis_tsuchigae_glmm01)
cor.test(predict(Squalidus_chankaensis_tsuchigae_glmm01, type = "response", allow.new.levels=TRUE), Squalidus_chankaensis_tsuchigae_dat$Squalidus_chankaensis_tsuchigae_PA)

#'- Squalidus_chankaensis_tsuchigae depth
Squalidus_chankaensis_tsuchigae_glmm01_dep_fit <- effect_plot(Squalidus_chankaensis_tsuchigae_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Squalidus_chankaensis_tsuchigae_glmm01_dep_fit_PA <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Squalidus_chankaensis_tsuchigae_PA"]])
Squalidus_chankaensis_tsuchigae_glmm01_dep_fit_dep <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Squalidus_chankaensis_tsuchigae_dep_res <- cbind(Squalidus_chankaensis_tsuchigae_glmm01_dep_fit_PA, Squalidus_chankaensis_tsuchigae_glmm01_dep_fit_dep)
colnames(Squalidus_chankaensis_tsuchigae_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Squalidus_chankaensis_tsuchigae_dep_res <- Squalidus_chankaensis_tsuchigae_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Squalidus_chankaensis_tsuchigae_depth_fig <- ggplot(Squalidus_chankaensis_tsuchigae_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Squalidus_chankaensis_tsuchigae depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Squalidus_chankaensis_tsuchigae_depth_fig

#'- Squalidus_chankaensis_tsuchigae velocity
Squalidus_chankaensis_tsuchigae_glmm01_vel_fit <- effect_plot(Squalidus_chankaensis_tsuchigae_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Squalidus_chankaensis_tsuchigae_glmm01_vel_fit_PA <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Squalidus_chankaensis_tsuchigae_PA"]])
Squalidus_chankaensis_tsuchigae_glmm01_vel_fit_vel <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Squalidus_chankaensis_tsuchigae_vel_res <- cbind(Squalidus_chankaensis_tsuchigae_glmm01_vel_fit_PA, Squalidus_chankaensis_tsuchigae_glmm01_vel_fit_vel)
colnames(Squalidus_chankaensis_tsuchigae_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Squalidus_chankaensis_tsuchigae_vel_res <- Squalidus_chankaensis_tsuchigae_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Squalidus_chankaensis_tsuchigae_velocity_fig <- ggplot(Squalidus_chankaensis_tsuchigae_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Squalidus_chankaensis_tsuchigae Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Squalidus_chankaensis_tsuchigae_velocity_fig

#'- Squalidus_chankaensis_tsuchigae substrate
Squalidus_chankaensis_tsuchigae_glmm01_sub_fit <- effect_plot(Squalidus_chankaensis_tsuchigae_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Squalidus_chankaensis_tsuchigae_glmm01_sub_fit_PA <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm01_sub_fit[["data"]][["Squalidus_chankaensis_tsuchigae_PA"]])
Squalidus_chankaensis_tsuchigae_glmm01_sub_fit_sub <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm01_sub_fit[["data"]][["substrate"]])
Squalidus_chankaensis_tsuchigae_sub_res <- cbind(Squalidus_chankaensis_tsuchigae_glmm01_sub_fit_PA, Squalidus_chankaensis_tsuchigae_glmm01_sub_fit_sub)
colnames(Squalidus_chankaensis_tsuchigae_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Squalidus_chankaensis_tsuchigae_sub_res <- Squalidus_chankaensis_tsuchigae_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Squalidus_chankaensis_tsuchigae_sub_res$sub <- factor(Squalidus_chankaensis_tsuchigae_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                              'pebble', 'cobble', 'boulder'))

#'- plot
Squalidus_chankaensis_tsuchigae_substrate_fig <- ggplot(Squalidus_chankaensis_tsuchigae_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Squalidus_chankaensis_tsuchigae Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Squalidus_chankaensis_tsuchigae_substrate_fig



#'## (1) Squalidus_chankaensis_tsuchigae model
##### Squalidus_chankaensis_tsuchigae group1 model #####

#'## select Group1

Squalidus_chankaensis_tsuchigae_dat02 <- Squalidus_chankaensis_tsuchigae_dat %>%
  filter(site_group == "GRP1")

#'##### glm
Squalidus_chankaensis_tsuchigae_glm02 <- glm(Squalidus_chankaensis_tsuchigae_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                             family=binomial(link = "logit"), data= Squalidus_chankaensis_tsuchigae_dat02, na.action ="na.fail")

summary(Squalidus_chankaensis_tsuchigae_glm02)

vif(Squalidus_chankaensis_tsuchigae_glm02)
#'##### glmm
#'- combination of site and season are used as a random effect ``
Squalidus_chankaensis_tsuchigae_glmm02 <- glmmTMB(Squalidus_chankaensis_tsuchigae_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                    (1|MBSNCD/site) + (1|season), family=binomial(link = "logit"), data= Squalidus_chankaensis_tsuchigae_dat02, na.action ="na.fail")

summary(Squalidus_chankaensis_tsuchigae_glmm02)



#'- Model validation
r.squaredGLMM(Squalidus_chankaensis_tsuchigae_glmm02)
cor.test(predict(Squalidus_chankaensis_tsuchigae_glmm02, type = "response", allow.new.levels=TRUE), Squalidus_chankaensis_tsuchigae_dat02$Squalidus_chankaensis_tsuchigae_PA)

#'- Squalidus_chankaensis_tsuchigae depth
Squalidus_chankaensis_tsuchigae_glmm02_dep_fit <- effect_plot(Squalidus_chankaensis_tsuchigae_glmm02, pred = depth, interval = TRUE, plot.points = FALSE)
Squalidus_chankaensis_tsuchigae_glmm02_dep_fit_PA <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["Squalidus_chankaensis_tsuchigae_PA"]])
Squalidus_chankaensis_tsuchigae_glmm02_dep_fit_depth <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Squalidus_chankaensis_tsuchigae_dep_res_grp1 <- cbind(Squalidus_chankaensis_tsuchigae_glmm02_dep_fit_PA, Squalidus_chankaensis_tsuchigae_glmm02_dep_fit_depth)
colnames(Squalidus_chankaensis_tsuchigae_dep_res_grp1) <- c("prob", "dep")

#'- normalize 0 to 1
Squalidus_chankaensis_tsuchigae_dep_res_grp1 <- Squalidus_chankaensis_tsuchigae_dep_res_grp1 %>%
  mutate(HSI_dep_grp1 = scale_values(prob))

#'- plot
Squalidus_chankaensis_tsuchigae_depth_fig1 <- ggplot(Squalidus_chankaensis_tsuchigae_dep_res_grp1, aes(x = dep, y = HSI_dep_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Squalidus_chankaensis_tsuchigae depth HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Squalidus_chankaensis_tsuchigae_depth_fig1

#'- Squalidus_chankaensis_tsuchigae velocity
Squalidus_chankaensis_tsuchigae_glmm02_vel_fit <- effect_plot(Squalidus_chankaensis_tsuchigae_glmm02, pred = velocity, interval = TRUE, plot.points = FALSE)
Squalidus_chankaensis_tsuchigae_glmm02_vel_fit_PA <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["Squalidus_chankaensis_tsuchigae_PA"]])
Squalidus_chankaensis_tsuchigae_glmm02_vel_fit_vel <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Squalidus_chankaensis_tsuchigae_vel_res_grp1 <- cbind(Squalidus_chankaensis_tsuchigae_glmm02_vel_fit_PA, Squalidus_chankaensis_tsuchigae_glmm02_vel_fit_vel)
colnames(Squalidus_chankaensis_tsuchigae_vel_res_grp1) <- c("prob", "vel")

#'- normalize 0 to 1
Squalidus_chankaensis_tsuchigae_vel_res_grp1 <- Squalidus_chankaensis_tsuchigae_vel_res_grp1 %>%
  mutate(HSI_vel_grp1 = scale_values(prob))

#'- plot
Squalidus_chankaensis_tsuchigae_velocity_fig1 <- ggplot(Squalidus_chankaensis_tsuchigae_vel_res_grp1, aes(x = vel, y = HSI_vel_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Squalidus_chankaensis_tsuchigae Velocity HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Squalidus_chankaensis_tsuchigae_velocity_fig1

#'- Squalidus_chankaensis_tsuchigae substrate
Squalidus_chankaensis_tsuchigae_glmm02_sub_fit <- effect_plot(Squalidus_chankaensis_tsuchigae_glmm02, pred = substrate, interval = FALSE, plot.points = FALSE)
Squalidus_chankaensis_tsuchigae_glmm02_sub_fit_PA <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm02_sub_fit[["data"]][["Squalidus_chankaensis_tsuchigae_PA"]])
Squalidus_chankaensis_tsuchigae_glmm02_sub_fit_sub <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm02_sub_fit[["data"]][["substrate"]])
Squalidus_chankaensis_tsuchigae_sub_res_grp1 <- cbind(Squalidus_chankaensis_tsuchigae_glmm02_sub_fit_PA, Squalidus_chankaensis_tsuchigae_glmm02_sub_fit_sub)
colnames(Squalidus_chankaensis_tsuchigae_sub_res_grp1) <- c("prob", "sub")

#'- normalize 0 to 1
Squalidus_chankaensis_tsuchigae_sub_res_grp1 <- Squalidus_chankaensis_tsuchigae_sub_res_grp1 %>%
  mutate(HSI_sub_grp1 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Squalidus_chankaensis_tsuchigae_sub_res_grp1$sub <- factor(Squalidus_chankaensis_tsuchigae_sub_res_grp1$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                                        'pebble', 'cobble', 'boulder'))

#'- plot
Squalidus_chankaensis_tsuchigae_substrate_fig1 <- ggplot(Squalidus_chankaensis_tsuchigae_sub_res_grp1, aes(x = sub, y = HSI_sub_grp1)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Squalidus_chankaensis_tsuchigae Substrate HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Squalidus_chankaensis_tsuchigae_substrate_fig1



##### Squalidus_chankaensis_tsuchigae group2 model #####

#'## select Group2
Squalidus_chankaensis_tsuchigae_dat03 <- Squalidus_chankaensis_tsuchigae_dat %>%
  filter(site_group == "GRP2")
sum( Squalidus_chankaensis_tsuchigae_dat03$Squalidus_chankaensis_tsuchigae_PA=="1") # Presence data(n = 21)
#'##### glm
Squalidus_chankaensis_tsuchigae_glm03 <- glm(Squalidus_chankaensis_tsuchigae_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate,
                                             family=binomial(link = "logit"), data= Squalidus_chankaensis_tsuchigae_dat03, na.action ="na.fail")

summary(Squalidus_chankaensis_tsuchigae_glm03)

vif(Squalidus_chankaensis_tsuchigae_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect
Squalidus_chankaensis_tsuchigae_glmm03 <- glmmTMB(Squalidus_chankaensis_tsuchigae_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate +
                                                    (1|MBSNCD/site) + (1|season),family=binomial(link = "logit"), data= Squalidus_chankaensis_tsuchigae_dat03, na.action ="na.fail")

summary(Squalidus_chankaensis_tsuchigae_glmm03)

#'- Model validation
r.squaredGLMM(Squalidus_chankaensis_tsuchigae_glmm03)
cor.test(predict(Squalidus_chankaensis_tsuchigae_glmm03, type = "response", allow.new.levels=TRUE), Squalidus_chankaensis_tsuchigae_dat03$Squalidus_chankaensis_tsuchigae_PA)

#'- Squalidus_chankaensis_tsuchigae depth
Squalidus_chankaensis_tsuchigae_glmm03_dep_fit <- effect_plot(Squalidus_chankaensis_tsuchigae_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Squalidus_chankaensis_tsuchigae_glmm03_dep_fit_PA <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Squalidus_chankaensis_tsuchigae_PA"]])
Squalidus_chankaensis_tsuchigae_glmm03_dep_fit_depth <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Squalidus_chankaensis_tsuchigae_dep_res_grp2 <- cbind(Squalidus_chankaensis_tsuchigae_glmm03_dep_fit_PA, Squalidus_chankaensis_tsuchigae_glmm03_dep_fit_depth)
colnames(Squalidus_chankaensis_tsuchigae_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Squalidus_chankaensis_tsuchigae_dep_res_grp2 <- Squalidus_chankaensis_tsuchigae_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Squalidus_chankaensis_tsuchigae_depth_fig2 <- ggplot(Squalidus_chankaensis_tsuchigae_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) +
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +
  labs(title= "Squalidus_chankaensis_tsuchigae depth HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Squalidus_chankaensis_tsuchigae_depth_fig2

#'- Squalidus_chankaensis_tsuchigae velocity
Squalidus_chankaensis_tsuchigae_glmm03_vel_fit <- effect_plot(Squalidus_chankaensis_tsuchigae_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Squalidus_chankaensis_tsuchigae_glmm03_vel_fit_PA <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Squalidus_chankaensis_tsuchigae_PA"]])
Squalidus_chankaensis_tsuchigae_glmm03_vel_fit_vel <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Squalidus_chankaensis_tsuchigae_vel_res_grp2 <- cbind(Squalidus_chankaensis_tsuchigae_glmm03_vel_fit_PA, Squalidus_chankaensis_tsuchigae_glmm03_vel_fit_vel)
colnames(Squalidus_chankaensis_tsuchigae_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Squalidus_chankaensis_tsuchigae_vel_res_grp2 <- Squalidus_chankaensis_tsuchigae_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Squalidus_chankaensis_tsuchigae_velocity_fig2 <- ggplot(Squalidus_chankaensis_tsuchigae_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) +
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +
  labs(title= "Squalidus_chankaensis_tsuchigae Velocity HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Squalidus_chankaensis_tsuchigae_velocity_fig2

#'- Squalidus_chankaensis_tsuchigae substrate
Squalidus_chankaensis_tsuchigae_glmm03_sub_fit <- effect_plot(Squalidus_chankaensis_tsuchigae_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Squalidus_chankaensis_tsuchigae_glmm03_sub_fit_PA <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm03_sub_fit[["data"]][["Squalidus_chankaensis_tsuchigae_PA"]])
Squalidus_chankaensis_tsuchigae_glmm03_sub_fit_sub <- as.data.frame(Squalidus_chankaensis_tsuchigae_glmm03_sub_fit[["data"]][["substrate"]])
Squalidus_chankaensis_tsuchigae_sub_res_grp2 <- cbind(Squalidus_chankaensis_tsuchigae_glmm03_sub_fit_PA, Squalidus_chankaensis_tsuchigae_glmm03_sub_fit_sub)
colnames(Squalidus_chankaensis_tsuchigae_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Squalidus_chankaensis_tsuchigae_sub_res_grp2 <- Squalidus_chankaensis_tsuchigae_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Squalidus_chankaensis_tsuchigae_sub_res_grp2$sub <- factor(Squalidus_chankaensis_tsuchigae_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                                        'pebble', 'cobble', 'boulder'))

#'- plot
Squalidus_chankaensis_tsuchigae_substrate_fig2 <- ggplot(Squalidus_chankaensis_tsuchigae_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +
  labs(title= "Squalidus_chankaensis_tsuchigae Substrate HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Squalidus_chankaensis_tsuchigae_substrate_fig2




#'## Combined figures

#'- merge mulitple effect plot
Squalidus_chankaensis_tsuchigae_dep_res_1 <- Squalidus_chankaensis_tsuchigae_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)

Squalidus_chankaensis_tsuchigae_dep_res_grp1_1 <- Squalidus_chankaensis_tsuchigae_dep_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp1)

Squalidus_chankaensis_tsuchigae_dep_res_grp2_1 <- Squalidus_chankaensis_tsuchigae_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)

Squalidus_chankaensis_tsuchigae_dep_res_comb <- rbind(Squalidus_chankaensis_tsuchigae_dep_res_1, Squalidus_chankaensis_tsuchigae_dep_res_grp1_1,
                                                      Squalidus_chankaensis_tsuchigae_dep_res_grp2_1)


#'- depth plot
Squalidus_chankaensis_tsuchigae_depth_fig4 <- ggplot(Squalidus_chankaensis_tsuchigae_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Squalidus_chankaensis_tsuchigae Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Squalidus_chankaensis_tsuchigae_depth_fig4 <-Squalidus_chankaensis_tsuchigae_depth_fig4+scale_color_manual(values=c("black","blue3", "red", "green3"))
Squalidus_chankaensis_tsuchigae_depth_fig4




#'- merge mulitple effect plot
Squalidus_chankaensis_tsuchigae_vel_res_1 <- Squalidus_chankaensis_tsuchigae_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)

Squalidus_chankaensis_tsuchigae_vel_res_grp1_1 <- Squalidus_chankaensis_tsuchigae_vel_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp1)

Squalidus_chankaensis_tsuchigae_vel_res_grp2_1 <- Squalidus_chankaensis_tsuchigae_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)

Squalidus_chankaensis_tsuchigae_vel_res_comb <- rbind(Squalidus_chankaensis_tsuchigae_vel_res_1, Squalidus_chankaensis_tsuchigae_vel_res_grp1_1,
                                                      Squalidus_chankaensis_tsuchigae_vel_res_grp2_1)


#'- velocity plot
Squalidus_chankaensis_tsuchigae_velocity_fig4 <- ggplot(Squalidus_chankaensis_tsuchigae_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Squalidus_chankaensis_tsuchigae Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Squalidus_chankaensis_tsuchigae_velocity_fig4 <- Squalidus_chankaensis_tsuchigae_velocity_fig4 + scale_color_manual(values=c("black","blue3", "red", "green3"))
Squalidus_chankaensis_tsuchigae_velocity_fig4



#'- merge mulitple effect plot
Squalidus_chankaensis_tsuchigae_sub_res_1 <- Squalidus_chankaensis_tsuchigae_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)

Squalidus_chankaensis_tsuchigae_sub_res_grp1_1 <- Squalidus_chankaensis_tsuchigae_sub_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp1)

Squalidus_chankaensis_tsuchigae_sub_res_grp2_1 <- Squalidus_chankaensis_tsuchigae_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)


Squalidus_chankaensis_tsuchigae_sub_res_comb <- rbind(Squalidus_chankaensis_tsuchigae_sub_res_1, Squalidus_chankaensis_tsuchigae_sub_res_grp1_1,
                                                      Squalidus_chankaensis_tsuchigae_sub_res_grp2_1)

#'- substrate plot
Squalidus_chankaensis_tsuchigae_substrate_fig4 <- ggplot(Squalidus_chankaensis_tsuchigae_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Squalidus_chankaensis_tsuchigae Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Squalidus_chankaensis_tsuchigae_substrate_fig4 <- Squalidus_chankaensis_tsuchigae_substrate_fig4 + scale_fill_manual(values=c("black","blue3", "red", "green3"))
Squalidus_chankaensis_tsuchigae_substrate_fig4


#'## combine figures

ggarrange(Squalidus_chankaensis_tsuchigae_depth_fig4, Squalidus_chankaensis_tsuchigae_velocity_fig4, Squalidus_chankaensis_tsuchigae_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 



##### Microphysogobio_longidorsalis model #####

#'## (2) Microphysogobio_longidorsalis model

#'### data management
Microphysogobio_longidorsalis_dat <- HSI_DB2 %>%
  mutate(Microphysogobio_longidorsalis_PA = ifelse(species == "Microphysogobio_longidorsalis", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Microphysogobio_longidorsalis_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Microphysogobio_longidorsalis_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Microphysogobio_longidorsalis_PA, depth, velocity, substrate, MBSNCD) 

head(Microphysogobio_longidorsalis_dat)


#'## remove site where species did not occur
Microphysogobio_longidorsalis_Pre <- HSI_DB2 %>%
  mutate(Microphysogobio_longidorsalis_PA = ifelse(species == "Microphysogobio_longidorsalis", 1, 0)) %>%
  aggregate(Microphysogobio_longidorsalis_PA ~ site, "sum") %>%
  mutate(Microphysogobio_longidorsalis_Pre = if_else(Microphysogobio_longidorsalis_PA > 0, "Yes", "No")) %>%
  select(site, Microphysogobio_longidorsalis_Pre)

Microphysogobio_longidorsalis_dat <- Microphysogobio_longidorsalis_dat %>%
  left_join(Microphysogobio_longidorsalis_Pre, by = "site") %>%
  filter(Microphysogobio_longidorsalis_Pre == "Yes") %>%
  as_tibble()

# Microphysogobio_longidorsalis_dat$Microphysogobio_longidorsalis_PA[Microphysogobio_longidorsalis_dat$substrate == "silt"] <- 0


#'## (1) Microphysogobio_longidorsalis model

#'### data management
Microphysogobio_longidorsalis_dat <- HSI_DB2 %>%
  mutate(Microphysogobio_longidorsalis_PA = ifelse(species == "Microphysogobio_longidorsalis", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Microphysogobio_longidorsalis_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Microphysogobio_longidorsalis_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section - to remove duplicae in presence data
  select(basin, site_group, site_by_season, site, date, season, sample_number, Microphysogobio_longidorsalis_PA, depth, velocity, substrate, MBSNCD)%>%
  mutate(Microphysogobio_longidorsalis_PA = ifelse(substrate == "silt", 0, Microphysogobio_longidorsalis_PA))













Microphysogobio_longidorsalis_Pre <- HSI_DB2 %>%
  mutate(Microphysogobio_longidorsalis_PA = ifelse(species == "Microphysogobio_longidorsalis", 1, 0)) %>%
  aggregate(Microphysogobio_longidorsalis_PA ~ site, "sum") %>%
  mutate(Microphysogobio_longidorsalis_Pre = if_else(Microphysogobio_longidorsalis_PA > 0, "Yes", "No")) %>%
  select(site, Microphysogobio_longidorsalis_Pre)

Microphysogobio_longidorsalis_dat <- Microphysogobio_longidorsalis_dat %>%
  left_join(Microphysogobio_longidorsalis_Pre, by = "site") %>%
  filter(Microphysogobio_longidorsalis_Pre == "Yes")

#check_species_site_type
list(unique(Microphysogobio_longidorsalis_dat$site_group))
sum(Microphysogobio_longidorsalis_dat$site_group == "GRP1")
sum(Microphysogobio_longidorsalis_dat$site_group == "GRP2")
sum(Microphysogobio_longidorsalis_dat$site_group == "GRP3")



#'### Global Model

#'##### glm
Microphysogobio_longidorsalis_glm01 <- glm(Microphysogobio_longidorsalis_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                           family=binomial(link = "logit"), data= Microphysogobio_longidorsalis_dat, na.action ="na.fail")

summary(Microphysogobio_longidorsalis_glm01)

vif(Microphysogobio_longidorsalis_glm01)



#'##### glmm
#'- combination of site and season are used as a random effect 
Microphysogobio_longidorsalis_glmm01 <- glmmTMB(Microphysogobio_longidorsalis_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                  (1|MBSNCD/site) + (1|season), family=binomial(link = "logit"), data= Microphysogobio_longidorsalis_dat, na.action ="na.fail")

summary(Microphysogobio_longidorsalis_glmm01)










#'- Model validation
r.squaredGLMM(Microphysogobio_longidorsalis_glmm01)
cor.test(predict(Microphysogobio_longidorsalis_glmm01, type = "response", allow.new.levels=TRUE), Microphysogobio_longidorsalis_dat$Microphysogobio_longidorsalis_PA)

#'- Microphysogobio_longidorsalis depth
Microphysogobio_longidorsalis_glmm01_dep_fit <- effect_plot(Microphysogobio_longidorsalis_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Microphysogobio_longidorsalis_glmm01_dep_fit_PA <- as.data.frame(Microphysogobio_longidorsalis_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Microphysogobio_longidorsalis_PA"]])
Microphysogobio_longidorsalis_glmm01_dep_fit_dep <- as.data.frame(Microphysogobio_longidorsalis_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Microphysogobio_longidorsalis_dep_res <- cbind(Microphysogobio_longidorsalis_glmm01_dep_fit_PA, Microphysogobio_longidorsalis_glmm01_dep_fit_dep)
colnames(Microphysogobio_longidorsalis_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Microphysogobio_longidorsalis_dep_res <- Microphysogobio_longidorsalis_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Microphysogobio_longidorsalis_depth_fig <- ggplot(Microphysogobio_longidorsalis_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Microphysogobio_longidorsalis depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Microphysogobio_longidorsalis_depth_fig

#'- Microphysogobio_longidorsalis velocity
Microphysogobio_longidorsalis_glmm01_vel_fit <- effect_plot(Microphysogobio_longidorsalis_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Microphysogobio_longidorsalis_glmm01_vel_fit_PA <- as.data.frame(Microphysogobio_longidorsalis_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Microphysogobio_longidorsalis_PA"]])
Microphysogobio_longidorsalis_glmm01_vel_fit_vel <- as.data.frame(Microphysogobio_longidorsalis_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Microphysogobio_longidorsalis_vel_res <- cbind(Microphysogobio_longidorsalis_glmm01_vel_fit_PA, Microphysogobio_longidorsalis_glmm01_vel_fit_vel)
colnames(Microphysogobio_longidorsalis_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Microphysogobio_longidorsalis_vel_res <- Microphysogobio_longidorsalis_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Microphysogobio_longidorsalis_velocity_fig <- ggplot(Microphysogobio_longidorsalis_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Microphysogobio_longidorsalis Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Microphysogobio_longidorsalis_velocity_fig

#'- Microphysogobio_longidorsalis substrate
Microphysogobio_longidorsalis_glmm01_sub_fit <- effect_plot(Microphysogobio_longidorsalis_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Microphysogobio_longidorsalis_glmm01_sub_fit_PA <- as.data.frame(Microphysogobio_longidorsalis_glmm01_sub_fit[["data"]][["Microphysogobio_longidorsalis_PA"]])
Microphysogobio_longidorsalis_glmm01_sub_fit_sub <- as.data.frame(Microphysogobio_longidorsalis_glmm01_sub_fit[["data"]][["substrate"]])
Microphysogobio_longidorsalis_sub_res <- cbind(Microphysogobio_longidorsalis_glmm01_sub_fit_PA, Microphysogobio_longidorsalis_glmm01_sub_fit_sub)
colnames(Microphysogobio_longidorsalis_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Microphysogobio_longidorsalis_sub_res <- Microphysogobio_longidorsalis_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Microphysogobio_longidorsalis_sub_res$sub <- factor(Microphysogobio_longidorsalis_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                          'pebble', 'cobble', 'boulder'))

#'- plot
Microphysogobio_longidorsalis_substrate_fig <- ggplot(Microphysogobio_longidorsalis_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Microphysogobio_longidorsalis Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Microphysogobio_longidorsalis_substrate_fig


##### Microphysogobio_longidorsalis group2 model #####

#'## select Group2
Microphysogobio_longidorsalis_dat03 <- Microphysogobio_longidorsalis_dat %>%
  filter(site_group == "GRP2")
sum( Microphysogobio_longidorsalis_dat03$Microphysogobio_longidorsalis_PA=="1") # Presence data(n = 21)
#'##### glm
Microphysogobio_longidorsalis_glm03 <- glm(Microphysogobio_longidorsalis_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate,
                                           family=binomial(link = "logit"), data= Microphysogobio_longidorsalis_dat03, na.action ="na.fail")

summary(Microphysogobio_longidorsalis_glm03)

vif(Microphysogobio_longidorsalis_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect
Microphysogobio_longidorsalis_glmm03 <- glmmTMB(Microphysogobio_longidorsalis_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate +
                                                  (1|MBSNCD/site) + (1|season),family=binomial(link = "logit"), data= Microphysogobio_longidorsalis_dat03, na.action ="na.fail")

summary(Microphysogobio_longidorsalis_glmm03)

#'- Model validation
r.squaredGLMM(Microphysogobio_longidorsalis_glmm03)
cor.test(predict(Microphysogobio_longidorsalis_glmm03, type = "response", allow.new.levels=TRUE), Microphysogobio_longidorsalis_dat03$Microphysogobio_longidorsalis_PA)

#'- Microphysogobio_longidorsalis depth
Microphysogobio_longidorsalis_glmm03_dep_fit <- effect_plot(Microphysogobio_longidorsalis_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Microphysogobio_longidorsalis_glmm03_dep_fit_PA <- as.data.frame(Microphysogobio_longidorsalis_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Microphysogobio_longidorsalis_PA"]])
Microphysogobio_longidorsalis_glmm03_dep_fit_depth <- as.data.frame(Microphysogobio_longidorsalis_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Microphysogobio_longidorsalis_dep_res_grp2 <- cbind(Microphysogobio_longidorsalis_glmm03_dep_fit_PA, Microphysogobio_longidorsalis_glmm03_dep_fit_depth)
colnames(Microphysogobio_longidorsalis_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Microphysogobio_longidorsalis_dep_res_grp2 <- Microphysogobio_longidorsalis_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Microphysogobio_longidorsalis_depth_fig2 <- ggplot(Microphysogobio_longidorsalis_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) +
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +
  labs(title= "Microphysogobio_longidorsalis depth HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Microphysogobio_longidorsalis_depth_fig2

#'- Microphysogobio_longidorsalis velocity
Microphysogobio_longidorsalis_glmm03_vel_fit <- effect_plot(Microphysogobio_longidorsalis_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Microphysogobio_longidorsalis_glmm03_vel_fit_PA <- as.data.frame(Microphysogobio_longidorsalis_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Microphysogobio_longidorsalis_PA"]])
Microphysogobio_longidorsalis_glmm03_vel_fit_vel <- as.data.frame(Microphysogobio_longidorsalis_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Microphysogobio_longidorsalis_vel_res_grp2 <- cbind(Microphysogobio_longidorsalis_glmm03_vel_fit_PA, Microphysogobio_longidorsalis_glmm03_vel_fit_vel)
colnames(Microphysogobio_longidorsalis_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Microphysogobio_longidorsalis_vel_res_grp2 <- Microphysogobio_longidorsalis_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Microphysogobio_longidorsalis_velocity_fig2 <- ggplot(Microphysogobio_longidorsalis_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) +
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +
  labs(title= "Microphysogobio_longidorsalis Velocity HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Microphysogobio_longidorsalis_velocity_fig2

#'- Microphysogobio_longidorsalis substrate
Microphysogobio_longidorsalis_glmm03_sub_fit <- effect_plot(Microphysogobio_longidorsalis_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Microphysogobio_longidorsalis_glmm03_sub_fit_PA <- as.data.frame(Microphysogobio_longidorsalis_glmm03_sub_fit[["data"]][["Microphysogobio_longidorsalis_PA"]])
Microphysogobio_longidorsalis_glmm03_sub_fit_sub <- as.data.frame(Microphysogobio_longidorsalis_glmm03_sub_fit[["data"]][["substrate"]])
Microphysogobio_longidorsalis_sub_res_grp2 <- cbind(Microphysogobio_longidorsalis_glmm03_sub_fit_PA, Microphysogobio_longidorsalis_glmm03_sub_fit_sub)
colnames(Microphysogobio_longidorsalis_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Microphysogobio_longidorsalis_sub_res_grp2 <- Microphysogobio_longidorsalis_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Microphysogobio_longidorsalis_sub_res_grp2$sub <- factor(Microphysogobio_longidorsalis_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                                    'pebble', 'cobble', 'boulder'))

#'- plot
Microphysogobio_longidorsalis_substrate_fig2 <- ggplot(Microphysogobio_longidorsalis_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +
  labs(title= "Microphysogobio_longidorsalis Substrate HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Microphysogobio_longidorsalis_substrate_fig2





#'## Combined figures

#'- merge mulitple effect plot
Microphysogobio_longidorsalis_dep_res_1 <- Microphysogobio_longidorsalis_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)


Microphysogobio_longidorsalis_dep_res_grp2_1 <- Microphysogobio_longidorsalis_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)

Microphysogobio_longidorsalis_dep_res_comb <- rbind(Microphysogobio_longidorsalis_dep_res_1, 
                                                    Microphysogobio_longidorsalis_dep_res_grp2_1)


#'- depth plot
Microphysogobio_longidorsalis_depth_fig4 <- ggplot(Microphysogobio_longidorsalis_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)", limits = c(0,2) , breaks = seq(0, 2, 0.5)) +     
  labs(title= "Microphysogobio longidorsalis Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Microphysogobio_longidorsalis_depth_fig4 <-Microphysogobio_longidorsalis_depth_fig4+scale_color_manual(values=c("black", "red"))
Microphysogobio_longidorsalis_depth_fig4




#'- merge mulitple effect plot
Microphysogobio_longidorsalis_vel_res_1 <- Microphysogobio_longidorsalis_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)


Microphysogobio_longidorsalis_vel_res_grp2_1 <- Microphysogobio_longidorsalis_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)

Microphysogobio_longidorsalis_vel_res_comb <- rbind(Microphysogobio_longidorsalis_vel_res_1, 
                                                    Microphysogobio_longidorsalis_vel_res_grp2_1)


#'- velocity plot
Microphysogobio_longidorsalis_velocity_fig4 <- ggplot(Microphysogobio_longidorsalis_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)", limits = c(0,3.1), breaks = seq(0, 3, 0.5)) +     
  labs(title= "Microphysogobio_longidorsalis Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Microphysogobio_longidorsalis_velocity_fig4 <- Microphysogobio_longidorsalis_velocity_fig4 + scale_color_manual(values=c("black", "red"))
Microphysogobio_longidorsalis_velocity_fig4



#'- merge mulitple effect plot
Microphysogobio_longidorsalis_sub_res_1 <- Microphysogobio_longidorsalis_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)

Microphysogobio_longidorsalis_sub_res_grp2_1 <- Microphysogobio_longidorsalis_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)



Microphysogobio_longidorsalis_sub_res_comb <- rbind(Microphysogobio_longidorsalis_sub_res_1, 
                                                    Microphysogobio_longidorsalis_sub_res_grp2_1)

#'- substrate plot
Microphysogobio_longidorsalis_substrate_fig4 <- ggplot(Microphysogobio_longidorsalis_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Microphysogobio_longidorsalis Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Microphysogobio_longidorsalis_substrate_fig4 <- Microphysogobio_longidorsalis_substrate_fig4 + scale_fill_manual(values=c("black","red"))
Microphysogobio_longidorsalis_substrate_fig4


#'## combine figures

ggarrange(Microphysogobio_longidorsalis_depth_fig4, Microphysogobio_longidorsalis_velocity_fig4, Microphysogobio_longidorsalis_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 






##### Hemibarbus_longirostris model #####

#'## (2) Hemibarbus_longirostris model

#'### data management
Hemibarbus_longirostris_dat <- HSI_DB2 %>%
  mutate(Hemibarbus_longirostris_PA = ifelse(species == "Hemibarbus_longirostris", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Hemibarbus_longirostris_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Hemibarbus_longirostris_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Hemibarbus_longirostris_PA, depth, velocity, substrate, MBSNCD) 

head(Hemibarbus_longirostris_dat)


#'## remove site where species did not occur
Hemibarbus_longirostris_Pre <- HSI_DB2 %>%
  mutate(Hemibarbus_longirostris_PA = ifelse(species == "Hemibarbus_longirostris", 1, 0)) %>%
  aggregate(Hemibarbus_longirostris_PA ~ site, "sum") %>%
  mutate(Hemibarbus_longirostris_Pre = if_else(Hemibarbus_longirostris_PA > 0, "Yes", "No")) %>%
  select(site, Hemibarbus_longirostris_Pre)

Hemibarbus_longirostris_dat <- Hemibarbus_longirostris_dat %>%
  left_join(Hemibarbus_longirostris_Pre, by = "site") %>%
  filter(Hemibarbus_longirostris_Pre == "Yes") %>%
  as_tibble()



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

#check_species_site_type
list(unique(Hemibarbus_longirostris_dat$site_group))
sum(Hemibarbus_longirostris_dat$site_group == "GRP1")
sum(Hemibarbus_longirostris_dat$site_group == "GRP2")
sum(Hemibarbus_longirostris_dat$site_group == "GRP3")



#'### Global Model

#'##### glm
Hemibarbus_longirostris_glm01 <- glm(Hemibarbus_longirostris_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                     family=binomial(link = "logit"), data= Hemibarbus_longirostris_dat, na.action ="na.fail")

summary(Hemibarbus_longirostris_glm01)

vif(Hemibarbus_longirostris_glm01)



#'##### glmm
#'- combination of site and season are used as a random effect 
Hemibarbus_longirostris_glmm01 <- glmmTMB(Hemibarbus_longirostris_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                            (1|MBSNCD/site) + (1|season), family=binomial(link = "logit"), data= Hemibarbus_longirostris_dat, na.action ="na.fail")

summary(Hemibarbus_longirostris_glmm01)










#'- Model validation
r.squaredGLMM(Hemibarbus_longirostris_glmm01)
cor.test(predict(Hemibarbus_longirostris_glmm01, type = "response", allow.new.levels=TRUE), Hemibarbus_longirostris_dat$Hemibarbus_longirostris_PA)

#'- Hemibarbus_longirostris depth
Hemibarbus_longirostris_glmm01_dep_fit <- effect_plot(Hemibarbus_longirostris_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Hemibarbus_longirostris_glmm01_dep_fit_PA <- as.data.frame(Hemibarbus_longirostris_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_longirostris_PA"]])
Hemibarbus_longirostris_glmm01_dep_fit_dep <- as.data.frame(Hemibarbus_longirostris_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Hemibarbus_longirostris_dep_res <- cbind(Hemibarbus_longirostris_glmm01_dep_fit_PA, Hemibarbus_longirostris_glmm01_dep_fit_dep)
colnames(Hemibarbus_longirostris_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Hemibarbus_longirostris_dep_res <- Hemibarbus_longirostris_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Hemibarbus_longirostris_depth_fig <- ggplot(Hemibarbus_longirostris_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_longirostris depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_longirostris_depth_fig

#'- Hemibarbus_longirostris velocity
Hemibarbus_longirostris_glmm01_vel_fit <- effect_plot(Hemibarbus_longirostris_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Hemibarbus_longirostris_glmm01_vel_fit_PA <- as.data.frame(Hemibarbus_longirostris_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_longirostris_PA"]])
Hemibarbus_longirostris_glmm01_vel_fit_vel <- as.data.frame(Hemibarbus_longirostris_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Hemibarbus_longirostris_vel_res <- cbind(Hemibarbus_longirostris_glmm01_vel_fit_PA, Hemibarbus_longirostris_glmm01_vel_fit_vel)
colnames(Hemibarbus_longirostris_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Hemibarbus_longirostris_vel_res <- Hemibarbus_longirostris_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Hemibarbus_longirostris_velocity_fig <- ggplot(Hemibarbus_longirostris_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_longirostris Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_longirostris_velocity_fig

#'- Hemibarbus_longirostris substrate
Hemibarbus_longirostris_glmm01_sub_fit <- effect_plot(Hemibarbus_longirostris_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Hemibarbus_longirostris_glmm01_sub_fit_PA <- as.data.frame(Hemibarbus_longirostris_glmm01_sub_fit[["data"]][["Hemibarbus_longirostris_PA"]])
Hemibarbus_longirostris_glmm01_sub_fit_sub <- as.data.frame(Hemibarbus_longirostris_glmm01_sub_fit[["data"]][["substrate"]])
Hemibarbus_longirostris_sub_res <- cbind(Hemibarbus_longirostris_glmm01_sub_fit_PA, Hemibarbus_longirostris_glmm01_sub_fit_sub)
colnames(Hemibarbus_longirostris_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Hemibarbus_longirostris_sub_res <- Hemibarbus_longirostris_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Hemibarbus_longirostris_sub_res$sub <- factor(Hemibarbus_longirostris_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                              'pebble', 'cobble', 'boulder'))

#'- plot
Hemibarbus_longirostris_substrate_fig <- ggplot(Hemibarbus_longirostris_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Hemibarbus_longirostris Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_longirostris_substrate_fig



#'## (1) Hemibarbus_longirostris model
##### Hemibarbus_longirostris group1 model #####

#'## select Group1

Hemibarbus_longirostris_dat02 <- Hemibarbus_longirostris_dat %>%
  filter(site_group == "GRP1")

#'##### glm
Hemibarbus_longirostris_glm02 <- glm(Hemibarbus_longirostris_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                     family=binomial(link = "logit"), data= Hemibarbus_longirostris_dat02, na.action ="na.fail")

summary(Hemibarbus_longirostris_glm02)

vif(Hemibarbus_longirostris_glm02)
#'##### glmm
#'- combination of site and season are used as a random effect ``
Hemibarbus_longirostris_glmm02 <- glmmTMB(Hemibarbus_longirostris_PA ~ poly(depth, 1) + poly(velocity, 2) + substrate + 
                                            (1|site) + (1|season), family=binomial(link = "logit"), data= Hemibarbus_longirostris_dat02, na.action ="na.fail")

summary(Hemibarbus_longirostris_glmm02)



#'- Model validation
r.squaredGLMM(Hemibarbus_longirostris_glmm02)
cor.test(predict(Hemibarbus_longirostris_glmm02, type = "response", allow.new.levels=TRUE), Hemibarbus_longirostris_dat02$Hemibarbus_longirostris_PA)

#'- Hemibarbus_longirostris depth
Hemibarbus_longirostris_glmm02_dep_fit <- effect_plot(Hemibarbus_longirostris_glmm02, pred = depth, interval = TRUE, plot.points = FALSE)
Hemibarbus_longirostris_glmm02_dep_fit_PA <- as.data.frame(Hemibarbus_longirostris_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_longirostris_PA"]])
Hemibarbus_longirostris_glmm02_dep_fit_depth <- as.data.frame(Hemibarbus_longirostris_glmm02_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Hemibarbus_longirostris_dep_res_grp1 <- cbind(Hemibarbus_longirostris_glmm02_dep_fit_PA, Hemibarbus_longirostris_glmm02_dep_fit_depth)
colnames(Hemibarbus_longirostris_dep_res_grp1) <- c("prob", "dep")

#'- normalize 0 to 1
Hemibarbus_longirostris_dep_res_grp1 <- Hemibarbus_longirostris_dep_res_grp1 %>%
  mutate(HSI_dep_grp1 = scale_values(prob))

#'- plot
Hemibarbus_longirostris_depth_fig1 <- ggplot(Hemibarbus_longirostris_dep_res_grp1, aes(x = dep, y = HSI_dep_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_longirostris depth HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_longirostris_depth_fig1

#'- Hemibarbus_longirostris velocity
Hemibarbus_longirostris_glmm02_vel_fit <- effect_plot(Hemibarbus_longirostris_glmm02, pred = velocity, interval = TRUE, plot.points = FALSE)
Hemibarbus_longirostris_glmm02_vel_fit_PA <- as.data.frame(Hemibarbus_longirostris_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_longirostris_PA"]])
Hemibarbus_longirostris_glmm02_vel_fit_vel <- as.data.frame(Hemibarbus_longirostris_glmm02_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Hemibarbus_longirostris_vel_res_grp1 <- cbind(Hemibarbus_longirostris_glmm02_vel_fit_PA, Hemibarbus_longirostris_glmm02_vel_fit_vel)
colnames(Hemibarbus_longirostris_vel_res_grp1) <- c("prob", "vel")

#'- normalize 0 to 1
Hemibarbus_longirostris_vel_res_grp1 <- Hemibarbus_longirostris_vel_res_grp1 %>%
  mutate(HSI_vel_grp1 = scale_values(prob))

#'- plot
Hemibarbus_longirostris_velocity_fig1 <- ggplot(Hemibarbus_longirostris_vel_res_grp1, aes(x = vel, y = HSI_vel_grp1)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_longirostris Velocity HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_longirostris_velocity_fig1

#'- Hemibarbus_longirostris substrate
Hemibarbus_longirostris_glmm02_sub_fit <- effect_plot(Hemibarbus_longirostris_glmm02, pred = substrate, interval = FALSE, plot.points = FALSE)
Hemibarbus_longirostris_glmm02_sub_fit_PA <- as.data.frame(Hemibarbus_longirostris_glmm02_sub_fit[["data"]][["Hemibarbus_longirostris_PA"]])
Hemibarbus_longirostris_glmm02_sub_fit_sub <- as.data.frame(Hemibarbus_longirostris_glmm02_sub_fit[["data"]][["substrate"]])
Hemibarbus_longirostris_sub_res_grp1 <- cbind(Hemibarbus_longirostris_glmm02_sub_fit_PA, Hemibarbus_longirostris_glmm02_sub_fit_sub)
colnames(Hemibarbus_longirostris_sub_res_grp1) <- c("prob", "sub")

#'- normalize 0 to 1
Hemibarbus_longirostris_sub_res_grp1 <- Hemibarbus_longirostris_sub_res_grp1 %>%
  mutate(HSI_sub_grp1 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Hemibarbus_longirostris_sub_res_grp1$sub <- factor(Hemibarbus_longirostris_sub_res_grp1$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                        'pebble', 'cobble', 'boulder'))

#'- plot
Hemibarbus_longirostris_substrate_fig1 <- ggplot(Hemibarbus_longirostris_sub_res_grp1, aes(x = sub, y = HSI_sub_grp1)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Hemibarbus_longirostris Substrate HSI Group1") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Hemibarbus_longirostris_substrate_fig1



##### Hemibarbus_longirostris group2 model #####

#'## select Group2
Hemibarbus_longirostris_dat03 <- Hemibarbus_longirostris_dat %>%
  filter(site_group == "GRP2")
sum( Hemibarbus_longirostris_dat03$Hemibarbus_longirostris_PA=="1") # Presence data(n = 21)

#'##### glm
Hemibarbus_longirostris_glm03 <- glm(Hemibarbus_longirostris_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate,
                                     family=binomial(link = "logit"), data= Hemibarbus_longirostris_dat03, na.action ="na.fail")

summary(Hemibarbus_longirostris_glm03)

vif(Hemibarbus_labeo_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect
Hemibarbus_longirostris_glmm03 <- glmmTMB(Hemibarbus_longirostris_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate +
                                            (1|MBSNCD/site) + (1|season),family=binomial(link = "logit"), data= Hemibarbus_longirostris_dat03, na.action ="na.fail")

summary(Hemibarbus_longirostris_glmm03)

#'- Model validation
r.squaredGLMM(Hemibarbus_longirostris_glmm03)
cor.test(predict(Hemibarbus_longirostris_glmm03, type = "response", allow.new.levels=TRUE), Hemibarbus_longirostris_dat03$Hemibarbus_longirostris_PA)

#'- Hemibarbus_longirostris depth
Hemibarbus_longirostris_glmm03_dep_fit <- effect_plot(Hemibarbus_longirostris_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Hemibarbus_longirostris_glmm03_dep_fit_PA <- as.data.frame(Hemibarbus_longirostris_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_longirostris_PA"]])
Hemibarbus_longirostris_glmm03_dep_fit_depth <- as.data.frame(Hemibarbus_longirostris_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Hemibarbus_longirostris_dep_res_grp2 <- cbind(Hemibarbus_longirostris_glmm03_dep_fit_PA, Hemibarbus_longirostris_glmm03_dep_fit_depth)
colnames(Hemibarbus_longirostris_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Hemibarbus_longirostris_dep_res_grp2 <- Hemibarbus_longirostris_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Hemibarbus_longirostris_depth_fig2 <- ggplot(Hemibarbus_longirostris_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) +
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +
  labs(title= "Hemibarbus_longirostris depth HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Hemibarbus_longirostris_depth_fig2

#'- Hemibarbus_longirostris velocity
Hemibarbus_longirostris_glmm03_vel_fit <- effect_plot(Hemibarbus_longirostris_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Hemibarbus_longirostris_glmm03_vel_fit_PA <- as.data.frame(Hemibarbus_longirostris_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Hemibarbus_longirostris_PA"]])
Hemibarbus_longirostris_glmm03_vel_fit_vel <- as.data.frame(Hemibarbus_longirostris_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Hemibarbus_longirostris_vel_res_grp2 <- cbind(Hemibarbus_longirostris_glmm03_vel_fit_PA, Hemibarbus_longirostris_glmm03_vel_fit_vel)
colnames(Hemibarbus_longirostris_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Hemibarbus_longirostris_vel_res_grp2 <- Hemibarbus_longirostris_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Hemibarbus_longirostris_velocity_fig2 <- ggplot(Hemibarbus_longirostris_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) +
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +
  labs(title= "Hemibarbus_longirostris Velocity HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Hemibarbus_longirostris_velocity_fig2

#'- Hemibarbus_longirostris substrate
Hemibarbus_longirostris_glmm03_sub_fit <- effect_plot(Hemibarbus_longirostris_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Hemibarbus_longirostris_glmm03_sub_fit_PA <- as.data.frame(Hemibarbus_longirostris_glmm03_sub_fit[["data"]][["Hemibarbus_longirostris_PA"]])
Hemibarbus_longirostris_glmm03_sub_fit_sub <- as.data.frame(Hemibarbus_longirostris_glmm03_sub_fit[["data"]][["substrate"]])
Hemibarbus_longirostris_sub_res_grp2 <- cbind(Hemibarbus_longirostris_glmm03_sub_fit_PA, Hemibarbus_longirostris_glmm03_sub_fit_sub)
colnames(Hemibarbus_longirostris_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Hemibarbus_longirostris_sub_res_grp2 <- Hemibarbus_longirostris_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Hemibarbus_longirostris_sub_res_grp2$sub <- factor(Hemibarbus_longirostris_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                        'pebble', 'cobble', 'boulder'))

#'- plot
Hemibarbus_longirostris_substrate_fig2 <- ggplot(Hemibarbus_longirostris_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +
  labs(title= "Hemibarbus_longirostris Substrate HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Hemibarbus_longirostris_substrate_fig2





#'## Combined figures

#'- merge mulitple effect plot
Hemibarbus_longirostris_dep_res_1 <- Hemibarbus_longirostris_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)

Hemibarbus_longirostris_dep_res_grp1_1 <- Hemibarbus_longirostris_dep_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp1)

Hemibarbus_longirostris_dep_res_grp2_1 <- Hemibarbus_longirostris_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)




Hemibarbus_longirostris_dep_res_comb <- rbind(Hemibarbus_longirostris_dep_res_1, Hemibarbus_longirostris_dep_res_grp1_1,
                                              Hemibarbus_longirostris_dep_res_grp2_1)


#'- depth plot
Hemibarbus_longirostris_depth_fig4 <- ggplot(Hemibarbus_longirostris_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_longirostris Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Hemibarbus_longirostris_depth_fig4 <-Hemibarbus_longirostris_depth_fig4+scale_color_manual(values=c("black","blue3", "red", "green3"))
Hemibarbus_longirostris_depth_fig4




#'- merge mulitple effect plot
Hemibarbus_longirostris_vel_res_1 <- Hemibarbus_longirostris_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)

Hemibarbus_longirostris_vel_res_grp1_1 <- Hemibarbus_longirostris_vel_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp1)

Hemibarbus_longirostris_vel_res_grp2_1 <- Hemibarbus_longirostris_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)



Hemibarbus_longirostris_vel_res_comb <- rbind(Hemibarbus_longirostris_vel_res_1, Hemibarbus_longirostris_vel_res_grp1_1,
                                              Hemibarbus_longirostris_vel_res_grp2_1)


#'- velocity plot
Hemibarbus_longirostris_velocity_fig4 <- ggplot(Hemibarbus_longirostris_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Hemibarbus_longirostris Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Hemibarbus_longirostris_velocity_fig4 <- Hemibarbus_longirostris_velocity_fig4 + scale_color_manual(values=c("black","blue3", "red", "green3"))
Hemibarbus_longirostris_velocity_fig4



#'- merge mulitple effect plot
Hemibarbus_longirostris_sub_res_1 <- Hemibarbus_longirostris_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)

Hemibarbus_longirostris_sub_res_grp1_1 <- Hemibarbus_longirostris_sub_res_grp1 %>%
  mutate(group = "Group1") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp1)

Hemibarbus_longirostris_sub_res_grp2_1 <- Hemibarbus_longirostris_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)


Hemibarbus_longirostris_sub_res_comb <- rbind(Hemibarbus_longirostris_sub_res_1, Hemibarbus_longirostris_sub_res_grp1_1,
                                              Hemibarbus_longirostris_sub_res_grp2_1)

#'- substrate plot
Hemibarbus_longirostris_substrate_fig4 <- ggplot(Hemibarbus_longirostris_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Hemibarbus_longirostris Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Hemibarbus_longirostris_substrate_fig4 <- Hemibarbus_longirostris_substrate_fig4 + scale_fill_manual(values=c("black","blue3", "red", "green3"))
Hemibarbus_longirostris_substrate_fig4


#'## combine figures

ggarrange(Hemibarbus_longirostris_depth_fig4, Hemibarbus_longirostris_velocity_fig4, Hemibarbus_longirostris_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv")





##### Sarcocheilichthys_variegatus_wakiyae model #####

#'## (2) Sarcocheilichthys_variegatus_wakiyae model

#'### data management
Sarcocheilichthys_variegatus_wakiyae_dat <- HSI_DB2 %>%
  mutate(Sarcocheilichthys_variegatus_wakiyae_PA = ifelse(species == "Sarcocheilichthys_variegatus_wakiyae", 1, 0)) %>% # create column based on P/A of target species
  distinct(site_by_season, sample_number, Sarcocheilichthys_variegatus_wakiyae_PA, .keep_all = TRUE) %>% # remove duplicate based on target species binary data
  arrange(basin, site_by_season, sample_number, desc(Sarcocheilichthys_variegatus_wakiyae_PA)) %>% # arrange order
  group_by(site_by_season, sample_number) %>% 
  slice(1) %>% # keep only first row by occasion and section
  select(basin, site_group, site_by_season, site, date, season, sample_number, Sarcocheilichthys_variegatus_wakiyae_PA, depth, velocity, substrate, MBSNCD) 

head(Sarcocheilichthys_variegatus_wakiyae_dat)


#'## remove site where species did not occur
Sarcocheilichthys_variegatus_wakiyae_Pre <- HSI_DB2 %>%
  mutate(Sarcocheilichthys_variegatus_wakiyae_PA = ifelse(species == "Sarcocheilichthys_variegatus_wakiyae", 1, 0)) %>%
  aggregate(Sarcocheilichthys_variegatus_wakiyae_PA ~ site, "sum") %>%
  mutate(Sarcocheilichthys_variegatus_wakiyae_Pre = if_else(Sarcocheilichthys_variegatus_wakiyae_PA > 0, "Yes", "No")) %>%
  select(site, Sarcocheilichthys_variegatus_wakiyae_Pre)

Sarcocheilichthys_variegatus_wakiyae_dat <- Sarcocheilichthys_variegatus_wakiyae_dat %>%
  left_join(Sarcocheilichthys_variegatus_wakiyae_Pre, by = "site") %>%
  filter(Sarcocheilichthys_variegatus_wakiyae_Pre == "Yes") %>%
  as_tibble()



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

#check_species_site_type
list(unique(Sarcocheilichthys_variegatus_wakiyae_dat$site_group))
sum(Sarcocheilichthys_variegatus_wakiyae_dat$site_group == "GRP1")
sum(Sarcocheilichthys_variegatus_wakiyae_dat$site_group == "GRP2")
sum(Sarcocheilichthys_variegatus_wakiyae_dat$site_group == "GRP3")



#'### Global Model

#'##### glm
Sarcocheilichthys_variegatus_wakiyae_glm01 <- glm(Sarcocheilichthys_variegatus_wakiyae_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate, 
                                                  family=binomial(link = "logit"), data= Sarcocheilichthys_variegatus_wakiyae_dat, na.action ="na.fail")

summary(Sarcocheilichthys_variegatus_wakiyae_glm01)

vif(Sarcocheilichthys_variegatus_wakiyae_glm01)



#'##### glmm
#'- combination of site and season are used as a random effect 
Sarcocheilichthys_variegatus_wakiyae_glmm01 <- glmmTMB(Sarcocheilichthys_variegatus_wakiyae_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate + 
                                                         (1|MBSNCD/site) + (1|season), family=binomial(link = "logit"), data= Sarcocheilichthys_variegatus_wakiyae_dat, na.action ="na.fail")

summary(Sarcocheilichthys_variegatus_wakiyae_glmm01)










#'- Model validation
r.squaredGLMM(Sarcocheilichthys_variegatus_wakiyae_glmm01)
cor.test(predict(Sarcocheilichthys_variegatus_wakiyae_glmm01, type = "response", allow.new.levels=TRUE), Sarcocheilichthys_variegatus_wakiyae_dat$Sarcocheilichthys_variegatus_wakiyae_PA)

#'- Sarcocheilichthys_variegatus_wakiyae depth
Sarcocheilichthys_variegatus_wakiyae_glmm01_dep_fit <- effect_plot(Sarcocheilichthys_variegatus_wakiyae_glmm01, pred = depth, interval = TRUE, plot.points = FALSE)
Sarcocheilichthys_variegatus_wakiyae_glmm01_dep_fit_PA <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["Sarcocheilichthys_variegatus_wakiyae_PA"]])
Sarcocheilichthys_variegatus_wakiyae_glmm01_dep_fit_dep <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm01_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Sarcocheilichthys_variegatus_wakiyae_dep_res <- cbind(Sarcocheilichthys_variegatus_wakiyae_glmm01_dep_fit_PA, Sarcocheilichthys_variegatus_wakiyae_glmm01_dep_fit_dep)
colnames(Sarcocheilichthys_variegatus_wakiyae_dep_res) <- c("prob", "dep")

#'- normalize 0 to 1
Sarcocheilichthys_variegatus_wakiyae_dep_res <- Sarcocheilichthys_variegatus_wakiyae_dep_res %>%
  mutate(HSI_dep = scale_values(prob))

#'- plot
Sarcocheilichthys_variegatus_wakiyae_depth_fig <- ggplot(Sarcocheilichthys_variegatus_wakiyae_dep_res, aes(x = dep, y = HSI_dep)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Sarcocheilichthys_variegatus_wakiyae depth HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Sarcocheilichthys_variegatus_wakiyae_depth_fig

#'- Sarcocheilichthys_variegatus_wakiyae velocity
Sarcocheilichthys_variegatus_wakiyae_glmm01_vel_fit <- effect_plot(Sarcocheilichthys_variegatus_wakiyae_glmm01, pred = velocity, interval = TRUE, plot.points = FALSE)
Sarcocheilichthys_variegatus_wakiyae_glmm01_vel_fit_PA <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["Sarcocheilichthys_variegatus_wakiyae_PA"]])
Sarcocheilichthys_variegatus_wakiyae_glmm01_vel_fit_vel <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm01_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Sarcocheilichthys_variegatus_wakiyae_vel_res <- cbind(Sarcocheilichthys_variegatus_wakiyae_glmm01_vel_fit_PA, Sarcocheilichthys_variegatus_wakiyae_glmm01_vel_fit_vel)
colnames(Sarcocheilichthys_variegatus_wakiyae_vel_res) <- c("prob", "vel")

#'- normalize 0 to 1
Sarcocheilichthys_variegatus_wakiyae_vel_res <- Sarcocheilichthys_variegatus_wakiyae_vel_res %>%
  mutate(HSI_vel = scale_values(prob))

#'- plot
Sarcocheilichthys_variegatus_wakiyae_velocity_fig <- ggplot(Sarcocheilichthys_variegatus_wakiyae_vel_res, aes(x = vel, y = HSI_vel)) +
  geom_line(size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Sarcocheilichthys_variegatus_wakiyae Velocity HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Sarcocheilichthys_variegatus_wakiyae_velocity_fig

#'- Sarcocheilichthys_variegatus_wakiyae substrate
Sarcocheilichthys_variegatus_wakiyae_glmm01_sub_fit <- effect_plot(Sarcocheilichthys_variegatus_wakiyae_glmm01, pred = substrate, interval = FALSE, plot.points = FALSE)
Sarcocheilichthys_variegatus_wakiyae_glmm01_sub_fit_PA <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm01_sub_fit[["data"]][["Sarcocheilichthys_variegatus_wakiyae_PA"]])
Sarcocheilichthys_variegatus_wakiyae_glmm01_sub_fit_sub <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm01_sub_fit[["data"]][["substrate"]])
Sarcocheilichthys_variegatus_wakiyae_sub_res <- cbind(Sarcocheilichthys_variegatus_wakiyae_glmm01_sub_fit_PA, Sarcocheilichthys_variegatus_wakiyae_glmm01_sub_fit_sub)
colnames(Sarcocheilichthys_variegatus_wakiyae_sub_res) <- c("prob", "sub")

#'- normalize 0 to 1
Sarcocheilichthys_variegatus_wakiyae_sub_res <- Sarcocheilichthys_variegatus_wakiyae_sub_res %>%
  mutate(HSI_sub = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Sarcocheilichthys_variegatus_wakiyae_sub_res$sub <- factor(Sarcocheilichthys_variegatus_wakiyae_sub_res$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                                        'pebble', 'cobble', 'boulder'))

#'- plot
Sarcocheilichthys_variegatus_wakiyae_substrate_fig <- ggplot(Sarcocheilichthys_variegatus_wakiyae_sub_res, aes(x = sub, y = HSI_sub)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Sarcocheilichthys_variegatus_wakiyae Substrate HSI") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 

Sarcocheilichthys_variegatus_wakiyae_substrate_fig



##### Sarcocheilichthys_variegatus_wakiyae group2 model #####

#'## select Group2
Sarcocheilichthys_variegatus_wakiyae_dat03 <- Sarcocheilichthys_variegatus_wakiyae_dat %>%
  filter(site_group == "GRP2")
sum( Sarcocheilichthys_variegatus_wakiyae_dat03$Sarcocheilichthys_variegatus_wakiyae_PA=="1") # Presence data(n = 21)
#'##### glm
Sarcocheilichthys_variegatus_wakiyae_glm03 <- glm(Sarcocheilichthys_variegatus_wakiyae_PA ~ scale(poly(depth, 2)) + scale(poly(velocity, 2)) + substrate,
                                                  family=binomial(link = "logit"), data= Sarcocheilichthys_variegatus_wakiyae_dat03, na.action ="na.fail")

summary(Sarcocheilichthys_variegatus_wakiyae_glm03)

vif(Sarcocheilichthys_variegatus_wakiyae_glm03)
#'##### glmm
#'- combination of site and season are used as a random effect
Sarcocheilichthys_variegatus_wakiyae_glmm03 <- glmmTMB(Sarcocheilichthys_variegatus_wakiyae_PA ~ poly(depth, 2) + poly(velocity, 2) + substrate +
                                                         (1|site) + (1|season),family=binomial(link = "logit"), data= Sarcocheilichthys_variegatus_wakiyae_dat03, na.action ="na.fail")

summary(Sarcocheilichthys_variegatus_wakiyae_glmm03)

#'- Model validation
r.squaredGLMM(Sarcocheilichthys_variegatus_wakiyae_glmm03)
cor.test(predict(Sarcocheilichthys_variegatus_wakiyae_glmm03, type = "response", allow.new.levels=TRUE), Sarcocheilichthys_variegatus_wakiyae_dat03$Sarcocheilichthys_variegatus_wakiyae_PA)

#'- Sarcocheilichthys_variegatus_wakiyae depth
Sarcocheilichthys_variegatus_wakiyae_glmm03_dep_fit <- effect_plot(Sarcocheilichthys_variegatus_wakiyae_glmm03, pred = depth, interval = TRUE, plot.points = FALSE)
Sarcocheilichthys_variegatus_wakiyae_glmm03_dep_fit_PA <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["Sarcocheilichthys_variegatus_wakiyae_PA"]])
Sarcocheilichthys_variegatus_wakiyae_glmm03_dep_fit_depth <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm03_dep_fit[["plot_env"]][["p"]][["data"]][["depth"]])
Sarcocheilichthys_variegatus_wakiyae_dep_res_grp2 <- cbind(Sarcocheilichthys_variegatus_wakiyae_glmm03_dep_fit_PA, Sarcocheilichthys_variegatus_wakiyae_glmm03_dep_fit_depth)
colnames(Sarcocheilichthys_variegatus_wakiyae_dep_res_grp2) <- c("prob", "dep")

#'- normalize 0 to 1
Sarcocheilichthys_variegatus_wakiyae_dep_res_grp2 <- Sarcocheilichthys_variegatus_wakiyae_dep_res_grp2 %>%
  mutate(HSI_dep_grp2 = scale_values(prob))

#'- plot
Sarcocheilichthys_variegatus_wakiyae_depth_fig2 <- ggplot(Sarcocheilichthys_variegatus_wakiyae_dep_res_grp2, aes(x = dep, y = HSI_dep_grp2)) +
  geom_line(size = 1) +
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +
  labs(title= "Sarcocheilichthys_variegatus_wakiyae depth HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Sarcocheilichthys_variegatus_wakiyae_depth_fig2

#'- Sarcocheilichthys_variegatus_wakiyae velocity
Sarcocheilichthys_variegatus_wakiyae_glmm03_vel_fit <- effect_plot(Sarcocheilichthys_variegatus_wakiyae_glmm03, pred = velocity, interval = TRUE, plot.points = FALSE)
Sarcocheilichthys_variegatus_wakiyae_glmm03_vel_fit_PA <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["Sarcocheilichthys_variegatus_wakiyae_PA"]])
Sarcocheilichthys_variegatus_wakiyae_glmm03_vel_fit_vel <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm03_vel_fit[["plot_env"]][["p"]][["data"]][["velocity"]])
Sarcocheilichthys_variegatus_wakiyae_vel_res_grp2 <- cbind(Sarcocheilichthys_variegatus_wakiyae_glmm03_vel_fit_PA, Sarcocheilichthys_variegatus_wakiyae_glmm03_vel_fit_vel)
colnames(Sarcocheilichthys_variegatus_wakiyae_vel_res_grp2) <- c("prob", "vel")

#'- normalize 0 to 1
Sarcocheilichthys_variegatus_wakiyae_vel_res_grp2 <- Sarcocheilichthys_variegatus_wakiyae_vel_res_grp2 %>%
  mutate(HSI_vel_grp2 = scale_values(prob))

#'- plot
Sarcocheilichthys_variegatus_wakiyae_velocity_fig2 <- ggplot(Sarcocheilichthys_variegatus_wakiyae_vel_res_grp2, aes(x = vel, y = HSI_vel_grp2)) +
  geom_line(size = 1) +
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +
  labs(title= "Sarcocheilichthys_variegatus_wakiyae Velocity HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Sarcocheilichthys_variegatus_wakiyae_velocity_fig2

#'- Sarcocheilichthys_variegatus_wakiyae substrate
Sarcocheilichthys_variegatus_wakiyae_glmm03_sub_fit <- effect_plot(Sarcocheilichthys_variegatus_wakiyae_glmm03, pred = substrate, interval = FALSE, plot.points = FALSE)
Sarcocheilichthys_variegatus_wakiyae_glmm03_sub_fit_PA <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm03_sub_fit[["data"]][["Sarcocheilichthys_variegatus_wakiyae_PA"]])
Sarcocheilichthys_variegatus_wakiyae_glmm03_sub_fit_sub <- as.data.frame(Sarcocheilichthys_variegatus_wakiyae_glmm03_sub_fit[["data"]][["substrate"]])
Sarcocheilichthys_variegatus_wakiyae_sub_res_grp2 <- cbind(Sarcocheilichthys_variegatus_wakiyae_glmm03_sub_fit_PA, Sarcocheilichthys_variegatus_wakiyae_glmm03_sub_fit_sub)
colnames(Sarcocheilichthys_variegatus_wakiyae_sub_res_grp2) <- c("prob", "sub")

#'- normalize 0 to 1
Sarcocheilichthys_variegatus_wakiyae_sub_res_grp2 <- Sarcocheilichthys_variegatus_wakiyae_sub_res_grp2 %>%
  mutate(HSI_sub_grp2 = scale_values(prob))

#'- specify factor level order
# this is to re-order the category of substrate class from silt to boulder
Sarcocheilichthys_variegatus_wakiyae_sub_res_grp2$sub <- factor(Sarcocheilichthys_variegatus_wakiyae_sub_res_grp2$sub, levels = c('silt', 'sand', 'gravel',
                                                                                                                                  'pebble', 'cobble', 'boulder'))

#'- plot
Sarcocheilichthys_variegatus_wakiyae_substrate_fig2 <- ggplot(Sarcocheilichthys_variegatus_wakiyae_sub_res_grp2, aes(x = sub, y = HSI_sub_grp2)) +
  geom_bar(stat = "identity", colour = "black") + # bar plot for substrate category
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +
  labs(title= "Sarcocheilichthys_variegatus_wakiyae Substrate HSI Group2") +
  theme_bw() +
  theme(text=element_text(face="bold", size=12),
        legend.position = "none",                     # remove legend panel
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white"))

Sarcocheilichthys_variegatus_wakiyae_substrate_fig2






#'## Combined figures

#'- merge mulitple effect plot
Sarcocheilichthys_variegatus_wakiyae_dep_res_1 <- Sarcocheilichthys_variegatus_wakiyae_dep_res %>%
  mutate(group = "Global") %>%
  select(group, dep ,HSI_dep)


Sarcocheilichthys_variegatus_wakiyae_dep_res_grp2_1 <- Sarcocheilichthys_variegatus_wakiyae_dep_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, dep ,HSI_dep = HSI_dep_grp2)


Sarcocheilichthys_variegatus_wakiyae_dep_res_comb <- rbind(Sarcocheilichthys_variegatus_wakiyae_dep_res_1, 
                                                           Sarcocheilichthys_variegatus_wakiyae_dep_res_grp2_1)


#'- depth plot
Sarcocheilichthys_variegatus_wakiyae_depth_fig4 <- ggplot(Sarcocheilichthys_variegatus_wakiyae_dep_res_comb) +
  geom_line(aes(x = dep, y = HSI_dep, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Depth (m)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Sarcocheilichthys_variegatus_wakiyae Depth") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Sarcocheilichthys_variegatus_wakiyae_depth_fig4 <-Sarcocheilichthys_variegatus_wakiyae_depth_fig4+scale_color_manual(values=c("black", "red"))
Sarcocheilichthys_variegatus_wakiyae_depth_fig4




#'- merge mulitple effect plot
Sarcocheilichthys_variegatus_wakiyae_vel_res_1 <- Sarcocheilichthys_variegatus_wakiyae_vel_res %>%
  mutate(group = "Global") %>%
  select(group, vel ,HSI_vel)


Sarcocheilichthys_variegatus_wakiyae_vel_res_grp2_1 <- Sarcocheilichthys_variegatus_wakiyae_vel_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, vel ,HSI_vel = HSI_vel_grp2)



Sarcocheilichthys_variegatus_wakiyae_vel_res_comb <- rbind(Sarcocheilichthys_variegatus_wakiyae_vel_res_1,
                                                           Sarcocheilichthys_variegatus_wakiyae_vel_res_grp2_1)


#'- velocity plot
Sarcocheilichthys_variegatus_wakiyae_velocity_fig4 <- ggplot(Sarcocheilichthys_variegatus_wakiyae_vel_res_comb) +
  geom_line(aes(x = vel, y = HSI_vel, color = group), size = 1) + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Velocity (m/s)",  breaks = seq(0, 2, 0.5)) +     
  labs(title= "Sarcocheilichthys_variegatus_wakiyae Velocity") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Sarcocheilichthys_variegatus_wakiyae_velocity_fig4 <- Sarcocheilichthys_variegatus_wakiyae_velocity_fig4 + scale_color_manual(values=c("black","red"))
Sarcocheilichthys_variegatus_wakiyae_velocity_fig4



#'- merge mulitple effect plot
Sarcocheilichthys_variegatus_wakiyae_sub_res_1 <- Sarcocheilichthys_variegatus_wakiyae_sub_res %>%
  mutate(group = "Global") %>%
  select(group, sub ,HSI_sub)


Sarcocheilichthys_variegatus_wakiyae_sub_res_grp2_1 <- Sarcocheilichthys_variegatus_wakiyae_sub_res_grp2 %>%
  mutate(group = "Group2") %>%
  select(group, sub ,HSI_sub = HSI_sub_grp2)


Sarcocheilichthys_variegatus_wakiyae_sub_res_comb <- rbind(Sarcocheilichthys_variegatus_wakiyae_sub_res_1, 
                                                           Sarcocheilichthys_variegatus_wakiyae_sub_res_grp2_1)

#'- substrate plot
Sarcocheilichthys_variegatus_wakiyae_substrate_fig4 <- ggplot(Sarcocheilichthys_variegatus_wakiyae_sub_res_comb) +
  geom_bar(stat = "identity", position = "dodge", aes(x = sub, y = HSI_sub, color = group, fill = group), size = 1, colour = "black") + 
  scale_y_continuous("HSI", limits=c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("Substrate") +     
  labs(title= "Sarcocheilichthys_variegatus_wakiyae Substrate") +     
  theme_bw() +
  theme(text=element_text(face="bold", size=12),  
        legend.position = "top",  
        panel.border=element_rect(colour='black'),
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)))  # remove grid
theme(panel.background = element_rect(fill = "white")) 
Sarcocheilichthys_variegatus_wakiyae_substrate_fig4 <- Sarcocheilichthys_variegatus_wakiyae_substrate_fig4 + scale_fill_manual(values=c("black","red"))
Sarcocheilichthys_variegatus_wakiyae_substrate_fig4


#'## combine figures

ggarrange(Sarcocheilichthys_variegatus_wakiyae_depth_fig4, Sarcocheilichthys_variegatus_wakiyae_velocity_fig4, Sarcocheilichthys_variegatus_wakiyae_substrate_fig4,
          ncol = 3, nrow = 1, align = "hv") 




