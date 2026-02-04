

#' # Session information
sessionInfo()

#' # Load library
rm(list=ls(all.names=T)) # clear all data
pacman::p_load(tidyverse) # load package
#' # Change language
Sys.getlocale() # check encoding
Sys.setlocale("LC_ALL","C") # remove language
Sys.setlocale("LC_ALL","Korean") # set as Korean


#'# Read data 
fish_data_2007_2010 <- read.csv("./Raw_data/fish_data_2007_2010.csv", header=TRUE)
fish_data_2007_2010

fish_data_hongcheon_2012 <- read.csv("./Raw_data/fish_data_hongcheon_2012.csv", header=TRUE, strip.white=TRUE)
fish_data_hongcheon_2012

fish_data_2013 <- read.csv("./Raw_data/fish_data_2013.csv", header=TRUE, strip.white=TRUE)
fish_data_2013

fish_data_2017 <- read.csv("./Raw_data/fish_data_2017.csv", header=TRUE, strip.white=TRUE)
fish_data_2017

fish_data_envflow_project_2018 <- read.csv("./Raw_data/fish_data_envflow_project_2018.csv", header=TRUE, strip.white=TRUE)
fish_data_envflow_project_2018

fish_data_envflow_project_2019 <- read.csv("./Raw_data/fish_data_envflow_project_2019.csv", header=TRUE, strip.white=TRUE)
fish_data_envflow_project_2019

fish_data_envflow_project_2020 <- read.csv("./Raw_data/fish_data_envflow_project_2020.csv", header=TRUE, strip.white=TRUE)
fish_data_envflow_project_2020

fish_data_envflow_project_2021 <- read.csv("./Raw_data/fish_data_envflow_project_2021.csv", header=TRUE, strip.white=TRUE)
fish_data_envflow_project_2021

fish_data_envflow_project_2022 <- read.csv("./Raw_data/fish_data_envflow_project_2022.csv", header=TRUE, strip.white=TRUE)
fish_data_envflow_project_2022

fish_data_s_y_basin_2020 <- read.csv("./Raw_data/fish_data_s_y_basin_2020.csv", header=TRUE, strip.white=TRUE)
fish_data_s_y_basin_2020

fish_data_s_y_basin_2021 <- read.csv("./Raw_data/fish_data_s_y_basin_2021.csv", header=TRUE, strip.white=TRUE)
fish_data_s_y_basin_2021

fish_data_s_y_basin_2022 <- read.csv("./Raw_data/fish_data_s_y_basin_2022.csv", header=TRUE, strip.white=TRUE)
fish_data_s_y_basin_2022


fish_data_Dopyeong_Stream <- read.csv("./Raw_data/fish_data_Dopyeong_Stream.csv", header=TRUE, strip.white=TRUE)
fish_data_Dopyeong_Stream

fish_data_Dutayeon <- read.csv("./Raw_data/Dutayeon_Data_2023_0908.csv", header=TRUE, strip.white=TRUE)
fish_data_Dutayeon

fish_data_mancheon <- read.csv("./Raw_data/fish_data_mancheon.csv", header=TRUE, strip.white=TRUE)
fish_data_mancheon



donggang <- read.csv("./Raw_data/donggang.csv", header=TRUE,strip.white=TRUE)
donggang

gotangyo <-  read.csv("./Raw_data/gotangyo.csv", header=TRUE)
gotangyo

inje <- read.csv("./Raw_data/inje.csv", header=TRUE)
inje

naeseongcheon_jungryu <- read.csv("./Raw_data/naeseongcheon_jungryu.csv", header=TRUE)
naeseongcheon_jungryu

yongdong_jeosuji <- read.csv("./Raw_data/yongdong_jeosuji.csv", header=TRUE)
yongdong_jeosuji

youngju_seocheon <- read.csv("./Raw_data/youngju_seocheon.csv", header=TRUE)%>%
  drop_na()
youngju_seocheon


daljeoncheon <- read.csv("./Raw_data/daljeoncheon.csv", header= TRUE)
daljeoncheon

geumdangcheon <- read.csv("./Raw_data/geumdangcheon.csv", header= TRUE)
geumdangcheon

yulmuncheon <- read.csv("./Raw_data/yulmuncheon.csv", header= TRUE)
yulmuncheon


palmicheon_haryu <- read.csv("./Raw_data/palmicheon_haryu.csv", header= TRUE)
palmicheon_haryu

bokhacheon <- read.csv("./Raw_data/bokhacheon.csv", header= TRUE)
bokhacheon

iljuk_cheongmicheon <- read.csv("./Raw_data/iljuk_cheongmicheon.csv", header= TRUE)
iljuk_cheongmicheon

sobaek_geumgye <- read.csv("./Raw_data/sobaek_geumgye.csv", header= TRUE)
sobaek_geumgye

whitebeerplant <- read.csv("./Raw_data/whitebeerplant.csv", header= TRUE)
whitebeerplant

munmakgyo <- read.csv("./Raw_data/munmakgyo.csv", header= TRUE)
munmakgyo

gapyeongcheon <- read.csv("./Raw_data/gapyeong_stream.csv", header= TRUE)
gapyeongcheon

jojongcheon <- read.csv("./Raw_data/jojong_stream.csv", header= TRUE)
jojongcheon



#'# Data management

#'## Change data format
fish_data_2007_2010_2 <- fish_data_2007_2010 %>%
  mutate(vegetation = "NA") %>%
  mutate(habitat_type = ifelse(habitat_type == "", NA, habitat_type)) %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type)


fish_data_hongcheon_2012_2 <- fish_data_hongcheon_2012 %>%
  mutate(vegetation = "NA") %>%
  mutate(substrate = "NA") %>%
  mutate(habitat_type = "NA") %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type)

fish_data_2013_2 <- fish_data_2013 %>%
  mutate(vegetation = "NA") %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site = site_name, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type)

fish_data_2017_2 <- fish_data_2017 %>%
  mutate(vegetation = "NA") %>%
  mutate(habitat_type = "NA") %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type)

fish_data_envflow_project_2018_2 <- fish_data_envflow_project_2018 %>%
  mutate(vegetation = "NA") %>%
  mutate(habitat_type = "NA") %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type)

fish_data_envflow_project_2019_2 <- fish_data_envflow_project_2019 %>%
  mutate(vegetation = "NA") %>%
  mutate(habitat_type = "NA") %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type)

fish_data_envflow_project_2020_2 <- fish_data_envflow_project_2020 %>%
  mutate(vegetation = "NA") %>%
  mutate(habitat_type = "NA") %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type)

fish_data_envflow_project_2021_2 <- fish_data_envflow_project_2021 %>%
  mutate(vegetation = "NA") %>%
  mutate(habitat_type = "NA") %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site = Site_Name, date, year, month, sample_number, 
         species = Species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type)

fish_data_envflow_project_2022_2 <- fish_data_envflow_project_2022 %>%
  mutate(vegetation = "NA") %>%
  mutate(habitat_type = "NA") %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type)

fish_data_s_y_basin_2020_2 <- fish_data_s_y_basin_2020 %>%
  mutate(habitat_type = "NA") %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type)

fish_data_s_y_basin_2021_2 <- fish_data_s_y_basin_2021 %>%
  mutate(habitat_type = "NA") %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number = individaul_number, 
         depth, velocity, substrate, vegetation, habitat_type)

fish_data_s_y_basin_2022_2 <- fish_data_s_y_basin_2022 %>%
  mutate(habitat_type = "NA") %>%
  mutate(date2 = date) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Substrate_num, vegetation, habitat_type)

fish_data_Dopyeong_Stream_2 <- fish_data_Dopyeong_Stream %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(velocity = velocity/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type = Habitat_Type)

fish_data_Dutayeon_2 <- fish_data_Dutayeon %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  mutate(
    individual_number = ifelse(duplicated(species), n_distinct(species), 1),
    individual_number = ifelse(is.na(species), NA, individual_number)
  )  %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, vegetation, habitat_type = Habitat_Type)



fish_data_mancheon_2 <- fish_data_mancheon %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  unite(site, c(site, plot), sep = "_", remove = TRUE) %>%
  group_by(site, sample_number) %>% 
  mutate(
    individual_number = ifelse(duplicated(species), n_distinct(species), 1),
    individual_number = ifelse(is.na(species), NA, individual_number)
  )  %>%
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate , vegetation, habitat_type) 


donggang_2 <- donggang %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()




gotangyo_2 <- gotangyo %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()



inje_2 <- inje %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()


naeseongcheon_jungryu_2 <- naeseongcheon_jungryu %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()


yongdong_jeosuji_2 <- yongdong_jeosuji %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin = basin_eng, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()


youngju_seocheon_2 <- youngju_seocheon %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100)%>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin , site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()


geumdangcheon_2 <- geumdangcheon %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()

palmicheon_haryu_2 <- palmicheon_haryu %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()


daljeoncheon_2 <- daljeoncheon %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()


yulmuncheon_2 <- yulmuncheon %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()


munmakgyo_2 <- munmakgyo %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()


sobaek_geumgye_2 <- sobaek_geumgye %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()

iljuk_cheongmicheon_2 <- iljuk_cheongmicheon %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()


bokhacheon_2 <- bokhacheon %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()



whitebeerplant_2 <- whitebeerplant %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()



gapyeongcheon_2 <- gapyeongcheon %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()




jojongcheonn_2 <- jojongcheon %>%
  mutate(length = "NA") %>%
  mutate(date2 = date) %>%
  mutate(depth = depth/100) %>%
  separate(date2, c("year", "month", "day"), sep = "-") %>% # split year, month, day
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = as.numeric(day)) %>%
  group_by(sample_number) %>%  
  select(basin, site, date, year, month, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate = Dominant_substrate, vegetation, habitat_type = Habitat_Type) %>% as.tibble()




#'## merge fish data
fish_df <- rbind(fish_data_2007_2010_2, fish_data_2013_2, fish_data_2017_2, fish_data_envflow_project_2018_2,
                 fish_data_envflow_project_2019_2, fish_data_envflow_project_2020_2, fish_data_envflow_project_2021_2,
                 fish_data_envflow_project_2022_2, fish_data_s_y_basin_2020_2, fish_data_s_y_basin_2021_2,
                 fish_data_s_y_basin_2022_2, fish_data_Dopyeong_Stream_2,fish_data_hongcheon_2012_2,fish_data_Dutayeon_2,fish_data_mancheon_2
                 ,donggang_2, gotangyo_2, inje_2, naeseongcheon_jungryu_2, yongdong_jeosuji_2, youngju_seocheon_2, geumdangcheon_2,
                 yulmuncheon_2, daljeoncheon_2, palmicheon_haryu_2,bokhacheon_2, iljuk_cheongmicheon_2, sobaek_geumgye_2, munmakgyo_2,whitebeerplant_2
                 ,gapyeongcheon_2, jojongcheonn_2 ) %>%
  mutate(basin = case_when(basin == "Geum_River" ~ "Guem_River",
                           basin == "Nakdon_River" ~ "Nakdong_River",
                           basin == "han_river" ~ "Han_River",
                           basin == "geum_river" ~ "Guem_River",
                           basin == "nakdong_River" ~ "Nakdong_River", 
                           basin == "nakdong_river" ~ "Nakdong_River",
                           TRUE ~ as.character(basin))) %>% # fix typo
  mutate(habitat_type = case_when(habitat_type == "Riflle" ~ "Riffle",
                                  habitat_type == "Ruiffle" ~ "Riffle",
                                  habitat_type == "Riiffle" ~ "Riffle",
                                  TRUE ~ as.character(habitat_type))) %>% # fix typo
  mutate(substrate2 = case_when(substrate == "1" ~ "silt",
                                substrate == "2" ~ "sand",
                                substrate == "3" ~ "gravel",
                                substrate == "4" ~ "pebble",
                                substrate == "5" ~ "cobble",
                                substrate == "6" ~ "boulder",
                                substrate == "7" ~ "boulder",
                                TRUE ~ as.character(substrate))) %>% # change number to substrate name
  mutate(season = case_when(month == "1" ~ "NA",
                            month == "3" ~ "spring",
                            month == "4" ~ "spring",
                            month == "5" ~ "spring",
                            month == "6" ~ "summer",
                            month == "7" ~ "summer",
                            month == "8" ~ "summer",
                            month == "9" ~ "fall",
                            month == "10" ~ "fall",
                            month == "11" ~ "fall")) %>% # change month to season
  select(basin, site, date, year, month, season, sample_number, 
         species, length, individual_number, 
         depth, velocity, substrate, substrate2, vegetation, habitat_type) %>%
  mutate(site = str_to_upper(site))%>%# change site ID to capital 
  filter(basin != "")

#'- check data property 
length(unique(fish_df$basin)) # n =  6basins
length(unique(fish_df$site)) # n =  211sites 
length(unique(fish_df$species)) # n =  125species
sum(fish_df$individual_number, na.rm = TRUE) # n = 95134  individuals
range(fish_df$year) # sample from 2007 to 2024
list(unique(fish_df$basin))
list(unique(fish_df$site))
list(unique(fish_df$date))
#'# Export final data into data_outcome folder
filename <- paste0("./data_outcome/fish_sampling_DB.csv")
write.csv(fish_df, filename)  



problematic_rows <-fish_df[c(51:58, 70:89), ]
print(problematic_rows)

