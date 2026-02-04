


#' # Session information
sessionInfo()

#' # Load library
rm(list=ls(all.names=T)) # clear all data
pacman::p_load(tidyverse,
               foreign) # load package


#' # Read data

#'## site GIS data

site_env <- read.dbf("./raw_data/epsg4326_fish_site1.dbf")%>%
  mutate(slope_cv = slp_mean/slp_stdev*100) %>%
  mutate(slope= (elev_max - elev_min)/length*100) %>%
  select(site, Latitude, Longitude, Channel_ID= str_ID, stream_order= storder_ma, segment_length= length,
         site_elevation = elev_mean, site_slope = slp_mean,general_slope= slope, slope_cv,area_mean, area_min, area_max,
         urban= urban_mean, forest= fore_mean, agriculture= agri_mean, 
         ws_elev= ws_de_mean , ws_agriculture= ws_ag_mean, ws_forest= ws_fo_mean , ws_urban= ws_ur_mean, weir_num= NUMPOINTS, c(29:47),
         SBSNCD, BBSNCD, MBSNCD) %>%
  mutate(site = str_to_upper(site))%>%
  as_tibble()

site_env_add <- read.dbf("./raw_data/epsg4326_fish_site_add1.dbf")%>%
  mutate(slope_cv = slp_mean/slp_stdev*100) %>%
  mutate(slope= (elev_max - elev_min)/length*100) %>%
  select(site, Latitude = lattitude, Longitude = longitude, Channel_ID= str_ID, stream_order= storder_ma,segment_length= length,
         site_elevation = elev_mean,  site_slope = slp_mean,general_slope= slope, slope_cv,area_mean, area_min, area_max,  
         urban= urban_mean, forest= fore_mean, agriculture= agri_mean, 
         ws_elev= ws_de_mean, ws_agriculture= ws_ag_mean, ws_forest= ws_fo_mean , ws_urban= ws_ur_mean, weir_num= NUMPOINTS, c(28:46),
         SBSNCD, BBSNCD, MBSNCD ) %>%
  mutate(site = str_to_upper(site))%>%
  as_tibble()



site_env_0905 <- read.dbf("./raw_data/epsg4326_0905_site_final.dbf")%>%
  mutate(slope_cv = slp_mean/slp_stdev*100) %>%
  mutate(slope= (elev_max - elev_min)/length*100) %>%
  select(site, Latitude = lattitude, Longitude = longitude, Channel_ID= str_ID, stream_order= storder_ma,segment_length= length,
         site_elevation = elev_mean,  site_slope = slp_mean,general_slope= slope, slope_cv,area_mean, area_min, area_max,  
         urban= urban_mean, forest= fore_mean, agriculture= agri_mean, 
         ws_elev= dem_mean, ws_agriculture= agr_mean, ws_forest= for_mean , ws_urban= urb_mean, weir_num= NUMPOINTS, c(29:47),
         SBSNCD, BBSNCD, MBSNCD) %>%
  mutate(site = str_to_upper(site))%>%
  as_tibble()


site_env_0906 <- read.dbf("./raw_data/epsg4326_0906_site_final.dbf")%>%
  mutate(slope_cv = slp_mean/slp_stdev*100) %>%
  mutate(slope= (elev_max - elev_min)/length*100) %>%
  select(site, Latitude = lattitude, Longitude = longitude ,Channel_ID= str_ID, stream_order= storder_ma,segment_length= length,
         site_elevation = elev_mean,  site_slope = slp_mean,general_slope= slope, slope_cv,area_mean, area_min, area_max,  
         urban= urban_mean, forest= fore_mean, agriculture= agri_mean, 
         ws_elev= dem_mean, ws_agriculture= agri_mean, ws_forest= for_mean , ws_urban= urb_mean, weir_num= NUMPOINTS, c(27:45),
         SBSNCD, BBSNCD, MBSNCD) %>%
  mutate(site = str_to_upper(site))%>%
  as_tibble()

site_env_dg <- read.dbf("./raw_data/epsg4326_dg_site_finalshp.dbf")%>%
  mutate(slope_cv = slp_mean/slp_stdev*100) %>%
  mutate(slope= (elev_max - elev_min)/length*100) %>%
  select(site,  Latitude = lattitude, Longitude = longitude, Channel_ID= str_ID, stream_order= storder_ma, segment_length= length,
         site_elevation = elev_mean,  site_slope = slp_mean,general_slope= slope, slope_cv,area_mean, area_min, area_max,  
         urban= urban_mean, forest= fore_mean, agriculture= agri_mean, 
         ws_elev= dem_mean, ws_agriculture= agri_mean, ws_forest= for_mean , ws_urban= urb_mean, weir_num= NUMPOINTS, c(26:44),
         SBSNCD, BBSNCD, MBSNCD) %>%
  mutate(site = str_to_upper(site))%>%
  as_tibble()


site_env_gj<- read.dbf("./raw_data/epsg4326_site_gapyeongjojong.dbf")%>%
  mutate(slope_cv = slp_mean/slp_stdev*100) %>%
  mutate(slope= (elev_max - elev_min)/length*100) %>%
  select(site, Latitude = lattitude, Longitude = longitude, Channel_ID= str_ID, stream_order= storder_ma, segment_length= length,
         site_elevation = elev_mean,  site_slope = slp_mean,general_slope= slope, slope_cv,area_mean, area_min, area_max,  
         urban= urban_mean, forest= fore_mean, agriculture= agri_mean, 
         ws_elev= dem_mean, ws_agriculture= agri_mean, ws_forest= for_mean , ws_urban= urb_mean, weir_num= NUMPOINTS, c(30:49),
         SBSNCD, BBSNCD, MBSNCD) %>%
  mutate(site = str_to_upper(site))%>%
  as_tibble()

#'merge dbf files

site_env1 <- site_env%>%
  rbind(site_env_add)%>%
  rbind(site_env_0905) %>%
  rbind(site_env_0906)%>%
  rbind(site_env_dg)%>%
  rbind(site_env_gj)


#' check data property
length(unique(site_env1$site) ) # n = 205
list(unique(site_env1$site))

#'# Export final data into data_outcome folder
filename <- paste0("./data_outcome/HSI_site_env_dat.csv")
write.csv(site_env1, filename,row.names=FALSE) 

