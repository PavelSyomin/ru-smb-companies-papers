library(dplyr)
library(forcats)
library(ggplot2)
library(readr)
library(scales)
library(sf)
library(tidyr)

# Data
smb_data <- read_csv("assets/smb.csv")
revexp_data <- read_csv("assets/revexp.csv")
empl_data <- read_csv("assets/empl.csv")

# Maps
regions_boundaries <- st_read("assets/ru.geojson")
ru_crs <- st_crs("+proj=aea +lat_0=0 +lon_0=100 +lat_1=68 +lat_2=44 +x_0=0 +y_0=0 +ellps=krass +towgs84=28,-130,-95,0,0,0,0 +units=m +no_defs")
regions <- read_csv("assets/regions.csv")
regions <- regions_boundaries %>% 
  left_join(regions, by = c("shapeISO" = "iso_code")) %>% 
  select(name, name_en = shapeName)
ru_svr <- st_read("assets/ru_svr.geojson")

# Regional distribution
regions_map <- smb_data %>%
  filter(
    start_date <= "2021-12-31",
    end_date >= "2021-12-31") %>% 
  count(region) %>% 
  left_join(regions, by = c("region" = "name")) %>% 
  st_as_sf() %>% 
  ggplot() +
  geom_sf(aes(fill = n), size = .1) +
  coord_sf(crs = ru_crs) +
  scale_fill_distiller(
    name = ifelse(exists("RU"), "Число юридических фирм", "Count of law firms"),
    palette = "YlGn",
    direction = 1) +
  theme_void(base_size = 11, base_family = "Times New Roman") +
  theme(
    legend.title.position = "top",
    legend.position = c(.05, 0), 
    legend.justification = c(0, 0),
    legend.direction = "horizontal")
regions_map

settlements_df <- smb_data %>%
  filter(
    start_date <= "2021-12-31",
    end_date >= "2021-12-31"
  ) %>% 
  count(region, area, settlement, lat, lon) %>% 
  drop_na(lat, lon) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

settlements_map <- regions %>% ggplot() +
  geom_sf(fill = "white", color = "gray50") +
  geom_sf(data = settlements_df, aes(color = n), size = .5, shape = 21) +
  coord_sf(crs = ru_crs) +
  scale_color_distiller(
    name = ifelse(exists("RU"), "Число юридических фирм", "Count of law firms"),
    palette = "YlGn",
    transform = "log10",
    direction = 1) +
  theme_void(base_size = 11, base_family = "Times New Roman") +
  theme(
    legend.title.position = "top",
    legend.position = c(.05, 0), 
    legend.justification = c(0, 0),
    legend.direction = "horizontal")
settlements_map

settlements_df_svr <- smb_data %>%
  filter(
    region == "Свердловская область",
    start_date <= "2021-12-31",
    end_date >= "2021-12-31"
  ) %>% 
  count(region, area, settlement, lat, lon) %>% 
  drop_na(lat, lon) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) 

settlements_map_svr <- ru_svr %>% 
  ggplot() +
  geom_sf(fill = "white", color = "gray50") +
  geom_sf(data = settlements_df_svr, aes(fill = n), size = 2, shape = 21) +
  coord_sf(crs = 4326) +
  scale_size_continuous() +
  scale_fill_distiller(
    name = ifelse(exists("RU"), "Число юридических фирм", "Count of law firms"),
    transform = "log10",
    palette = "YlGn",
    direction = 1) +
  theme_void(base_size = 11, base_family = "Times New Roman")
settlements_map_svr

profit_empl_df <- smb_data %>%
  filter(
    start_date <= "2021-12-31",
    end_date >= "2021-12-31"
  ) %>% 
  left_join(filter(revexp_data, year == 2021), by = "tin") %>% 
  left_join(filter(empl_data, year == 2021), by = "tin") %>% 
  mutate(profit = revenue - expenditure) %>% 
  select(tin, region, area, settlement, lat, lon, profit, empl = employees_count) %>% 
  group_by(region, area, settlement, lat, lon) %>% 
  summarise(
    profit = sum(profit, na.rm = TRUE) ,
    empl = sum(empl, na.rm = TRUE)
  ) %>% 
  filter(profit > 0, empl > 0) %>% 
  drop_na(lat, lon) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = "EPSG:4326")

settlements_map_profit <- profit_empl_df %>% 
  ggplot() +
  geom_sf(data = regions, fill = "white", color = "gray50") +
  geom_sf(aes(color = profit), size = .5, shape = 21) +
  coord_sf(crs = ru_crs) +
  scale_color_distiller(
    name = ifelse(exists("RU"), "Прибыль юридических фирм", "Profit of law firms"),
    transform = "log10",
    palette = "YlGn",
    direction = 1) +
  theme_void(base_size = 11, base_family = "Times New Roman") +
  theme(
    legend.title.position = "top",
    legend.position = c(.05, 0), 
    legend.justification = c(0, 0),
    legend.direction = "horizontal")
settlements_map_profit

settlements_map_empl <- profit_empl_df %>% 
  ggplot() +
  geom_sf(data = regions, fill = "white", color = "gray50") +
  geom_sf(aes(color = empl), size = .5, shape = 21) +
  coord_sf(crs = ru_crs) +
  scale_color_distiller(
    name = ifelse(exists("RU"), "Число работников\nюридических фирм", "Count of employees\nat law firms"),
    transform = "log10",
    palette = "YlGn",
    direction = 1) +
  theme_void(base_size = 11, base_family = "Times New Roman") +
  theme(
    legend.title.position = "top",
    legend.position = c(.05, 0), 
    legend.justification = c(0, 0),
    legend.direction = "horizontal")
settlements_map_empl
