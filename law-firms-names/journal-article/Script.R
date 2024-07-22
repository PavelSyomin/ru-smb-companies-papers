library(dplyr)
library(ggplot2)
library(glue)
library(here)
library(readr)
library(readxl)
library(sf)
library(stringr)
library(tidyr)
library(units)

sf_use_s2(FALSE)

data <- read_csv(here("../../ru-smb-companies/legal/panel.csv"))
er <- read_csv(here("journal-article", "economic-regions.csv"))
tiles <- read_csv(here("common", "russia-tiles.csv"))
ru <- st_read(here("common", "ru.geojson"))
ru_crs <- st_crs("+proj=aea +lat_0=0 +lon_0=100 +lat_1=68 +lat_2=44 +x_0=0 +y_0=0 +ellps=krass +towgs84=28,-130,-95,0,0,0,0 +units=m +no_defs")
org_forms <- c(
  "общество с ограниченной ответственностью",
  "общество с дополнительной ответственностью",
  "акционерное общество",
  "открытое акционерное общество",
  "закрытое акционерное общество",
  "публичное акционерное общество",
  "непубличное акционерное общество",
  "хозяйственное партнерство",
  "полное товарищество",
  "товарищество на вере"
)

firms <- data %>% 
  filter(
    kind == 1, # companies only
    year <= 2021, # 2016–2021
  ) %>% 
  mutate(
    name_id = row_number(),
    settlement = case_when(
      region == "Москва" ~ "Москва",
      region == "Санкт-Петербург" ~ "Санкт-Петербург",
      TRUE ~ settlement)
  ) %>% 
  select(
    name_id,
    tin,
    name = org_name, 
    region, 
    settlement, 
    lat,
    lon, 
    year
  ) %>% 
  mutate(
    name = str_to_lower(name),
    name = str_remove(name, paste0(org_forms, collapse = "|")),
    name = str_remove_all(name, fixed("\"")),
    name = str_trim(name)
  )
unique_names <- distinct(firms, name)

# Save unique names to a file for processing with YandexGPT API
write_csv(
  unique_names,
  glue("names-{strftime(Sys.time(), '%Y-%m-%d-%H-%M-%S')}.csv")
)

# Load vectors obtained from YandexGPT API
vectors <- read_csv(here("common", "names-vectors.csv"))

# Clustering
## Find optimal number of clusters
n_clusters <- 2:100
scores <- sapply(
  n_clusters, 
  function(x) kmeans(
    slice_sample(select(vectors, dim_0:dim_255), n = 1000),
    centers = x)$tot.withins
)
ggplot(tibble(x = n_clusters, y = scores), aes(x = x, y = y)) +
  geom_line()

# Cluster names and save a sample for manual analysis
set.seed(42)
fit <- kmeans(select(vectors, dim_0:dim_255), centers = 50)
vectors$cluster <- fit$cluster
vectors %>% 
  group_by(cluster) %>% 
  slice_sample(n = 20) %>% 
  select(name) %>% 
  write_csv(
    glue("names-sample-{strftime(Sys.time(), '%Y-%m-%d-%H-%M-%S')}.csv")
  )

# Load the results of manual analysis and join them with names
cluster_labels <- read_excel(here("journal-article", "cluster-labels.xlsx"))
firms_lifetime <- count(firms, tin, name = "lifetime")
clustered <- vectors %>% 
  left_join(cluster_labels) %>% 
  right_join(firms) %>% 
  right_join(firms_lifetime) %>% 
  filter(name != "", lifetime >= 3) %>% 
  select(
    name_id, tin, name, region, settlement,
    lat, lon, year, 
    cluster, keywords, tag, group, strategy, western) %>% 
  drop_na(name, tin, region, settlement, lat, lon, year) %>% 
  distinct(name, tin, .keep_all = TRUE)

# Distance between regions
region_vectors <- vectors %>% 
  right_join(firms) %>% 
  right_join(firms_lifetime) %>% 
  filter(name != "", lifetime >= 3) %>% 
  distinct(name, tin, .keep_all = TRUE) %>% 
  drop_na(region) %>% 
  group_by(region) %>% 
  summarise(across(dim_0:dim_255, mean), cnt = n()) %>% 
  filter(cnt > quantile(cnt, .05)) %>% 
  select(-cnt)

neighboring_regions <- as_tibble(cbind(
  iso = ru$shapeISO,
  as.data.frame(st_intersects(ru, ru, sparse = FALSE, remove_self = TRUE))
))
colnames(neighboring_regions) <- c("iso", ru$shapeISO)
neighboring_regions <- pivot_longer(
  neighboring_regions,
  -iso,
  names_to = "iso_2",
  values_to = "is_neighbor"
)

centroids <- st_centroid(ru) %>% select(iso = shapeISO)
distances <- st_distance(centroids)
units(distances) <- "km"
colnames(distances) <- centroids$iso
geo_distances <- pivot_longer(
  cbind(iso = centroids$iso, as.data.frame(distances)),
  cols = -iso, 
  names_to = "iso_2", 
  values_to = "geo_distance"
)

region_names_distances <- as_tibble(cbind(
  select(region_vectors, region),
  as.matrix(dist(select(region_vectors, -region), diag = FALSE))
))
region_names_distances[upper.tri(region_names_distances, diag = FALSE)] <- NA 
colnames(region_names_distances) <- c("region", pull(select(region_names_distances, region)))
region_names_distances <- region_names_distances %>% 
  pivot_longer(-region, names_to = "region_2", values_to = "namedist") %>% 
  left_join(er) %>% 
  left_join(select(
    er, 
    region_2 = region, 
    economic_region_2 = economic_region,
    iso_code_2 = iso_code
  )) %>% 
  drop_na(namedist) %>% 
  mutate(within = economic_region == economic_region_2) %>% 
  left_join(neighboring_regions, by = c("iso_code" = "iso", "iso_code_2" = "iso_2")) %>% 
  left_join(geo_distances, by = c("iso_code" = "iso", "iso_code_2" = "iso_2"))

region_distances_stat <- region_distances %>% 
  group_by(economic_region, within) %>% 
  summarise(dist = median(dist))

economic_regions_distances <- region_distances %>% 
  ggplot(aes(x = within, y = dist)) +
  geom_violin(draw_quantiles = c(.5)) +
  geom_jitter(size = 1, shape = 21, alpha = .3) +
  geom_text(
    aes(
      y = stage(dist, after_stat = 0),
      label = after_stat(paste("M == ", round(middle, 2)))
    ),
    stat = "boxplot",
    vjust = -.25,
    parse = TRUE,
    size = 3,
    family = "Segoe UI Semilight",
  ) +
  scale_x_discrete(
    name = "Distance type",
    labels = c("Inside–outside", "Within")
  ) +
  facet_wrap(vars(economic_region), ncol = 3) +
  labs(
    y = "Distance between regions"
  ) +
  theme_bw(base_family = "Segoe UI Semilight", base_size = 9)
economic_regions_distances

arrange(region_distances, -dist)
region_distances %>% 
  ggplot(aes(x = dist)) +
  geom_freqpoly()

neighbor_plot <- region_names_distances %>% 
  ggplot(aes(x = namedist, y = is_neighbor)) +
  geom_boxplot() +
  geom_text(
    aes(
      label = after_stat(paste("M == ", round(xmiddle, 2))),
      x = stage(namedist, after_stat = xmiddle)),
    stat = "boxplot",
    vjust = -.75,
    parse = TRUE,
    size = 3,
    family = "Segoe UI Semilight",
    angle = 90
  ) + 
  scale_y_discrete(labels = c("Other\nregions", "Neighboring\nregions")) + 
  labs(
    x = "Euclidean distance between region vectors",
    y = ""
  ) +
  theme_bw(base_family = "Segoe UI Semilight", base_size = 9)
neighbor_plot

t.test(namedist ~ is_neighbor, data = region_names_distances)
wilcox.test(namedist ~ is_neighbor, data = region_names_distances)
region_names_distances %>% group_by(is_neighbor) %>% summarise(m = median(namedist))

distance_plot <- region_names_distances %>% 
  ggplot(aes(x = geo_distance, y = namedist)) +
  geom_point(color = "grey50", size = 1, shape = 21) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate(
    "label", 
    x = as_units(Inf, "km"),
    y = Inf, 
    label = paste(
      "r[Pearson] == ",
      round(cor.test(~namedist+geo_distance, data = region_names_distances)$estimate, 2)
    ),
    hjust = 1,
    vjust = 1,
    parse = TRUE,
    size = 3,
    family = "Segoe UI Semilight",
    label.size = 0,
  ) + 
  labs(
    x = "Geographical distance between regions",
    y = "Distance between region vectors"
  ) +
  theme_bw(base_family = "Segoe UI Semilight", base_size = 9)
distance_plot

cor.test(~namedist+geo_distance, data = filter(region_names_distances))$estimate

# Look at distribution by various classes
count(clustered, tag, sort = TRUE)
count(clustered, group, sort = TRUE)
count(clustered, strategy, sort = TRUE)
count(clustered, western, sort = TRUE)

# Tile grid map of regions by naming group
count(clustered, region, group) %>% 
  filter(group != "Misc") %>% 
  group_by(region) %>% 
  summarise(group = group, share = n / sum(n)) %>% 
  right_join(tiles, by = c("region" = "name")) %>% 
  drop_na(share) %>% 
  ggplot(aes(x = col, y = -row)) +
  geom_raster(aes(fill = share)) +
  geom_text(aes(label = code_en), size = 3) +
  coord_fixed() +
  facet_wrap(vars(group), ncol = 2)

# Field-based vs service-based names
naming_strategies <- clustered %>% 
  count(region, group) %>% 
  filter(group != "Misc") %>% 
  group_by(region) %>% 
  arrange(-n) %>% 
  summarise(
    group = group, 
    share = round(100 * n / sum(n)),
    main_group = first(group)
  ) %>% 
  pivot_wider(
    id_cols = c("region", "main_group"), 
    names_from = group,
    values_from = share
  ) %>% 
  right_join(tiles, by = c("region" = "name")) %>% 
  right_join(er) %>% 
  drop_na(main_group) %>% 
  ggplot(aes(x = col, y = -row)) +
  geom_tile(aes(fill = main_group, width = 1, height = 1)) +
  geom_text(
    aes(x = col - .45, y = -row + .45, label = code_en),
    size = 3, hjust = 0, vjust = 1, family = "Segoe UI Semilight"
  ) +
  geom_text(
    aes(x = col - .25, y = -row, label = Law),
    size = 3, family = "Segoe UI Semilight"
  ) +
  geom_text(
    aes(x = col - .25, y = -row - .35, label = Service),
    size = 3, family = "Segoe UI Semilight"
  ) +
  geom_text(
    aes(x = col + .25, y = -row - .25, label = Form),
    size = 3, family = "Segoe UI Semilight"
  ) +
  scale_color_brewer(palette = "Set3") +
  coord_fixed() +
  theme_void(base_family = "Segoe UI Semilight", base_size = 9)
naming_strategies

naming_strategies <- clustered %>% 
  count(region, group) %>% 
  filter(group != "Misc") %>% 
  group_by(region) %>% 
  arrange(-n) %>% 
  summarise(
    group = group, 
    share = round(100 * n / sum(n)),
    main_group = first(group)
  ) %>% 
  pivot_wider(
    id_cols = c("region", "main_group"), 
    names_from = group,
    values_from = share
  ) %>% 
  mutate(ratio = (Law - Service) / (Law + Service)) %>% 
  right_join(tiles, by = c("region" = "name")) %>%
  drop_na(main_group) %>% 
  ggplot(aes(x = col, y = -row)) +
  geom_tile(aes(fill = main_group, width = 1, height = 1)) +
  geom_text(
    aes(x = col - .45, y = -row + .45, label = code_en),
    size = 3, hjust = 0, vjust = 1, family = "Segoe UI Semilight"
  ) +
  coord_fixed() +
  theme_void(base_family = "Segoe UI Semilight", base_size = 9) +
  theme(
    legend.position = c(1, 0),
    legend.direction = "horizontal",
    legend.justification = c(1, 0),
    legend.title.position = "top"
  )
naming_strategies

clustered %>% 
  filter(group != "Misc") %>%
  count(region, group, sort = TRUE) %>% 
  group_by(region) %>% 
  slice_max(n) %>% 
  right_join(er) %>% 
  ungroup() %>% 
  count(economic_region, group)
