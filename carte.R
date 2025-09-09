library(tidyverse)
library(sf)
library(cartography)

shp_path <- "./data/shape_2022/region_cam_2022.shp"

#Read shape file
mtq <- st_read(shp_path)
plot(st_geometry(mtq))

df_region <- data.frame(code_rs=c("ADA","CEN","EST","EXT","LIT","NOR","NOW","SUD","SUW","OUE"),
                        prevalence = c(1.65,3.85,2.75,1.1,0,4.4,1.65,1.1,1.1,0),
                        region=c("Adamawa", "Centre", "East", "Far North", "Littoral", "North", "North-west", "South", "South-west", "West"))

View(df_region)

mtq <- mtq %>% 
  left_join(df_region, by = c("Code_RS" = "code_rs"))

#Plot in map
choroLayer(x = mtq, var = "prevalence",
           method = "quantile", nclass = 4,
           border = "grey40",
           legend.pos = "topright", legend.values.rnd = 1,
           legend.title.txt = "Prevalence (%) per region")

mtq <- mtq %>% 
  mutate(labels_map = paste(paste(region, prevalence, sep = "; "), "%", sep = ""))

labelLayer(x = mtq, txt = "labels_map",
           halo = TRUE, overlap = FALSE)

north(pos = "topleft")

barscale(size = 5)
