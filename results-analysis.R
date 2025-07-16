#Load libraries
library(conflicted)
library(tidyverse)
library(janitor)
library(rio)
library(lubridate)
library(sf)
library(cartography)
library(tidyplots)
library(ggalluvial)

#1- Sociodemographic data analysis

#Setup Working Directory
setwd("~/Documents/projects/Klebgen Analysis")

#Metadata path
data_file <-"/Users/macbookpro/Library/CloudStorage/OneDrive-Personnel/LNSP/Dossier Partage LNSP/KlebGen/MetaData/Klebsiella pneumoniae sequencing tracking Africa CDC_Cameroon.xlsx"

#Read sociodemographics  data sheet and rename columns 
socio.demo.data <- import(data_file, which = "LNSP-sociodemo") %>%
  clean_names()

#Data Summary 
skimr::skim(socio.demo.data)

#Rename columns
socio.demo.data <- socio.demo.data %>%
  rename(country = geo_loc_name_country,
         region = geo_loc_name_admin1_e_g_region,
         laboratory_name = geo_loc_name_admin1_e_g_location,
         collection_date = sample_collection_date_yyyy_mm_dd,
         patient_type = patient_in_out)

#Columns selection and data filtering 
main_sociodemo_data <- socio.demo.data %>%
  select(origine_sample_id, collection_date, date_arrivee_lnsp, region, laboratory_name, 
         class_age, gender, patient_type, origine_lab_organism) %>%
  dplyr::filter(origine_lab_organism!="AUTRE") %>%
  mutate(gender = case_match(gender,
                             "M" ~ "Male",
                             "F" ~ "Female"),
         patient_type = str_to_lower(patient_type),
         region = str_to_lower(region),
         laboratory_name = case_match(laboratory_name,
                                 "SW-06-RHL" ~ "RHL",
                                 "LT-HGOPED" ~ "HGOPED",
                                 "HRBAF" ~ "HR BAFOUSSAM",
                                 "HRBUEA" ~ "HR BUEA",
                                 "HRGAROUA"~ "HR GAROUA",
                                 .default = laboratory_name),
         region = case_match(region,
                             "extreme-nord" ~ "extreme nord",
                             "south-west"~ "sud ouest",
                             .default = region))


main_sociodemo_data %>% select(collection_date, region) %>%
  drop_na()%>%
  mutate(collection_date = as.Date(collection_date),
         collection_month = paste(year(collection_date), month(collection_date), sep = "-"),
         collection_month = case_match(collection_month,
                                       "2026-3" ~ "2025-3",
                                       "2823-8"~"2023-8",
                                       .default = collection_month)) %>%
  ggplot(aes(x=collection_month)) +
  geom_bar(aes(fill=region), color="black")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "collection month", y = "Total samples collected", fill="Region")

#Shape file path
shp_path <- "./data/shape_2022/region_cam_2022.shp"

#Read shape file
mtq <- st_read(shp_path)
plot(st_geometry(mtq))

#Total samples collected per region
df_region <- main_sociodemo_data %>%
  group_by(region) %>%
  summarise(total_samples = n())%>%
  mutate(region = str_to_title(region))

df_region$code_rs <- c("CEN", "EXT", "LIT", "NOR", "OUE", "SUW")


#Add in map
mtq <- mtq %>% 
  left_join(df_region, by = c("Code_RS" = "code_rs"))

#Plot in map
choroLayer(x = mtq, var = "total_samples",
           method = "quantile", nclass = 4,
           border = "grey40",
           legend.pos = "topright", legend.values.rnd = 1,
           legend.title.txt = "Total samples collected")

labelLayer(x = mtq, txt = "Code_RS",
           halo = TRUE, overlap = FALSE)

#Distribution of strain
main_sociodemo_data %>%
  select(origine_lab_organism) %>%
  drop_na() %>%
  tidyplot(color=origine_lab_organism) %>%
  add_pie() %>%
  adjust_legend_title("Strain") %>%
  adjust_size(width = 150, height = 150) 

#Distribution of strain per health facilities 
lab_list <- main_sociodemo_data %>% 
  select(laboratory_name) %>%
  group_by(laboratory_name) %>%
  summarise(total=n()) %>%
  arrange(desc(total))

main_sociodemo_data %>% 
  select(laboratory_name, origine_lab_organism) %>%
  mutate(laboratory_name = factor(laboratory_name, levels = lab_list$laboratory_name)) %>%
  ggplot() +
  geom_bar(aes(x=laboratory_name, fill=origine_lab_organism), color="black")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "Health facility", y = "Total samples collected", fill="Strain")

#Distribution per age group and strain
main_sociodemo_data %>% 
  select(origine_lab_organism, class_age) %>%
  mutate(class_age = factor(class_age, levels=c("<= 30 days", "1 - 12 Months", "1 - 25 years", ">25 years")))%>%
  ggplot() +
  geom_bar(aes(x=origine_lab_organism, fill=class_age), color="black", position = "fill")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "Strain", y = "Proportion of samples collected", fill="Age group")


#Distribution per NPHL arrival date 
main_sociodemo_data %>% 
  select(date_arrivee_lnsp) %>%
  drop_na()%>%
  mutate(date_arrivee_lnsp = as.Date(date_arrivee_lnsp),
         collection_month = paste(year(date_arrivee_lnsp), month(date_arrivee_lnsp), sep = "-")) %>%
  ggplot() +
  geom_bar(aes(x=collection_month), color="black", fill="grey")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "NPHL arrival month", y = "Total samples collected")

#Distribution per NPHL arrival date and region
main_sociodemo_data %>% 
  select(date_arrivee_lnsp, region) %>%
  drop_na()%>%
  mutate(date_arrivee_lnsp = as.Date(date_arrivee_lnsp),
         collection_month = paste(year(date_arrivee_lnsp), month(date_arrivee_lnsp), sep = "-")) %>%
  ggplot() +
  geom_bar(aes(x=collection_month, fill=region), color="black", position = "fill")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "NPHL arrival month", y = "Total samples collected", fill = "Region")

#Distribution per NPHL arrival date and labs
main_sociodemo_data %>% 
  select(date_arrivee_lnsp, laboratory_name) %>%
  drop_na()%>%
  mutate(date_arrivee_lnsp = as.Date(date_arrivee_lnsp),
         collection_month = paste(year(date_arrivee_lnsp), month(date_arrivee_lnsp), sep = "-")) %>%
  ggplot() +
  geom_bar(aes(x=collection_month, fill=laboratory_name), color="black", position = "fill")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "NPHL arrival month", y = "Total samples collected", fill = "Health facility")


#Alluvium distribution per NPHL arrival date, gender, age group and health facilities 
main_sociodemo_data %>% 
  select(date_arrivee_lnsp, laboratory_name, region, gender, class_age) %>%
  drop_na()%>%
  mutate(date_arrivee_lnsp = as.Date(date_arrivee_lnsp),
         collection_month = paste(year(date_arrivee_lnsp), month(date_arrivee_lnsp), sep = "-")) %>%
  ggplot(aes(axis1=collection_month,axis2=laboratory_name, axis3=class_age, axis4=gender)) +
  geom_alluvium(color="darkgrey")+
  geom_stratum(width = 0.4) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("collection month", "Health facility", "Gender", "Age group"))

#Distribution per NPHL arrival date and gender
main_sociodemo_data %>% 
  select(date_arrivee_lnsp, gender) %>%
  drop_na()%>%
  mutate(date_arrivee_lnsp = as.Date(date_arrivee_lnsp),
         collection_month = paste(year(date_arrivee_lnsp), month(date_arrivee_lnsp), sep = "-")) %>%
  ggplot() +
  geom_bar(aes(x=collection_month, fill=gender), color="black", position = "fill")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "NPHL arrival month", y = "Total samples collected", fill = "Gender")


#Distribution per NPHL arrival date and age group
main_sociodemo_data %>% 
  select(date_arrivee_lnsp, class_age) %>%
  drop_na()%>%
  mutate(date_arrivee_lnsp = as.Date(date_arrivee_lnsp),
         collection_month = paste(year(date_arrivee_lnsp), month(date_arrivee_lnsp), sep = "-")) %>%
  ggplot() +
  geom_bar(aes(x=collection_month, fill=class_age), color="black", position = "fill")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "NPHL arrival month", y = "Total samples collected", fill = "Age group")



#Import metadata file
metadata_path <- "data/Terra/KlebGen_checker - Cameroon.tsv"
meta_data <- import(metadata_path)%>%
  clean_names() %>%
  rename(sample_id=specimen_collector_sample_id,
         collection_date = sample_collection_date_yyyy_mm_dd, 
         region=geo_loc_name_admin1_e_g_region,
         site = admin2) %>%
  select(sample_id, collection_date, year:gender, region, site , specimen_type)

head(meta_data)

#2- KPSC Analysis


#3- E. coli Analysis



