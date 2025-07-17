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
setwd("~/Documents/projects/Klebgen_Analysis")

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

gender_tab <- table(main_sociodemo_data$gender)
gender_tab[1]/gender_tab[2]

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

#Distribution of species per health facilities 
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
  labs(x = "Health facility", y = "Total samples collected", fill="Species")

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

#2- TheiaProk results 

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

theiaprok_results_path <- "data/Terra/CMR_All_Terra_results.tsv" 

theiaprok_results <- import(theiaprok_results_path) %>%
  clean_names() %>% 
  rename(sample_id = entity_cmr_all_id, 
         microbiology_results = organism_detected_phenotypically,
         gambit_results = gambit_predicted_taxon,
         plasmids=plasmidfinder_plasmids,
         num_reads_clean = fastq_scan_num_reads_clean_pairs,
         num_reads_raw=fastq_scan_num_reads_raw_pairs,
         ) %>%
  mutate(microbiology_results = case_match(microbiology_results, "E.coli" ~"Escherichia coli", 
                                           "K.pneumoniae" ~ "Klebsiella Pneumoniae",
                                           .default = microbiology_results))

theiaprok_results %>% 
  select(microbiology_results, gambit_results) %>%
  dplyr::filter(gambit_results!="") %>% 
  ggplot(aes(fill=microbiology_results, y=gambit_results)) +
  geom_bar(position = "fill") +
  labs(fill="Microbiology taxon", y="Gambit predicted taxon")

nb_kleb <- theiaprok_results %>% 
  select(microbiology_results, gambit_results) %>%
  dplyr::filter(gambit_results %in% c("Klebsiella quasipneumoniae", "Klebsiella pneumoniae",
                                      "Klebsiella", "Klebsiella variicola") &
                microbiology_results == "Klebsiella Pneumoniae")%>%
  count()

nb_ecoli <- theiaprok_results %>% 
  select(microbiology_results, gambit_results) %>%
  dplyr::filter(gambit_results %in% c("Escherichia coli/Shigella", "Escherichia coli",
                                      "Shigella sonnei", "Shigella flexneri", "Shigella boydii") &
                  microbiology_results=="Escherichia coli")%>%
  count()

nb_aerogenes <- theiaprok_results %>% 
  select(microbiology_results, gambit_results) %>%
  dplyr::filter(gambit_results == "Klebsiella aerogenes" &
                  microbiology_results=="K.aerogenes")%>%
  count()

#Percentage of match results
100*(nb_aerogenes + nb_ecoli + nb_kleb)/count(theiaprok_results)

#Plot N50 Klebsiella
theiaprok_results %>% 
  select(n50_value, gambit_results) %>%
  dplyr::filter(gambit_results %in% c("Klebsiella quasipneumoniae", "Klebsiella pneumoniae",
                                      "Klebsiella", "Klebsiella variicola")) %>%
  ggplot(aes(x=n50_value)) +
  geom_histogram(fill="grey", color = "black") +
  geom_vline(xintercept = mean(theiaprok_results$n50_value, na.rm = TRUE),
             linetype = "dashed", color = "red") +
  geom_vline(xintercept = median(theiaprok_results$n50_value, na.rm = TRUE),
             linetype = "dashed", color = "blue") +
  labs(x="N50")+
  theme_minimal()

#Plot N50 Ecoli
theiaprok_results %>% 
    select(n50_value, gambit_results) %>%
    dplyr::filter(gambit_results %in% c("Escherichia coli/Shigella", "Escherichia coli",
                                        "Shigella sonnei", "Shigella flexneri", "Shigella boydii")) %>%
    ggplot(aes(x=n50_value)) +
    geom_histogram(fill="grey", color = "black") +
    geom_vline(xintercept = mean(theiaprok_results$n50_value, na.rm = TRUE),
               linetype = "dashed", color = "red") +
    geom_vline(xintercept = median(theiaprok_results$n50_value, na.rm = TRUE),
               linetype = "dashed", color = "blue") +
  labs(x="N50")+
    theme_minimal()

#Plot reads number Klebsiella
theiaprok_results %>% 
  select(num_reads_clean, gambit_results) %>%
  dplyr::filter(gambit_results %in% c("Klebsiella quasipneumoniae", "Klebsiella pneumoniae",
                                      "Klebsiella", "Klebsiella variicola")) %>%
  ggplot(aes(x=num_reads_clean))+
  geom_histogram(fill="grey", color = "black") +
  geom_vline(xintercept = mean(theiaprok_results$num_reads_clean, na.rm = TRUE),
             linetype = "dashed", color = "red") +
  geom_vline(xintercept = median(theiaprok_results$num_reads_clean, na.rm = TRUE),
             linetype = "dashed", color = "blue") +
  labs(x="Number of Reads")+
  theme_minimal()

#Plot reads number Ecoli
theiaprok_results %>% 
  select(num_reads_clean, gambit_results) %>%
  dplyr::filter(gambit_results %in% c("Escherichia coli/Shigella", "Escherichia coli",
                                      "Shigella sonnei", "Shigella flexneri", "Shigella boydii")) %>%
  ggplot(aes(x=num_reads_clean))+
  geom_histogram(fill="grey", color = "black") +
  geom_vline(xintercept = mean(theiaprok_results$num_reads_clean, na.rm = TRUE),
             linetype = "dashed", color = "red") +
  geom_vline(xintercept = median(theiaprok_results$num_reads_clean, na.rm = TRUE),
             linetype = "dashed", color = "blue") +
  labs(x="Number of Reads")+
  theme_minimal()

#3- KPSC Analysis

#Import kleborate results for Klebsiella pneumoniae
kleborate_kpsc <- import("data/kleborate/KPSC/kleborate_results/klebsiella_pneumo_complex_output.txt") %>% 
  clean_names() %>%
  mutate(strain = sub("_.*", "", strain))
  

kleborate_kpsc <- kleborate_kpsc %>% select(strain, species, strain, st, yb_st, yersiniabactin, colibactin, aerobactin, 
                                            salmochelin, rmp_adc,virulence_score, k_type, o_type,a_gly_acquired:resistance_score)

glimpse(kleborate_kpsc)

merged_kpsc_data <- left_join(kleborate_kpsc, meta_data, by = c("strain" = "sample_id"))

glimpse(merged_kpsc_data)

#Distribution des ST
merged_kpsc_data %>% 
  count(st, sort = TRUE) %>%
  mutate(st = factor(st, levels = st)) %>% 
  ggplot(aes(x = st, y = n)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(y = "Total isolates", x = "ST")

#Distribution des ST per region
prevalent_st <- merged_kpsc_data %>% 
  count(st, sort = TRUE) %>%
  mutate(st = factor(st, levels = st))


merged_kpsc_data %>% 
  count(st, region ,sort = TRUE) %>%
  mutate(st = factor(st, levels = prevalent_st$st)) %>% 
  ggplot(aes(x = st, y = n, fill=region)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(y = "Total isolates", x = "ST", fill= "Region")


#Distribution des ST per site
merged_kpsc_data %>% 
  count(st, site ,sort = TRUE) %>%
  mutate(st = factor(st, levels = prevalent_st$st)) %>% 
  ggplot(aes(x = st, y = n, fill=site)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(y = "Total isolates", x = "ST", fill= "Site")

#Distribution des O Type par region
prevalent_o <- merged_kpsc_data %>% 
  count(o_type, sort = TRUE)

merged_kpsc_data %>% 
  mutate(o_type = factor(o_type, levels = prevalent_o$o_type)) %>% 
  ggplot(aes(y=o_type, , fill=region))+
  geom_bar(color="black")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(y="O Type", x="Total isolates")

#Distribution des K Type par region
prevalent_k <- merged_kpsc_data %>% 
  count(k_type, sort = TRUE)

merged_kpsc_data %>% 
  mutate(k_type = factor(k_type, levels = prevalent_k$k_type)) %>% 
  ggplot(aes(x=k_type, fill=region))+
  geom_bar(color="black")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x="K Type", y="Total isolates")

#K-type associated with most prevalent ST
prop_data <- merged_kpsc_data %>%
  dplyr::filter(st %in% prevalent_st[1:4, ]$st) %>%
  count(st, k_type) %>%
  group_by(st) %>%
  mutate(prop = n / sum(n))

ggplot(prop_data, aes(x = st, y = prop, fill = k_type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = scales::percent(prop, accuracy = 1),
                group = k_type),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(y = "Proportion of isolates", x = "ST", fill = "K Type")

#O-type associated with most provalent ST
prop_data2 <- merged_kpsc_data %>%
  dplyr::filter(st %in% prevalent_st[1:4, ]$st) %>%
  count(st, o_type) %>%
  group_by(st) %>%
  mutate(prop = n / sum(n))

ggplot(prop_data2, aes(x = st, y = prop, fill = o_type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = scales::percent(prop, accuracy = 1),
                group = o_type),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_fill_brewer(palette = "Set1") +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(y = "Proportion of isolates", x = "ST", fill = "O Type")

#Resistance score per ST
merged_kpsc_data %>% 
  mutate(resistance_score=as.character(resistance_score),
         st = factor(st, levels = prevalent_st$st))%>%
  ggplot(aes(y = resistance_score, x=st)) +
  geom_bin2d() + 
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "ST", y = "Resistance score", fill="No. Isolates")

#Virulence score per ST
merged_kpsc_data %>% 
  mutate(virulence_score=as.character(virulence_score),
         st = factor(st, levels = prevalent_st$st))%>%
  ggplot(aes(y = virulence_score, x=st)) +
  geom_bin2d() + 
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "ST", y = "Virulence score", fill="No. Isolates")

gene_matrix <- merged_kpsc_data %>%
  select(strain, st, region,a_gly_acquired:bla_carb_acquired) %>%           
  pivot_longer(cols = a_gly_acquired:bla_carb_acquired,              
               names_to = "amr_class", values_to = "genes") %>%
  dplyr::filter(!is.na(genes)) %>%                           
  separate_rows(genes, sep = ";") %>%                               
  mutate(
    gene = str_remove_all(genes, "\\^|\\*|\\.v\\d+"),                 
    present = 1
  ) %>%
  distinct(strain, st, gene, .keep_all = TRUE) %>%
  pivot_wider(names_from = gene, values_from = present, values_fill = 0) %>%
  dplyr::filter(genes != '-') %>%
  select(strain:amr_class, aadA5:dfrA25)

glimpse(gene_matrix)

#Distribution of acquired AMR per ST
gene_matrix %>% 
  mutate(amr_class = sub("_.*", "", amr_class),
         amr_class = case_match(amr_class,
                                "a" ~ "a_gly",
                                .default = amr_class)) %>%
  ggplot(aes(x=st, y=amr_class))+
  geom_bin2d()+
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "ST", y = "Acquired Drug Class Resistance", fill="No. Isolates")

annotation_data <- gene_matrix %>%
  select(strain, st, region) %>%
  unique() %>%
  column_to_rownames(var = "strain") %>% 
  arrange(region) %>%
  select(region, st)%>%
  dplyr::filter(st %in% prevalent_st$st[1:30])

amr_genes_matrix <- gene_matrix %>%
  select(strain, aadA5:dfrA25) %>%
  group_by(strain) %>%
  summarize(across(everything(), ~ as.integer(any(. == 1))), .groups = "drop") %>%
  mutate(strain = make.unique(strain)) %>%
  column_to_rownames(var="strain") 

amr_genes_matrix <- amr_genes_matrix[rownames(annotation_data), ]

pheatmap::pheatmap(
  amr_genes_matrix,
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = annotation_data,
  display_numbers = F,
  show_rownames = FALSE,
  color = colorRampPalette(c("white", "red"))(100),
)

plasmids_matrix <- theiaprok_results %>% 
  select(sample_id, plasmids) %>%
  dplyr::filter(plasmids!="") %>%
  separate_rows(plasmids, sep=",") %>%
  mutate(present = 1) %>%
  drop_na() %>%
  distinct(sample_id, plasmids, .keep_all = TRUE) %>%
  pivot_wider(names_from = plasmids, values_from = present, values_fill = 0)%>%
  column_to_rownames(var="sample_id") 

kleb_plasmids_matrix <- plasmids_matrix[rownames(annotation_data), ]

pheatmap::pheatmap(
  kleb_plasmids_matrix,
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = annotation_data,
  display_numbers = F,
  show_rownames = FALSE,
  color = colorRampPalette(c("white", "red"))(100),
)

merged_kpsc_data %>%
  select(strain, year, month, day, site, st, region, k_type, o_type) %>% 
  mutate(Country = "Cameroon") %>%
  rename(id=strain, Year=year, Month=month, Day=day, Site=site, ST=st, Region= region, O_type=o_type, K_type=k_type) %>%
  select(id, Year, Month, Day, Country, Site, ST, Region, O_type, K_type) %>%
  export(file = "data/transmission/transmission_klebsiella.tsv", format = "tsv")

#4- E. coli Analysis

#Import kleborate results for Ecoli/Shigella samples
kleborate_ecoli <- import("data/kleborate/Ecoli/kleborate_results/escherichia_output.txt") %>%
  clean_names() %>%
  mutate(strain = sub("_.*", "", strain)) %>%
  select(strain, species, st, clonal_complex, pathotype, lee_st, lee_lineage, clermont_type, clermont_profile, o_type:serotype, aminoglycoside:other_classes)

glimpse(kleborate_ecoli)

merged_ecoli_data <- left_join(kleborate_ecoli, meta_data, by = c("strain" = "sample_id"))

glimpse(merged_ecoli_data)

#Plot distribution of STs
merged_ecoli_data %>% 
  count(st, sort = T) %>%
  mutate(st = factor(st, levels = st)) %>% 
  ggplot(aes(x = st, y = n)) +
  geom_bar(stat = "identity", fill="grey", color="black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(y = "Total isolates", x = "ST")

prevalent_st <- merged_ecoli_data %>% 
  count(st, sort = T) %>%
  mutate(st = factor(st, levels = st))

merged_ecoli_data %>% 
  count(st, region ,sort = TRUE) %>%
  mutate(st = factor(st, levels = prevalent_st$st)) %>% 
  ggplot(aes(x = st, y = n, fill=region)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(y = "Total isolates", x = "ST", fill= "Region")

#Distribution des O Type par region
prevalent_o <- merged_ecoli_data %>% 
  count(o_type, sort = TRUE)

merged_ecoli_data %>% 
  mutate(o_type = factor(o_type, levels = prevalent_o$o_type)) %>% 
  ggplot(aes(x=o_type, , fill=region))+
  geom_bar(color="black")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x="O Type", y="Total isolates")

#Distribution des H Type par region
prevalent_h <- merged_ecoli_data %>% 
  count(h_type, sort = TRUE)

merged_ecoli_data %>% 
  mutate(h_type = factor(h_type, levels = prevalent_h$h_type)) %>% 
  ggplot(aes(x=h_type, fill=region))+
  geom_bar(color="black")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x="H Type", y="Total isolates")

#O-type associated with most prevalent ST
prop_data3 <- merged_ecoli_data %>%
  dplyr::filter(st %in% prevalent_st[1:5, ]$st, o_type!="-") %>%
  count(st, o_type) %>%
  group_by(st) %>%
  mutate(prop = n / sum(n)) 

ggplot(prop_data3, aes(x = st, y = prop, fill = o_type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = scales::percent(prop, accuracy = 1),
                group = o_type),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(y = "Proportion of isolates", x = "ST", fill = "O Type")

#H-type associated with most prevalent ST
prop_data4 <- merged_ecoli_data %>%
  dplyr::filter(st %in% prevalent_st[1:5, ]$st, h_type!="-") %>%
  count(st, h_type) %>%
  group_by(st) %>%
  mutate(prop = n / sum(n)) 

ggplot(prop_data4, aes(x = st, y = prop, fill = h_type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = scales::percent(prop, accuracy = 1),
                group = h_type),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(y = "Proportion of isolates", x = "ST", fill = "H Type")


gene_ecoli_matrix <- merged_ecoli_data %>%
  select(strain, st, region,aminoglycoside:betalactam) %>%           
  pivot_longer(cols = aminoglycoside:betalactam,              
               names_to = "amr_class", values_to = "genes") %>%
  dplyr::filter(!is.na(genes)) %>%                           
  separate_rows(genes, sep = ";") %>%                               
  mutate(
    gene = str_remove_all(genes, "\\^|\\*|\\.v\\d+"),                 
    present = 1
  ) %>%
  distinct(strain, st, gene, .keep_all = TRUE) %>%
  pivot_wider(names_from = gene, values_from = present, values_fill = 0) %>%
  dplyr::filter(genes != '-')

gene_ecoli_matrix %>%
  ggplot(aes(x=st, y=amr_class))+
  geom_bin2d()+
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "ST", y = "Acquired Drug Class Resistance", fill="No. Isolates")

annotation_data <- gene_ecoli_matrix %>%
  select(strain, st, region) %>%
  unique() %>%
  column_to_rownames(var = "strain") %>% 
  arrange(region) %>%
  select(region, st)%>%
  dplyr::filter(st %in% prevalent_st$st[1:20])

amr_genes_matrix <- gene_ecoli_matrix %>%
  select(strain, aadA1:cat86) %>%
  group_by(strain) %>%
  summarize(across(everything(), ~ as.integer(any(. == 1))), .groups = "drop") %>%
  mutate(strain = make.unique(strain)) %>%
  column_to_rownames(var="strain")
amr_genes_matrix <- amr_genes_matrix[rownames(annotation_data), ]

pheatmap::pheatmap(
  amr_genes_matrix,
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = annotation_data,
  display_numbers = F,
  show_rownames = FALSE,
  color = colorRampPalette(c("white", "red"))(100),
)

ecoli_plasmids_matrix <- plasmids_matrix[rownames(annotation_data), ]

pheatmap::pheatmap(
  ecoli_plasmids_matrix,
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = annotation_data,
  display_numbers = F,
  show_rownames = FALSE,
  color = colorRampPalette(c("white", "red"))(100),
)

merged_ecoli_data %>%
  select(strain, year, month, day, site, st, region, o_type, h_type) %>% 
  mutate(Country = "Cameroon") %>%
  rename(id=strain, Year=year, Month=month, Day=day, Site=site, ST=st, Region= region) %>%
  select(id, Year, Month, Day, Country, Site, ST, Region, o_type, h_type) %>%
  export(file = "data/transmission/transmission_ecoli.tsv", format = "tsv")

