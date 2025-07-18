pacman::p_load(
  rio,          
  here,            
  tidyverse,      
  ape,            
  ggtree,          
  treeio,          
  ggnewscale,
  janitor,
  tidytree) 
library(scales)


#Setup Working Directory
setwd("~/Documents/projects/Klebgen_Analysis")

#Import metadata file
metadata_path <- "data/Terra/KlebGen_checker - Cameroon.tsv"
meta_data <- import(metadata_path)%>%
  clean_names() %>%
  rename(sample_id=specimen_collector_sample_id,
         collection_date = sample_collection_date_yyyy_mm_dd, 
         region=geo_loc_name_admin1_e_g_region,
         site = admin2) %>%
  select(sample_id, collection_date, year:gender, region, site , specimen_type)


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
                                           .default = microbiology_results)) %>%
  select(sample_id, plasmids)


#Import kleborate results for Klebsiella pneumoniae
kleborate_kpsc <- import("data/kleborate/KPSC/kleborate_results/klebsiella_pneumo_complex_output.txt") %>% 
  clean_names() %>%
  mutate(strain = sub("_.*", "", strain))


kleborate_kpsc <- kleborate_kpsc %>% select(strain, species, strain, st, yb_st, yersiniabactin, colibactin, aerobactin, 
                                            salmochelin, rmp_adc,virulence_score, k_type, o_type,a_gly_acquired:resistance_score)


merged_kpsc_data <- left_join(kleborate_kpsc, meta_data, by = c("strain" = "sample_id"))


#Import kleborate results for Ecoli/Shigella samples
kleborate_ecoli <- import("data/kleborate/Ecoli/kleborate_results/escherichia_output.txt") %>%
  clean_names() %>%
  mutate(strain = sub("_.*", "", strain)) %>%
  select(strain, species, st, clonal_complex, pathotype, lee_st, lee_lineage, clermont_type, clermont_profile, o_type:serotype, aminoglycoside:other_classes)


merged_ecoli_data <- left_join(kleborate_ecoli, meta_data, by = c("strain" = "sample_id"))

#1- KPSC

kleb_tree <- treeio::read.nexus("data/kleborate/KPSC/ska_output/rooted_bestTree.kleb_tree.nwk")

head(kleb_tree$tip.label)

merged_kpsc_data$strain %in% kleb_tree$tip.label

kleb_tree$tip.label %in% merged_kpsc_data$strain

ggtree(kleb_tree)

ggtree(kleb_tree,) %<+% merged_kpsc_data +
  geom_tiplab(size = 1.5) +                
  geom_text2(
    mapping = aes(subset = !isTip,
                  label = node),
    size = 2,
    color = "darkred",
    hjust = 1,
    vjust = 1) 


p <- ggtree(kleb_tree,) %<+% merged_kpsc_data +
  geom_tiplab(size = 1.5) 

p

p + geom_hilight(  
  node = 164,
  fill = "steelblue",
  extend = 0.0017) +  
  geom_hilight(           
    node = 305,
    fill = "yellow",
    extend = 0.0017) 

viewClade(p, node = 164)

viewClade(p, node = 305)

clade1_kleb_tree <- tree_subset(kleb_tree, node = 164, levels_back = 0)


clade2_kleb_tree <- tree_subset(
  kleb_tree,
  node = 305, 
  levels_back = 0)


ggtree(clade1_kleb_tree) 

ggtree(clade2_kleb_tree) 


clade1_kpsc_data <- merged_kpsc_data %>%
  dplyr::filter(merged_kpsc_data$strain %in% clade1_kleb_tree$tip.label) %>%
  tibble::column_to_rownames("strain")

clade2_kpsc_data <- merged_kpsc_data %>%
  dplyr::filter(merged_kpsc_data$strain %in% clade2_kleb_tree$tip.label) %>%
  tibble::column_to_rownames("strain")

#Clade 1 plot

p1 <- ggtree(clade1_kleb_tree) %<+% clade1_kpsc_data +
  geom_tiplab(size = 1.5) 

p1

region_df <- clade1_kpsc_data %>%
  select(region)

otype_df <- clade1_kpsc_data %>%
  select(o_type)

ktype_df <- clade1_kpsc_data %>%
  select(k_type)

virulence_df <- clade1_kpsc_data %>%
  select(virulence_score)

resistance_df <- clade1_kpsc_data %>%
  select(resistance_score)

h1 <- gheatmap(p1, region_df, 
               offset = 0.0002,
               width = 0.05, 
               color = "black", 
               colnames = FALSE) +
  scale_fill_manual(name = "Region",
                    values = c("Center" = "#00d1b1",
                               "South-West" = "grey",
                               "Littoral" = "purple",
                               "North" = "pink",
                               "West" = "red")) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.box = "vertical", 
        legend.direction = "vertical",
        legend.margin = margin())

h1

h2 <- h1 + new_scale_fill()

h3 <- gheatmap(h2, otype_df,   
               offset = 0.0004, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "O Type",
                    values = c(
                      "O2β" = "#1f77b4",         
                      "O13" = "#ff7f0e",          
                      "O3αβ" = "#2ca02c",        
                      "O1αβ,2α" = "#d62728",     
                      "O4" = "#9467bd",          
                      "O1αβ,2β" = "#8c564b",     
                      "O5" = "#e377c2",           
                      "O2α" = "#7f7f7f",         
                      "O3γ" = "#bcbd22",          
                      "O10" = "#17becf",          
                      "unknown (OL15)" = "#000000"  
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())

h3

h4 <- h3 + new_scale_fill()

k_types <- c(
  "K27", "unknown (KL183)", "K55", "K25", "K17", "unknown (KL149)", "unknown (KL108)", "K9",
  "unknown (KL102)", "K62", "Capsule null", "K2", "unknown (KL151)", "K23", "K7", "K12",
  "K3", "K54", "unknown (KL106)", "K64", "K48", "K28", "K51", "K22", "K43", "unknown (KL171)",
  "unknown (KL126)", "unknown (KL114)", "K15", "unknown (KL112)", "unknown (KL125)", "K14",
  "unknown (KL122)", "K11", "K52", "K39", "unknown (KL157)", "K37", "K24", "unknown (KL132)",
  "unknown (KL155)", "K31", "unknown (KL146)", "K81", "K38", "K47", "K53", "K16",
  "unknown (KL116)", "unknown (KL105)", "K42", "K13", "K10", "K74"
)

k_type_colors <- setNames(hue_pal()(length(k_types)), k_types)

h5 <- gheatmap(h4, ktype_df,   
               offset = 0.0006, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "k Type",
                    values = k_type_colors) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())


h5

h6 <- h5 + new_scale_fill()

h7 <- gheatmap(h6, virulence_df, 
               offset =  0.0008,  
               width = 0.05,
               color = "black", 
               colnames = FALSE)+
  scale_fill_continuous(name = "Virulence score",
                        low = "green", high = "red",
                        breaks = c(0,1.00, 3.00),
                        na.value = "white")+
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1))+
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.direction = "vertical",
        legend.box = "horizontal", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 2)))

h7

#Clade 2 plot

p1 <- ggtree(clade2_kleb_tree) %<+% clade1_kpsc_data +
  geom_tiplab(size = 1.5) 

p1

region_df <- clade2_kpsc_data %>%
  select(region)

otype_df <- clade2_kpsc_data %>%
  select(o_type)

ktype_df <- clade2_kpsc_data %>%
  select(k_type)

virulence_df <- clade2_kpsc_data %>%
  select(virulence_score)

resistance_df <- clade2_kpsc_data %>%
  select(resistance_score)

h1 <- gheatmap(p1, region_df, 
               offset = 0.0006,
               width = 0.05, 
               color = "black", 
               colnames = T) +
  scale_fill_manual(name = "Region",
                    values = c("Center" = "#00d1b1",
                               "South-West" = "grey",
                               "Littoral" = "purple",
                               "North" = "pink",
                               "West" = "red")) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.box = "vertical", 
        legend.direction = "vertical",
        legend.margin = margin())

h1

h2 <- h1 + new_scale_fill()

h3 <- gheatmap(h2, otype_df,   
               offset = 0.0013, 
               width = 0.05,
               color = "black",
               colnames = T)+
  scale_fill_manual(name = "O Type",
                    values = c(
                      "O2β" = "#1f77b4",         
                      "O13" = "#ff7f0e",          
                      "O3αβ" = "#2ca02c",        
                      "O1αβ,2α" = "#d62728",     
                      "O4" = "#9467bd",          
                      "O1αβ,2β" = "#8c564b",     
                      "O5" = "#e377c2",           
                      "O2α" = "#7f7f7f",         
                      "O3γ" = "#bcbd22",          
                      "O10" = "#17becf",          
                      "unknown (OL15)" = "#000000"  
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())

h3

h4 <- h3 + new_scale_fill()

k_types <- c(
  "K27", "unknown (KL183)", "K55", "K25", "K17", "unknown (KL149)", "unknown (KL108)", "K9",
  "unknown (KL102)", "K62", "Capsule null", "K2", "unknown (KL151)", "K23", "K7", "K12",
  "K3", "K54", "unknown (KL106)", "K64", "K48", "K28", "K51", "K22", "K43", "unknown (KL171)",
  "unknown (KL126)", "unknown (KL114)", "K15", "unknown (KL112)", "unknown (KL125)", "K14",
  "unknown (KL122)", "K11", "K52", "K39", "unknown (KL157)", "K37", "K24", "unknown (KL132)",
  "unknown (KL155)", "K31", "unknown (KL146)", "K81", "K38", "K47", "K53", "K16",
  "unknown (KL116)", "unknown (KL105)", "K42", "K13", "K10", "K74"
)

k_type_colors <- setNames(hue_pal()(length(k_types)), k_types)

h5 <- gheatmap(h4, ktype_df,   
               offset = 0.002, 
               width = 0.05,
               color = "black",
               colnames = T)+
  scale_fill_manual(name = "k Type",
                    values = k_type_colors) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())


h5

h6 <- h5 + new_scale_fill()

h7 <- gheatmap(h6, resistance_df, 
               offset =  0.003,  
               width = 0.05,
               color = "black", 
               colnames = T,
               colnames_angle = 90)+
  scale_fill_continuous(name = "Resistance score",
                        low = "green", high = "red",
                        breaks = c(0,1.00, 3.00),
                        na.value = "white")+
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1))+
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.direction = "vertical",
        legend.box = "horizontal", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 2)))

h7

#2- Ecoli

ecoli_tree <- treeio::read.nexus("data/kleborate/Ecoli/ska_output/rooted_bestTree.ecoli_tree.nwk")


merged_ecoli_data$strain %in% ecoli_tree$tip.label

ecoli_tree$tip.label %in% merged_ecoli_data$strain

ggtree(ecoli_tree,) %<+% merged_ecoli_data +
  geom_tiplab(size = 1.5) +                
  geom_text2(
    mapping = aes(subset = !isTip,
                  label = node),
    size = 2,
    color = "darkred",
    hjust = 1,
    vjust = 1) 

p <- ggtree(ecoli_tree,) %<+% merged_ecoli_data +
  geom_tiplab(size = 1.5) 

p

p + geom_hilight(  
  node = 131,
  fill = "steelblue",
  extend = 0.0017) +  
  geom_hilight(           
    node = 188,
    fill = "yellow",
    extend = 0.0017) +  
  geom_hilight(           
    node = 217,
    fill = "#d62728",
    extend = 0.0017)+
  geom_hilight(           
    node = 237,
    fill = "#9467bd",
    extend = 0.0017)


clade1_ecoli_tree <- tree_subset(ecoli_tree, node = 131, levels_back = 0)


clade2_ecoli_tree <- tree_subset(
  ecoli_tree,
  node = 188, 
  levels_back = 0)


clade3_ecoli_tree <- tree_subset(
  ecoli_tree,
  node = 217, 
  levels_back = 0)


clade4_ecoli_tree <- tree_subset(
  ecoli_tree,
  node = 237, 
  levels_back = 0)

ggtree(clade4_ecoli_tree)


clade1_ecoli_data <- merged_ecoli_data %>%
  dplyr::filter(merged_ecoli_data$strain %in% clade1_ecoli_tree$tip.label) %>%
  tibble::column_to_rownames("strain")

clade2_ecoli_data <- merged_ecoli_data %>%
  dplyr::filter(merged_ecoli_data$strain %in% clade2_ecoli_tree$tip.label) %>%
  tibble::column_to_rownames("strain")

clade3_ecoli_data <- merged_ecoli_data %>%
  dplyr::filter(merged_ecoli_data$strain %in% clade3_ecoli_tree$tip.label) %>%
  tibble::column_to_rownames("strain")

clade4_ecoli_data <- merged_ecoli_data %>%
  dplyr::filter(merged_ecoli_data$strain %in% clade4_ecoli_tree$tip.label) %>%
  tibble::column_to_rownames("strain")


#Plot Tree Clade 1


p1 <- ggtree(clade1_ecoli_tree) %<+% clade1_ecoli_data +
  geom_tiplab(size = 1.5) 

p1

st_df <- clade1_ecoli_data %>%
  select(st)

lab_df <- clade1_ecoli_data %>%
  select(site)

otype_df <- clade1_ecoli_data %>%
  select(o_type)

htype_df <- clade1_ecoli_data %>%
  select(h_type)


h1 <- gheatmap(p1, st_df, 
               offset = 0.0002,
               width = 0.05, 
               color = "black", 
               colnames = FALSE) +
  scale_fill_manual(name = "ST",
                    values = c("ST12"        = "#1b9e77",
                               "ST131"       = "#d95f02",
                               "ST1193"      = "#7570b3",
                               "ST73-1LV"    = "#e7298a",
                               "ST10106"     = "#66a61e",
                               "ST136-1LV"   = "#e6ab02",
                               "ST73"        = "#a6761d",
                               "ST636"       = "#666666",
                               "ST13793"     = "#1f78b4",
                               "ST706-1LV"   = "#b2df8a",
                               "ST95"        = "#33a02c",
                               "ST131-1LV"   = "#fb9a99",
                               "ST136"       = "#a6cee3",
                               "ST636-1LV"   = "#cab2d6",
                               "ST92-1LV"    = "#ff7f00",
                               "ST127"       = "#fdbf6f",
                               "ST355-1LV"   = "#b15928",
                               "ST1193-1LV"  = "#8dd3c7")) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.box = "vertical", 
        legend.direction = "vertical",
        legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h1

h2 <- h1 + new_scale_fill()

h3 <- gheatmap(h2, otype_df,   
               offset = 0.0004, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "O Type",
                    values = c(
                      "O4"   = "#1f77b4",  # bleu
                      "O16"  = "#ff7f0e",  # orange
                      "O25"  = "#2ca02c",  # vert
                      "O75"  = "#d62728",  # rouge
                      "O6"   = "#9467bd",  # violet
                      "O112" = "#8c564b",  # brun
                      "O8"   = "#e377c2",  # rose
                      "O21"  = "#7f7f7f",  # gris
                      "O18"  = "#bcbd22",  # vert olive
                      "O51"  = "#17becf",  # bleu clair
                      "O83"  = "#aec7e8",  # bleu pâle
                      "O2"   = "#ffbb78"   # orange pâle  
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h3

h4 <- h3 + new_scale_fill()

h5 <- gheatmap(h4, htype_df,   
               offset = 0.0006, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "H Type",
                    values = c(
                      "H1"  = "#1f77b4",  # bleu
                      "H5"  = "#ff7f0e",  # orange
                      "H4"  = "#2ca02c",  # vert
                      "H14" = "#d62728",  # rouge
                      "H7"  = "#9467bd",  # violet
                      "-"   = "#7f7f7f",  # gris pour type inconnu
                      "H10" = "#8c564b",  # brun
                      "H31" = "#17becf"   # bleu clair 
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h5


h6 <- h5 + new_scale_fill()

h7 <- gheatmap(h6, lab_df,   
               offset = 0.0008, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "Site",
                    values = c(
                      "CHE-CNPS"  = "#1f78b4",  # bleu
                      "HRGAROUA"  = "#33a02c",  # vert
                      "SW-06-RHL" = "#e31a1c",  # rouge
                      "HCY"       = "#ff7f00",  # orange
                      "CHUY"      = "#6a3d9a",  # violet foncé
                      "HRBAF"     = "#a6cee3",  # bleu clair
                      "HMR-1"     = "#fb9a99",  # rose clair
                      "HLD"       = "#b2df8a",  # vert clair
                      "HGOPY"     = "#cab2d6"   # mauve clair 
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h7

#Plot Tree Clade 2


p1 <- ggtree(clade2_ecoli_tree) %<+% clade2_ecoli_data +
  geom_tiplab(size = 1.5) 

p1

st_df <- clade2_ecoli_data %>%
  select(st)

lab_df <- clade2_ecoli_data %>%
  select(site)

otype_df <- clade2_ecoli_data %>%
  select(o_type)

htype_df <- clade2_ecoli_data %>%
  select(h_type)

h1 <- gheatmap(p1, st_df, 
               offset = 0.0002,
               width = 0.05, 
               color = "black", 
               colnames = FALSE) +
  scale_fill_manual(name = "ST",
                    values = c("ST410"      = "#1f77b4",  # bleu
                               "ST156"      = "#ff7f0e",  # orange
                               "ST40"       = "#2ca02c",  # vert
                               "ST90"       = "#d62728",  # rouge
                               "ST448"      = "#9467bd",  # violet
                               "ST1196"     = "#8c564b",  # marron
                               "ST424"      = "#e377c2",  # rose
                               "ST410-1LV"  = "#7f7f7f",  # gris
                               "ST155"      = "#bcbd22",  # vert lime
                               "ST295-1LV"  = "#17becf",  # bleu clair
                               "ST224-2LV"  = "#aec7e8",  # bleu pastel
                               "ST29"       = "#ffbb78",  # orange clair
                               "ST5614"     = "#98df8a",  # vert clair
                               "ST156-1LV"  = "#ff9896",  # rouge clair
                               "ST94-1LV"   = "#c5b0d5",  # violet clair
                               "ST58"       = "#c49c94",  # rose terre
                               "ST10955"    = "#f7b6d2",  # rose pâle
                               "ST5981-4LV" = "#dbdb8d"   # jaune pâle
                               )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.box = "vertical", 
        legend.direction = "vertical",
        legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h1

h2 <- h1 + new_scale_fill()

h3 <- gheatmap(h2, otype_df,   
               offset = 0.0004, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "O Type",
                    values = c(
                      "O8"      = "#1f77b4",  # bleu
                      "O100"    = "#ff7f0e",  # orange
                      "O21"     = "#2ca02c",  # vert
                      "O148"    = "#d62728",  # rouge
                      "O76"     = "#9467bd",  # violet
                      "-"       = "#7f7f7f",  # gris (type non spécifié)
                      "O131"    = "#8c564b",  # marron
                      "O64"     = "#e377c2",  # rose
                      "O27"     = "#bcbd22",  # lime
                      "O9/O166" = "#17becf",  # bleu clair
                      "O82"     = "#ff9896",  # rouge clair
                      "O174"    = "#98df8a",  # vert clair
                      "O141"    = "#c5b0d5"  
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h3

h4 <- h3 + new_scale_fill()

h5 <- gheatmap(h4, htype_df,   
               offset = 0.0006, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "H Type",
                    values = c(
                      "H9"  = "#1f77b4",   # bleu
                      "H45" = "#ff7f0e",   # orange
                      "H21" = "#2ca02c",   # vert
                      "H28" = "#d62728",   # rouge
                      "H10" = "#9467bd",   # violet
                      "H12" = "#8c564b",   # marron
                      "H27" = "#e377c2",   # rose
                      "H30" = "#7f7f7f",   # gris
                      "H11" = "#bcbd22",   # lime
                      "H14" = "#17becf",   # turquoise
                      "H8"  = "#aec7e8",   # bleu clair
                      "H19" = "#ffbb78",   # orange clair
                      "H5"  = "#98df8a",   # vert clair
                      "H25" = "#c5b0d5"    # violet clair
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h5


h6 <- h5 + new_scale_fill()

h7 <- gheatmap(h6, lab_df,   
               offset = 0.0008, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "Site",
                    values = c(
                      "CHUY"      = "#1f77b4",   # bleu
                      "HRBUEA"    = "#ff7f0e",   # orange
                      "HRGAROUA"  = "#2ca02c",   # vert
                      "HGOPY"     = "#d62728",   # rouge
                      "SW-06-RHL" = "#9467bd",   # violet
                      "HCY"       = "#8c564b",   # marron
                      "HMR-1"     = "#e377c2",   # rose
                      "HR MAROUA" = "#7f7f7f",   # gris
                      "HRBAF"     = "#bcbd22",   # lime
                      "HLD"       = "#17becf",   # turquoise
                      "CHE-CNPS"  = "#aec7e8"    # bleu clair 
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h7

#Plot Tree Clade 3


p1 <- ggtree(clade3_ecoli_tree) %<+% clade3_ecoli_data +
  geom_tiplab(size = 1.5) 

p1

st_df <- clade3_ecoli_data %>%
  select(st)

lab_df <- clade3_ecoli_data %>%
  select(site)

otype_df <- clade3_ecoli_data %>%
  select(o_type)

htype_df <- clade3_ecoli_data %>%
  select(h_type)


h1 <- gheatmap(p1, st_df, 
               offset = 0.0002,
               width = 0.05, 
               color = "black", 
               colnames = FALSE) +
  scale_fill_manual(name = "ST",
                    values = c("NA"          = "#999999",
                               "ST540"       = "#e41a1c",
                               "ST8763"      = "#377eb8",
                               "ST617"       = "#4daf4a",
                               "ST167"       = "#984ea3",
                               "ST617-2LV"   = "#ff7f00",
                               "ST167-1LV"   = "#a65628",
                               "ST10"        = "#f781bf",
                               "ST1139"      = "#ffff33",
                               "ST226"       = "#a6cee3",
                               "ST206"       = "#1f78b4",
                               "ST617-1LV"   = "#b2df8a",
                               "ST450"       = "#33a02c",
                               "ST46"        = "#fb9a99"
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.box = "vertical", 
        legend.direction = "vertical",
        legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h1

h2 <- h1 + new_scale_fill()

h3 <- gheatmap(h2, otype_df,   
               offset = 0.0004, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "O Type",
                    values = c(
                      "O141"       = "#e41a1c",
                      "O9"         = "#377eb8",
                      "O51"        = "#4daf4a",
                      "-"          = "#999999", 
                      "O8"         = "#984ea3",
                      "O92"        = "#ff7f00",
                      "O153/O178"  = "#a65628",
                      "O54"        = "#f781bf",
                      "O45"        = "#9999ff",
                      "O12"        = "#66c2a5"
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h3

h4 <- h3 + new_scale_fill()

h5 <- gheatmap(h4, htype_df,   
               offset = 0.0006, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "H Type",
                    values = c(
                      "H5"  = "#1b9e77",
                      "H10" = "#d95f02",
                      "H9"  = "#7570b3",
                      "H21" = "#e7298a",
                      "H32" = "#66a61e",
                      "H16" = "#e6ab02"
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h5


h6 <- h5 + new_scale_fill()

h7 <- gheatmap(h6, lab_df,   
               offset = 0.0008, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "Site",
                    values = c(
                      "CHUY"      = "#1f77b4",   # bleu
                      "HRBUEA"    = "#ff7f0e",   # orange
                      "HRGAROUA"  = "#2ca02c",   # vert
                      "HGOPY"     = "#d62728",   # rouge
                      "SW-06-RHL" = "#9467bd",   # violet
                      "HCY"       = "#8c564b",   # marron
                      "HMR-1"     = "#e377c2",   # rose
                      "HR MAROUA" = "#7f7f7f",   # gris
                      "HRBAF"     = "#bcbd22",   # lime
                      "HLD"       = "#17becf",   # turquoise
                      "CHE-CNPS"  = "#aec7e8"    # bleu clair 
                    )) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h7



#Plot Tree Clade 4


p1 <- ggtree(clade4_ecoli_tree) %<+% clade4_ecoli_data +
  geom_tiplab(size = 1.5) 

p1

st_df <- clade4_ecoli_data %>%
  select(st)

lab_df <- clade4_ecoli_data %>%
  select(site)

otype_df <- clade4_ecoli_data %>%
  select(o_type)

htype_df <- clade4_ecoli_data %>%
  select(h_type)

st_colors <- setNames(hue_pal()(length(unique(st_df$st))), unique(st_df$st))

lab_colors <- setNames(hue_pal()(length(unique(lab_df$site))), unique(lab_df$site))

otype_colors <- setNames(hue_pal()(length(unique(otype_df$o_type))), unique(otype_df$o_type))

htype_colors <- setNames(hue_pal()(length(unique(htype_df$h_type))), unique(htype_df$h_type))

h1 <- gheatmap(p1, st_df, 
               offset = 0.0002,
               width = 0.05, 
               color = "black", 
               colnames = FALSE) +
  scale_fill_manual(name = "ST",
                    values = st_colors) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.box = "vertical", 
        legend.direction = "vertical",
        legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h1

h2 <- h1 + new_scale_fill()

h3 <- gheatmap(h2, otype_df,   
               offset = 0.0004, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "O Type",
                    values = otype_colors) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h3

h4 <- h3 + new_scale_fill()

h5 <- gheatmap(h4, htype_df,   
               offset = 0.0006, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "H Type",
                    values = htype_colors) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h5


h6 <- h5 + new_scale_fill()

h7 <- gheatmap(h6, lab_df,   
               offset = 0.0008, 
               width = 0.05,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "Site",
                    values = lab_colors) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(ncol = 2,byrow = TRUE))

h7

