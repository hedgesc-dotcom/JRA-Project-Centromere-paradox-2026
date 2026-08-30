
#loading all packages to build the foundations of this project 
library(tidyverse)
library(stringr)
library(janitor)

#loading in major files

tsv_dir <- "C:/Users/hedge/OneDrive/Documents/R_JRA_projects/90 new files"

master_table <- read_csv(
  "C:/Users/hedge/OneDrive/Documents/R_JRA_projects/TE_results/LTR_lineage_master_table.csv"
)

#===========================================================
# GET ALL CLS FILES
#===========================================================

cls_files <- list.files(
  tsv_dir,
  pattern = "\\.cls\\.tsv$",
  full.names = TRUE
)

#===========================================================
# FUNCTION TO PROCESS ONE FILE
#===========================================================

process_cls_file <- function(file){
  
  message("Processing: ", basename(file))
  
  df <- read_tsv(
    file,
    show_col_types = FALSE
  )
  
  df <- df %>%
    rename(
      TE = `#TE`
    )
  #-----------------------------------
  # Species name from filename
  #-----------------------------------
  
  species_id <- basename(file) |>
    str_remove("_edta.*")
  
  df <- df %>%
    mutate(
      Species_ID = species_id,
      
      status = case_when(
        str_detect(TE, "intact") ~ "intact",
        str_detect(TE, "fragment") ~ "fragment",
        TRUE ~ "other"
      )
    )
  
  #-----------------------------------
  # Coordinates
  #-----------------------------------
  
  coords <- str_match(
    df$TE,
    ":(\\d+)\\.\\.(\\d+)"
  )
  
  df$Start <- as.numeric(coords[,2])
  df$End <- as.numeric(coords[,3])
  
  df <- df %>%
    mutate(
      Length_bp = End - Start + 1
    )
  
  #-----------------------------------
  # Domain presence
  #-----------------------------------
  
  df <- df %>%
    mutate(
      GAG  = if_else(str_detect(Domains, "GAG\\|"), 1, 0),
      PROT = if_else(str_detect(Domains, "PROT\\|"), 1, 0),
      RT   = if_else(str_detect(Domains, "RT\\|"), 1, 0),
      RH   = if_else(str_detect(Domains, "RH\\|"), 1, 0),
      INT  = if_else(str_detect(Domains, "INT\\|"), 1, 0)
    )
  
  #-----------------------------------
  # Coding capacity
  #-----------------------------------
  
  df <- df %>%
    mutate(
      Coding_capacity =
        GAG +
        PROT +
        RT +
        RH +
        INT
    )
  
  return(df)
}

#===========================================================
# PROCESS ALL FILES
#===========================================================

all_elements <- map_dfr(
  cls_files,
  process_cls_file
)

#===========================================================
# SPECIES METADATA
#===========================================================


taxonomy <- read_csv(
  "C:/Users/hedge/Downloads/species taxonomy.csv"
)

species_metadata <- taxonomy %>%
  mutate(
    Genome = str_remove(fasta, "\\.fa$"),
    
    # Fix mismatched genome names
    Genome = case_when(
      Genome == "daTanVulg1.hap1.1" ~ "daTanVulg1.1",
      Genome == "rosCan_S27_v1.fasta.F2B.ChrOnly" ~
        "rosCan_S27_v1",
      TRUE ~ Genome
    )
  ) %>%
  transmute(
    Species_ID = Genome,
    Major_group = taxa1,
    Taxonomic_Order = taxa2,
    Genus = genus,
    Species = species
  )
#===========================================================
# JOIN METADATA
#===========================================================


all_elements <- all_elements %>%
  left_join(
    species_metadata,
    by = "Species_ID"
  )



#===========================================================
# CLEAN TABLE
#===========================================================

element_master_table <- all_elements %>%
  rename(
    TE_Order = Order
  ) %>%
  select(
    Major_group,
    Taxonomic_Order,
    Genus,
    Species,
    Species_ID,
    TE,
    TE_Order,
    Superfamily,
    Clade,
    status,
    Complete,
    Start,
    End,
    Length_bp,
    GAG,
    PROT,
    RT,
    RH,
    INT,
    Coding_capacity,
    Domains
  )
#===========================================================
# EXPORT
#===========================================================

write_csv(
  element_master_table,
  "TE_element_master_table.csv"
)


saveRDS(
  element_master_table,
  "TE_element_master_table.rds"
)






files <- list.files(
  "C:/Users/hedge/OneDrive/Documents/R_JRA_projects/90 new files",
  pattern = "\\.tsv$",
  full.names = TRUE
)
filtered_list <- list()

output_dir <- "C:/Users/hedge/OneDrive/Documents/R_JRA_projects/TE_results"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

master_list <- list()
for(file in files){
  #for(file in files[1:2]){
  
  tryCatch({
    
    
    
    cat(
      "\n============================\n",
      "Starting:",
      basename(file),
      "\n============================\n"
    )
    
    
    species <- sub(
      "_edta_filtered.*$",
      "",
      basename(file)
    )
    
    species_dir <- file.path(
      output_dir,
      species
    )
    
    dir.create(
      species_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    cat(
      "\nProcessing:",
      species,
      "\n"
    )
    
    # ==========================================================
    # LOAD DATA
    # ==========================================================
    
    te.data <- read_tsv(
      file,
      show_col_types = FALSE
    )
    
    colnames(te.data) <- gsub(
      "^#",
      "",
      colnames(te.data)
    )
    
    
    
    # ==========================================================
    # EVERYTHING ELSE FROM YOUR PIPELINE
    # ==========================================================
    
    
    
    te.data  <- te.data  %>%
      clean_names()
    
    # ============================================================
    # Superfamily assignment
    # ============================================================
    
    df_superfamily <- te.data  %>%
      mutate(
        
        edta_superfamily = case_when(
          
          str_detect(te, regex("gypsy", TRUE)) ~ "Gypsy",
          str_detect(te, regex("copia", TRUE)) ~ "Copia",
          
          str_detect(te, regex("LTR_retrotransposon", TRUE)) ~ "LTR_other",
          
          str_detect(te, regex("mutator", TRUE)) ~ "Mutator",
          str_detect(te, regex("cacta", TRUE)) ~ "CACTA",
          str_detect(te, regex("harbinger|pong", TRUE)) ~ "Harbinger",
          str_detect(te, regex("hat", TRUE)) ~ "hAT",
          str_detect(te, regex("mariner|tc1", TRUE)) ~ "Tc1_Mariner",
          
          str_detect(te, regex("merlin", TRUE)) ~ "Merlin",
          str_detect(te, regex("sola1", TRUE)) ~ "Sola1",
          
          str_detect(te, regex("helitron", TRUE)) ~ "Helitron",
          
          str_detect(te, regex("LINE|L1_LINE", TRUE)) ~ "LINE",
          str_detect(te, regex("SINE", TRUE)) ~ "SINE",
          
          str_detect(te, regex("DNA_transposon", TRUE)) ~ "DNA_other",
          
          TRUE ~ "unknown"
        ),
        
        tesorter_superfamily = case_when(
          superfamily == "MuDR_Mutator" ~ "Mutator",
          superfamily == "EnSpm_CACTA" ~ "CACTA",
          superfamily == "PIF_Harbinger" ~ "Harbinger",
          is.na(superfamily) ~ "unknown",
          TRUE ~ superfamily
        )
      )
    
    # ============================================================
    # Filter dataset
    # ============================================================
    
    df_filtered <- df_superfamily %>%
      filter(
        edta_superfamily %in% c(
          "LTR_other",
          "Gypsy",
          "Copia"
        ) |
          (
            edta_superfamily == tesorter_superfamily &
              !edta_superfamily %in% c(
                "LTR_other",
                "Gypsy",
                "Copia",
                "unknown"
              )
          )
      ) %>%
      mutate(
        final_superfamily = tesorter_superfamily,
        
        start = as.numeric(
          str_extract(te, "(?<=:)[0-9]+")
        ),
        
        end = as.numeric(
          str_extract(te, "(?<=\\.\\.)[0-9]+")
        ),
        
        length_bp = end - start + 1,
        
        status = case_when(
          str_detect(te, "intact") ~ "intact",
          str_detect(te, "fragment") ~ "fragment",
          TRUE ~ "other"
        )
      )
    
    
    filtered_list[[species]] <- df_filtered %>%
      mutate(
        Genome = species
        
      )
    
    
    
    # ============================================================
    # LTR lineage summary
    # ============================================================
    
    ltr_summary <- df_filtered %>%
      filter(
        final_superfamily %in% c(
          "Gypsy",
          "Copia"
        )
      ) %>%
      group_by(
        final_superfamily,
        clade
      ) %>%
      summarise(
        Number_of_elements = n(),
        
        Intact_elements =
          sum(status == "intact"),
        Total_DNA_bp =
          sum(length_bp, na.rm = TRUE),
        Intact_DNA_bp =
          sum(
            ifelse(
              status == "intact",
              length_bp,
              0
            ),
            na.rm = TRUE
          ),
        .groups = "drop"
      )
    
    total_ltr_bp <-
      sum(ltr_summary$Total_DNA_bp)
    
    total_intact_ltr_bp <-
      sum(ltr_summary$Intact_DNA_bp)
    
    ltr_summary <- ltr_summary %>%
      mutate(
        Percent_of_all_LTR_DNA =
          round(
            100 * Total_DNA_bp /
              total_ltr_bp,
            1
          ),
        Percent_of_all_intact_LTR_DNA =
          round(
            100 * Intact_DNA_bp /
              total_intact_ltr_bp,
            1
          )
      )
    
    
    
    
    
    
    
    
    
    
    # ============================================================
    # Add summary rows
    # ============================================================
    
    ltr_summary <- ltr_summary %>%
      
      rename(
        
        Superfamily = final_superfamily,
        
        Lineage = clade
        
      )
    
    
    # ============================================================
    # MASTER TABLE
    # ============================================================
    
    
    master_species <- ltr_summary %>%
      
      filter(
        !Lineage %in% c(
          "unknown",
          "mixture"
        )
      ) %>%
      
      mutate(
        
        Genome = as.character(species)
      )
    
    master_list[[species]] <- master_species
    
  }, error = function(e){
    
    cat(
      "\nERROR:",
      species,
      "\n",
      conditionMessage(e),
      "\n"
    )
    
  })
  
}






#===========================================================
# COMBINE ALL SPECIES
#===========================================================

master_table <- bind_rows(
  master_list
)


taxonomy <- read_csv(
  "C:/Users/hedge/Downloads/species taxonomy.csv",
  show_col_types = FALSE
) %>%
  
  rename(
    Genome = fasta,
    Major_group = taxa1,
    Order = taxa2,
    Genus = genus,
    Species = species
  ) %>%
  
  mutate(
    
    Genome = sub(
      "\\.fa$",
      "",
      Genome
    ),
    
    Genome = sub(
      "\\.fasta\\.F2B\\.ChrOnly$",
      "",
      Genome
    ),
    
    Genome = ifelse(
      Genome == "daTanVulg1.hap1.1",
      "daTanVulg1.1",
      Genome
    )
    
  )

master_table <- master_table %>%
  
  left_join(
    taxonomy,
    by = "Genome"
  )





length(master_list)
names(master_list)

master_table <- master_table %>%
  
  select(
    
    Major_group,
    Order,
    Genus,
    Species,
    
    Genome,
    
    Superfamily,
    Lineage,
    
    Number_of_elements,
    
    Intact_elements,
    
    Total_DNA_bp,
    Intact_DNA_bp,
    
    
    
    Percent_of_all_LTR_DNA,
    
    
    Percent_of_all_intact_LTR_DNA
    
  ) %>%
  
  arrange(
    
    Major_group,
    Order,
    Genus,
    Species,
    
    Superfamily,
    
    desc(Percent_of_all_intact_LTR_DNA)
    
  )




#===========================================================
# COMBINE ALL FILTERED TE RECORDS
#===========================================================

all_filtered <- bind_rows(
  filtered_list
)

write_csv(
  all_filtered,
  file.path(
    output_dir,
    "All_filtered_TE_records.csv"
  )
)
#===========================================================
# CHECK TOP CHROMO-OUTGROUP SPECIES
#===========================================================

top_chromo_species <- master_table %>%
  filter(
    Lineage == "chromo-outgroup"
  ) %>%
  arrange(
    desc(Percent_of_all_LTR_DNA)
  ) %>%
  select(
    Genome,
    Major_group,
    Order,
    Genus,
    Species,
    Number_of_elements,
    Intact_elements,
    Total_DNA_bp,
    Intact_DNA_bp,
    Percent_of_all_LTR_DNA,
    Percent_of_all_intact_LTR_DNA
  ) %>%
  slice_head(n = 20)

print(
  top_chromo_species,
  n = Inf
)







#===========================================================
# PRETTY VERSION FOR HUMAN READING
#===========================================================

pretty_table <- master_table

# blank repeated taxonomy / genome fields

pretty_table <- pretty_table %>%
  
  group_by(
    Major_group,
    Order,
    Genus,
    Species,
    Genome
  ) %>%
  
  mutate(
    
    Major_group = ifelse(
      row_number() == 1,
      Major_group,
      ""
    ),
    
    Order = ifelse(
      row_number() == 1,
      Order,
      ""
    ),
    
    Genus = ifelse(
      row_number() == 1,
      Genus,
      ""
    ),
    
    Species = ifelse(
      row_number() == 1,
      Species,
      ""
    ),
    
    Genome = ifelse(
      row_number() == 1,
      Genome,
      ""
    )
    
  ) %>%
  
  ungroup()

# blank repeated superfamily names

pretty_table <- pretty_table %>%
  
  mutate(
    Superfamily = ifelse(
      duplicated(
        paste(
          Major_group,
          Order,
          Genus,
          Species,
          Genome,
          Superfamily
        )
      ),
      "",
      Superfamily
    )
  )


lineage_order <- master_table %>%
  filter(Lineage != "Grand Total") %>%
  group_by(Lineage) %>%
  summarise(
    median_occ = median(
      Percent_of_all_LTR_DNA,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(median_occ)) %>%
  pull(Lineage)

master_table <- master_table %>%
  mutate(
    Lineage = factor(
      Lineage,
      levels = rev(lineage_order)
    )
  )

master_table <- master_table %>%
  filter(Genome != "All")


#===========================================================
# ARCHITECTURES BY LINEAGE AND PLANT GROUP
#===========================================================



lineages_keep <- master_table %>%
  
  filter(
    Genus != "All",
    Species != "All"
  ) %>%
  
  group_by(Lineage) %>%
  
  summarise(
    
    n_species_1pct =
      sum(
        Percent_of_all_LTR_DNA >= 1,
        na.rm = TRUE
      ),
    
    max_occupancy =
      max(
        Percent_of_all_LTR_DNA,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  ) %>%
  
  filter(
    n_species_1pct > 1 |
      max_occupancy > 4
  ) %>%
  
  pull(Lineage)


write_csv(
  master_table,
  file.path(
    output_dir,
    "LTR_lineage_master_table.csv"
  )
)

write_csv(
  pretty_table,
  file.path(
    output_dir,
    "LTR_lineage_master_table_pretty.csv"
  )
)




cat(
  "\nSaved: LTR_lineage_master_table.csv\n"
)

setdiff(
  unique(master_table$Genome),
  unique(taxonomy$Genome)
)



master_table %>%
  filter(is.na(Species))






#Boxplots

#===========================================================
# OCCUPANCY BOXPLOTS
#===========================================================

library(patchwork)

#-----------------------------------------------------------
# Fill missing lineage observations with zeros
#-----------------------------------------------------------
all_lineages <- master_table %>%
  filter(Lineage != "Grand Total") %>%
  distinct(
    Lineage,
    Superfamily
  )

all_genomes <- master_table %>%
  filter(Genome != "All") %>%
  distinct(Genome)

master_table_complete <- all_genomes %>%
  crossing(all_lineages) %>%
  left_join(
    master_table,
    by = c(
      "Genome",
      "Lineage",
      "Superfamily"
    )
    
  ) %>%
  mutate(
    Number_of_elements =
      replace_na(Number_of_elements, 0),
    
    Intact_elements =
      replace_na(Intact_elements, 0),
    
    Total_DNA_bp =
      replace_na(Total_DNA_bp, 0),
    
    Intact_DNA_bp =
      replace_na(Intact_DNA_bp, 0),
    
    Percent_of_all_LTR_DNA =
      replace_na(Percent_of_all_LTR_DNA, 0),
    
    Percent_of_all_intact_LTR_DNA =
      replace_na(Percent_of_all_intact_LTR_DNA, 0)
  )

#-----------------------------------------------------------
# Number of species containing each lineage
#-----------------------------------------------------------

lineage_counts <- master_table %>%
  filter(
    Lineage != "Grand Total",
    Number_of_elements > 0
  ) %>%
  group_by(Lineage) %>%
  summarise(
    species_count = n_distinct(Genome),
    .groups = "drop"
  )

lineage_labels <- setNames(
  paste0(
    lineage_counts$Lineage,
    " (",
    lineage_counts$species_count,
    " spp)"
  ),
  lineage_counts$Lineage
)

#-----------------------------------------------------------
# Total occupancy plot
#-----------------------------------------------------------

occupancy_total_plot <- ggplot(
  master_table_complete,
  aes(
    x = Lineage,
    y = Percent_of_all_LTR_DNA,
    fill = Superfamily
  )
) +
  
  geom_boxplot(
    alpha = 0.8,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.5,
    size = 1.2,
    colour = "black"
  ) +
  
  coord_flip() +
  
  scale_x_discrete(
    labels = lineage_labels
  ) +
  
  scale_y_continuous(
    limits = c(0, 100)
  ) +
  
  scale_fill_manual(
    values = c(
      Copia = "#08519C",
      Gypsy = "#D94801"
    )
  ) +
  
  labs(
    title = "Total LTR DNA Space",
    x = "LTR Lineage (Number of Species)",
    y = "Occupancy of Total LTR DNA Space (%)",
    caption = "Points represent individual species. Species lacking a lineage were included as zero values."
  ) +
  
  theme_bw(
    base_size = 14
  ) +
  
  theme(
    plot.title =
      element_text(face = "bold")
  )

#-----------------------------------------------------------
# Intact occupancy plot
#-----------------------------------------------------------

occupancy_intact_plot <- ggplot(
  master_table_complete,
  aes(
    x = Lineage,
    y = Percent_of_all_intact_LTR_DNA,
    fill = Superfamily
  )
) +
  
  geom_boxplot(
    alpha = 0.8,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.5,
    size = 1.2,
    colour = "black"
  ) +
  
  coord_flip() +
  
  scale_x_discrete(
    labels = lineage_labels
  ) +
  
  scale_y_continuous(
    limits = c(0, 100)
  ) +
  
  scale_fill_manual(
    values = c(
      Copia = "#08519C",
      Gypsy = "#D94801"
    )
  ) +
  
  labs(
    title = "Intact LTR DNA Space",
    x = "LTR Lineage (Total Elements)",
    y = "Occupancy of Intact LTR DNA Space (%)"
  ) +
  
  theme_bw(
    base_size = 14
  ) +
  
  theme(
    plot.title =
      element_text(face = "bold")
  )

#-----------------------------------------------------------
# Combine plots
#-----------------------------------------------------------

combined_plot <-
  occupancy_total_plot +
  occupancy_intact_plot +
  plot_layout(ncol = 2)

#-----------------------------------------------------------
# Save
#-----------------------------------------------------------

ggsave(
  file.path(
    output_dir,
    "LTR_lineage_DNA_space_occupancy_boxplots.png"
  ),
  combined_plot,
  width = 22,
  height = 12,
  dpi = 600,
  bg = "white"
)




#chromo groups show a high occupancy which is unusual


all_filtered %>%
  filter(
    clade %in% c(
      "chromo-outgroup",
      "chromo-unclass",
      "non-chromo-outgroup"
    )
  ) %>%
  count(
    Genome,
    clade,
    domains,
    sort = TRUE
  ) %>%
  print(n = 100)


chromo_summary <- all_filtered %>%
  filter(
    clade %in% c(
      "chromo-outgroup",
      "chromo-unclass",
      "non-chromo-outgroup"
    )
  ) %>%
  count(
    Genome,
    clade,
    domains,
    sort = TRUE
  )

write_csv(
  chromo_summary,
  file.path(
    output_dir,
    "Chromo_group_domain_summary.csv"
  )
)





#chromooutgroups are when te sorter identifies it as a chromogroup but cannot assign it to a known chromovirus lineage
#the values for these unknown chromo groups are pretty high. Several bryophyte genomes contain unexpectedly large amounts 
#of chromovirus-related Gypsy DNA, and much of this diversity appears to fall outside the currently recognised chromovirus 
#lineages, resulting in TEsorter assigning many elements to chromo-outgroup and chromo-unclass categories.











#Building the Alignment tree

species_list <- master_table %>%
  distinct(Genus, Species) %>%
  unite("Species_name", Genus, Species, sep = " ")




species_for_timetree <- master_table %>%
  distinct(Genus, Species) %>%
  arrange(Genus, Species) %>%
  mutate(
    Species_name = paste(Genus, Species)
  ) %>%
  pull(Species_name)

writeLines(
  species_for_timetree,
  file.path(
    output_dir,
    "TimeTree_species_list.txt"
  )
)

#statistical tests

#do occupancy differ from each other
kruskal.test(
  Percent_of_all_LTR_DNA ~ Lineage,
  data = master_table
)

#which one differs the most

library(FSA)

dunn_results <- dunnTest(
  Percent_of_all_LTR_DNA ~ Lineage,
  data = master_table,
  method = "bh"
)

dunn_results$res


plot_data <- element_master_table %>%
  
  filter(
    status == "intact"
  ) %>%
  
  rename(
    Lineage = Clade
  ) %>%
  
  filter(
    Lineage %in% lineages_keep
  ) %>%
  
  filter(
    !Lineage %in% c(
      "unknown",
      "chromo-outgroup",
      "chromo-unclass",
      "non-chromo-outgroup"
    )
  ) %>%
  
  mutate(
    
    Group_panel =
      case_when(
        
        Major_group == "Monocots" ~
          "Monocots",
        
        Major_group == "Dicots" ~
          "Dicots",
        
        Major_group %in%
          c(
            "Algae",
            "Bryophyta"
          ) ~
          "Algae + Bryophytes"
        
      )
    
  )

p <- ggplot(
  plot_data,
  aes(
    x = Group_panel,
    y = Length_bp,
    colour = Group_panel
  )
) +
  
  ggbeeswarm::geom_quasirandom(
    width = 0.25,
    alpha = 0.03,
    size = 0.15
  ) +
  
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    colour = "black",
    linewidth = 0.4
  ) +
  
  facet_wrap(
    ~ Lineage,
    ncol = 5
  ) +
  
  scale_colour_manual(
    values = c(
      "Monocots" = "#8DA0CB",
      "Dicots" = "#E78AC3",
      "Algae + Bryophytes" = "#66C2A5"
    )
  ) +
  
  labs(
    title =
      "Lengths of Intact Elements from Successful LTR-Retrotransposon Lineages",
    
    subtitle =
      "Points represent individual intact elements; horizontal lines show medians",
    
    x =
      "",
    
    y =
      "Length (bp)"
  ) +
  
  theme_bw()

ggsave(
  "Lengths_by_lineage_and_group.png",
  p,
  width = 20,
  height = 14,
  dpi = 600,
  bg = "white"
)


#-----------------------------------------------------------
# ELEMENT DATA
#-----------------------------------------------------------

plot_df <- element_master_table %>%
  
  filter(
    status == "intact",
    Superfamily %in% c(
      "Gypsy",
      "Copia"
    )
  ) %>%
  
  rename(
    Lineage = Clade
  ) %>%
  
  filter(
    Lineage %in% lineages_keep
  ) %>%
  
  filter(
    !Lineage %in% c(
      "unknown",
      "chromo-outgroup",
      "chromo-unclass",
      "non-chromo-outgroup"
    )
  ) %>%
  
  mutate(
    
    Architecture =
      str_extract_all(
        Domains,
        "GAG|PROT|RT|RH|INT|CHDCR"
      ) %>%
      sapply(
        function(x)
          paste(
            unique(x),
            collapse = " + "
          )
      ),
    
    Architecture =
      ifelse(
        Architecture == "",
        "No domains",
        Architecture
      ),
    
    Group_panel =
      case_when(
        
        Major_group == "Monocots" ~
          "Monocots",
        
        Major_group == "Dicots" ~
          "Dicots",
        
        Major_group %in%
          c(
            "Algae",
            "Bryophyta"
          ) ~
          "Algae + Bryophytes"
        
      )
    
  )

#-----------------------------------------------------------
# ADD "ALL"
#-----------------------------------------------------------

all_panel <- plot_df %>%
  
  mutate(
    Group_panel = "All"
  )

plot_df <- bind_rows(
  plot_df,
  all_panel
)

#-----------------------------------------------------------
# TOP 5 ARCHITECTURES
#-----------------------------------------------------------

top_arch <- plot_df %>%
  
  count(
    Lineage,
    Group_panel,
    Architecture
  ) %>%
  
  group_by(
    Lineage,
    Group_panel
  ) %>%
  
  slice_max(
    n,
    n = 5,
    with_ties = FALSE
  ) %>%
  
  ungroup()

plot_df <- plot_df %>%
  
  semi_join(
    top_arch,
    by = c(
      "Lineage",
      "Group_panel",
      "Architecture"
    )
  )

#-----------------------------------------------------------
# SUMMARISE
#-----------------------------------------------------------

plot_summary <- plot_df %>%
  
  count(
    Lineage,
    Group_panel,
    Architecture
  )

#-----------------------------------------------------------
# ORDER LINEAGES
#-----------------------------------------------------------

lineage_order <- plot_df %>%
  
  count(
    Lineage,
    name = "Total"
  ) %>%
  
  arrange(
    desc(Total)
  ) %>%
  
  pull(Lineage)

plot_summary <- plot_summary %>%
  
  mutate(
    Lineage =
      factor(
        Lineage,
        levels = lineage_order
      ),
    
    Group_panel =
      factor(
        Group_panel,
        levels = c(
          "All",
          "Monocots",
          "Dicots",
          "Algae + Bryophytes"
        )
      )
  )

#-----------------------------------------------------------
# PLOT
#-----------------------------------------------------------

p <- ggplot(
  plot_summary,
  aes(
    x = n,
    y = Architecture
  )
) +
  
  geom_col(
    aes(fill = Group_panel),
    colour = "black",
    linewidth = 0.15
  ) +
  
  scale_fill_manual(
    values = c(
      "All" = "grey70",
      "Monocots" = "#8DA0CB",
      "Dicots" = "#E78AC3",
      "Algae + Bryophytes" = "#66C2A5"
    )
  ) +
  
  facet_grid(
    Lineage ~ Group_panel,
    scales = "free_y"
  ) +
  
  scale_x_continuous(
    labels = scales::comma
  ) +
  
  labs(
    title =
      "Domain Architectures of Successful LTR-Retrotransposon Lineages",
    
    subtitle =
      "Top architectures within each lineage and plant group",
    
    x =
      "Number of Intact Elements",
    
    y =
      "Architecture"
  ) +
  
  theme_bw(
    base_size = 12
  )

ggsave(
  "Architecture_by_lineage_and_group.png",
  p,
  width = 18,
  height = 22,
  dpi = 600,
  bg = "white"
)

library(tidyverse)

plot_df <- element_master_table %>%
  
  filter(
    status == "intact"
  ) %>%
  
  rename(
    Lineage = Clade
  ) %>%
  
  filter(
    Lineage %in% lineages_keep
  ) %>%
  
  mutate(
    
    Major_group =
      factor(
        Major_group,
        levels = c(
          "Algae",
          "Bryophyta",
          "Monocots",
          "Dicots"
        )
      )
    
  )

p <- ggplot(
  plot_df,
  aes(
    x = Length_bp,
    fill = Major_group
  )
) +
  
  geom_histogram(
    binwidth = 500,
    position = "identity",
    alpha = 0.4
  ) +
  
  facet_wrap(
    ~ Lineage,
    scales = "free_y",
    ncol = 5
  ) +
  
  scale_fill_manual(
    values = c(
      "Algae"="#66C2A5",
      "Bryophyta"="#FC8D62",
      "Monocots"="#8DA0CB",
      "Dicots"="#E78AC3"
    )
  ) +
  
  labs(
    title =
      "Length Distributions of Intact TE Elements",
    
    x =
      "Length (bp)",
    
    y =
      "Frequency",
    
    fill =
      "Major Group"
  ) +
  
  theme_bw()

ggsave(
  "Length_histograms_by_lineage.png",
  p,
  width = 20,
  height = 14,
  dpi = 600
)
















library(tidyverse)

bin_size <- 100000

plot_df <- element_master_table %>%
  
  filter(
    status == "intact"
  ) %>%
  
  rename(
    Lineage = Clade
  ) %>%
  
  filter(
    Lineage %in% lineages_keep
  )

occupancy_df <- plot_df %>%
  
  rowwise() %>%
  
  mutate(
    
    Bin =
      list(
        seq(
          floor(Start / bin_size),
          floor(End / bin_size)
        )
      )
    
  ) %>%
  
  unnest(Bin)

hist_df <- occupancy_df %>%
  
  count(
    Lineage,
    Major_group,
    Bin
  ) %>%
  
  mutate(
    
    Position_bp =
      Bin * bin_size
    
  )

p <- ggplot(
  hist_df,
  aes(
    x = Position_bp,
    y = n,
    fill = Major_group
  )
) +
  
  geom_col(
    width = bin_size
  ) +
  
  facet_wrap(
    ~ Lineage,
    scales = "free_y",
    ncol = 5
  ) +
  
  scale_fill_manual(
    values = c(
      "Algae"="#66C2A5",
      "Bryophyta"="#FC8D62",
      "Monocots"="#8DA0CB",
      "Dicots"="#E78AC3"
    )
  ) +
  
  scale_x_continuous(
    labels = scales::comma
  ) +
  
  labs(
    title =
      "Genomic Distribution of Intact TE Elements",
    
    x =
      "Genome Position (bp)",
    
    y =
      "Frequency",
    
    fill =
      "Major Group"
  ) +
  
  theme_bw()

ggsave(
  "Genome_position_histograms_by_lineage.png",
  p,
  width = 20,
  height = 14,
  dpi = 600
)

length_plot_data <- element_master_table %>%
  
  filter(
    status == "intact"
  ) %>%
  
  select(
    Major_group,
    Clade,
    Length_bp
  )

write_csv(
  length_plot_data,
  "TE_lengths_for_python.csv"
)


nrow(length_plot_data)


set.seed(123)

length_plot_sample <- length_plot_data %>%
  
  sample_n(
    100000
  )

write_csv(
  length_plot_sample,
  "TE_lengths_sample.csv"
)


# work out number of total elements vs number of intact elements. 

element_master_table %>%
  
  summarise(
    
    Total_elements = n(),
    
    Intact_elements =
      sum(
        status == "intact",
        na.rm = TRUE
      )
    
  )