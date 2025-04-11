rm(list = ls())
setwd("C:/Users/chris/OneDrive - Stellenbosch University/Documents/sANNIC_workflow")
# rstudioapi::writeRStudioPreference("data_viewer_max_columns", 100L)

## Load packages ###############################################################

#install.packages("remotes")
#install.pachages("devtools")
#devtools::install_github("matildabrown/rWCVP")
#remotes::install_github('matildabrown/rWCVPdata')


library('tidyverse')
library('rWCVP')
#install.packages("rWCVPdata", repos=c("https://matildabrown.github.io/drat", "https://cloud.r-project.org"))
library(rWCVPdata)

wcvp_check_version()
wcvp_version()


## Define functions ############################################################

apply_special_cases <- function (df) {
  # This function will input the full dataset, (native and non-native),
  # add a column named "special_case" to all observations in the input csv,
  # and change the taxon name to the special_cases csv input names
  special_cases_applied<-
    df%>%
    mutate(special_case = case_when(taxon_species_name %in% filter(Special_cases,Rank == "Species")$Name ~ T,
                                    taxon_genus_name %in% filter(Special_cases,Rank == "Genus")$Name ~ T,
                                    taxon_subfamily_name %in% filter(Special_cases,Rank == "Subfamily")$Name  ~ T,
                                    taxon_family_name %in% filter(Special_cases,Rank == "Family")$Name  ~ T,
                                    taxon_id %in% filter(Special_cases,Rank == "Special")$TaxonID ~ T,
                                    TRUE ~ F))%>%
    mutate(taxon = case_when(special_case == T &
                               taxon_species_name %in% filter(Special_cases,Rank == "Species")$Name ~ paste(taxon_species_name),
                             special_case == T &
                               taxon_genus_name %in% filter(Special_cases,Rank == "Genus")$Name ~ paste(taxon_genus_name),
                             special_case == T &
                               taxon_subfamily_name %in% filter(Special_cases,Rank == "Subfamily")$Name  ~ paste(taxon_subfamily_name),
                             special_case == T &
                               taxon_family_name %in% filter(Special_cases,Rank == "Family")$Name  ~ paste(taxon_family_name),
                             special_case == T &
                               taxon_id %in% filter(Special_cases,Rank == "Special",Name == "sect Oenothera")$TaxonID ~ "sect Oenothera",
                             TRUE ~ taxon))
  return(special_cases_applied)
}

multi_correct <- function (df){
  # Rules of the function - multi_correct
  # If no multiple matches, keep
  # If multiple matches and one is accepted, just keep accepted
  # If multiple matches, none are accepted and just one is synonym, keep synonym
  
  #synonym_codes <- c("Synonym", "Orthographic", "Artificial Hybrid", "Unplaced")
  
  synonym_codes <- c("Synonym")
  
  accepted_matches <-
    df %>%
    group_by(taxon)%>%
    mutate(n_matches = n())%>%
    mutate(multi_correct = case_when(multiple_matches == FALSE ~ "keeper",
                                     multiple_matches == TRUE & wcvp_status == "Accepted" ~ "keeper",
                                     TRUE ~ "remove"))%>%
    mutate(multi_correct = case_when(multiple_matches == TRUE & 
                                       wcvp_status%in%synonym_codes == TRUE & 
                                       "keeper" %in% multi_correct == FALSE ~ "keeper",
                                     TRUE ~ multi_correct))%>%
    ungroup()
  
  #Shows species for which no legitimate name was found
  no_legit_match<-
    df%>%
    distinct(taxon)%>%
    mutate(keeper_found = (df%>%distinct(taxon))$taxon %in% (accepted_matches%>%filter(multi_correct == "keeper"))$taxon)%>%
    filter(keeper_found == F)
  
  #Shows species for which multiple legitimate names were found, `dups` is the number of multiples found
  multi_legit_match<-
    accepted_matches%>%
    filter(multi_correct == "keeper")%>%
    group_by(taxon)%>%
    summarise(dups = length(taxon))%>%
    filter(dups>1)
  
  accepted_matches2<-
    accepted_matches%>%
    mutate(multi_correct = ifelse(taxon%in%no_legit_match$taxon | taxon%in%multi_legit_match$taxon,"manual_check_needed",multi_correct))
  
  return(accepted_matches2)
  
}

load_special_cases <- function (file_path = " ",csv2 = T){  
  if (file.exists(file_path)) {
    if (csv2) {
      Special_cases<-read_csv2(file_path)
      Special_cases<-
        Special_cases%>%
        mutate(Name = str_trim(Name),
               Rank = str_trim(Rank),
               TaxonID = str_trim(TaxonID))%>%
        separate_rows(TaxonID,sep = ",")
    } else {
      Special_cases<-read_csv(file_path)
      Special_cases<-
        Special_cases%>%
        mutate(Name = str_trim(Name),
               Rank = str_trim(Rank),
               TaxonID = str_trim(TaxonID))%>%
        separate_rows(TaxonID,sep = ",")}
  } else
  {Special_cases<-data.frame(Name = character(),
                             Rank = character(),
                             TaxonID_incl = character())
  
  message("Couldn't find file path. Assuming no special cases are applied.")
  }
  return(Special_cases)
}

load_native_taxa <- function (file_path = " ",csv2 = T){
  if (file.exists(file_path)) {
    if (csv2) {
      manual_native_input<-read_csv2(file_path)
    } else {
      manual_native_input<-read_csv(file_path)}
  } else
  {manual_native_input<-data.frame(Name = character())
  message("Couldn't find file path. Assuming no manual native input is applied.")
  }
  return(manual_native_input)
}

load_nonnative_taxa <- function (file_path = " ",csv2 = T){
  if (file.exists(file_path)) {
    if (csv2) {
      manual_nonnative_input<-read_csv2(file_path)
    } else {
      manual_nonnative_input<-read_csv(file_path)}
  } else
  {manual_nonnative_input<-data.frame(Name = character())
  message("Couldn't find file path. Assuming no manual non-native input is applied.")
  }
  return(manual_nonnative_input)
}

load_manual_checked_data <- function (file_path = " ",csv2 = T){
  if (file.exists(file_path)) {
    if (csv2) {
      manual_checked_data<-read_csv2(file_path)%>%
        select(-any_of("X"))
    } else {
      manual_checked_data<-read_csv(file_path)%>%
        select(-any_of("X"))}
  } else{
    manual_checked_data <- data.frame(taxon = character(),
                                      checked = logical(),
                                      match_correct = logical(),
                                      match_type = character(),
                                      wcvp_id_full = numeric()
    )
    message("Couldn't find file path. All manual_check_needed entries will be discarded.")
  }
  return(manual_checked_data)
}

## Input #######################################################################


inat_input<-read_csv(c("Input/observations-524244.csv","Input/observations-524250.csv","Input/observations-524265.csv"))

Special_cases<-load_special_cases("Input/Problem_taxa.csv",csv2 = T)

manual_native_input<-load_native_taxa("Input/Manual native input.csv",csv2 = T)

manual_nonnative_input<-load_nonnative_taxa("Input/Manual non-native input.csv",csv2 = T)

get_wgsrpd3_codes("South Africa")

region_of_interest<-"CPP"

wcvp_distributions <-rWCVPdata::wcvp_distributions
wcvp_names <- rWCVPdata::wcvp_names



## 2a Apply problematic taxa ######################################################

# All species information is in one column and a rank column is added
inat_input_2<-
  inat_input%>%
  mutate(taxon = case_when(taxon_hybrid_name != "" ~ paste(taxon_hybrid_name),
                           taxon_species_name != "" ~ paste(taxon_species_name),
                           taxon_genus_name != "" ~ paste(taxon_genus_name),
                           taxon_family_name != "" ~ paste(taxon_family_name),
                           TRUE ~ NA))%>%
  mutate(rank = case_when(taxon == taxon_genus_name ~ "Genus",
                          taxon == taxon_species_name ~ "species",
                          taxon == taxon_family_name ~ "Family",
                          taxon == taxon_hybrid_name ~ "hybrid",
                          taxon == taxon_subfamily_name ~ "subfamily",
                          TRUE ~ NA))%>%
  filter(!is.na(rank))

# Apply problematic taxa and filtering for research grade
inat_input_3<-
  inat_input_2%>%
  apply_special_cases()%>%
  filter(quality_grade == "research" | special_case == T)


## 2b Standardise names ########################################################

# Summarise inat data into list
inat_list<-
  inat_input_3%>%
  filter(rank != "Family")%>%
  group_by(taxon,special_case)%>%
  summarise(num_obs = length(taxon))%>%
  ungroup()

# Match inat list with WCVP
# note: this can take a couple of minutes (~16min on my machine)
#       make `fuzzy=FALSE` for faster run time
inat_list_match <- 
  wcvp_match_names(inat_list,
                   name_col="taxon",
                   author_col=NULL,
                   fuzzy=TRUE,
                   progress_bar=TRUE)

# Resolve multiple matches and flag poor matches and no distributions
inat_list_resolved<-
  inat_list_match%>%
  mutate(wcvp_id_full = ifelse(!is.na(wcvp_accepted_id),wcvp_accepted_id,wcvp_id),
         distribution_data_available = wcvp_id_full %in% wcvp_distributions$plant_name_id)%>%
  multi_correct()%>%
  filter(multi_correct != "remove")%>%
  mutate(need_manual_check = case_when(special_case ~ FALSE,
                                       multi_correct == "manual_check_needed" ~ TRUE,
                                       (match_similarity != 1.000 | is.na(match_similarity)) & 
                                         (match_type != "Fuzzy (phonetic)" & match_similarity < 0.9) ~ TRUE,
                                       distribution_data_available == F & multi_correct != "remove" ~ TRUE,
                                       TRUE ~ FALSE))

# Write a csv file with all taxa that need to be manually checked
write_csv2(inat_list_resolved%>%
             filter(need_manual_check == TRUE)%>%
             add_column(match_correct = "",
                        nativeness = "",
                        checked = ""),
           "Output/manual_check_needed.csv")


## 2c Obtain native ranges #####################################################

#Obtaining native ranges for all species listed in inat_list_resolved
native_ranges<-
  inat_list_resolved%>%
  left_join(wcvp_names%>%
              select(plant_name_id,family,taxon_name,lifeform_description,climate_description)%>%
              rename(wcvp_id_full = plant_name_id))%>%
  mutate(wcvp_id_full = ifelse(!is.na(wcvp_accepted_id),wcvp_accepted_id,wcvp_id))%>%
  select(wcvp_id_full,taxon,family)%>%
  left_join(wcvp_distributions%>%
              mutate(presence = ifelse(introduced == 1,2,1))%>%
              select(area_code_l3,
                     plant_name_id,
                     presence)%>%
              rename(wcvp_id_full = plant_name_id)%>%
              distinct()%>%
              mutate()%>%
              pivot_wider(names_from = area_code_l3,
                          values_from = presence,
                          values_fill = 0))

# Add nativeness column to the list
inat_list_resolved_2<-
  inat_list_resolved%>%
  filter(need_manual_check == FALSE)%>%
  mutate(nativeness = case_when(special_case ~ "nonnative",
                                wcvp_id_full %in% filter(native_ranges,if_any(all_of(region_of_interest),~. == 1))$wcvp_id_full ~ "native",
                                wcvp_id_full %in% filter(native_ranges,if_all(all_of(region_of_interest),~. == 2))$wcvp_id_full ~ "nonnative",
                                wcvp_id_full %in% filter(native_ranges,if_all(all_of(region_of_interest),~. == 0))$wcvp_id_full ~ "nonnative",
                                TRUE ~ NA))

## 2d Manual Checking ##########################################################


# Read back in the csv if changes were made
manual_check_needed_checked<-
  load_manual_checked_data("Output/manual_check_needed_checked.csv")%>%
  filter(checked == "checked")%>%
  mutate(across(match_type:wcvp_id_full, ~ ifelse(match_correct == F,NA,.)))%>%
  select(-checked,-match_correct,-match_type)

# Bind the manually checked records to the rest of the list
inat_list_resolved_3<-
  inat_list_resolved_2%>%
  bind_rows(manual_check_needed_checked)


## 2e Apply manual native / non-native #########################################

inat_list_resolved_4<-
  inat_list_resolved_3%>%
  mutate(nativeness = case_when(taxon %in% manual_native_input$Name ~ "native",
                                taxon %in% manual_nonnative_input$Name ~ "nonnative",
                                TRUE ~ nativeness))

## Output ######################################################################

nonnative_inventory<-
  inat_list_resolved_4%>%
  filter(nativeness == "nonnative",
         wcvp_rank != "Genus" | special_case == T)

# Export as csv

write_csv2(nonnative_inventory,"Output/Western Cape urban non-native plant inventory 1 February 2025.csv")
write_csv2(inat_input_3%>%
             filter(taxon%in%nonnative_inventory$taxon),"Output/Western Cape urban non-native records 1 February 2025.csv")


# ------------------------------ the end ---------------------------------------


## Show native ranges of selected taxon ##

r<-wcvp_distribution(taxon = "Quercus nigra",
                     introduced = F)
wcvp_distribution_map(r)
#view(r)

