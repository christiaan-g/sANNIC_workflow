rm(list = ls())
setwd("C:/Users/chris/OneDrive - Stellenbosch University/Documents/sANNIC_workflow")
# rstudioapi::writeRStudioPreference("data_viewer_max_columns", 100L)
set.seed(1234)

## Unpack packages ############################################################

library(tidyverse)
library(iNEXT)
library(sf)
library(units)
library(iNEXT.4steps)


## Import datasets ###############################################################

# Importing Western Cape urban shapefiles
WC_urban_shapefiles<-
  st_read("Input/WC_urban_areas.kml")%>%
  rename(Town_name = Name)

# Importing the non-native plant inventory
nonnative_inventory<-
  read_csv2("Output/Western Cape urban non-native plant inventory 1 February 2025.csv")

# Importing records of non-native plants
nonnative_records<-
  read_csv2("Output/Western Cape urban non-native records 1 February 2025.csv")


## Preparing data ###########################################################

# Joining non-native plant records with urban centre shapefiles
WC_records_points<-
  st_as_sf(nonnative_records, coords = c("longitude","latitude"),crs = 4326)%>%
  st_transform(st_crs(WC_urban_shapefiles))%>%
  st_join(WC_urban_shapefiles)

# Dropping geometry, main dataset
WC_alien_plants<-
  WC_records_points%>%
  st_drop_geometry()%>%
  filter(!is.na(Town_name))

# Calculating the urban area (m^2) of each urban centre
Urban_areas<-
  WC_urban_shapefiles%>%
  mutate(area = st_area(.))%>%
  st_drop_geometry()

## Creating matrices for urban stats ###########################################

# Converting main dataset into matrix format with one columns, data coarsening applied
total_alien_matrix_coarse<-
  WC_alien_plants%>%
  group_by(taxon,year(observed_on),month(observed_on),day(observed_on),user_id)%>%
  summarise(obs = length(taxon))%>%
  mutate(num_obs = 1)%>%
  group_by(taxon)%>%
  summarise(num_obs = length(taxon))%>%
  column_to_rownames(var = "taxon")%>%
  as.matrix()

# Converting main dataset into a site-by-species matrix, data coarsening applied
site_by_species_coarse<-
  WC_alien_plants%>%
  group_by(Town_name,taxon,user_id,year(observed_on),month(observed_on),day(observed_on))%>%
  summarise(obs = length(taxon))%>%
  mutate(obs = 1)%>%
  group_by(Town_name,taxon)%>%
  summarise(num_obs =length(taxon))%>%
  ungroup()%>%
  pivot_wider(names_from = Town_name,
              values_from = num_obs)%>%
  select(where(~ !all(is.na(.))))%>%
  column_to_rownames(var = "taxon")%>%
  mutate(across(everything(), ~ replace_na(., 0)))%>%
  as.matrix()%>%t()

# Transposing site-by-species into species-by-site
species_by_site_coarse<-
  t(site_by_species_coarse)


## iNEXT 4 steps ##############################################################

# Calculating a sample completeness profile of the total urban areas
total_completeness_coarse<-
  Completeness(
    total_alien_matrix_coarse,
    q = seq(0, 2, 1),
    datatype = "abundance",
    nboot = 100,
    conf = 0.95,
    nT = NULL
  )
total_completeness_coarse

# Calculating sample completeness profiles for each urban area separately
completeness_coarse<-
  Completeness(
    species_by_site_coarse,
    #species_by_site,
    q = seq(0, 2, 1),
    datatype = "abundance",
    nboot = 100,
    conf = 0.95,
    nT = NULL
  )

# Removing sites with falling profile and isolating q = 1
urban_completeness<-
  completeness_coarse%>%
  group_by(Assemblage)%>%
  filter(all(diff(Estimate.SC)>0))%>%
  filter(Order.q == 1)%>%
  select(Assemblage,Estimate.SC,SC.LCL,SC.UCL
         )%>%
  rename(Site = Assemblage)


# Produces summary table of each urban centre,including urban area (m^2), 
#   the number of observations, number of species, 
#   and sample coverage (q = 1) estimates
Urban_summary_table<-
  Urban_areas%>%
  rename(Site = Town_name)%>%
  left_join(WC_records_points%>%
              group_by(Town_name)%>%
              summarise(obs = length(taxon),
                        spp = length(unique(taxon)))%>%
              st_drop_geometry()%>%
              rename(Site = Town_name))%>%
  left_join(urban_completeness)%>%
  drop_units()

view(Urban_summary_table)

## Prevalence index ############################################################

# Produces a site-by-species data frame with species ranks
a<-
  site_by_species_coarse%>%
  as.data.frame()%>%
  rownames_to_column("Site")%>%
  left_join(select(urban_completeness,Site,Estimate.SC))%>%
  pivot_longer(names_to = "Species",cols = -Site)%>%
  group_by(Site)%>%
  dplyr::mutate(rank = ifelse(Species == "Estimate.SC",value,rank(-value)))%>%
  dplyr::mutate(rank = ifelse(value == 0 & Species != "Estimate.SC",Inf,rank))%>%
  select(-value)%>%
  pivot_wider(names_from = Site,values_from = rank)

# Produces a site-by-species presence-absence data frame
b<-site_by_species_coarse%>%
  as.data.frame()%>%
  rownames_to_column("Site")%>%
  pivot_longer(names_to = "Species",cols = -Site)%>%
  group_by(Site)%>%
  dplyr::mutate(pres = ifelse(value == 0,0,1))%>%
  select(-value)%>%
  pivot_wider(names_from = Site,values_from = pres)


# Calculates Prevalence Indices of all species (see Gildenhuys et al. (unpubl.) for details) 
d1<-
  a%>%
  mutate(across(everything(),~replace_na(.x,0)))%>%
  rowwise()%>%
  column_to_rownames("Species")%>%
  mutate(across(everything(), ~ last(.x)/ sqrt(.x) ))%>%
  rownames_to_column("Species")%>%
  filter(Species != "Estimate.SC")%>%
  rowwise()%>%
  mutate(PI = sum(c_across(-Species)))%>%
  select(Species,PI)

# Calculates number of urban centres where each species is found
d2<-
  b%>%
  rowwise()%>%
  mutate(num_areas = sum((c_across(-Species))))%>%
  select(Species, num_areas)

# Calculates the total number of observations of each species
d3<-
  nonnative_inventory%>%
  select(taxon,num_obs)%>%
  rename(Species = taxon,
         total_obs = num_obs)%>%view()

# Combines all the previous calculations into a species summary table
e<-full_join(d1,d2)%>%
  full_join(d3)%>%
  dplyr::arrange(desc(PI))

view(e)


## Prevalence index and urban summary output ###################################

write_csv2(Urban_summary_table,"Output/Urban summary table.csv")
write_csv2(e,"Output/Nonnative inventory prevalence index scores.csv")

## Completeness profile visuals ################################################

# Produces Appendix C
AppendixC<-
  completeness_coarse%>%
  mutate(label = ifelse(is.na(Estimate.SC)|Estimate.SC == 1|
                          Estimate.SC == -Inf|Estimate.SC == "NaN", "NA", NA),
         Estimate.SC = ifelse(Estimate.SC == 1,NA,Estimate.SC))%>%
  mutate(SC.UCL = ifelse(is.na(Estimate.SC),NA,SC.UCL),
         SC.LCL = ifelse(is.na(Estimate.SC),NA,SC.LCL))%>%
  ggplot(aes(Order.q,Estimate.SC))+
  geom_line(linewidth = 0.3,colour = "black")+
  geom_errorbar(aes(ymin = SC.LCL,ymax = SC.UCL),colour = "blue",linewidth = 0.3)+
  facet_wrap(~Assemblage)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 9))+
  ylab("Estimate SC")+
  xlab("Order q")+
  geom_text(aes(label = label),size = 2,y = 0.1)

tiff('Output/AppendixC.tiff', units="cm", width=25, height=16, res=400)
AppendixC
dev.off()


# Produces figure 4
figure4<-
completeness_coarse%>%
  filter(Assemblage == "Cape Town"|Assemblage =="Paarl & Wellington"|
           Assemblage =="Yzterfontein"|Assemblage =="Beaufort West")%>%
  mutate(Assemblage = fct_relevel(Assemblage,c("Cape Town" ,"Paarl & Wellington",
                                               "Yzterfontein" ,"Beaufort West")))%>%
  mutate(Assemblage = fct_recode(Assemblage,"a)  Cape Town" = "Cape Town" ,
                                 "b)  Paarl & Wellington" = "Paarl & Wellington",
                                 "c)  Yzerfontein" = "Yzterfontein",
                                 "d)  Beaufort West" = "Beaufort West"))%>%
  ggplot(aes(Order.q,Estimate.SC))+
  geom_line(colour = "black",linewidth = 0.4)+
  geom_errorbar(aes(ymin = SC.LCL,ymax = SC.UCL),colour = "blue",linewidth = 0.4)+
  facet_wrap(~Assemblage)+
  theme_bw()+
  theme(#panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        text = element_text(size=12, colour = "black"))+
  ylab("Estimate SC")+
  xlab("Order q")

tiff('Output/figure4.tiff', units="cm", width=13, height=10, res=300)
figure4
dev.off()

#----------- the end -------------