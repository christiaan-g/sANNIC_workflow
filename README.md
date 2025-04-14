# sANNIC workflow
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15210704.svg)](https://doi.org/10.5281/zenodo.15210704)

This repository includes the code and data related to the execution of the semi-Automated Non-Native Inventory Compilation (sANNIC) workflow. This workflow is designed to facilitate the compilation of non-native vascular plant inventories from iNaturalist data. It follows four steps: 1) acquire data from iNaturalist, a community science platform which allows users to contribute records of organisms, 2) compile a non-native inventory referencing the World Checklist of Vascular Plants (WCVP), 3) create and curate an iNaturalist project of the non-native species inventory, and 4) conduct a completeness assessment of the inventory generated.

This workflow was demonstrated using the non-native plants of Western Cape urban areas case study, South Africa https://www.inaturalist.org/projects/non-native-plants-of-western-cape-urban-areas.

## Visual
![Updated workflow 2025](https://github.com/user-attachments/assets/67471b4a-2916-4e58-b870-045ed3c3003b)

## Instructions

### Step 1: Acquire data

Export records for any vascular plant group of interest in any region of interest from iNaturalist. The taxonomic group of interest and region of interest can be specified in the iNaturalist search filters. For region of interest, either an existing 'place' or a manually specified 'place' can be used https://www.inaturalist.org/places.

Columns required for this workflow are:

-   id 

-   observed_on

-   quality_grade

-   captive_cultivated

-   taxon_id

-   taxon_family_name

-   taxon_subfamily_name

-   taxon_genus_name

-   taxon_species_name

-   taxon_hybrid_name

-   taxon_subspecies_name

### Step 2: Run ANNIC R script

Run the "Step2 ANNIC R script.R" to produce a non-native inventory. It requires the R packages tidyverse, rWCVP, and rWCVPdata.

#### Input

Input an iNaturalist csv (comma separated value) export with at least the columns specified in Step 1. Any additional columns such as latitude and longitude may be included. Optional (but recommended) inputs include a spreadsheet of problematic taxa (2a) and spreadsheets for manual native/non-native species (2e). See the Metadata.xlsx and the case study example on the required columns. These can be inported as csv or csv2 (semi-colon separated values). When importing, specify the format using the csv2 argument in the 'load_special_cases', 'load_native_taxa', and 'load_nonnative_taxa' functions csv2 is the default.

#### 2a) Apply problematic taxa (optional)

Allows for taxa that are difficult to identify to be collapsed to a higher taxonomic rank. Requires a list of problematic taxa as an input, including the taxonomic name they will be collapsed to (as it appears in iNaturalist, e.g. _Avena_) and the rank of that taxon (e.g. genus). It also allows for the flexible manual specification of taxa to be collapsed based on their taxon_id specified in iNaturalist. taxon_id for any taxon can be obtained from the url of the taxon page e.g. https://www.inaturalist.org/taxa/53178-Plantago-lanceolata.

#### 2b) Standardise names

This section standardises the taxonomic names of the vascular plants using the World Checklist of Vascular Plants (WCVP). All names not considered problematic taxa will be collapsed to species level (i.e. subspecific names will be ignored). Names will be matched to the WCVP using the rWCVP and rWCVPdata packages (Brown et al. 2023, Govaerts et al. 2021). By default fuzzy matching is enabled. Multiple matches are resolved using the 'multi_correct' function.

#### 2c) Obtain native ranges

Native ranges of all taxa are acquired from the matched names in the WCVP. Taxa that have their native ranges entirely outside the region of interest or are considered introduced in the region of interest will be flagged as "nonnative", and taxa with native ranges within the region of interest will be flagged as "native".

#### 2d) Manual checking (optional)

Involves manual checking of taxa which had poor matches during step 2b or for which native range information was not available in step 2c. Edit the output "Output/manual_check_needed.csv" and fill in the columns *match_correct* (TRUE or FALSE whether the match suggestion is correct), *nativeness* ("native" or "nonnative"), and *checked* ("checked" or "remove" to confirm whether the row should be kept or discarded). After edits, save file as "Output/manual_check_needed_checked.csv" and continue with the workflow. If this this step is ignored the taxa with manual checks needed will be ignored and not included in the final inventory.

#### 2e) Apply manual check native/nonnative (optional)

Allows for the manual specification of a species as native or non-native using input files. This is useful to correct any errors in the WCVP. These files are imported using the "load_native_taxa" and "load_nonnative_taxa" functions. If an input file is not found the script will assume that no native/non-native taxa are considered for manual checking and the automated flagging is accepted.

#### Output

Will output two spreadsheets (csv2 by default). A list of non-native species and associated information, and a dataset of all the observations of those non-native species.

### Step 3: Curate iNaturalist project

The next step is bringing the compiled list back into iNaturalist. Create an iNaturalist project. Through the project settings, manually add all the species from the non-native list and specify the place(s) of interest. This creates a focussed collection of all the potential non-native plants in the region of interest. If there are too many species they can be added at the genus level (if all species in the genus are non-native to the region of interest) or can be added to multiple parallel projects (for example grouped alphabetically by family).

Creating the project is just the first part of step 3. The crucial part is to then go through all observations within that project to check their accuracy. This involves verifying identifications, identifying records as accurately as possible, and making sure species are correctly flagged as cultivated if needed. This step allows for the improvement of data quality and better understanding of the non-native flora in the region of interest. The level of curation required depends on the taxon and region of interest. Refer to Richardson and Potgieter (2024) for details.

### Step 4: Conduct completeness assessment

Conduct a completeness assessment of the data to assess how well the taxa of interest are sampled within the region of interest. This is all about understanding how comprehensive the dataset is. iNaturalist data is opportunistic so sampling intensity is not even across the spatial scales. The concept of a sample completeness profile helps us estimate what proportion of the total non-native plant diversity we actually managed to observe. It takes into account not just the overall number of species but also how common or rare they are in our dataset. We used the *Completeness* function in the iNEXT.4steps package to calculate the completeness profile in our case study. Uncertainty is obtained using bootstrapping. Refer to Chao et al. (2020) for more information.

In our urban Western Cape case study we calculated the profile for the dataset in its entirety and for each urban centre separately. For each species we also calculated a Prevalence Index based on the number of urban centres where each species was found and the abundance ranking of those species within each urban centre, weighted by the Sample Coverage calculated for that urban centre. All code we used for this analysis can be found at "Output/Step4 Completeness assessment.R".

## References

Brown MJM, Walker BE, Black N, Govaerts R, Ondo I, Turner R, Nic Lughadha E (2023) rWCVP: A companion R package to the World Checklist of Vascular Plants. New Phytologist 240: 1355-1365. https://doi.org/10.1111/nph.18919

Chao A, Kubota Y, Zelený D, Chiu CH, Li CF, Kusumoto B, Yasuhara M, Thorn S, Wei CL, Costello MJ, Colwell RK (2020) Quantifying sample completeness and comparing diversities among assemblages. Ecological Research 35: 292–314. https://doi.org/10.1111/1440-1703.12102

Govaerts R, Nic Lughadha E, Black N, Turner R, Paton A (2021) The World Checklist of Vascular Plants, a continuously updated resource for exploring global plant diversity. Scientific Data 8: 215. https://doi.org/10.1038/s41597-021-00997-6

Richardson DM, Potgieter LJ (2024) A living inventory of planted trees in South Africa derived from iNaturalist. South African Journal of Botany 173: 365–379. https://doi.org/10.1016/j.sajb.2024.08.006
