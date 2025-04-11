# sANNIC workflow

This repository includes the code and data related to the execution of the semi-Automated Non-Native Inventory Compilation (sANNIC) workflow. This workflow is specifically designed to facilitate the compilation of non-native vascular plant inventories from iNaturalist data. It follows four steps: 1) acquire data from iNaturalist, a community science platform which allows users to contribute records of organisms, 2) compile a non-native inventory referencing the World Checklist of Vascular Plants (WCVP), 3) create and curate an iNaturalist project of the non-native species inventory, and 4) conduct a completeness assessment of the inventory generated.

## Visual

## Instructions

### Step 1: Acquire data

Export records for any vascular plant group of interest in any location of interest from iNaturalist.

Columns required for this workflow are:

-   id observed_on

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

Run the "Step2 ANNIC R script.R" to produce a non-native inventory. It requires the R packages, tidyverse, rWCVP, and rWCVPdata.

#### Input

Inputs an iNaturalist csv (comma separated value) export with at least the columns specified in Step 1. Any additional columns such as latitude and longitude may be included. Optional (but recommended) inputs include a spreadsheet of problematic taxa (2a) and spreadsheets for manual native/non-native species (2e). These can be uploaded as csv or csv2 (semi-colon separated values). csv2 is the default.

#### 2a) Apply problematic taxa (optional)

Allows for taxa that are difficult to identify to be collapsed to a higher (or lower) taxonomic rank. Requires a list of problematic taxa as an input, including the taxonomic name they will be collapsed to (as it appears in iNaturalist) and the rank of that taxon (e.g. name: Avena; rank: genus). It also allows for the flexible manual specification of taxa to be collapsed based in their taxon_id specified in iNaturalist (this can be useful for subgeneric groups or other non-standard groupings not accounted for in iNaturalist).

#### 2b) Standardise names

This section standardises the taxonomic names of the vascular plants using the World Checklist of Vascular Plants (WCVP). All names not considered problematic taxa will be collapsed to species level (i.e. subspecific names will be ignored).

#### 2c) Obtain native ranges

Native ranges of all taxa are acquired from the matched names in the WCVP. Taxa that have their native ranges entirely outside the region of interest or are considered introduced in the region of interest will be flagged as "nonnative" and taxa with native ranges within the region of interest will be flagged as "native".

#### 2d) Manual checking (optional)

Involves manual checking of taxa which had poor matches during step 2b or for which native range information was not available in step 2c. Edit the output "Output/manual_check_needed.csv" and fill in the columns *match_correct* (TRUE or FALSE whether the match suggestion is correct), *nativeness* ("native" or "nonnative"), and *checked* ("checked" or "remove" to confirm whether the row should be kept or discarded). After edits, safe file as "Output/manual_check_needed_checked.csv". If this this step is ignored the taxa with manual checks needed will be ignored and not included in the final inventory.

#### 2e) Apply manual check native/nonnative (optional)

Allows for the manual specification of a species as native or non-native using input files. This is useful to correct any errors in the WCVP. If no input file is found the script will assume that no native or non-native taxa are considered for manual checking and the automated flagging is accepted.

#### Output

Will output two spreadsheets (csv2 by default). A list of non-native species and associated information, and a dataset of all the observations of those non-native species.

### Step 3: Curate iNaturalist project

The next step is bringing the compiled list back into iNaturalist. Create an iNaturalist project with the region of interest. Through the project settings, manually add all the species from the non-native list, creating a focussed collection of all the potential non-native plants in the region of interest. If there are too many species they can be added at the genus level (if all species in the genus are non-native to the region of interest) or can be added to multiple parallel projects (for example grouped alphabetically).

Creating the project is just the first part of step 3. The crucial part is to then go through all of the observations within that project and check their accuracy. This involves verifying identifications, identifying records as accurately as possible, and making sure species are correctly flagged as cultivated if needed. This step allows for the improvement of data quality and better understanding of the non-native flora in the region of interest. The intensity of curation required depends on the taxon and region of interest. Refer to Gildenhuys et al. (unpubl.) for details.

### Step 4: Conduct completeness assessment

Conduct a completeness assessment of the data to assess how well sampled the taxa of interest are in the region of interest. This is all about understanding how comprehensive the dataset is. iNaturalist data is opportunistic so sampling intensity is not even across spatial scales. The concept of a sample completeness profile helps us estimate what proportion of the total non-native plant diversity we actually managed to observe. It takes into account not just the overall number of species but also how common or rare they are in our dataset. We used the *Completeness* function in the iNEXT.4steps package to calculate the completeness profile in our case study. Uncertainty is obtained using bootstrapping.

In our urban Western Cape case study we calculated the profile for the dataset in its entirety and for each urban centre separately. We then calculated a Prevalence Index of all species in our dataset. All code we used for this analysis can be found at "Output/Step4 Completeness assessment.R".
