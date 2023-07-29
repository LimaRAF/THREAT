#########################################H
#### PREPARING THE SUPPLEMENTARY DATA ####
#########################################H
rm(list = ls())

# READING NECESSARY FILES ------------------------------------------------------

## Taxonomy
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 3, 4)]
full.tax <- readRDS("data/threat_final_taxonomy.rds")[, c(1:8)]
full.tax$species.correct2 <- full.tax$species
full.info <- dplyr::left_join(full.tax, tax)

## Synonyms
synms <- read.csv("data/sis_connect/synonyms_threat.csv", encoding = "UTF-8")
synms <- aggregate(synms$name, list(synms$internal_taxon_id), paste0, collapse = "|")
names(synms) <- c("internal_taxon_id", "Synonyms")
full.info <- dplyr::left_join(full.info, synms)

## Common names
common <- read.csv("data/sis_connect/commonnames_threat.csv", encoding = "UTF-8")
common <- aggregate(common$name, list(common$internal_taxon_id), paste0, collapse = "|")
names(common) <- c("internal_taxon_id", "CommonNames")
full.info <- dplyr::left_join(full.info, common)

## Countries
countries0 <- readRDS("data/countries_per_species.rds")
countries0 <- countries0[!is.na(countries0$iso2), ]
countries <- aggregate(countries0$iso2, list(countries0$species.correct2), 
                       function(x) paste0(unique(x), collapse = "|"))
names(countries) <- c("species.correct2", "CountryCode")
full.info <- dplyr::left_join(full.info, countries)

## Habitat, growth form and other species info
hab1 <- read.csv("data/threat_habitats_preliminar1.csv", encoding = "UTF-8")
hab2 <- read.csv("data/threat_habitats.csv", encoding = "UTF-8")
hab <- dplyr::left_join(hab2, hab1)
hab <- hab[,c(1,2,5,6,13,20,32)]
hab$GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName <- 
  gsub("\\|NA*", "", hab$GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName)
full.info <- dplyr::left_join(full.info, hab)

## Species uses
usos <- read.csv("data/sis_connect/usetrade_threat.csv", encoding = "UTF-8")
usos <- aggregate(usos$UTEndUse.UTEndUseSubfield.UTEndUseName, list(usos$internal_taxon_id), 
                       function(x) paste0(unique(x), collapse = "|"))
names(usos) <- c("internal_taxon_id", "UsesIUCN")
full.info <- dplyr::left_join(full.info, usos)

## Occurrence data 
oc.data <- readRDS("data/threat_occ_data_final.rds")
sources <- as.data.frame.matrix(table(oc.data$tax, oc.data$source))
sources$species.correct2 <- row.names(sources)

## Level of confidence in Taxonomic identifications
tax.level <- as.data.frame.matrix(table(oc.data$tax, oc.data$tax.check.final))
tax.level$medium <- round(tax.level$medium/2, 0) 
tax.level$TaxConfLevel <- round(100 * (apply(tax.level[,c("high", "medium")], 1, sum)/
                            apply(tax.level, 1, sum)), 2)
tax.level$species.correct2 <- row.names(tax.level)
occ.data <- dplyr::left_join(sources, tax.level[,c("TaxConfLevel", "species.correct2")])
full.info <- dplyr::left_join(full.info, occ.data)


## COMBINED INFO FOR ALL SPECIES ##
all.spp <- readRDS("data/all_spp_crits_tax_cites.rds")
# all.spp1 <- all.spp[c(1:9, 13:24, 37:42, 44:50, 52:58, 62, 71, 73, 77:80, 88, 89, 91)]
colunas0 <- c("internal_taxon_id","species",
              "assessment.period", "reduction_A12", "reduction_A12_obs", "basis_d",
              "A2", "category_A", "category_A_code",
              "reduction_A12.10ys", "reduction_A12.15ys", "reduction_A12.20ys",
              "reduction_A12.25ys", "reduction_A12.30ys", "reduction_A12.35ys",
              "reduction_A12.40ys", "reduction_A12.45ys", "reduction_A12.50ys",
              "EOO", "AOO", "Nbe_occs", "Nbe_subPop", "nbe_loc_total", "protected", "declineB", 
              "sever.frag", "B1", "B2", "category_B", "category_B_code", "any.decline", "cont.decline", 
              #"C1", 
              "reduction_3gen", "reduction_2gen", "reduction_1gen",
              "C2", "category_C", "category_C_code", "pop.size", "pop.size.low", "D", "D.AOO", 
              "D2.Loc", "category_D", "category_D_code", "category", "main.criteria", "aux.criteria",
              "last.record", "prop.coly.na", "prop.high.tax", "only.type", "endemic", 
              "prop.EOO.in.StrictUCs", "prop.EOO.forest", "status.reflora", "redlistCategory", 
              "redlistCriteria", "yearPublished", "category.regional", "cat.reg.clean", "appendix")
all.spp1 <- all.spp[, c(colunas0)]
full.info <- dplyr::left_join(full.info, all.spp1)

## Other IUCN related info
#Assessments
assess <- read.csv("data/sis_connect/assessments_threat.csv", encoding = "UTF-8")
cols <- c("internal_taxon_id", "BiogeographicRealm.realm", "PopulationTrend.value")
assess1 <- assess[, cols]
full.info <- dplyr::left_join(full.info, assess1)

#All fields
all <- read.csv("data/sis_connect/allfields_threat.csv", encoding = "UTF-8")
cols <- c("internal_taxon_id", "HabitatContinuingDecline.isDeclining", 
          "PopulationContinuingDecline.isDeclining",
          # "PopulationDeclineGenerations1.range", 
          # "PopulationDeclineGenerations2.range", 
          # "PopulationDeclineGenerations3.range", 
          "PopulationReductionPast.direction",
          "GenerationLength.range",
          "PopulationReductionPastBasis.value", "SevereFragmentation.isFragmented")
all1 <- all[, cols]
full.info <- dplyr::left_join(full.info, all1)


# FILTERING, RENAMING AND SAVING -----------------------------------------------
colunas <- readxl::read_xlsx("./text/renameSuppData.xlsx")
rpl.cols <- colunas$new
names(rpl.cols) <- colunas$original
rpl.cols1 <- rpl.cols[colunas$final.order > 0]

#Filtering
full.info1 <- full.info[ , names(full.info) %in% names(rpl.cols1)]

#Replacing names
rpl.cols2 <- rpl.cols1[order(nchar(names(rpl.cols1)), decreasing = TRUE)]
names(full.info1) <- stringr::str_replace_all(names(full.info1), rpl.cols2)  

#Reordering
tmp <- colunas[colunas$final.order > 0, ]
tmp$new[!tmp$new %in% names(full.info1)]
names(full.info1)[!names(full.info1) %in% tmp$new]
final.order <- tmp$final.order[match(tmp$new, names(full.info1), nomatch = 0)]
full.info2 <- full.info1[, order(final.order)]

#Final edits
full.info3 <- full.info2[ , -which(names(full.info2) %in% "Genus")] 
full.info3$GeneralHabitatsIUCN <- 
  gsub("Subtropical/Tropical", "(Sub)Tropical", full.info3$GeneralHabitatsIUCN, fixed = TRUE)

full.info3$PropOccsInProtectedAreas <- round(full.info3$PropOccsInProtectedAreas, 1)
full.info3$PropEOOInProtectedAreas <- round(full.info3$PropEOOInProtectedAreas, 1)
full.info3$PropEOOIsHabitat <- round(full.info3$PropEOOIsHabitat, 1)

full.info3$HabitatContinuingDecline[full.info3$HabitatContinuingDecline %in% 1] <-
  "Decreasing"
full.info3$HabitatContinuingDecline[full.info3$HabitatContinuingDecline %in% 0] <-
  "Increasing"
full.info3$HabitatContinuingDecline[full.info3$HabitatContinuingDecline %in% "U"] <-
  "Unknown"

full.info3$PopulationContinuingDecline[full.info3$PopulationContinuingDecline %in% 1] <-
  "Decreasing"
full.info3$PopulationContinuingDecline[full.info3$PopulationContinuingDecline %in% 0] <-
  "Increasing"
full.info3$PopulationContinuingDecline[full.info3$PopulationContinuingDecline %in% "U"] <-
  "Unknown"

names(full.info3)[names(full.info3) %in% c("EOO", "AOO")] <- 
  paste0(names(full.info3)[names(full.info3) %in% c("EOO", "AOO")], "_km2")

full.info3$EcologicalGroup <- stringr::str_to_title(gsub("_", " ", full.info3$EcologicalGroup)) 


## Not sure why ordering is not working
colunas1 <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Species", "SpeciesAuthority",
              "Synonyms", "CommonNames", "CountryCodesISO2", "BiogeographicRealm",
              "GeneralHabitatsIUCN", "GenerationLength", "GrowthFormIUCN",
              "EcologicalGroup", "MaximumHeight", "UsesIUCN", 
              "HerbariumOccurrences", "ExtraInventoryOccurrences", 
              "NonDuplicatedOccurrences", "PropTaxConfidenceHigh", 
              "AssessmentPeriod", 
              #"AssessA1", 
              "AssessA2", "PopulationTrend", 
              "HabitatContinuingDecline", "PopulationContinuingDecline", 
              #"PopulationReductionBasis",
              "PopReduction",
              "Category_A", "Category_A_code", "Exploitation_basis_d", "PopReductionObservations",
              "PopReduction10ysGL", "PopReduction15ysGL", "PopReduction20ysGL",
              "PopReduction25ysGL", "PopReduction30ysGL", "PopReduction35ysGL", 
              "PopReduction40ysGL", "PopReduction45ysGL", "PopReduction50ysGL",
              "EOO_km2", "AOO_km2", "NumberSubPop", 
              "NumberLocations", "SevereFragmentation", "AssessB1", "AssessB2",
              "CategoryB", "CategoryB_code", 
              "PopulationDecline1Generation", "PopulationDecline2Generations", 
              "PopulationDecline3Generations", "PopulationSize", "PopulationSizeLowEstimate",
              #"AssessC1", 
              "AssessC2", "CategoryC", "CategoryC_code", "CategoryD", 
              "CategoryD_code", "IUCN_CategoryOriginal", "IUCN_CategoryRegional",
              "IUCN_CategoryFinal", "MainCriteria", "AuxiliaryCriteria",
              "YearLastSeen", "OnlyKnownFromType", "EndemismCategory", 
              "PropOccsInProtectedAreas", "PropEOOInProtectedAreas", "PropEOOIsHabitat",
              "PreviousCNCFloraCategory", "PreviousRedListCategory", 
              "PreviousRedListCriteria", "CITES_appendix")
colunas1[!colunas1 %in% names(full.info3)]
full.info4 <- full.info3[ , colunas1]

## Replacing NA by empty characters
for (i in seq_len(dim(full.info4)[2])) {
  check_these <- full.info4[,i] %in% c("", NA, NA_character_)
  if (any(check_these))
    full.info4[check_these, i] <- ""
}

##Saving
writexl::write_xlsx(full.info4, "./text/Science/DataS1_new.xlsx",
                    format_headers = FALSE)
