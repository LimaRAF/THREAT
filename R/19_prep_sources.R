#### PREPARING THE APPENDICES ####

## Occurrence data
path <- "C:/Users/renat/Documents/raflima/Pos Doc/Manuscritos/Artigo AF checklist/data analysis/AppendixA_records.csv" 
herbs <- read.csv(path)
herbs <- herbs[order(herbs$citation), ]
herbs1 <- herbs[!herbs$citation %in% "GBIF.org (11 December 2019) GBIF Occurrence Download https://doi.org/10.15468/dl.mzmat2",]

citacao <- herbs1$citation
# JABOT exclusive
replace_these <-  grepl("^JBRJ ", citacao, perl = TRUE)  
citacao[replace_these] <-
  paste(herbs1$collectionName[replace_these], citacao[replace_these], sep= ". ")

# speciesLink exclusive
replace_these <- !is.na(herbs1$GBIF.datasetKey) & 
                  !grepl("accessed via GBIF", citacao) &
                    !grepl("[0-9][0-9][0-9][0-9]\\)", citacao, perl = TRUE) &
                      !grepl("Occurrence dataset", citacao, perl = TRUE)
                        !is.na(herbs1$collectionName)
citacao[replace_these] <-
  paste0(herbs1$collectionName[replace_these], 
        ". speciesLink network. INCT - Herbário Virtual da Flora e dos Fungos. Available at: http://specieslink.net/search. Accessed in ",
        herbs1$publication.or.download.date[replace_these],".")

# Removing duplicates and saving
citacao1 <- citacao[order(citacao)]
#citacao1 <- gsub(" accessed via GBIF.org on 2019-12-11.",".", citacao1)
citacao1 <- gsub(" $", ".", gsub("<a href=", "", citacao1))

#Saving
writeLines(citacao1, "text/SupplementaryText1.txt")


# -------------------------------------------------------------------------
## Occurrence data

# forest density
path <- "C:/Users/renat/Documents/raflima/Pos Doc/Manuscritos/Artigo Hyperdominance/selected_sites_DA.csv" 
DA <- read.csv(path, as.is=TRUE) # all types of forests inside and around the AF
#Getting SiteCode list for all AF sites and those inside the 100 km  buffer, with more common dbh cutoffs and far IFFSC samples
ids = DA$select %in% c("TRUE", "TRUE_buff50km","TRUE_undeclared") | (DA$select %in% "TRUE_buff100km" & !grepl("Caatinga|Cerrado", DA$forest_type))
DA1 = DA[ids,]
DA1 = DA1[DA1$dbh %in% c("DBH>=4.8-5.0cm","DBH>=4.4-4.5cm"),]
DA1 = DA1[DA1$DA>450,]
DA1 = DA1[DA1$DA<3400,]

# species abundances
path1 <- "C:/Users/renat/Documents/raflima/Pos Doc/Manuscritos/Artigo Hyperdominance/plot_metadata.csv" 
sites <- read.csv(path1, as.is=TRUE, na.strings = c(""," ",NA))

refIDs <- sort(unique(DA1$refID, sites$refID))

# Getting TreeCo references
path2 <- "C:/Users/renat/Documents/raflima/Pos Doc/Databases/References/references.xlsx" 
refs <- readxl::read_xlsx(path2)
refs1 <- refs[refs$refID %in% refIDs,]
refs2 <- refs1[!duplicated(refs1$refID),]

# Removing duplicates and saving
citacoes <- paste0(refs2$Reference,
                   " (refID= ", 
                   refs2$refID,")")
citacoes <- citacoes[order(citacoes)]
citacoes <- citacoes[!duplicated(citacoes)]
#citacao1 <- gsub(" accessed via GBIF.org on 2019-12-11.",".", citacao1)
#citacoes1 <- gsub(" $", ".", gsub("<a href=", "", citacoes1))

#Saving
writeLines(citacoes, "text/SupplementaryText2.txt")
