#########################################################H
#########################################################H
#### PREPARING THE PREVIOUS IUCN RED LIST ASSESSMENTS ####
#########################################################H
#########################################################H

# The data you have downloaded falls under the IUCN Red List Terms and
# Conditions of Use which can be viewed on the IUCN Red List website at
# www.iucnredlist.org/terms/terms-of-use.

# Depending on the type of non-spatial download you requested, the zip file that
# you receive will contain multiple data tables in Comma Delimited, CSV text
# file format. They are best viewed within applications to create
# spreadsheets (e.g. MS-Excel) and databases (e.g. MS-Access).

# To view multiple csv tables properly, they should be imported into a database
# application and the different tables joined using the "internalTaxonId" when
# viewing taxon information, or the "assessmentID" field to see assessment
# information, especially in cases where downloads include multiple assesments
# for a taxon (e.g. regional and global assessments).

# IMPORTANT ISSUES TO NOTE:
# 1- When importing the csv files you need to specifiy that the files are in UTF-8
#    format otherwise the extended characters will not import and display correctly.

# 2 - Many fields (especially text fields) will contain html coding to render the
#    text properly on The IUCN Red List website; it is imposible to strip this out
#    of the exports. eliminate this by pasting into an intermediate text editor.

# CONTENTS OF THE ZIP FILES:
  
# simple_summary.csv [a flat table containing the taxonomy for each taxon, Red
# List Category, Red List Criteria, version of criteria used, and current
# population trend].

# additional csv tables [these contain much more information as specified by the
# type of download selected and by the settings selected under your user
# profile. These tables need to be imported into an appropriate Application (see
# below) and joined using the "internalTaxonId" and/or "assessmentId" fields
# (see above)].


# TO IMPORT A TABLE INTO USE:
# IMPORT OPTION 		VALUE
# Data Type 		Comma Delimited
# File Origin 		Unicode (UTF-8)
# Text Qualifier 		{"}
# Field Delimiter 	{,}
# Multi-value Delimiter	{|}
# First Row 		Field Names

rm(list=ls())
gc()

#### LOADING PACKAGES  ####
# require(data.table)


#### LOADING THE FILES  ####
path <- "data/data-raw/redlist_species_data_27cbaeb0-4d95-4128-a993-a11673f1e581" # folder containing all files from version 2022-2
files <- list.files(path, pattern = "csv", full.names = TRUE)
files <- files[!grepl("dois.csv", files)]

#### GETTING ONLY THE TAXON IDS IN THREAT ####
iucn2 <- readRDS("data/threat_iucn_taxonomy.rds")
taxon.id <- "internalTaxonId"
colar <- function(x) paste(x, collapse = "|")
colar.un <- function(x) paste(unique(x), collapse = "_|_")
rep.na <- function(x) {
  x[is.na(x)] <- ""
  x
}  

all.files <- vector("list", length(files))
for (i in seq_along(all.files)) {
  tmp <- read.csv(files[i], sep = ",", encoding = "UTF-8", as.is = TRUE, na.string = c(NA,""," "))
  tmp <- tmp[tmp[[taxon.id]] %in% as.double(iucn2$Redlist_id), ]
  
  if (dim(tmp)[1] > 0) {
    path.to.save <- gsub("\\.csv", "_clean.rds", files[i])
    saveRDS(tmp, path.to.save)  
    
    if (any(duplicated(tmp$internalTaxonId))) {
      cat(i,"\n")
      tmp0 <- tmp[,sapply(tmp, function(x) !all(is.na(x)))]
      for (j in 1:dim(tmp0)[2]) tmp0[[j]] <- rep.na(tmp0[[j]])
      
      tmp1 <- aggregate(. ~ internalTaxonId, tmp0, colar.un)
      tmp1 <- tmp1[match(tmp1[[taxon.id]], tmp0[[taxon.id]], nomatch = 0),]
      tmp1 <- unique(tmp1)
      tmp1[] <- lapply(tmp1, gsub, pattern = "^_\\|_|_\\|_$", replacement = "")
      
      path.to.save1 <- gsub("\\.csv", "_clean_edited.rds", files[i])
      saveRDS(tmp1, path.to.save1)  
    } 
  }  
}


#### MERGING THE AVAILABLE INFO IN LESS FILES ####
# _clean.rds : previous assessments info filtered for THREAT species
# _clean_edited.rds : previous assessments info filtered for THREAT species and
#                     aggregated for each taxon ID (countries, habitats, conservation and research needs, etc...)
ed.files <- list.files(path, pattern = "rds", full.names = TRUE)

## For the comparison between assessments
ass.files <- c("assessments", "taxonomy", "all_other_fields", "simple_summary")

iucn_assess <- NULL
for (i in seq_along(ass.files)) {
  path.i <- ed.files[grepl(ass.files[i], ed.files)]
  tmp.i <- readRDS(path.i)
  if (is.null(iucn_assess)) iucn_assess <- tmp.i else iucn_assess <- dplyr::left_join(iucn_assess, tmp.i, by = taxon.id) 
}  
iucn_assess <- iucn_assess[, !grepl("\\.y$", names(iucn_assess))]
iucn_assess <- iucn_assess[, !grepl("\\.x\\.x$", names(iucn_assess))]
names(iucn_assess) <- gsub("\\.x$", "", names(iucn_assess)) 
iucn_assess <- iucn_assess[, !duplicated(names(iucn_assess))]
saveRDS(iucn_assess, "data/IUCN_2022_v2_assessments_THREAT.rds")


## For the construction of the SIS Connect Files
sis.files <- c("common_names", "conservation_needed", "countries", "credits", "habitats",
               "plant_specific", "references", "research_needed", "synonyms", "threats",
               "usetrade")

iucn_sis <- NULL
for (i in seq_along(sis.files)) {
  path.i <- ed.files[grepl(sis.files[i], ed.files)]
  path.i <- path.i[!grepl("_edited", path.i)]
  tmp.i <- readRDS(path.i)
  tmp.i$sic.connect.file <- sis.files[i]
  
  if ("code" %in% names(tmp.i))
    tmp.i$code <- as.character(tmp.i$code)
  
  if (is.null(iucn_sis)) iucn_sis <- tmp.i else iucn_sis <- dplyr::bind_rows(iucn_sis, tmp.i) 
}
saveRDS(iucn_sis, "data/IUCN_2022_v2_sis_connect_THREAT.rds")

# #### SOME PRELIMINARY EXPLORATION FOR PREVIOUS SPECIES  ####
# iucn_assess <- readRDS("data/IUCN_2022_v2_assessments_THREAT.rds")
# iucn_sis <- readRDS("data/IUCN_2022_v2_sis_connect_THREAT.rds")
# 
# ## Generation length
# 100*table(iucn_assess$GenerationLength.range, useNA = "always")/dim(iucn_assess)[1]
# iucn_gl <- iucn_assess[!is.na(iucn_assess$GenerationLength.range),
#                        c("scientificName","redlistCategory","redlistCriteria","yearPublished","GenerationLength.range")]
# 100*dim(iucn_gl)[1]/dim(iucn_assess)[1]
# iucn_gl$mean <- sapply(strsplit(iucn_gl$GenerationLength.range, "-"), function(x) mean(as.double(x)))
# iucn_gl.myrcia <- iucn_gl[order(as.double(iucn_gl$mean)),][grepl("Myrcia", iucn_gl[order(as.double(iucn_gl$mean)),"scientificName"]),]
# table(iucn_gl.myrcia$GenerationLength.range)
# 
# 
# hab1 <- readRDS("data/threat_habitats_preliminar1.rds")
# hab1.myrcia <- hab1[hab1$Name_submitted %in% iucn_gl.myrcia$scientificName,]
# hab1.myrcia <- hab1.myrcia[match(iucn_gl.myrcia$scientificName, hab1.myrcia$Name_submitted),]
# plot(jitter(iucn_gl.myrcia$mean), jitter(hab1.myrcia$GL), 
#      xlim=c(10,80), ylim = c(10,80), xlab = "GL IUCN", ylab = "GL THREAT",
#      main = "GL for Myrcia spp"); abline(0,1)
