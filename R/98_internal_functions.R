#' 
#' @title Remove Unwanted Spaces
#' 
#' @param x a character or vector
#'
#' @return the character `x` without trailing or double spaces
#'  
#' @keywords internal
#'
#' @noRd
#' 
.squish <- function (x) {
  x <- gsub("\\s\\s+", " ", as.character(x), perl = TRUE)
  # x <- gsub("  ", " ", x, perl = TRUE)
  x <- gsub("^ | $", "", x, perl = TRUE)
  return(x)
}
#' 
#' @title Build Organism Name
#' 
#' @description Combine diffent table columns with species name information (i.e.
#'   genus, epiteth, infra-epiteth) into a single organism name
#' 
#' @param x the data frame with the taxonomic information
#' @param col.names the name of the columns containing the information to be
#'   combined in the desired order
#'
#' @return a vector with the combined information
#' 
#'  
#' @keywords internal
#' 
#' @noRd
#' 
.build.name <- function(x, col.names = c("Genus_original", "Species_original"))
  {
  
  if (any(!col.names %in% colnames(x)))
    stop("One or more names in 'col.names' were not found in 'x'")
  
  # cols <- names(x)[names(x) %in% col.names]
  cols <- names(x)[match(col.names, names(x), nomatch = 0)]

  if (length(cols) > 1) {
    
    #organismName <- apply(x[, cols], 1, paste, collapse = " ")
    organismName <- do.call(paste, x[, cols])
    organismName[organismName %in% c("NA NA", "NA NA NA",
                                     "NULL NULL", "NULL NULL NULL")] <- NA_character_
    organismName <- gsub(" NA NA$| NULL NULL$", "", organismName, perl = TRUE)
    organismName <- gsub(" NA$| NULL$", "", organismName, perl = TRUE)
    organismName <- .squish(organismName)
    return(organismName)
  } else {
    warning("Less than two columns found; skipping...")
  }    
}
#' 
#' @title Complete SIS Connect Files
#' 
#' @description Combine the output SIS Connect files from THREAT assessments 
#'   with those from the available IUCN assessments
#' 
#' @param x the target data frame with the SIS connect file to be completed
#' @param type the type of SIS connect file 
#' @param sis the data frame with the SIS connect file from the IUCN Red List
#' @param x.spp vector with the full species names (binomials) in `x` for
#'   merging
#' @param iucn.spp vector with the full species names (binomials) in `sis` for
#'   merging
#' @param replace should the different info found in `sis` be merged, replaced
#'   or added to `x`? Defaults to TRUE.
#' 
#' @return the input data frame x with the extra columns resulting from the
#'   merge bewteen input and IUCN SIS connect files
#' 
#'  
#' @keywords internal
#' 
#' @noRd
#' 
.merge_sis_connect <- function(x = NULL, type = NULL,
                               sis = NULL, 
                               x.spp = NULL,
                               tax.table = NULL, 
                               replace = TRUE) {
  
  # x <- allfields
  # type <- "allfields"
  # sis <- prev.assess
  # x.spp <- prev.assess$species.correct2
  # tax.table <- tax
  # replace <- TRUE
  
  # Format standardizations
  sis[] <- lapply(sis, gsub, pattern = "&amp;", replacement = "&", fixed = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "\t$", replacement = "", perl = TRUE)
  sis[] <- lapply(sis, gsub, pattern = " \n$", replacement = "", perl = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "<em>", replacement = "", fixed = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "</em>", replacement = "", fixed = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "<i>", replacement = "", fixed = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "</i>", replacement = "", fixed = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "<i>", replacement = "", fixed = TRUE)
  sis[] <- lapply(sis, .squish)
  
  if (type == "taxonomy") {
    x$scientificName <- x.spp
    
    iucn.spp.col <- "scientificName"
    
    colunas <- c("kingdomName", "phylumName", "className", "orderName", 
                      "familyName", "genusName", "speciesName", "infraType", 
                      "infraName", "infraAuthority", "authority","taxonomicNotes")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- c("kingdom", "phylum", "classname", "ordername", "family", 
                   "genus", "species", "infraType", "infra_name", 
                   "infra_authority", "taxonomicAuthority", 
                   "TaxonomicNotes.value")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    x.temp <- dplyr::left_join(x, sis.prev, suffix = c("", ".iucn_rl"), 
                               by = iucn.spp.col)
    rep.col <- paste0("merged.", type)
    x.temp[[rep.col]] <- FALSE
      
    if (replace) {
      
      vector.true <- rep(TRUE, length(nomes))
      names(vector.true) <- nomes
      no.rep <- c("infraType", "infra_authority", "infra_name")
      vector.true[names(vector.true) %in% no.rep] <- FALSE
      rep.tax <- nomes[vector.true]

      vector.action <- rep("replace", length(nomes))
      names(vector.action) <- nomes
      concatenar <- c("TaxonomicNotes.value")
      vector.action[names(vector.action) %in% concatenar] <- "concatenate"
      ignore <- c("taxonomicAuthority")
      vector.action[names(vector.action) %in% ignore] <- ""
      
      for (i in seq_along(rep.tax)) {
        col.i <- rep.tax[i]
        col.rep.i <- 
          colnames(x.temp)[which(grepl(col.i, colnames(x.temp), perl = TRUE))]
        flag.i <- !is.na(x.temp[[col.rep.i[2]]]) &
                    tolower(gsub(" ", "", x.temp[[col.rep.i[1]]], perl = TRUE)) != tolower(gsub(" ", "", x.temp[[col.rep.i[2]]], perl = TRUE))
        if (any(flag.i)) {
          vector.action.i <- vector.action[names(vector.action) %in% col.i]
          new.col.i <- paste0("merged.", col.i)
          
          if (vector.action.i == "replace") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            x.temp[flag.i, new.col.i] <- x.temp[flag.i, col.rep.i[2]]
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
          
          if (vector.action.i == "concatenate") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            
            f1 <- function(x) paste(unique(x), collapse = "|")
            x.temp[flag.i, new.col.i] <- 
              apply(x.temp[flag.i, col.rep.i], 1, f1)
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
        }
      }  
    }
    
    new.columns <- colnames(x.temp)[!colnames(x.temp) %in% colnames(x)]
    x1 <- cbind.data.frame(x, x.temp[ , new.columns])
    
    x2 <- x1[order(as.double(gsub("sp", "", x1$internal_taxon_id))), ]
    
    return(x2)
  }

  if (type == "commonnames") { 
    
    iucn.spp.col <- "scientificName"
    
    colunas <- c("name", "language")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- c("name", "language")
    names(sis.prev) <- c("validName", nomes)
    sis.prev$combo <- 
      paste(sis.prev[["validName"]], sis.prev[[colunas[1]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[["validName"]], x[[colunas[1]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- "validName"
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = "validName")
    x1 <- x[ , c("validName", nomes, "internal_taxon_id")]
    x1$source.info <- "Brazilian_Flora_Online_v.393.285"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    # sis.prev.miss2 <- sis.prev.miss1[,match(names(x), names(sis.prev.miss1), nomatch = 0)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    
    return(x3)  
  }
  
  if (type == "countries") {
    
    iucn.spp.col <- "scientificName"
    
    colunas <- c("name", "code", "presence", "origin", "seasonality", "formerlyBred")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- paste0("CountryOccurrence.CountryOccurrenceSubfield.", colunas)
    nomes[1:2] <- c("CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceName", 
                    "CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceLookup")
    names(sis.prev) <- c("species.correct2", nomes)
    
    check_these <- !sis.prev[[nomes[1]]] %in% x$CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceName
    diff <- sis.prev[[nomes[1]]][check_these]
    sort(table(diff))
    diff1 <- plantR::prepCountry(countrycode::countrycode(diff, "country.name", "iso3c"), 
                                 to.lower = FALSE, special.char = FALSE, rm.abbrev = FALSE)
    diff1 <- gsub("u\\.s\\. virgin islands", "United States Virgin Islands", 
                  diff1, perl = TRUE, ignore.case = TRUE)
    diff1 <- gsub("St. Barthelemy", "St. Barthélemy", 
                  diff1, perl = TRUE, ignore.case = TRUE)
    sis.prev[[nomes[1]]][check_these] <- diff1

    sis.prev$combo <- 
      paste(sis.prev[["species.correct2"]], sis.prev[[nomes[1]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[["species.correct2"]], x[[nomes[1]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    tmp1 <- dplyr::left_join(as.data.frame(x), sis.prev.common, 
                             by = "combo", suffix = c("", ".y"))
    rep_these <- !is.na(tmp1[[paste0(nomes[4], ".y")]]) &
                  !tmp1[[paste0(nomes[4], ".y")]] %in% "Native" &
                  tmp1[[nomes[4]]] %in% "Native" 
    x[rep_these, nomes[c(4,6)]] <- tmp1[rep_these, paste0(nomes[c(4,6)], ".y")]
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    sis.prev.miss$combo <- NULL
    x$combo <- NULL
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = "species.correct2")
    
    x1 <- x[ , c("species.correct2", nomes, "internal_taxon_id")]
    x1$source.info <- "Brazilian_Flora_Online_v.393.285"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    
    ## Taging countries with uncertain occurrences (probable undeclared cultivated individuals)
    # not doing it (need to know for which new countries may be uncertain e.g. Araucaria and which are not e.g. Cabralea canjerana)
    # for (i in seq_along(unique(x1$internal_taxon_id))) {
    #   id.i <- unique(x1$species.correct2)[i]
    #   x.i <- x1[x1$species.correct2 %in% id.i, ]
    #   if (id.i %in% sis.prev$species.correct2) {
    #     iucn.i <- sis.prev[sis.prev$species.correct2 %in% id.i, ]
    #     new.countries.i <- unique(x.i[[2]][!x.i[[2]] %in% iucn.i[[2]]])
    #     
    #     if (length(new.countries.i) > 0)
    #       x1[x1$species.correct2 %in% id.i &
    #            x1$CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceName %in% new.countries.i, 
    #          "CountryOccurrence.CountryOccurrenceSubfield.origin"] <- ""
    #   }
    # }
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    
    return(x3)  
    
  }
  
  if (type == "habitats") {
    
    iucn.spp.col <- "scientificName"
    
    colunas <- c("name", "code", "majorImportance", "season", "suitability")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- c("GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName", 
               "GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup",
               "GeneralHabitats.GeneralHabitatsSubfield.majorImportance",
               "GeneralHabitats.GeneralHabitatsSubfield.season",
               "GeneralHabitats.GeneralHabitatsSubfield.suitability")
    names(sis.prev) <- c("validName", nomes)
    
    sis.prev <- sis.prev[!sis.prev[[nomes[2]]] %in% 
                           c(paste0("9.", 1:10), 
                             paste0("10.", 1:4),
                             paste0("11.", 1:6),
                             paste0("12.", 1:7),
                             paste0("13.", 1:5),
                             paste0("15.", 1:13)),]
    
    # diferenca <- sis.prev[[nomes[2]]][!sis.prev[[nomes[2]]] %in%
    #                                             unique(x[[nomes[2]]])]
    # sort(table(diferenca))

    sis.prev$combo <- 
      paste(sis.prev[["validName"]], sis.prev[[nomes[1]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    x[["validName"]] <- x.spp
    combo.x <- paste(x[["validName"]], x[[nomes[1]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- gsub(" – ", " ", combo.x, perl = TRUE)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
      
    sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    tmp1 <- dplyr::left_join(x, sis.prev.common, 
                             by = "combo", suffix = c("", ".y"))
    rep_these <- !is.na(tmp1[[paste0(nomes[3], ".y")]]) &
                  tmp1[[nomes[3]]] %in% "" 
    x[rep_these, nomes[c(3,5)]] <- tmp1[rep_these, paste0(nomes[c(3,5)], ".y")]

    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- "validName"
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = "validName")
    x1 <- x[ , c("validName", nomes, "internal_taxon_id")]
    x1$source.info <- "Brazilian_Flora_Online_v.393.285"
    x$combo <- NULL
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss1$combo <- NULL
    sis.prev.miss2 <- sis.prev.miss1[, names(x1)]
  
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    x3$validName <- NULL
    
    return(x3)
  }
  
  if (type == "plant_specific") {
    
    x$scientificName <- x.spp
    
    iucn.spp.col <- "scientificName"
    
    colunas <- c("name", "code")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- c("PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName", 
               "PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsLookup")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    # x.temp <- dplyr::left_join(x, sis.prev, suffix = c("", ".iucn_rl"), 
    #                            by = iucn.spp.col)
    # rep.col <- paste0("merged.", type)
    # x.temp[[rep.col]] <- FALSE
    
    # Removing vines, epiphytes and herbs
    sis.prev <- sis.prev[grepl(" - ",iucn.sis.plant.spec$name),]
    
    sis.prev$combo <- 
      paste(sis.prev[[iucn.spp.col]], sis.prev[[nomes[2]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    x[[iucn.spp.col]] <- x.spp
    combo.x <- paste(x[[iucn.spp.col]], x[[nomes[2]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- gsub(" – ", " ", combo.x, perl = TRUE)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- iucn.spp.col
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = iucn.spp.col)
    sis.prev.miss2 <- sis.prev.miss1[,match(names(x), names(sis.prev.miss1), nomatch = 0)]
    x1 <- x[ , c(iucn.spp.col, nomes, "internal_taxon_id")]
    x1$source.info <- "Brazilian_Flora_Online_v.393.285"
    sis.prev.miss2$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss2$combo <- NULL
    sis.prev.miss3 <- sis.prev.miss2[,match(names(x1), names(sis.prev.miss2), nomatch = 0)]
    
    x2 <- rbind.data.frame(x1, sis.prev.miss3)
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    x3[[iucn.spp.col]] <- NULL
    
    ## Removing duplicated info for the same species with lower resolution (i.e. "T" and "S")
    dup.IDs <- x3$internal_taxon_id[duplicated(x3$internal_taxon_id)]
    remove_these <- rep(TRUE, dim(x3)[1])

    for (i in seq_along(unique(dup.IDs))) {
      id.i <- unique(dup.IDs)[i]
      check_these <- x3$internal_taxon_id %in% id.i
      code.i <- 
        x3$PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsLookup[check_these]
      rep_these <- nchar(code.i) == 1 
      if (any(rep_these) & !all(rep_these)) {
        remove_these[check_these][rep_these] <- FALSE
      }
    }
    
    x4 <- x3[remove_these, ]
    return(x4)
  }
  
  if (type == "synonyms") {
   
    iucn.spp.col <- "scientificName"
    
    colunas <- c("name", "speciesName", "speciesAuthor", "infraType", 
                 "infrarankAuthor")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- c("name", "speciesName", "speciesAuthor", "infraType", 
               "infrarankAuthor")
    names(sis.prev) <- c("species.correct2" , nomes)
    sis.prev$combo <- 
      paste(sis.prev[["species.correct2"]], sis.prev[[colunas[1]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE))
    sis.prev$combo <- .squish(tolower(gsub(" [orth. error]", "", sis.prev$combo, perl = TRUE)))
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    x1 <- dplyr::left_join(x, tax.table1, by = "internal_taxon_id")
    
    combo.x <- paste(x1[["species.correct2"]], x[[colunas[1]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    sis.prev.miss$combo <- NULL
    
    x2 <- x1[ , c("species.correct2", nomes, "infrarankName", "internal_taxon_id")]
    
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = "species.correct2")
    
    sis.prev.miss1$infrarankName <- NA
    sis.prev.miss2 <- sis.prev.miss1[,match(names(x2), names(sis.prev.miss1), nomatch = 0)]

    x2$source.info <- "Brazilian_Flora_Online_v.393.285"
    sis.prev.miss2$source.info <- "IUCN_Red_List_v.2022-2"
    x3 <- rbind.data.frame(x2, sis.prev.miss2)
    
    return(x3)
  }

  if (type == "conservation_needed") {
    
    iucn.spp.col <- "scientificName"
    
    x[[iucn.spp.col]] <- x.spp
    x$ConservationActions.ConservationActionsSubfield.note.iucn_rl <- ""
    
    colunas <- c("name", "code", "note")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- paste0("ConservationActions.ConservationActionsSubfield.", colunas)
    nomes[1:3] <- c("ConservationActions.ConservationActionsSubfield.ConservationActionsName", 
                    "ConservationActions.ConservationActionsSubfield.ConservationActionsLookup",
                    "ConservationActions.ConservationActionsSubfield.note.iucn_rl")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    sis.prev$ConservationActions.ConservationActionsSubfield.note <- ""
    
    sis.prev <-
      sis.prev[!sis.prev[[nomes[2]]] %in% c("4.1"), ]
    # sis.prev[[nomes[1]]] <- gsub("", "", sis.prev[[nomes[1]]], perl = TRUE)
    
    # check_these <- !sis.prev[[nomes[2]]] %in% x[[nomes[2]]]
    # diff <- sis.prev[[nomes[2]]][check_these]
    # table(diff)
    
    sis.prev$combo <- 
      paste(sis.prev[[iucn.spp.col]], sis.prev[[nomes[2]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[[iucn.spp.col]], x[[nomes[2]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    tmp1 <- dplyr::left_join(x, sis.prev.common, 
                             by = "combo", suffix = c("", ".y"))

    # Adding previous notes (not doing it - too difficult, using separate columns)
    # rep_these <- !is.na(tmp1[paste0(nomes[3], ".y")])
    # if (any(rep_these))
    #   x[rep_these, nomes[3]] <- tmp1[rep_these, paste0(nomes[3], ".y")]
    
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% x$combo, ]
    sis.prev.miss$combo <- NULL
    x$combo <- NULL
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- iucn.spp.col
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = iucn.spp.col)
    
    x1 <- x[ , c(iucn.spp.col, nomes, 
                 "ConservationActions.ConservationActionsSubfield.note", "internal_taxon_id")]
    x1$source.info <- "THREAT_project"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    x2$remove_these <- TRUE
    
    ## Replacing the categories Scale Unknown by categories with known threat scale
    # not doing it: ask CNCFlora what they prefer
    unknown <- c("3.4")
    for (i in seq_along(unique(x1$internal_taxon_id))) {
      id.i <- unique(x1$scientificName)[i]
      x.i <- x1[x1$scientificName %in% id.i, ]
      if (id.i %in% sis.prev$scientificName) {
        iucn.cats.i <- sis.prev[sis.prev$scientificName %in% id.i, ]
        new.cats.i <- unique(x.i[[3]][!x.i[[3]] %in% iucn.cats.i[[3]]])

        unknown.cats.i <- new.cats.i[new.cats.i %in% unknown]
        if (length(unknown.cats.i) > 0) {
          for (j in seq_along(unknown.cats.i)) {
            if (any(grepl(unknown.cats.i[j], iucn.cats.i[[3]])))
              x2[x2$scientificName %in% id.i &
                   x2[[3]] %in% unknown.cats.i[j],
                 "remove_these"] <- FALSE
          }
        }
      }
    }
    x3 <- x2[x2$remove_these, ]
    x3$remove_these <- NULL

    x3 <- x3[order(as.double(gsub("sp", "", x3$internal_taxon_id))), ]
    x3 <- x3[ , c("ConservationActions.ConservationActionsSubfield.ConservationActionsLookup",
                  "ConservationActions.ConservationActionsSubfield.ConservationActionsName",
                  "ConservationActions.ConservationActionsSubfield.note", 
                  "internal_taxon_id",
                  "ConservationActions.ConservationActionsSubfield.note.iucn_rl")]
    return(x3)
  }
  
  if (type == "research_needed") {
    
    iucn.spp.col <- "scientificName"
    
    x[[iucn.spp.col]] <- x.spp
    x$Research.ResearchSubfield.note.iucn_rl <- ""
    
    colunas <- c("name", "code", "note")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- paste0("Research.ResearchSubfield.", colunas)
    nomes[1:3] <- c("Research.ResearchSubfield.ResearchName", 
                    "Research.ResearchSubfield.ResearchLookup",
                    "Research.ResearchSubfield.note.iucn_rl")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    sis.prev$Research.ResearchSubfield.note <- ""
    
    sis.prev <-
      sis.prev[!sis.prev[[nomes[2]]] %in% c("4"), ]
    # sis.prev[[nomes[1]]] <- gsub("", "", sis.prev[[nomes[1]]], perl = TRUE)
    
    # check_these <- !sis.prev[[nomes[2]]] %in% x[[nomes[2]]]
    # diff <- sis.prev[[nomes[2]]][check_these]
    # table(diff)
    
    sis.prev$combo <- 
      paste(sis.prev[[iucn.spp.col]], sis.prev[[nomes[2]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[[iucn.spp.col]], x[[nomes[2]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    # sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    # tmp1 <- dplyr::left_join(x, sis.prev.common, 
    #                          by = "combo", suffix = c("", ".y"))
    # 
    # Adding previous notes (not doing it - too difficult, using separate columns)
    # rep_these <- !is.na(tmp1[paste0(nomes[3], ".y")]) &
    #                 !tmp1[paste0(nomes[3], ".y")] %in% "Population size, distribution & trends"
    # if (any(rep_these))
    #   x[rep_these, nomes[3]] <- tmp1[rep_these, paste0(nomes[3], ".y")]
    
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% x$combo, ]
    sis.prev.miss$combo <- NULL
    x$combo <- NULL
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- iucn.spp.col
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = iucn.spp.col)
    
    x1 <- x[ , c(iucn.spp.col, nomes, 
                 "Research.ResearchSubfield.note", "internal_taxon_id")]
    x1$source.info <- "THREAT_project"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    x3 <- x3[ , c("Research.ResearchSubfield.ResearchLookup",
                  "Research.ResearchSubfield.ResearchName",
                  "Research.ResearchSubfield.note", 
                  "internal_taxon_id",
                  "Research.ResearchSubfield.note.iucn_rl")]
    return(x3)
    
    
  }

  if (type == "threats") {
    
    iucn.spp.col <- "scientificName"
    
    x[[iucn.spp.col]] <- x.spp
    
    colunas <- c("name", "code", "stressCode", "stressName", "text",
                 "ancestry", "ias", "internationalTrade",
                 "scope", "severity", "timing", "virus")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- paste0("Threats.ThreatsSubfield.", colunas)
    nomes[1:2] <- c("Threats.ThreatsSubfield.ThreatsName", 
                    "Threats.ThreatsSubfield.ThreatsLookup")
    nomes[3:4] <- c("Threats.ThreatsSubfield.StressesSubfield.stress", 
                    "Threats.ThreatsSubfield.StressesSubfield.stressdesc")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    
    sis.prev <-
      sis.prev[!sis.prev$Threats.ThreatsSubfield.ThreatsLookup %in% c("1", "12.1"), ]
    
    # check_these <- !sis.prev[[nomes[2]]] %in% x$Threats.ThreatsSubfield.ThreatsLookup
    # diff <- sis.prev[[nomes[2]]][check_these]
    # table(diff)
    
    sis.prev$combo <- 
      paste(sis.prev[[iucn.spp.col]], sis.prev[[nomes[2]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[[iucn.spp.col]], x[[nomes[2]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    tmp1 <- dplyr::left_join(x, sis.prev.common, 
                             by = "combo", suffix = c("", ".y"))
    tmp1[[paste0(nomes[4], ".y")]] <- 
      stringr::str_to_title(tmp1[[paste0(nomes[4], ".y")]])
    tmp1[[paste0(nomes[4], ".y")]] <- 
      gsub(" Of ", " of ",tmp1[[paste0(nomes[4], ".y")]], perl = TRUE)
    
    
    # Adding missing stress codes
    rep_these <- !is.na(tmp1[[paste0(nomes[3], ".y")]])
    if (any(rep_these)) {
      prev.i <- strsplit(tmp1[[paste0(nomes[3], ".y")]][rep_these], "|", fixed = TRUE)
      prev.names.i <- strsplit(tmp1[[paste0(nomes[4], ".y")]][rep_these], "|", fixed = TRUE)
      prev.ii <- mapply(function(x, y) { 
                        names(x) <- y 
                        return(x) }, 
                        prev.i, prev.names.i)
      new.i <- strsplit(tmp1[[nomes[3]]][rep_these], "|", fixed = TRUE)  
      new.names.i <- strsplit(tmp1[[nomes[4]]][rep_these], "|", fixed = TRUE)
      new.ii <- mapply(function(x, y) { 
        names(x) <- y 
        return(x) }, 
        new.i, new.names.i)

      combo <- mapply(function(x, y) { 
                        vector <- c(unlist(x), unlist(y))
                        vector1 <- vector[!duplicated(tolower(names(vector)))]
                        ordem <- order(as.numeric(gsub("\\.3\\.",".3",vector1, perl = TRUE)))
                        vector2 <- vector1[ordem]
                        return(vector2)
                        }, new.ii, prev.ii)
      new.iii <- sapply(combo, paste0, collapse = "|")
      tmp1[[nomes[3]]][rep_these] <- new.iii
      
      new.names.iii <- sapply(combo, function(x) paste0(names(x), collapse = "|"))
      tmp1[[nomes[4]]][rep_these] <- new.names.iii
    }

    for (i in c(5:12)) {
      rep_these <- !is.na(tmp1[[paste0(nomes[i], ".y")]]) &
                    !tmp1[[paste0(nomes[i], ".y")]] %in% c("") &
                    tmp1[[nomes[i]]] %in% c("", "Unknown")
      if (any(rep_these))
        x[rep_these, nomes[i]] <- tmp1[rep_these, paste0(nomes[i], ".y")]
    }
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% x$combo, ]
    sis.prev.miss$combo <- NULL
    x$combo <- NULL
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- iucn.spp.col
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = iucn.spp.col)
    
    x1 <- x[ , c(iucn.spp.col, nomes, "internal_taxon_id")]
    x1$source.info <- "THREAT_project"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    # x2$remove_these <- rep(TRUE, dim(x2)[1])
    # 
    ## Replacing the categories Scale Unknown by categories with known threat scale
    # not doing it: ask CNCFlora what they prefer
    # unknown <- c("2.1.4", "2.3.4", "5.2.4", "5.3.5", "5.4.6")
    # for (i in seq_along(unique(x1$internal_taxon_id))) {
    #   id.i <- unique(x1$scientificName)[i]
    #   x.i <- x1[x1$scientificName %in% id.i, ]
    #   if (id.i %in% sis.prev$scientificName) {
    #     iucn.cats.i <- sis.prev[sis.prev$scientificName %in% id.i, ]
    #     new.cats.i <- unique(x.i[[3]][!x.i[[3]] %in% iucn.cats.i[[3]]])
    #     
    #     unknown.cats.i <- new.cats.i[new.cats.i %in% unknown]
    #     if (length(unknown.cats.i) > 0) {
    #       for (j in seq_along(unknown.cats.i)) {
    #         if (any(grepl(substr(unknown.cats.i[j],1,4), iucn.cats.i[[3]])))
    #           x2[x2$scientificName %in% id.i &
    #                x2[[3]] %in% unknown.cats.i[j],
    #              "remove_these"] <- FALSE
    #       }
    #     }
    #   }
    # }
    # x3 <- x2[x2$remove_these, ]
    # x3$remove_these <- NULL
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    
    return(x3)  
  }
  
  if (type == "usetrade") {
    
    iucn.spp.col <- "scientificName"
   
    x[[iucn.spp.col]] <- x.spp

    colunas <- c("name", "code", "international", "national", "other", "subsistence")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- paste0("UTEndUse.UTEndUseSubfield.", colunas)
    nomes[1:2] <- c("UTEndUse.UTEndUseSubfield.UTEndUseName", 
                    "UTEndUse.UTEndUseSubfield.UTEndUseLookup")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    
    sis.prev <-
      sis.prev[!sis.prev$UTEndUse.UTEndUseSubfield.UTEndUseLookup %in% c("15", "18"), ]
    
    # check_these <- !sis.prev[[nomes[2]]] %in% x$UTEndUse.UTEndUseSubfield.UTEndUseLookup
    # diff <- sis.prev[[nomes[2]]][check_these]
    # sort(table(diff))

    sis.prev$combo <- 
      paste(sis.prev[[iucn.spp.col]], sis.prev[[nomes[2]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[[iucn.spp.col]], x[[nomes[2]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    tmp1 <- dplyr::left_join(x, sis.prev.common, 
                             by = "combo", suffix = c("", ".y"))
    
    for (i in 3:6) {
      rep_these <- !is.na(tmp1[[paste0(nomes[i], ".y")]]) &
                    tmp1[[paste0(nomes[i], ".y")]] %in% "true" &
                    tmp1[[nomes[i]]] %in% "false" 
      x[rep_these, nomes[i]] <- tmp1[rep_these, paste0(nomes[i], ".y")]
    }
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    sis.prev.miss$combo <- NULL
    x$combo <- NULL
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- iucn.spp.col
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = iucn.spp.col)
    
    x1 <- x[ , c(iucn.spp.col, nomes, "internal_taxon_id")]
    x1$source.info <- "TreeCo_database_v.4.1"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    
    return(x3)  
  }

  if (type == "assessments") {
    
    iucn.spp.col <- "species.correct2"

    x[[iucn.spp.col]] <- x.spp
    
    colunas <- c("realm", 
                 "conservationActions", 
                 "habitat", 
                 #"MapStatus", 
                 "population",
                 "populationTrend",
                 "range",
                 "yearPublished",
                 "criteriaVersion",
                 #"RedListCriteria.dataDeficientReason",
                 #"RedListCriteria.isManual",
                 "redlistCategory",
                 "redlistCriteria",
                 #"RedListCriteria.manualCriteriaString",
                 "possiblyExtinct",
                 #"RedListCriteria.possiblyExtinctCandidate",
                 "yearLastSeen",
                 "rationale",
                 #"RedListReasonsForChange.catCritChanges",
                 #"RedListReasonsForChange.changeReasons",
                 #"RedListReasonsForChange.otherReason",
                 #"RedListReasonsForChange.timeframe",
                 #"RedListReasonsForChange.type",
                 "systems",
                 "threats",
                 "useTrade",
                 "language")
    sis.prev <- sis[, c("internal_taxon_id", iucn.spp.col, colunas)]
    
    nomes <- c("BiogeographicRealm.realm", 
                 "ConservationActionsDocumentation.narrative", 
                 "HabitatDocumentation.narrative", 
                 #"MapStatus.status", 
                 "PopulationDocumentation.narrative",
                 "PopulationTrend.value",
                 "RangeDocumentation.narrative",
                 "RedListAssessmentDate.value",
                 "RedListCriteria.critVersion",
                 #"RedListCriteria.dataDeficientReason",
                 #"RedListCriteria.isManual",
                 "RedListCriteria.manualCategory",
                 "RedListCriteria.manualCriteria",
                 #"RedListCriteria.manualCriteriaString",
                 "RedListCriteria.possiblyExtinct",
                 #"RedListCriteria.possiblyExtinctCandidate",
                 "RedListCriteria.yearLastSeen",
                 "RedListRationale.value",
                 #"RedListReasonsForChange.catCritChanges",
                 #"RedListReasonsForChange.changeReasons",
                 #"RedListReasonsForChange.otherReason",
                 #"RedListReasonsForChange.timeframe",
                 #"RedListReasonsForChange.type",
                 "System.value",
                 "ThreatsDocumentation.value",
                 "UseTradeDocumentation.value",
                 "language.value")
    
    names(sis.prev) <- c("internal_taxon_id", iucn.spp.col, nomes)

    x.temp <- dplyr::left_join(x, sis.prev, suffix = c("", ".iucn_rl"), 
                               by = c("internal_taxon_id", iucn.spp.col))
    rep.col <- paste0("merged.", type)
    x.temp[[rep.col]] <- FALSE
    
    if (replace) {
      
      vector.true <- rep(TRUE, length(nomes))
      names(vector.true) <- nomes
      no.rep <- c("RedListAssessmentDate.value", "System.value", "language.value")
      vector.true[names(vector.true) %in% no.rep] <- FALSE
      rep.tax <- nomes[vector.true]
      
      vector.action <- rep("replace", length(nomes))
      names(vector.action) <- nomes
      concatenar <- c("BiogeographicRealm.realm", 
                      "ConservationActionsDocumentation.narrative",
                      "HabitatDocumentation.narrative",
                      "PopulationDocumentation.narrative",
                      "RangeDocumentation.narrative",
                      "RedListRationale.value",
                      "ThreatsDocumentation.value",
                      "UseTradeDocumentation.value")
      vector.action[names(vector.action) %in% concatenar] <- "concatenate"
      ignore <- c("PopulationTrend.value", "RedListCriteria.critVersion", 
                  "RedListCriteria.manualCategory", 
                  "RedListCriteria.manualCriteria", 
                  "RedListCriteria.possiblyExtinct",
                  "RedListCriteria.yearLastSeen")
      vector.action[names(vector.action) %in% ignore] <- ""
      
      for (i in seq_along(rep.tax)) {
        col.i <- rep.tax[i]
        col.rep.i <- 
          colnames(x.temp)[which(grepl(col.i, colnames(x.temp), perl = TRUE))]
        flag.i <- !is.na(x.temp[[col.rep.i[2]]]) &
          tolower(gsub(" ", "", x.temp[[col.rep.i[1]]], perl = TRUE)) != tolower(gsub(" ", "", x.temp[[col.rep.i[2]]], perl = TRUE))
        if (any(flag.i)) {
          vector.action.i <- vector.action[names(vector.action) %in% col.i]
          new.col.i <- paste0("merged.", col.i)
          
          if (vector.action.i == "replace") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            x.temp[flag.i, new.col.i] <- x.temp[flag.i, col.rep.i[2]]
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
          
          if (vector.action.i == "concatenate") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            
            f1 <- function(x) paste(unique(x), collapse = "|")
            x.temp[flag.i, new.col.i] <- 
              apply(x.temp[flag.i, col.rep.i], 1, f1)
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
        }
      }  
    }
    
    new.columns <- colnames(x.temp)[!colnames(x.temp) %in% colnames(x)]
    x1 <- cbind.data.frame(x, x.temp[ , new.columns])
    
    x2 <- x1[order(as.double(gsub("sp", "", x1$internal_taxon_id))), ]
    
    return(x2)
    
  }

  if (type == "allfields") {
    
    iucn.spp.col <- "species.correct2"
    
    x[[iucn.spp.col]] <- x.spp
    
    colunas <- c("AOO.range",
                 "EOO.range",
                 "NoThreats.noThreats",
                 "ElevationLower.limit",
                 "ElevationUpper.limit",
                 "PopulationSize.range",
                 "ThreatsUnknown.value",
                 "LocationsNumber.range",
                 "GenerationLength.range",
                 "MovementPatterns.pattern",
                 "SubpopulationNumber.range",
                 "AreaRestricted.isRestricted",
                 "CropWildRelative.isRelative",
                 "YearOfPopulationEstimate.value",
                 "InPlaceEducationControlled.value",
                 "SevereFragmentation.isFragmented",
                 "InPlaceResearchRecoveryPlan.value",
                 "InPlaceLandWaterProtectionInPA.value",
                 "InPlaceSpeciesManagementExSitu.value",
                 "InPlaceResearchMonitoringScheme.value",
                 "InPlaceEducationSubjectToPrograms.value",
                 "InPlaceSpeciesManagementHarvestPlan.value",
                 "InPlaceLandWaterProtectionAreaPlanned.value",
                 "InPlaceEducationInternationalLegislation.value",
                 "InPlaceLandWaterProtectionInvasiveControl.value",
                 "InPlaceLandWaterProtectionSitesIdentified.value",
                 "InPlaceLandWaterProtectionPercentProtected.value")
    sis.prev <- sis[, c("internal_taxon_id", iucn.spp.col, colunas)]
    nomes <- colunas
    
    x.temp <- dplyr::left_join(x, sis.prev, suffix = c("", ".iucn_rl"), 
                               by = c("internal_taxon_id", iucn.spp.col))
    rep.col <- paste0("merged.", type)
    x.temp[[rep.col]] <- FALSE
    
    if (replace) {
      
      vector.true <- rep(TRUE, length(nomes))
      names(vector.true) <- nomes
      no.rep <- c("MovementPatterns.pattern")
      vector.true[names(vector.true) %in% no.rep] <- FALSE
      rep.tax <- nomes[vector.true]
      
      vector.action <- rep("replace", length(nomes))
      names(vector.action) <- nomes
      concatenar <- c("NoThreats.noThreats", 
                      "ElevationLower.limit",
                      "ElevationUpper.limit",
                      "PopulationSize.range",
                      "ThreatsUnknown.value",
                      "RedListRationale.value",
                      "ThreatsDocumentation.value",
                      "UseTradeDocumentation.value",
                      "AreaRestricted.isRestricted", 
                      "CropWildRelative.isRelative", 
                      "SevereFragmentation.isFragmented")
      vector.action[names(vector.action) %in% concatenar] <- "concatenate"
      ignore <- c("AOO.range", "EOO.range", "LocationsNumber.range",
                  "SubpopulationNumber.range", "GenerationLength.range",
                  "SubpopulationNumber.range", "YearOfPopulationEstimate.value",
                  "InPlaceEducationControlled.value", 
                  "InPlaceResearchRecoveryPlan.value",
                  "InPlaceLandWaterProtectionInPA.value",
                  "InPlaceSpeciesManagementExSitu.value",
                  "InPlaceResearchMonitoringScheme.value",
                  "InPlaceEducationSubjectToPrograms.value",
                  "InPlaceSpeciesManagementHarvestPlan.value",
                  "InPlaceLandWaterProtectionAreaPlanned.value",
                  "InPlaceEducationInternationalLegislation.value",
                  "InPlaceLandWaterProtectionInvasiveControl.value",
                  "InPlaceLandWaterProtectionSitesIdentified.value",
                  "InPlaceLandWaterProtectionPercentProtected.value")
      vector.action[names(vector.action) %in% ignore] <- ""
      
      for (i in seq_along(rep.tax)) {
        col.i <- rep.tax[i]
        col.rep.i <- 
          colnames(x.temp)[which(grepl(col.i, colnames(x.temp), perl = TRUE))]
        flag.i <- !is.na(x.temp[[col.rep.i[2]]]) &
          tolower(gsub(" ", "", x.temp[[col.rep.i[1]]], perl = TRUE)) != tolower(gsub(" ", "", x.temp[[col.rep.i[2]]], perl = TRUE))
        if (any(flag.i)) {
          vector.action.i <- vector.action[names(vector.action) %in% col.i]
          new.col.i <- paste0("merged.", col.i)
          
          if (vector.action.i == "replace") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            x.temp[flag.i, new.col.i] <- x.temp[flag.i, col.rep.i[2]]
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
          
          if (vector.action.i == "concatenate") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            
            f1 <- function(x) paste(unique(x), collapse = "|")
            x.temp[flag.i, new.col.i] <- 
              apply(x.temp[flag.i, col.rep.i], 1, f1)
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
        }
      }  
    }
    
    new.columns <- colnames(x.temp)[!colnames(x.temp) %in% colnames(x)]
    x1 <- cbind.data.frame(x, x.temp[ , new.columns])
    
    x2 <- x1[order(as.double(gsub("sp", "", x1$internal_taxon_id))), ]
    
    return(x2)
    
  }
  
  if (type == "references") {
    
  }
  
  if (type == "credits") {
    
  }
}
