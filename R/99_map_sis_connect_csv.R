################################################################
#### WHERE EACH IUCN SIS CONNECT FILES ARE BEING GENERATED? ####
################################################################

## CODES DONE, FILES IN FOLDER
# "commonnames_threat.csv" => 18_prep_sis_connect.R - Appendix ZZ
# "countries_threat.csv" => 18_prep_sis_connect.R - Appendix 
# "credits_threat.csv" => 18_prep_sis_connect.R - Appendix         
# "endemism_threat.rds" => 18_prep_sis_connect.R - Appendix 
# "habitats_threat.csv" => 18_prep_sis_connect.R - Appendix 
# "plantspecific_threat.csv" => 18_prep_sis_connect.R - Appendix CC
# "prev_assessments_threat.rds" => 18_prep_sis_connect.R - Appendix 
# "synonyms_threat.csv" => 18_prep_sis_connect.R - Appendix YY
# "taxonomy_threat.csv" => 18_prep_sis_connect.R - Appendix XX
# "threats_threat.csv" => 18_prep_sis_connect.R - Appendix 
# "usetrade_threat.csv" => 18_prep_sis_connect.R - Appendix 
# "conservationneeded_threat.csv" => 18_prep_sis_connect.R - Appendix
# "researchneeded_threat.csv" => 18_prep_sis_connect.R - Appendix
# "assessments_threat.csv" => 18_prep_sis_connect.R - Appendix
# "allfields_threat.csv" => 18_prep_sis_connect.R - Appendix
# "references_threat.csv" => 18_prep_sis_connect.R - Appendix

## Uploading all files and saving as a single zip
csv.files <- list.files("data/sis_connect/", full.names = TRUE)
csv.files <- csv.files[grepl("_threat.csv|_threat.rds", csv.files, perl = TRUE)]
caminho <- paste0(here::here(), "/data/sis_connect/THREAT")
dir.create(caminho)

for (i in seq_len(length(csv.files))) {
  if (grepl("csv$", csv.files[i], perl = TRUE)) { 
    arquivo <- read.csv(csv.files[i], encoding = "UTF-8")
  } else {
    
    if (grepl("zip$", csv.files[i], perl = TRUE)) {
      planilha <- gsub('data/sis_connect/', '', csv.files[i], fixed = TRUE)
      planilha1 <- gsub('_threat', '', planilha, perl = TRUE)
      file.copy(from = paste0(here::here(), "/data/sis_connect/", planilha),
                to = paste0(caminho,"/", planilha1))
      next
    } else {
      
      if (grepl("polygons", csv.files[i], perl = TRUE)) {
        planilha <- gsub('data/sis_connect/', '', csv.files[i], fixed = TRUE)
        planilha1 <- gsub('_threat', '', planilha, perl = TRUE)
        file.copy(from = paste0(here::here(), "/data/sis_connect/", planilha),
                  to = paste0(caminho,"/", planilha1))
        
      } else {
        arquivo <- readRDS(csv.files[i])
      }  
    }    
  }

  planilha <- gsub('data/sis_connect/', '', csv.files[i], fixed = TRUE)
  planilha <- gsub('_threat.csv|_threat.rds', '', planilha, perl = TRUE)
  
  write.csv(arquivo, paste0(caminho, "/", planilha ,".csv"), 
            row.names = FALSE, fileEncoding = "UTF-8")
}

files2zip <- dir(caminho, full.names = TRUE)
zip(zipfile = caminho, files = files2zip)
unlink(caminho, recursive = TRUE)
