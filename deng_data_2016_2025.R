#download and combine dengue virus files
rm(list = ls())

library(readr)
library(dplyr)
library(data.table)

# condensing municipality information rules
# if COMUNINF exists, use that
# otherwise, use residence
# source indicators will track what municipality ID was used: 
# infection_recorded --> COMUNINF exists, different from residence
# autochthonous_residence --> COMUNINF exists, same as residence
# residence_fallback --> COMUNINF does not exist, used residence instead

years <- 16:25
folder_path <- "/Users/allison/Desktop/Research_BentoLab/SP26/data"

# columns to keep when cleaning
req_col <- c("TP_NOT", "DT_NOTIFIC", "SEM_NOT", "NU_ANO", "ID_MN_RESI", "COMUNINF", "CLASSI_FIN", "CRITERIO", "CS_SEXO", "NU_IDADE_N")

# lists to store dataframes
file_list <- list.files(path = "/Users/allison/Desktop/Research_BentoLab/SP26/data") # list of raw datafiles
clean_list <- list() # list with cleaned dataframes

for (year in years) {
  name <- paste0("DENGBR", year, ".csv")
  path <- file.path(paste0(folder_path,"/dengue"), name)
  
  if (file.exists(path)) {
    df <- fread(path, select = req_col, na.strings = c("", "0"))
    
    # ensure all req_col exist
    missing_cols <- setdiff(req_col, names(df))
    if(length(missing_cols) > 0) df[, (missing_cols) := NA]
    
    df <- df[, ..req_col]  # order columns
    
    # create source_indicator and MUNI immediately
    df[, COMUNINF := fifelse(COMUNINF %in% c("", "0"), NA_character_, as.character(COMUNINF))]
    df[, source_indicator := fifelse(!is.na(COMUNINF) & COMUNINF != as.character(ID_MN_RESI),
                                     "infection_recorded",
                                     fifelse(!is.na(COMUNINF) & COMUNINF == as.character(ID_MN_RESI),
                                             "autochthonous_residence",
                                             "residence_fallback"))]
    df[, MUNI := fifelse(!is.na(COMUNINF), COMUNINF, as.character(ID_MN_RESI))]
    
    # reorder columns for output
    df <- df[, .(TP_NOT, DT_NOTIFIC, SEM_NOT, NU_ANO, MUNI,
                 CLASSI_FIN, CRITERIO, CS_SEXO, NU_IDADE_N, source_indicator)]
    
    clean_list[[length(clean_list)+1]] <- df
    cat("Processed", name, "\n")
  } else {
    cat(name, "not found\n")
  }
}

# combine all years
combined <- rbindlist(clean_list)

# fwrite
fwrite(combined, file.path(folder_path, "deng_data_2016_2025.csv"))