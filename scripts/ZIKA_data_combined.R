#download and combine zika virus files
rm(list = ls())

library(readr)
library(dplyr)

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
clean_list <- list() # list with cleaned dataframes

for (year in years) {
  name <- paste0("ZIKABR",year,".csv")
  path <- file.path(paste0(folder_path,"/zika"), name)
  
  if (file.exists(path)) {
    df <- read_csv(path, show_col_types = FALSE)
    
    remove <- setdiff(req_col, names(df))
    df[remove] <- NA
    df <- df[req_col]
    clean_list[[length(clean_list)+1]] <- df
    cat("Processed", name, "\n")
  } 
  else {
    cat(name, "not found\n")
  }
}  

combined <- bind_rows(clean_list)

combined <- combined %>%
  mutate(
    COMUNINF = na_if(as.character(COMUNINF), ""),
    COMUNINF = na_if(as.character(COMUNINF), "0")
  ) %>%
  mutate(
    source_indicator = case_when(
      !is.na(COMUNINF) & COMUNINF != ID_MN_RESI ~ "infection_recorded",
      !is.na(COMUNINF) & COMUNINF == ID_MN_RESI ~ "autochthonous_residence",
      is.na(COMUNINF) ~ "residence_fallback"
    ),
    MUNI = ifelse(!is.na(COMUNINF), COMUNINF, ID_MN_RESI)
  ) %>%
  select(TP_NOT,DT_NOTIFIC, SEM_NOT, NU_ANO, MUNI, CLASSI_FIN, CRITERIO, CS_SEXO, NU_IDADE_N, source_indicator)

combined <- combined %>%
  mutate(
    UF_code = substr(MUNI, 1, 2)  # 2-digit state code
  )

uf_lookup <- tibble(
  UF_code = c(
    "11","12","13","14","15","16","17",
    "21","22","23","24","25","26","27","28","29",
    "31","32","33","35",
    "41","42","43",
    "50","51","52","53"
  ),
  UF = c(
    "RO","AC","AM","RR","PA","AP","TO",
    "MA","PI","CE","RN","PB","PE","AL","SE","BA",
    "MG","ES","RJ","SP",
    "PR","SC","RS",
    "MS","MT","GO","DF"
  )
)
combined <- combined %>%
  mutate(UF_code = substr(MUNI, 1, 2)) %>%
  left_join(uf_lookup, by = "UF_code")

write_csv(combined, file.path(folder_path, "combined_zika_2016_2025.csv"))
  
