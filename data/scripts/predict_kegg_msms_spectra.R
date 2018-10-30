kegg_compounds <- read.csv("E:/project/KMet/data/kegg_compound.csv", stringsAsFactors=FALSE)
keep <- !is.na(kegg_compounds$smile)
kegg_compounds <- kegg_compounds[keep,]
keep <- sapply(kegg_compounds$Formula, function(f) nchar(f)) >= 3
kegg_compounds <- kegg_compounds[keep,]

cfm_id <- "D:/CFM-ID/cfm-id-2.4_win32/cfm-predict.exe"
para1 <- "D:/CFM-ID/cfm-id-2.4_win32/metab_se/param_output0.log"
para2 <- "D:/CFM-ID/cfm-id-2.4_win32/metab_se/param_config.txt"
output <- "E:/project/KMet/data/kegg_predicted_msms_spectra"

for (i in 1:nrow(kegg_compounds)){
  smiles <- kegg_compounds$smile[i]
  output_file <- paste(output, '/', kegg_compounds$ID[i], '.txt', sep='')
  string <- paste(cfm_id, smiles, 0.01, para1, para2, 0, output_file)
  shell(string)
}
