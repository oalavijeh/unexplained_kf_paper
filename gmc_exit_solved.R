# load labkey
library(Rlabkey)

# Set the baseURL

labkey.setDefaults(baseUrl = "https://labkey-embassy.gel.zone/labkey/")

# Write your SQL query here - if a new metric is in a non-listed table in the FROM section - needs defining
# SELECT are coloumns of interest.
# WHERE are the search function - disease, type of file (Array GenotypingGenomic/Repeat/Standard/Structural VCF or BAM),build(GRCh38 or GRCh37),participant (Proband, Relative),pick europeans with >0.9

query <- "SELECT rd.participant_id, rd.participant_type, par.consanguinity, par.participant_phenotypic_sex, par.year_of_birth, par.mother_affected, par.father_affected, par.full_brothers_affected, par.full_sisters_affected, rd.normalised_specific_disease, rd.plate_key, gc.acmg_classification, gc.case_solved_family, gc.additional_comments, gc.gene_name, gc.assembly, gc.chromosome, gc.position, gc.reference, gc.alternate
FROM rare_disease_analysis AS rd
LEFT JOIN gmc_exit_questionnaire gc
ON rd.participant_id = gc.participant_id
LEFT JOIN participant AS par
ON gc.participant_id = par.participant_id
WHERE rd.normalised_specific_disease = 'Cystic kidney disease'
OR rd.normalised_specific_disease = 'Congenital Anomaly of the Kidneys and Urinary Tract (CAKUT)'
OR rd.normalised_specific_disease = 'Proteinuric renal disease'
OR rd.normalised_specific_disease = 'Unexplained kidney failure in young people'
OR rd.normalised_specific_disease = 'Extreme early-onset hypertension'
OR rd.normalised_specific_disease = 'Familial IgA nephropathy and IgA vasculitis'
OR rd.normalised_specific_disease = 'Renal tract calcification (or Nephrolithiasis or nephrocalcinosis)'
OR rd.normalised_specific_disease = 'Renal tubular acidosis'
OR rd.normalised_specific_disease = 'Atypical haemolytic uraemic syndrome'
OR rd.normalised_specific_disease = 'Primary membranoproliferative glomerulonephritis'
OR rd.normalised_specific_disease = 'Familial haematuria'
AND rd.participant_type = 'Proband'"
         

mysql <- labkey.executeSql(
  schemaName="lists",                                                 # Do not change this
  colNameOpt = "rname",                                               # Do not change this
  maxRows = 100000000,                                                # Do not change this
  folderPath="/main-programme/main-programme_v16_2022-10-13",          # This can be changed to different main programme releases
  sql = query                                                         # This can be changed to your query of choice



)

y <-data.frame(mysql)
y<- y[!duplicated(y$participant_id),]
esrf <-subset(y, participant_type =='Proband' & normalised_specific_disease == "Cystic kidney disease")
esrf[is.na(esrf)] <- "Unknown"
solved_esrf <- subset(esrf, case_solved_family == "yes")
#for loop to replace
replace_list <- c("TSC2;PKD1.*" = "TSC2/PKD1","PKD1;.*"= "PKD1","PKD2;.*" = "PKD2",
                  "DNAJB11.*" = "DNAJB11","PKHD1.*" = "PKHD1","BBS1.*"="BBS1")
for (item in names(replace_list)) {
  solved_esrf$gene_name <- gsub(item, replace_list[item], solved_esrf$gene_name)
}
#get table of solved 
solved_esrf_genes <- as.data.frame(table(solved_esrf$gene_name))




#chatGPT solutions to find and replace
replace_col_contents_multiple_conditions <- function(df, col, pattern1, pattern2, replace_with) {
  rows_to_replace <- grep(pattern1, df[, col])
  rows_to_replace <- rows_to_replace[grepl(pattern2, df[rows_to_replace, col])]
  df[rows_to_replace, col] <- replace_with
  return(df)
}
solved_esrf_replaced <- replace_col_contents_multiple_conditions(solved_esrf, "gene_name", "PKD1", "TSC2", "PKD1/TSC2")

#partial match but include ; after PKD1  
replace_col_contents_partial <- function(df, col, pattern, replace_with) {
  rows_to_replace <- grep(pattern, df[, col], perl = T)
  df[rows_to_replace, col] <- replace_with
  return(df)
}
solved_esrf_replaced <- replace_col_contents_partial(solved_esrf, "gene_name", "PKD1;", "PKD1")
solved_esrf_replaced <- replace_col_contents_partial(solved_esrf, "gene_name", "PKD2", "PKD2")
