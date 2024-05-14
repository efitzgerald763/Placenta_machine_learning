library(fastDummies)


#--------------- Metabolic ---------------------
Metabolic <- fread("/data1/meaneylab/eamon/Placenta_WGCNA/Placenta-Eamon/clinical_data_placenta1N2.txt") %>%
  dplyr::rename(Subject.ID = SubjectID) %>%
  dplyr::select(Subject.ID, mother_gdm, ogtt_fasting, ogtt_2hour, ppBMI)
Metabolic$mother_gdm <- ifelse(Metabolic$mother_gdm == "No", 0, 1)

# Convert all columns to numeric except Subject.ID
Metabolic <- Metabolic %>% 
  mutate(across(-Subject.ID, as.numeric))

#--------------- Birth metrics ---------------------

BMI <- read_spss("/n0/gusto/DATA/BMI/GUSTO - ZBMI all 06Aug2019 - restructured - post_ageQC.sav")%>%
  dplyr::rename(Subject.ID = SubjectID) %>%
  dplyr::select(Subject.ID, GA,weight.day1,zbmi.day1)

# Convert all columns to numeric except Subject.ID
BMI <- BMI %>% 
  mutate(across(-Subject.ID, as.numeric))




#--------------- Delivery_info ---------------------
Delivery_info <- openxlsx::read.xlsx("/data1/meaneylab/eamon/FORMA_(332)_CS_28-DEC-2016_DeliveryGUSTOCaseReportForm.xlsx", sheet = 1)%>%
  dplyr :: select(1,11,30,31,32,34,40,47,54,59,62,117,118,126,130) %>%
  dplyr::rename(Subject.ID = 1) %>%
  dplyr::rename(Antepart_haemorrhage = "11") %>% # 1_No not_answered
  dplyr::rename (Ruptured_membranes = "30") %>% # 0_no not_answered
  dplyr::rename (Mat_pyrexia = "31") %>% # 0_no not_answered
  dplyr::rename (Chorio = "32") %>% # 0_no not_answered
  dplyr::rename (Antibiotics = "34") %>% # 0_no not_answered
  dplyr::rename (Tocolysis = "40")  %>% # 0_no not_answered not_applicable
  dplyr::rename(Mat_steroids = "47") %>% # 1_No 2_Incomplete_less_24hrs_bef_birth 3_Complete_1_course 4_Two_courses_or_more not_answered
  dplyr::rename (Analgesia = "54") %>% # 0_None not_answered
  dplyr::rename (Mode_of_deliv = "59") %>% # Categories not_answered
  dplyr::rename (Spon_induc_labour = "62") %>% #1_Spontaneous 2_Induced 3_na not_applicable
  dplyr::rename (APGAR_1_min = "117") %>% # 1-9 score 99_not_recorded
  dplyr::rename (APGAR_5_min = "118") %>% # 1-10 score 99_not_recorded
  dplyr::rename (Placenta_weight = "126") %>% # Continuous
  dplyr::rename (Neonatal_jaund = "130") # 0_No_significant_jaundice not_answered

Delivery_info = Delivery_info[-1,]

# Convert values to NA
na_values <- c("not_answered", "not_applicable", "99_not_recorded", "not recorded", "3_na")
Delivery_info <- Delivery_info %>%
  mutate(across(everything(), ~ifelse(. %in% na_values, NA, .)))

# Binary classifications
No_values <- c("1_No", "0_no", "0_None", "0_No_significant_jaundice","1_Spontaneous")
Delivery_info <- Delivery_info %>%
  mutate(across(everything(), ~ifelse(. %in% No_values, 0, .)))

# Convert non-0 and non-NA entries to 1 in binary columns
Delivery_info <- Delivery_info %>%
  mutate(across(c(Antepart_haemorrhage, Ruptured_membranes, Mat_pyrexia, Chorio,
                  Antibiotics, Tocolysis, Analgesia, Neonatal_jaund, Spon_induc_labour),
                ~ ifelse(. != "0" & !is.na(.), "1", .)))

# Ordinal encoding for Maternal steroids
levels_order <- c("0", 
                  "2_Incomplete_less_24hrs_bef_birth", 
                  "3_Complete_1_course", 
                  "4_Two_courses_or_more")

# Convert the column to an ordered factor
Delivery_info$Mat_steroids <- factor(Delivery_info$Mat_steroids, levels = levels_order, ordered = TRUE)

# OneHotEncoding for delivery mode
Delivery_info$Mode_of_deliv <- as.factor(Delivery_info$Mode_of_deliv)
Delivery_info <- dummy_cols(Delivery_info, select_columns = c("Mode_of_deliv"), ignore_na = TRUE)

# Remove unencoded column
Delivery_info$Mode_of_deliv <- NULL

# Continuous column
Delivery_info$Placenta_weight <- as.numeric(Delivery_info$Placenta_weight)

#---------------- Sex ----------------------------------
sex <- openxlsx::read.xlsx("/share/projects/gusto/DATA/FORMS/FORMA_(340)_CS_29-SEPT-2018UPDATE/compiled data/FormA340_20181013.xlsx", sheet = 1, check.names = T) %>%
  dplyr :: select(SubjectID, sex) %>%
  dplyr :: rename(Subject.ID = SubjectID)

sex$sex <- ifelse(sex$sex == "Male", 0, 1)

#---------------------- Cytokines ---------------------- 
# Remove adiponectin dilutions, leaving only the final data
ELISA <- openxlsx::read.xlsx("/data1/meaneylab/eamon/Placenta_study/Cord cytokines/ELISA/SIgN-ELISA_data.xlsx")%>%
  dplyr::select(-c(`Adipo_ngml.(2000X.dilution)`,`Adipo_ngml.(neat)`))

# Remove data with CV greater than 20% between plates
Luminex <- openxlsx::read.xlsx("/data1/meaneylab/eamon/Placenta_study/Cord cytokines/Luminex/SIgN-Luminex_data.xlsx")%>%
  dplyr::select(-c(`Eotaxin.(pg/ml)`,`Ghrelin.(pg/ml)`, `IL-2.(pg/ml)`, `IL-7.(pg/ml)`, `MIG.(pg/ml)`, `MIP-3alpha.(pg/ml)`,
                   `T4.(pg/ml)`,`Estradiol.(pg/ml)`, `Cortisol.(pg/ml)`, `IGFBP-4.(pg/ml)`, `PIGF-1.(pg/ml)`, `PAI-1.(pg/ml)`,
                   `TGF-beta1.(pg/ml)`, `hIgG4.(pg/ml)`, `hIgG3.(pg/ml)`, `hIgG2.(pg/ml)`, `hIgG1.(pg/ml)`, `IgM.(pg/ml)`,
                   `IgA.(pg/ml)`))

# Remove data with CV greater than 20% between plates
Quanterix <- openxlsx::read.xlsx("/data1/meaneylab/eamon/Placenta_study/Cord cytokines/Quanterix/Quanterix_cytokines_data.xlsx")%>%
  dplyr::rename(PSCID = SubjectID)%>%
  dplyr::select(-c("IFNg","IL-4"))

Cord_cyto <- full_join(ELISA, Luminex, by = "PSCID") %>%
  full_join(Quanterix,         by = "PSCID")

colnames(Cord_cyto) <- c('PSCID', 'Free_testosterone', 'Total_testosterone', 'Adiponectin', 'C_peptide',
                         'IL1RA', 'IP10', 'MCP1', 'MIP1a', 'MIP1b', 'VEGFA', 'GLP1', 'IL12p40', 'Leptin',
                         'Glucagon', 'Growth_hormone', 'IGFBP3', 'IGFBP7', 'Insulin', 'Prolactin', 'FSH',
                         'LH', 'TSH', 'IgE', 'CRP', 'IL10', 'IL6', 'TNFa')
Cord_cyto[, -c(1)] <- log10(Cord_cyto[, -c(1)])

Cord_cyto <- Cord_cyto %>%
  dplyr::rename(Subject.ID = PSCID)

#-------------------------- EPDS --------------------------------

STAI_EPDS <- readxl::read_excel ("/n0/gusto/DATA/FORMS/FORMA_(219)_CS_25-SEPT-2018UPDATE/clean dataset/FormA219_20180925 contains EPD, STAI.xlsx") %>%
  dplyr::rename(Subject.ID = SubjectID) %>%
  dplyr::select(Subject.ID, mother_age_delivery,mother_highest_education,mother_income,parity,EPDS_imputed_pw26)


# Ordinal encoding for mother_income
levels_order <- c("0_999", 
                  "1000_1999", 
                  "2000_3999", 
                  "4000_5999",
                  "more_than_6000")

# Convert the column to an ordered factor
STAI_EPDS$mother_income <- factor(STAI_EPDS$mother_income, levels = levels_order, ordered = TRUE)

# Ordinal encoding for mother_income
levels_order <- c("no_education", 
                  "primary", 
                  "secondary",
                  "gce",
                  "ite_ntc",
                  "university")

# Convert the column to an ordered factor
STAI_EPDS$mother_highest_education <- factor(STAI_EPDS$mother_highest_education, levels = levels_order, ordered = TRUE)

# Convert all columns to numeric except Subject.ID
STAI_EPDS <- STAI_EPDS %>% 
  mutate(across(-Subject.ID, as.numeric))

#-------------------- Maternal age @ delivery -------------------- 
Mat_age_at_deliv <- openxlsx::read.xlsx("/n6/n0_backup/n0/gusto/DATA/Demographics/GUSTO_DOB, Gender & Ethnicity_mat.age_20220919.xlsx") %>%
  dplyr::rename(Subject.ID = subjid) %>%
  dplyr::select(Subject.ID, Maternal_age_at_delivery) %>%
  na.omit()


#-------------------- PRS --------------------
names.PRS <- c("PRS_0.0001", "PRS_0.001", "PRS_0.01", "PRS_0.05", "PRS_0.1", "PRS_0.2", "PRS_0.5")

b249 <- fread('/n9/cortexjobs/b_prs/jobs/b249_eamon_crossdisorderpgc.gusto.kids_b249/result/b249_prs.score.csv', header = T, stringsAsFactors = F) %>% 
  mutate_at(names.PRS, ~(scale(.) %>% as.vector))%>% 
  mutate(ID1 = gsub("B", "", ID1)) 

b255 <- fread("/n9/cortexjobs/b_prs/jobs/b255_eamon_mddHoward19.gusto.kids_b255/result/b255_prs.score.csv", header = T, stringsAsFactors = F) %>% 
  mutate_at(names.PRS, ~(scale(.) %>% as.vector))%>% 
  mutate(ID1 = gsub("B", "", ID1))

colnames(b249) <- c('Subject.ID', paste0(colnames(b249)[-1], '_', 'crossdisorderpgc'))

colnames(b255) <- c('Subject.ID', paste0(colnames(b255)[-1], '_', 'depression'))

# Select PRS and subject columns before merging
df1_selected <- b249 %>%
  dplyr::select(Subject.ID, contains("PRS"))

df2_selected <- b255 %>%
  dplyr::select(Subject.ID, contains("PRS"))

PRS_df <- inner_join(df1_selected, df2_selected, by = "Subject.ID")

# Convert all columns to numeric except Subject.ID
PRS_df <- PRS_df %>% 
  mutate(across(-Subject.ID, as.numeric))


load("/data1/meaneylab/eamon/Placenta_WGCNA/rWGCNA/Reup/Two_batches/Final_image.RData")

allTraits <- rownames_to_column(allTraits, var = "Sample") %>%
  dplyr::rename(Subject.ID = SubjectID) %>%
  dplyr::select(Sample, Subject.ID)

MEs <- consMEsOnDatExpr$eigengenes
MEs <- rownames_to_column(MEs, var = "Sample")

ME_2 <- merge(MEs,allTraits,by="Sample")
ME_2$Sample <- NULL

ME_2 <- ME_2 %>% 
  mutate(across(-Subject.ID, as.numeric))




Full_data <- merge(ME_2, PRS_df, by = "Subject.ID", all = T) %>%
  merge(      Metabolic, by = "Subject.ID", all = T) %>%
  merge(      Cord_cyto, by = "Subject.ID", all = T) %>%
  merge(      BMI, by = "Subject.ID", all = T)  %>%
  merge(      Delivery_info, by = "Subject.ID", all = T) %>%
  merge(      Mat_age_at_deliv, by = "Subject.ID", all = T) %>%
  merge(      sex, by = "Subject.ID", all = T) %>%
  merge(      STAI_EPDS, by = "Subject.ID", all = T)

write.csv(Full_data,"/data1/meaneylab/eamon/Placenta_WGCNA/ML/ML_input.csv")
