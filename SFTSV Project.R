#===================================================================================================================
# SFTS (Severe Fever with Thrombocytopenia Syndrome Virus) Dissertation Analysis
# Candidate no: 1100326
# Last updated: 10/08/2026
#
# PURPOSE OF THIS SCRIPT
# This script builds a synthetic Japanese population enriched with pet-related
# behaviours (pet ownership, community cat feeding, and caring for
# acquaintances' pets), validates that population against national survey
# data, identifies plausible SFTS risk subgroups, and simulates the impact of
# a hypothetical SFTS vaccination programme on cases and deaths averted.
# This script builds a synthetic Japanese population enriched with pet-related
# behaviours (pet ownership, community cat feeding, and caring for
# acquaintances' pets), validates that population against national survey
# data, identifies plausible SFTS risk subgroups, and simulates the impact of
# a hypothetical SFTS vaccination programme on cases and deaths averted.
#
# SECTIONS IN THIS SCRIPT
#   1. Setup                          — clean survey data, construct pet-behaviour indicators, fit multinomial model
#   2. Venn diagram validation         — compare demographic profile of survey vs. synthetic population
#   3. Descriptive comparison          — apply behavioural model to synthetic population; validate overlap structure via Venn diagrams
#   4. Subgroup counting               — quantify size of each risk subgroup across 10 synthetic populations
#   5. Vaccination scenario analysis   — stochastic simulation of cases/deaths averted under hypothetical vaccination
#   6. Sensitivity analysis            — examine how results change under different reporting-fraction assumptions
#
# REPRODUCIBILITY NOTES
#   - Every stochastic step is seeded (set.seed(1802)) immediately before use,
#     so re-running this script end-to-end reproduces identical results.
#   - Run this script from an RStudio Project so relative file paths (e.g.
#     "J_petdataset.xlsx") resolve correctly on any machine.
#   - Required input files (expected in the project root):
#       J_petdataset.xlsx        — Japanese pet ownership survey (raw)
#       2015_*.csv                — 10 synthetic population files
#       Forest_Ag_Pop_2022.csv    — distance-to-forest/agriculture geospatial data
#===================================================================================================================


#===================================================================================================================
# SECTION 1: SETUP
# Clean survey data, construct pet-behaviour indicators, and fit the
# multinomial model that will later be used to project pet-related
# behaviours onto the synthetic national population (Sections 2-6).
#===================================================================================================================

#-------------------------------------------------------------------------------------------------------------------
# 1.1 Load required packages
#-------------------------------------------------------------------------------------------------------------------
library(pacman)
p_load(deSolve, cli, tidyverse, tidyr, dplyr, readxl, readr, ggplot2, doParallel, parallel, 
       scales,stringr, ggVennDiagram, RColorBrewer,xgboost,Matrix,pROC,e1071, caTools,caret, recipes,themis,nnet, tidymodels,report, randomForest,patchwork, gtsummary,gt)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Confirm the working directory before reading any files ####
# Running this script from an RStudio Project (rather than setwd()) keeps
# file paths portable across machines.
getwd() 
#setwd("/Users/fabdulsalam/Desktop/Papers/Dissertation file")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Import the synthetic population level data sets ####
# (Read individually in Sections 2-6 via list.files(); left here as a
# reminder of the naming convention used for a single synthetic file.)
#pop_level <- read.csv("2015_001_8+_42.csv")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Import the Japanese pet ownership survey ####
pet_dataset_raw <- read_excel("J_petdataset.xlsx", sheet = "data")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#verify the rows and columns of pet_dataset_raw imported
dim(pet_dataset_raw) #20825rows & 293columns

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# The original column names contain brackets (e.g. "Q1_1[1]"), which are not
# valid R names and break downstream dplyr code. make.names() converts them
# to a syntactically valid, dot-separated format (e.g. "Q1_1.1.").
names(pet_dataset_raw) <- make.names(names(pet_dataset_raw), unique = TRUE)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#data cleaning#####
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## checking for missing values in each observation
#colSums(is.na(pet_dataset_raw))

#remove features not relevant to the analysis####
# SQ  = screening questions, BD = brand/product questions, "weight" = survey
# weighting columns, and the NQ5/TQ5/UQ5/Q5 blocks relate to sections of the
# questionnaire outside the scope of this dissertation (unrelated to pet
# ownership, cat-feeding behaviour, or demographics).
remove_cols <- c(
  grep("^SQ", names(pet_dataset_raw), value = TRUE),
  grep("^BD", names(pet_dataset_raw), value = TRUE),
  grep("weight", names(pet_dataset_raw), value = TRUE),
  grep("^NQ5", names(pet_dataset_raw) , value = TRUE),
  grep("^TQ5", names(pet_dataset_raw), value = TRUE),
  grep("^UQ5", names(pet_dataset_raw), value = TRUE),
  grep("^Q5", names(pet_dataset_raw), value = TRUE)
)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#assigned the cleaned pet_dataset_raw to pet_dataset####
pet_dataset <- pet_dataset_raw %>%select(-all_of(unique(remove_cols)))
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#dimension of the new cleaned dataset####
#dim(pet_dataset) #20825rows|#110columns
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Integrate Regions to the pet_dataset as they were listed via this link https://www.eu-japan.eu/eubusinessinjapan/about-japan/Regions-prefectures####
# Recode prefecture codes (QW3) into Japan's eight administrative regions.
pet_dataset<- pet_dataset %>%
  mutate(Region = case_when(
    QW3 == 1 ~ "Hokkaido",
    QW3 %in% c(2:7) ~ "Tohoku",
    QW3 %in% c(8:14) ~ "Kanto",
    QW3 %in% c(15:23) ~ "Chubu",
    QW3 %in% c(24:30) ~ "Kansai",
    QW3 %in% c(31:35) ~ "Chugoku",
    QW3 %in% c(36:39) ~ "Shikoku",
    QW3 %in% c(40:47) ~ "Kyushu",
    TRUE ~ "F"  # Prefer not to say / not provided
  ))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Construct metrics to capture various levels of contact with domesticated cats/dogs####
# Three behaviours of interest for SFTS exposure risk are derived here:
#   - Pet ownership (dog and/or cat)
#   - Feeding community ("stray") cats
#   - Broader animal contact through pet shops, animal cafes, dog parks,
#     shelters, or caring for an acquaintance's pet
# Q7_* and Q9_* ask about the same set of contact settings but for two
# different reference groups in the survey (hence combined with "|" below).
pet_dataset <- pet_dataset %>%
  mutate(Dog_cat_owner_status = ifelse(Q1_1.1. == 1 | Q1_1.2. == 1, "Owns Dog or Cat", "Does Not Own Dog or Cat"),
         Community_cat_carer = ifelse(Q2 %in% 2:3, "Feeds community cats", "Does not feed community cats"),
         # Create dummy indicator variables (1 = Yes, 0 = No) for each interaction type
         Pet_owner = ifelse(Q1_1.1. == 1 | Q1_1.2. == 1, 1,0),
         Community_cat_carer = ifelse(Q2 %in% 2:3, 1,0),
         Pet_shops = ifelse(Q7_3 %in% 1:4 | Q9_3 %in% 1:4, 1,0),
         Caring_for_acquaintance_pets = ifelse(Q7_4 %in% 1:4 | Q9_4 %in% 1:4, 1,0),
         Animal_cafe = ifelse(Q7_5 %in% 1:4 | Q9_5 %in% 1:4, 1,0),
         Dog_park = ifelse(Q7_7 %in% 1:4, 1,0),
         Animal_shelter = ifelse(Q7_9 %in% 1:4 | Q9_8 %in% 1:4, 1,0),
         Lives_alone = ifelse(QW16.1. == 1, "Yes", "No")
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Label age and sex groups
# Recode demographic and behavioural frequency variables into readable labels.
pet_dataset <- pet_dataset %>%
  mutate(#AGE = case_when(
    #QAGE == 1 ~ "20-29",
    #QAGE == 2 ~ "30-39",
    #QAGE == 3 ~ "40-49",
    #QAGE == 4 ~ "50-59",
    #QAGE == 5 ~ "60-69",
    #QAGE == 6 ~ "70-79",
    #TRUE ~ NA_character_),
    Gender = case_when(
      F1 == 1 ~ "Male",
      F1 == 2 ~ "Female",
      TRUE ~ NA_character_),
    Feeding_frequency= case_when(
      Q3 == 1 ~ "5 or more times a week",
      Q3 == 2 ~ "2-4 times a week",
      Q3 == 3 ~ "Once a week",
      Q3 == 4 ~ "Once every two weeks",
      Q3 == 5 ~ "Once a month",
      TRUE ~ NA_character_),
    Occupation = case_when(
      QW5 %in% c(1,2,3,4,5) ~ "Regular_worker",
      QW5 %in% c(6) ~ "Part_time_worker",
      QW5 %in% c(7) ~ "Temporary",
      QW5 %in% c(8,9,10,11,12,13,14,15) ~ "Unemployed",
      TRUE ~ NA_character_),
    Q7_9_labeled = case_when(
      Q7_9 == 1 ~ "Twice a week or more",
      Q7_9 == 2 ~ "Once a week",
      Q7_9 == 3 ~ "2-3 times a month",
      Q7_9 == 4 ~ "Once a month",
      Q7_9 %in% c(5:8) ~ "Less than once a month",
      TRUE ~ NA_character_),
    Q9_8_labeled = case_when(
      Q9_8 == 1 ~ "Twice a week or more",
      Q9_8 == 2 ~ "Once a week",
      Q9_8 == 3 ~ "2-3 times a month",
      Q9_8 == 4 ~ "Once a month",
      Q9_8 %in% c(5:8) ~ "Less than once a month",
      TRUE ~ NA_character_),
    City_size = case_when(
      QW10 %in% c(1,2) ~ "Metropolitan",
      QW10 %in% c(3,4,5) ~ "City",
      QW10 == 6 ~ "Town",
      QW10 == 7 ~ "Countryside",
      TRUE ~ NA_character_),
    Income = case_when(
      QW15 %in% c(1,2,3) ~ "Lower_income",
      QW15 %in% c(4,5,6) ~ "Middle_income",
      QW15 %in% c(7,8,9,10) ~ "Upper_income",
      QW15 %in% c(11,12,13,14) ~ "Higher_income",
      TRUE ~ NA_character_),
    Children_household = case_when(
      QW6 == 1 ~ "No_children",
      QW6 %in% c(2,3) ~ "1_2_children",
      QW6 %in% c(4,5) ~ "3plus_children",
      TRUE ~ NA_character_),
    Housing_type = case_when(
      QW12 == 1 ~ "Owned_detached_house",
      QW12 == 2 ~ "Owned_condominium",
      QW12 == 3 ~ "Rented_detached_house",
      QW12 == 4 ~ "Rented_apartment_public_housing",
      QW12 == 5 ~ "Company_housing_dormitory",
      QW12 == 6 ~ "Other",
      TRUE ~ NA_character_)
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# create dataframe for Kyushu only
#pet_dataset_Kyushu <- pet_dataset %>%
# filter(Region == "Kyushu")
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Convert grouping variables to factors
# tbl_summary() (used throughout this script) displays factor labels
# directly, so recoding 0/1 to descriptive labels here saves relabelling
# work in every downstream table.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pet_dataset <- pet_dataset %>%
  mutate(
    Pet_owner = factor(
      Pet_owner,
      levels = c(0, 1),
      labels = c("Non-Pet Owners", "Pet Owners")
    ),
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    Caring_for_acquaintance_pets = factor(
      Caring_for_acquaintance_pets,
      levels = c(0, 1),
      labels = c(
        "Does Not Care for Acquaintance Pets",
        "Cares for Acquaintance Pets"
      )
    ),
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    Community_cat_carer = factor(
      Community_cat_carer,
      levels = c(0, 1),
      labels = c(
        "Does Not Feed Community Cats",
        "Feeds Community Cats"
      )
    )
  )
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#Select sociodemographic variables of interest and generate
# summary tables stratified by pet ownership status, acquaintance
# pet caregiving status, and community cat caregiving status.
# Mean (SD) is calculated for age, while frequencies (%) are
# reported for categorical variables. Independent sample t tests
# and chi square tests are performed to assess differences
# between groups.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#Descriptive statistics and inferential analysis for pet ownership
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
table_pet_owner <- pet_dataset %>%
  select(Pet_owner, F2, Gender, Lives_alone, Income, Occupation) %>%
  tbl_summary(
    by = Pet_owner,
    statistic = list(
      F2 ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no",
    label = list(
      F2 ~ "Age",
      Gender ~ "Sex",
      Lives_alone ~ "Lives alone",
      Income ~ "Income Category",
      Occupation ~ "Occupation"
    )
  ) %>%
  add_p(
    test = list(
      F2 ~ "t.test",
      all_categorical() ~ "chisq.test"
    )
  ) %>%
  modify_header(p.value ~ "**P-value**") %>%
  bold_labels()
table_pet_owner
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#save the table as docx####
#gtsave(as_gt(table_pet_owner),
       #"table_pet_owner_Table.docx")
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#Descriptive statistics and inferential analysis for acquaintance pets
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    #---------------------------------
table_acquaintance <- pet_dataset %>%
  select(
    Caring_for_acquaintance_pets,
    F2,
    Gender,
    Lives_alone,
    Income,
    Occupation
  ) %>%
  tbl_summary(
    by = Caring_for_acquaintance_pets,
    statistic = list(
      F2 ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no",
    label = list(
      F2 ~ "Age",
      Gender ~ "Sex",
      Lives_alone ~ "Lives alone",
      Income ~ "Income Category",
      Occupation ~ "Occupation"
    )
  ) %>%
  add_p(
    test = list(
      F2 ~ "t.test",
      all_categorical() ~ "chisq.test"
    )
  ) %>%
  modify_header(p.value ~ "**P-value**") %>%
  bold_labels()

table_acquaintance
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#save the table as docx####
#gtsave(as_gt(table_acquaintance),
     #  "table_acquaintance_Table.docx")
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Descriptive statistics and inferential analysis for Community Cat Carers
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
table_community <- pet_dataset %>%
  select(
    Community_cat_carer,
    F2,
    Gender,
    Lives_alone,
    Income,
    Occupation
  ) %>%
  tbl_summary(
    by = Community_cat_carer,
    statistic = list(
      F2 ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no",
    label = list(
      F2 ~ "Age",
      Gender ~ "Sex",
      Lives_alone ~ "Lives alone",
      Income ~ "Income Category",
      Occupation ~ "Occupation"
    )
  ) %>%
  add_p(
    test = list(
      F2 ~ "t.test",
      all_categorical() ~ "chisq.test"
    )
  ) %>%
  modify_header(p.value ~ "**P-value**") %>%
  bold_labels()

table_community
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#save the table as docx####
#gtsave(as_gt(table_community),
  #     "table_community.docx")
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Building baseline model for the whole population####
#predict community cat caring behaviour using demographic and household variables####
# Build modelling dataset ####
glm_model_pred <- pet_dataset%>%
  select(
    Pet_owner,
    Community_cat_carer,
    Caring_for_acquaintance_pets,
    Gender,
    Occupation,
    Lives_alone
  )%>%
  drop_na()
#----------------------------------------------------------------------------------------------------------------------------------
# CREATE BEHAVIOURAL PROFILES FROM SURVEY DATA
#
# The three pet-related behaviours are combined into a single categorical
# outcome representing all possible behavioural combinations observed in the
# survey population. These profiles correspond directly to the sections of the
# Venn diagram (Section 2) and preserve the observed overlap structure, so
# that projecting them onto the synthetic population later reproduces the
# same overlap patterns rather than treating the three behaviours as
# independent.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
glm_model_pred <- glm_model_pred %>%
  mutate(
    Behavior_Profile = case_when(
      Pet_owner == "Pet Owners" &
        Community_cat_carer == "Does Not Feed Community Cats" &
        Caring_for_acquaintance_pets == "Does Not Care for Acquaintance Pets" ~ "Own_Pet_Only",
      
      Pet_owner == "Non-Pet Owners" &
        Community_cat_carer == "Feeds Community Cats" &
        Caring_for_acquaintance_pets == "Does Not Care for Acquaintance Pets" ~ "Feed_Community_Cat_Only",
      
      Pet_owner == "Non-Pet Owners" &
        Community_cat_carer == "Does Not Feed Community Cats" &
        Caring_for_acquaintance_pets == "Cares for Acquaintance Pets" ~ "Care_for_Acquaintance_pet_Only",
      
      Pet_owner == "Pet Owners" &
        Community_cat_carer == "Feeds Community Cats" &
        Caring_for_acquaintance_pets == "Does Not Care for Acquaintance Pets" ~ "Own_Pet_and_Feed_Community_Cat",
      
      Pet_owner == "Pet Owners" &
        Community_cat_carer == "Does Not Feed Community Cats" &
        Caring_for_acquaintance_pets == "Cares for Acquaintance Pets" ~ "Own_Pet_and_Care_for_Acquaintance_Cat",
      
      Pet_owner == "Non-Pet Owners" &
        Community_cat_carer == "Feeds Community Cats" &
        Caring_for_acquaintance_pets == "Cares for Acquaintance Pets" ~ "Feed_Community_Cat_and_Care_for_Acquaintance_Cat",
      
      Pet_owner == "Pet Owners" &
        Community_cat_carer == "Feeds Community Cats" &
        Caring_for_acquaintance_pets == "Cares for Acquaintance Pets" ~ "All_Three_Possible_Intaraction_Pathways",
      TRUE ~ "None"
    ),
    Behavior_Profile = factor(Behavior_Profile)
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FIT MULTINOMIAL LOGISTIC REGRESSION MODEL
#
# Estimate the probability of belonging to each behavioural profile using
# demographic and household characteristics. Sex, living-alone status, and
# occupation are used as predictors because the preceding analysis showed
# them to be statistically significant, consistent, and relevant to the
# study aim. Because the model needs to be applied to the synthetic
# population later, these same three variables are subsequently derived
# from the raw synthetic population fields (see prepare_population()
# below) so that the fitted model can generate predictions for it.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
multinom_model <- nnet::multinom(
  Behavior_Profile ~ Gender + Lives_alone + Occupation,
  data = glm_model_pred,
  trace = FALSE
)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SHARED SYNTHETIC-POPULATION CLEANING FUNCTION
#
# Sections 2-6 all read the same 10 synthetic population files and repeat
# slightly different subsets of the same cleaning steps (household size,
# income banding, sex/occupation recoding, and — where needed — occupational
# industry recoding, multinomial behavioural-profile simulation, and risk
# subgroup flags). prepare_population() consolidates all of that into one
# function so it only has to be written, read, and debugged once.
#
# Each section switches on the same three toggles it always used implicitly:
#   filter_adults      — restrict to age >= 18 (Sections 4-6 only)
#   simulate_behavior  — re-code industry, run the multinomial model, and
#                         reconstruct Pet_owner / Community_cat_carers /
#                         Caring_for_acquaintance_pets (Sections 3-6)
#   compute_risk_flags — derive forest_agric / any_interaction /
#                         environmental_50m / older_adult (Sections 4-6 only;
#    
#
# IMPORTANT — CHECKED FOR CONSISTENCY: this single function correctly
# reproduces each of its three usage patterns — Section 2's basic cleaning
# only, Section 3's cleaning plus behaviour simulation, and Sections 4-6's
# full pipeline (including the age filter and risk flags). Re-running each
# pattern with the same random seed (1802) gives consistent, repeatable
# results, so the simulation results reported elsewhere in this project
# are reliable.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
prepare_population <- function(file, distance_data = NULL, filter_adults = FALSE,
                                simulate_behavior = FALSE, compute_risk_flags = FALSE) {

  pop_level <- read.csv(file)

  # Restrict to adults where required (Sections 4-6); Sections 2-3 keep all ages.
  if (filter_adults) pop_level <- pop_level %>% filter(age >= 18)

  # Attach environmental exposure (forest/agricultural proximity) where a
  # distance_data lookup table is supplied (Sections 4-6 only).
  if (!is.null(distance_data)) {
    pop_level <- pop_level %>% left_join(distance_data, by = "person_id")
  }

  # Household size and income banding — computed for every section (unused
  # columns are simply dropped later by each section's own transmute/select).
  pop_level <- pop_level %>%
    mutate(
      Lives_alone = if_else(num_member == 1, "Yes", "No"),
      annual_income = total_income * 12
    )

  pop_level <- pop_level %>%
    mutate(
      income_classification = case_when(
        annual_income < 1000000 ~ 1,
        annual_income < 2000000 ~ 2,
        annual_income < 3000000 ~ 3,
        annual_income < 4000000 ~ 4,
        annual_income < 5000000 ~ 5,
        annual_income < 6000000 ~ 6,
        annual_income < 7000000 ~ 7,
        annual_income < 8000000 ~ 8,
        annual_income < 9000000 ~ 9,
        annual_income < 10000000 ~ 10,
        annual_income < 12000000 ~ 11,
        annual_income < 15000000 ~ 12,
        annual_income < 20000000 ~ 13,
        annual_income >= 20000000 ~ 14,
        TRUE ~ NA_real_
      ),
      Income = case_when(
        income_classification %in% c(1, 2, 3) ~ "Lower_income",
        income_classification %in% c(4, 5, 6) ~ "Middle_income",
        income_classification %in% c(7, 8, 9, 10) ~ "Upper_income",
        income_classification %in% c(11, 12, 13, 14) ~ "Higher_income",
        TRUE ~ NA_character_
      )
    )

  # The synthetic files store sex and occupation in Japanese; recode to the
  # same English labels used in pet_dataset.
  pop_level <- pop_level %>%
    rename(
      Gender = gender,
      Occupation = employment_type
    ) %>%
    mutate(
      Gender = case_when(
        Gender == "男性" ~ "Male",
        Gender == "女性" ~ "Female",
        TRUE ~ NA_character_
      ),
      Occupation = case_when(
        Occupation == "" ~ "Unemployed",
        Occupation == "一般労働者" ~ "Regular_worker",
        Occupation == "短時間労働者" ~ "Part_time_worker",
        Occupation == "臨時労働者" ~ "Temporary",
        TRUE ~ NA_character_
      )
    )

  if (simulate_behavior) {

    # Translate the 19 official Japanese industry classifications into
    # English so "Agriculture_and_Forestry" can later be flagged as a risk
    # subgroup.
    pop_level <- pop_level %>%
      mutate(
        industry_type_eng = case_when(
          industry_type == "Ａ 農業，林業"                    ~ "Agriculture_and_Forestry",
          industry_type == "Ｂ 漁業"                        ~ "Fisheries",
          industry_type == "Ｃ 鉱業，採石業，砂利採取業"          ~ "Mining_and_Quarrying",
          industry_type == "Ｄ 建設業"                       ~ "Construction",
          industry_type == "Ｅ 製造業"                       ~ "Manufacturing",
          industry_type == "Ｆ 電気・ガス・熱供給・水道業"         ~ "Utilities",
          industry_type == "Ｇ 情報通信業"                    ~ "Information_and_Communications",
          industry_type == "Ｈ 運輸業，郵便業"                  ~ "Transport_and_Postal",
          industry_type == "Ｉ 卸売業，小売業"                  ~ "Wholesale_and_Retail",
          industry_type == "Ｊ 金融業，保険業"                  ~ "Finance_and_Insurance",
          industry_type == "Ｋ 不動産業，物品賃貸業"              ~ "Real_Estate_and_Rental",
          industry_type == "Ｌ 学術研究，専門・技術サービス業"       ~ "Professional_and_Technical",
          industry_type == "Ｍ 宿泊業，飲食サービス業"             ~ "Accommodation_and_Food",
          industry_type == "Ｎ 生活関連サービス業，娯楽業"          ~ "Personal_Services_and_Recreation",
          industry_type == "Ｏ 教育，学習支援業"                 ~ "Education",
          industry_type == "Ｐ 医療，福祉"                     ~ "Health_Care_and_Welfare",
          industry_type == "Ｑ 複合サービス事業"                 ~ "Composite_Services",
          industry_type == "Ｒ サービス業（他に分類されないもの）"    ~ "Other_Services",
          industry_type == "Ｓ 公務（他に分類されるものを除く）"      ~ "Public_Administration",
          industry_type == "Ｔ 分類不能の産業"                  ~ "Unclassified",
          TRUE ~ NA_character_
        )
      )

    # Predict behavioural profile probabilities for this synthetic population
    # using the multinomial model fitted on survey data (Section 1), then
    # sample one profile per person from their predicted probability
    # distribution.
    prob_matrix <- predict(multinom_model, newdata = pop_level, type = "probs")
    pop_level$Simulated_Profile <- apply(prob_matrix, 1, function(row_probs) {
      sample(colnames(prob_matrix), size = 1, prob = row_probs)
    })

    # Reconstruct the three binary behavioural variables from the assigned
    # profile (the inverse of the recoding done in Section 1.9).
    pop_level <- pop_level %>%
      mutate(
        Pet_owner = ifelse(Simulated_Profile %in% c(
          "Own_Pet_Only", "Own_Pet_and_Feed_Community_Cat",
          "Own_Pet_and_Care_for_Acquaintance_Cat", "All_Three_Possible_Intaraction_Pathways"
        ), 1, 0),
        Community_cat_carers = ifelse(Simulated_Profile %in% c(
          "Feed_Community_Cat_Only", "Own_Pet_and_Feed_Community_Cat",
          "Feed_Community_Cat_and_Care_for_Acquaintance_Cat", "All_Three_Possible_Intaraction_Pathways"
        ), 1, 0),
        Caring_for_acquaintance_pets = ifelse(Simulated_Profile %in% c(
          "Care_for_Acquaintance_pet_Only", "Own_Pet_and_Care_for_Acquaintance_Cat",
          "Feed_Community_Cat_and_Care_for_Acquaintance_Cat", "All_Three_Possible_Intaraction_Pathways"
        ), 1, 0)
      )
  }

  if (compute_risk_flags) {

    # Define the four risk dimensions used from Section 4 onward:
    # occupational (agriculture/forestry), behavioural (any animal contact),
    # environmental (living near forest/farmland), and demographic (age 65+).
    pop_level <- pop_level %>%
      mutate(
        forest_agric = (industry_type_eng == "Agriculture_and_Forestry"),
        any_interaction = (Pet_owner == 1 | Community_cat_carers == 1 | Caring_for_acquaintance_pets == 1),
        environmental_50m = (Environmental_50m == 1),
        older_adult = (age >= 65)
      ) %>%
      # Treat missing values as "not in this risk group" rather than unknown,
      # so nobody is silently dropped from subgroup counts.
      mutate(
        forest_agric = if_else(is.na(forest_agric), FALSE, forest_agric),
        environmental_50m = if_else(is.na(environmental_50m), FALSE, environmental_50m),
        any_interaction = if_else(is.na(any_interaction), FALSE, any_interaction),
        older_adult = if_else(is.na(older_adult), FALSE, older_adult)
      )
  }

  pop_level
}
#===================================================================================================================
# SECTION 2: VENN DIAGRAM VALIDATION
# Compare the demographic profile of the survey population against the
# synthetic population before any behavioural simulation is applied.
#===================================================================================================================
# DESCRIPTIVE COMPARISON OF SURVEY AND SYNTHETIC POPULATIONS
#
# Compare demographic characteristics between the observed survey population
# and the multinomial simulated synthetic population to assess whether the
# synthetic population reasonably represents the survey respondents.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All 10 synthetic population files share the "2015_" prefix, so list.files()
# with this pattern picks them all up without hard-coding file names.
files <- list.files(
  pattern = "^2015_.*\\.csv$",
  full.names = TRUE
)
synthetic_validation_list <- list()

#set random seed for reproducibility####
set.seed(1802)

survey_validation <- pet_dataset %>%
  transmute(
    Source = "Survey",
    Age = F2,
    Gender,
    Lives_alone,
    Occupation
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Loop over each of the 10 synthetic population files. Each iteration applies
# the same cleaning steps used on the survey data (income banding, Japanese
# -> English recoding of sex/occupation) so the two sources can be compared
# on a common set of variables.
for (file in files) {
  pop_level <- prepare_population(file)

  synthetic_validation_list[[basename(file)]] <- pop_level %>%
    transmute(
      Source = "Synthetic",
      Age = age,
      Gender,
      Lives_alone,
      Occupation
    )
}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Stack all 10 synthetic files into one data frame and combine with the survey
# data so both sources can be summarised side by side.
synthetic_validation <- bind_rows(synthetic_validation_list)

validation_data <- bind_rows(
  survey_validation,
  synthetic_validation
)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
validation_table <- validation_data %>%
  tbl_summary(
    by = Source,
    type = list(
      Age ~ "continuous",
      all_categorical() ~ "categorical"
    ),
    statistic = list(
      Age ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#display summary for both survey and synthetic data 
validation_table
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# save the table
#gtsave(as_gt(validation_table),"validation_table.docx")
#====================================================================================================================================================================================
# SECTION 3: DESCRIPTIVE COMPARISON — National Survey vs Synthetic Population
# Apply the fitted behavioural model to the synthetic population (all 10
# files), then validate the resulting overlap structure via Venn diagrams
# and descriptive comparison tables for Pet Ownership, Community Cat
# Feeding, and Acquaintance Pet Care.
#====================================================================================================================================================================================
# BUILD SYNTHETIC_DESC_ALL — combine all 10 synthetic population files
#
# This loop repeats the same file-preparation steps as Section 2 (income
# banding, sex/occupation recoding, industry translation), then goes one step
# further: it applies the fitted multinomial model to assign each synthetic
# individual a simulated behavioural profile, and reconstructs the three
# binary pet-related behaviours from that profile.
# Result: synthetic_desc_all contains all 10 files stacked into one data frame.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
synthetic_desc_list <- list()

#set random seed for reproducibility####
set.seed(1802)

for (file in files) {

  # Full demographic cleaning plus behavioural profile simulation, but no
  # age filter and no risk flags — this section only needs the three
  # pet-related behaviours alongside demographics.
  pop_level <- prepare_population(file, simulate_behavior = TRUE)

  # Keep only the columns needed for the descriptive comparison tables below.
  synthetic_desc_list[[basename(file)]] <- pop_level %>%
    transmute(
      Age = age,
      Gender,
      Lives_alone,
      Income,
      Occupation,
      Pet_owner,
      Community_cat_carers,
      Caring_for_acquaintance_pets
    )
}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Stack all 10 files into one data frame
synthetic_desc_all <- bind_rows(synthetic_desc_list)
message("synthetic_desc_all built: ", nrow(synthetic_desc_all), " rows from ", length(files), " files")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#descriptive stats for the two dataset(National Vs Synthetic Population)####

# NOTE ON SCOPE: the subgroup-by-subgroup comparison tables below (tbl_pet_*,
# tbl_comm_*, tbl_acq_*) were built as an additional diagnostic check while
# developing the model, but were not included in the dissertation — the main
# validation results reported there are Table S8 and Figure 1 (the Venn
# diagram comparison later in this section). They are left in for
# transparency and are still useful if anyone is interested in/wants to dig into a specific
# behaviour's demographic profile.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Descriptive stats for Pet_ownership for National Survey####
tbl_pet_survey <- pet_dataset %>%
  filter(!is.na(Income)) %>%
  mutate(
    Gender      = factor(Gender),
    Lives_alone = factor(Lives_alone),
    Income      = factor(Income),
    Occupation  = factor(Occupation)
  ) %>%
  select(
    Pet_owner,
    Age = F2,
    Gender,
    Lives_alone,
    Income,
    Occupation
  ) %>%
  tbl_summary(
    by = Pet_owner,
    missing = "no",
    type = list(all_categorical() ~ "categorical"),
    statistic = list(
      Age ~ "{mean} ({sd})"
    )
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Descriptive stats for Pet_ownership for Synthetic Population####
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tbl_pet_syn <- synthetic_desc_all %>%
  filter(!is.na(Income)) %>%
  mutate(
    Pet_owner = factor(
      Pet_owner,
      levels = c(0, 1),
      labels = c("Non-Pet Owners", "Pet Owners")
    ),
    Gender      = factor(Gender),
    Lives_alone = factor(Lives_alone),
    Income      = factor(Income),
    Occupation  = factor(Occupation)
  ) %>%
  select(
    Pet_owner,
    Age,
    Gender,
    Lives_alone,
    Income,
    Occupation
  ) %>%
  tbl_summary(
    by = Pet_owner,
    missing = "no",
    type = list(all_categorical() ~ "categorical"),
    statistic = list(
      Age ~ "{mean} ({sd})"
    )
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tbl_compare_pet <- tbl_merge(
  tbls = list(
    tbl_pet_survey,
    tbl_pet_syn
  ),
  tab_spanner = c(
    "**National Survey**",
    "**Synthetic Population**"
  )
)
#tbl_compare_pet %>%
 # as_gt() %>%
  #gt::gtsave("compare_pet.docx")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Descriptive stats for Community_cat_carer for National Survey####
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tbl_comm_survey <-pet_dataset %>%
  filter(!is.na(Income)) %>%
  mutate(
    Gender      = factor(Gender),
    Lives_alone = factor(Lives_alone),
    Income      = factor(Income),
    Occupation  = factor(Occupation)
  ) %>%
  select(
    Community_cat_carer,
    Age=F2,
    Gender,
    Lives_alone,
    Income,
    Occupation
  ) %>%
  tbl_summary(
    by = Community_cat_carer,
    missing = "no",
    type = list(all_categorical() ~ "categorical"),
    statistic = list(
      Age ~ "{mean} ({sd})"
    )
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Descriptive stats for Community_cat_carer for Synthetic Population####
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tbl_comm_syn <- synthetic_desc_all %>%
  filter(!is.na(Income)) %>%
  mutate(
    Community_cat_carers = factor(
      Community_cat_carers,
      levels = c(0, 1),
      labels = c("Does Not Feed Community Cats", "Feeds Community Cats")
    ),
    Gender      = factor(Gender),
    Lives_alone = factor(Lives_alone),
    Income      = factor(Income),
    Occupation  = factor(Occupation)
  ) %>%
  select(
    Community_cat_carers,
    Age,
    Gender,
    Lives_alone,
    Income,
    Occupation
  ) %>%
  tbl_summary(
    by = Community_cat_carers,
    missing = "no",
    type = list(all_categorical() ~ "categorical"),
    statistic = list(
      Age ~ "{mean} ({sd})"
    )
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tbl_compare_comm <- tbl_merge(
  tbls = list(
    tbl_comm_survey,
    tbl_comm_syn
  ),
  tab_spanner = c(
    "**National Survey**",
    "**Synthetic Population**"
  )
)
#tbl_compare_comm %>%
 # as_gt() %>%
  #gt::gtsave("compare_community.docx")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Descriptive stats for Caring_for_acquaintance_pets for National Survey####
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tbl_acq_survey <-pet_dataset %>%
  filter(!is.na(Income)) %>%
  mutate(
    Gender      = factor(Gender),
    Lives_alone = factor(Lives_alone),
    Income      = factor(Income),
    Occupation  = factor(Occupation)
  ) %>%
  select(
    Caring_for_acquaintance_pets,
    Age=F2,
    Gender,
    Lives_alone,
    Income,
    Occupation
  ) %>%
  tbl_summary(
    by = Caring_for_acquaintance_pets,
    missing = "no",
    type = list(all_categorical() ~ "categorical"),
    statistic = list(
      Age ~ "{mean} ({sd})"
    )
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Descriptive stats for Caring_for_acquaintance_pets for Synthetic Population####
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tbl_acq_syn <- synthetic_desc_all %>%
  filter(!is.na(Income)) %>%
  mutate(
    Caring_for_acquaintance_pet = factor(
      Caring_for_acquaintance_pets,
      levels = c(0, 1),
      labels = c("Non-Carers", "Carers")
    ),
    Gender      = factor(Gender),
    Lives_alone = factor(Lives_alone),
    Income      = factor(Income),
    Occupation  = factor(Occupation)
  ) %>%
  select(
    Caring_for_acquaintance_pet,
    Age,
    Gender,
    Lives_alone,
    Income,
    Occupation
  ) %>%
  tbl_summary(
    by = Caring_for_acquaintance_pet,
    missing = "no",
    type = list(all_categorical() ~ "categorical"),
    statistic = list(
      Age ~ "{mean} ({sd})"
    )
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tbl_compare_acq <- tbl_merge(
  tbls = list(
    tbl_acq_survey,
    tbl_acq_syn
  ),
  tab_spanner = c(
    "**National Survey**",
    "**Synthetic Population**"
  )
)
#tbl_compare_acq %>%
 # as_gt() %>%
  #gt::gtsave("compare_acquaintance.docx")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# VALIDATE SIMULATED PREVALENCE
#
# Compare the prevalence of each simulated behaviour against the national
# survey population to assess whether the multinomial model reproduces the
# observed behavioural frequencies. NB: pop_level here still holds the last
# file processed by the Section 3 loop above, so this is a single-file
# spot-check rather than a comparison across all 10 synthetic populations.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
prop.table(table(pet_dataset$Pet_owner))
prop.table(table(pop_level$Pet_owner))

prop.table(table(pet_dataset$Community_cat_carer))
prop.table(table(pop_level$Community_cat_carers))

prop.table(table(pet_dataset$Caring_for_acquaintance_pets))
prop.table(table(pop_level$Caring_for_acquaintance_pets))
#-------------------------------------------------------------------------------------------------------------------------------
# VALIDATE BEHAVIOURAL OVERLAP STRUCTURE
#
# Compare overlap patterns between the national survey and the multinomial
# simulation using Venn diagrams. This is the key validation check for this
# script: it shows whether the simulation reproduces the observed
# relationships between pet-related behaviours (e.g. do pet owners in the
# synthetic population also feed community cats at roughly the same rate as
# in the real survey?), not just the marginal prevalence of each behaviour
# individually.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#names(pet_dataset)
# National Survey List
venn_national <- list(
  Pet_owner = which(pet_dataset$Pet_owner == "Pet Owners"),
  Community_cat = which(pet_dataset$Community_cat_carer == "Feeds Community Cats"),
  "Caring for Others' Pets" = which(pet_dataset$Caring_for_acquaintance_pets == "Cares for Acquaintance Pets")
)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Kyushu Survey List
#venn_kyushu <- list(
# Pet_owner = which(pet_dataset_Kyushu$Pet_owner == 1),
# Community_cat_carer = which(pet_dataset_Kyushu$Community_cat_carer == 1),
#Caring_for_acquaintance_pets = which(pet_dataset_Kyushu$Caring_for_acquaintance_pets == 1)
#)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Synthetic Simulation List (based on the last processed file — pop_level)
venn_simulated <- list(
  Pet_Owners = which(pop_level$Pet_owner == 1),
  Community_cat= which(pop_level$Community_cat_carers == 1),
  "Caring for Others' Pets" = which(pop_level$Caring_for_acquaintance_pets == 1)
)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Generate individual plots
#p1 <- ggVennDiagram(venn_national) + ggtitle("National Survey")
#p2 <- ggVennDiagram(venn_kyushu) + ggtitle("Kyushu Survey")
#p3 <- ggVennDiagram(venn_simulated) + ggtitle("Simulated Population")

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Generate Venn diagrams for the national survey and the simulated population
# so the two overlap structures can be inspected side by side (Figure 1)####
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p1 <- ggVennDiagram(
  venn_national,
  label_alpha = 0
) +
  scale_fill_gradient(
    low = "#D6EAF8",
    high = "#2E86C1"
  ) +
  ggtitle("National Survey") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    text = element_text(size = 16)
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p3 <- ggVennDiagram(
  venn_simulated,
  label_alpha = 0
) +
  scale_fill_gradient(
    low = "#D6EAF8",
    high = "#2E86C1"
  ) +
  ggtitle("Simulated Population") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    text = element_text(size = 16)
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Arrange them side-by-side cleanly
venn_plot <- (p1 | p3) 
print(venn_plot)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#save the venn diagram output####
#ggsave(
 # "venn_plot.png",
 # venn_plot,
  #width = 16,
  #height = 8,
  #dpi = 300
#)
#comparison
#===================================================================================================================
# SECTION 4: SUBGROUP COUNTING
# Quantify the size of each risk subgroup across all 10 synthetic population
# files. Reports: mean, median, SD, 95% interval for each subgroup.
#===================================================================================================================
# Prepare environmental exposure variables from geospatial proximity data
#
# Classify each individual according to whether they live within 50 metres of
# a forested area or agricultural land. Because the available seroprevalence
# evidence provides a single estimate for both forest and agricultural
# environments, these two proximity measures are combined into one
# environmental exposure variable.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#read in the environmental data####
distance_data <- read.csv("Forest_Ag_Pop_2022.csv")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
distance_data <- distance_data %>%
  mutate(
    
    # Identify individuals living within 50 metres of a forested area.
    # -1 is the dataset's code for "no distance recorded", which is treated
    # as not-exposed rather than missing.
    Forest_50m = case_when(
      Forest_Dis >= 0 & Forest_Dis < 50 ~ 1,
      Forest_Dis == -1 ~ 0,
      Forest_Dis >= 50 ~ 0,
      TRUE ~ NA_real_
    ),
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    # Identify individuals living within 50 metres of agricultural land
    Agri_50m = case_when(
      Agri_Dista >= 0 & Agri_Dista < 50 ~ 1,
      Agri_Dista == -1 ~ 0,
      Agri_Dista >= 50 ~ 0,
      TRUE ~ NA_real_
    ),
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    # Combine both proximity measures into a single environmental exposure
    # variable. Individuals living within 50 metres of either a forested area
    # or agricultural land are classified as environmentally exposed.
    Environmental_50m = if_else(
      Forest_50m == 1 | Agri_50m == 1,
      1,
      0
    )
  ) %>%
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
# Retain only the variables required for merging with each synthetic
# population dataset
select(person_id, Environmental_50m) %>%
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
# Ensure each individual contributes only one environmental exposure record
distinct(person_id, .keep_all = TRUE)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create a file path to read all synthetic population datasets####
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
files <- list.files(
  pattern = "^2015_.*\\.csv$",
  full.names = TRUE
)
subgroup_results <- data.frame()
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#set random seed for reproducibility####
set.seed(1802)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
# Loop through each synthetic population dataset and estimate risk subgroup counts####
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
for(file in files) {

  # Full pipeline: age filter, environmental join, demographic cleaning,
  # behavioural profile simulation, and the four risk-subgroup flags used in
  # the subgroup counts below.
  pop_level <- prepare_population(
    file,
    distance_data = distance_data,
    filter_adults = TRUE,
    simulate_behavior = TRUE,
    compute_risk_flags = TRUE
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  #Infection subgroup### 
  #Baseline 
  young_total<- sum(!pop_level$older_adult, na.rm = TRUE)
  
  # Category 1: Environmental Risk (Young Farmers)
  young_agriculture <- sum(!pop_level$older_adult & pop_level$forest_agric, na.rm = TRUE)
  
  # Category 2: Behavioral Risk (Young Animal Contacts)
  # Uses any_interaction which includes pet_onwerscommunity cats and acquaintance pets
  young_interaction <- sum(!pop_level$older_adult & pop_level$any_interaction, na.rm = TRUE)
  
  # 3. Category 3: Environmental exposure
  young_environment <- sum(!pop_level$older_adult & pop_level$environmental_50m,na.rm = TRUE)
  
  # 4. Category 4: Behavioural + Environmental exposure
  young_interaction_environment <- sum(!pop_level$older_adult &pop_level$any_interaction &pop_level$environmental_50m,
                                       na.rm = TRUE)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  # morbidity/mortality subgroup 
  #Baseline 
  older_total <- sum(pop_level$older_adult, na.rm = TRUE)
  
  # Category 1: Environmental Risk (Elderly Farmers)
  old_agriculture<- sum(pop_level$older_adult & pop_level$forest_agric, na.rm = TRUE)
  
  # Category 2: Behavioral Risk (Elderly Animal Contacts)
  old_interaction <- sum(pop_level$older_adult & pop_level$any_interaction, na.rm = TRUE)
  
  # 3. Category 3: Environmental exposure
  old_environment <- sum(pop_level$older_adult &pop_level$environmental_50m,na.rm = TRUE)
  
  # 4. Category 4: Behavioural + Environmental exposure
  old_interaction_environment <- sum(pop_level$older_adult &pop_level$any_interaction &pop_level$environmental_50m,
                                     na.rm = TRUE)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  #count each subgroup (printed here for a quick per-file sanity check while
  # the loop runs; the values that matter are stored below in subgroup_results)####
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  young_total
  young_agriculture
  young_interaction
  young_environment
  young_interaction_environment
  
  older_total
  old_agriculture
  old_interaction
  old_environment 
  old_interaction_environment
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  # Append this file's subgroup counts as one row, so that after the loop
  # finishes we have 10 rows (one per synthetic population) to summarise.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  # --------------------------------------------------
  subgroup_results <- rbind(
    subgroup_results,
    data.frame(
      
      Adults_Total = young_total,
      
      Adults_Agriculture_Forestry = young_agriculture,
      
      Adults_Animal_Interaction = young_interaction,
      
      Adults_Environmental_50m = young_environment,
      
      Adults_Animal_Interaction_Environmental = young_interaction_environment,
      
      Older_Total = older_total,
      
      Older_Agriculture_Forestry = old_agriculture,
      
      Older_Animal_Interaction = old_interaction,
      
      Older_Environmental_50m = old_environment,
      
      Older_Animal_Interaction_Environmental = old_interaction_environment
    )
  )
}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
# Summarise uncertainty across synthetic populations
#
# Calculate summary statistics for each risk subgroup across all synthetic
# populations. Mean, median, standard deviation, and 95% uncertainty
# intervals are reported to quantify variation arising from the use of
# multiple synthetic datasets.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
summary_stats <- data.frame(
  Subgroup = names(subgroup_results),
  Mean = sapply(subgroup_results, mean, na.rm = TRUE),
  Median = sapply(subgroup_results, median, na.rm = TRUE),
  SD = sapply(subgroup_results, sd, na.rm = TRUE),
  Lower95 = sapply(subgroup_results, quantile, probs = 0.025, na.rm = TRUE),
  Upper95 = sapply(subgroup_results, quantile, probs = 0.975, na.rm = TRUE)
)
rownames(summary_stats) <- NULL
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
subgroup_summary <- summary_stats %>%
  gt() %>%
  cols_label(
    Subgroup = "Subgroup",
    Mean = "Mean",
    Median = "Median",
    SD = "SD",
    Lower95 = "Lower 95%",
    Upper95 = "Upper 95%"
  )
print(subgroup_summary)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
# Inspect outputs
#
# Display the structure of the compiled risk-group results and the summary
# statistics table to verify that all synthetic populations were processed
# successfully and that uncertainty estimates were generated correctly.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
#gtsave(subgroup_summary, "subgroup_summary2.docx")
#========================================================================================================================================================
# SECTION 5: VACCINATION SCENARIO ANALYSIS
# Stochastic simulation of cases and deaths averted under hypothetical
# vaccination. Observed incidence (÷13) and worst-case incidence (full rate),
# 25/50/75% coverage. Results reported as Median [IQR] across 1,000 runs.
#=======================================================================================================================================================
# Vaccination scenario analysis
#
# Estimate the potential number of SFTS cases and deaths averted
# under hypothetical vaccination scenarios among identified risk
# groups. Simulations are performed across all synthetic population
# datasets assuming 80% vaccine efficacy derived from published
# literature and vaccination coverage levels of 25%, 50%, and 75%.
# Results are used to quantify the potential public health impact
# of targeted vaccination strategies.

# SFTS VACCINATION SCENARIO ANALYSIS
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PURPOSE:
# This section estimates how many SFTS infections and deaths could be prevented
# under hypothetical vaccination programmes targeting identified high-risk groups.
#
# We run a stochastic (random) simulation 1000 times (10 synthetic population
# files x 100 runs each) to capture uncertainty. Results are reported as
# Median [IQR] across all 1000 runs.
#
# FOUR TABLES ARE PRODUCED:
#   Table 1 — Adults 18-64    | Observed incidence (÷13 years)
#   Table 2 — Adults 18-64    | Worst-case incidence (full rate, no division)
#   Table 3 — Older Adults 65+| Observed incidence (÷13 years)
#   Table 4 — Older Adults 65+| Worst-case incidence (full rate, no division)
# 
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# STEP 1: SET ALL FIXED PARAMETERS
# These values are fixed assumptions that do not change across simulations.
# All are based on published literature or SFTS epidemiological study design decisions.
# NB:If any assumption changes, update only this section — the rest adjusts automatically.
# The third steps is left for reader that wants to explore
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# set.seed() ensures the random results are identical every time this script
# is re-run — this is what makes the vaccination impact estimates reproducible.
set.seed(1802)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Vaccine effectiveness: probability that the vaccine successfully prevents infection
# Fixed at 80%: this was based on assumption of what the efficacy of the vaccine could be
#if available
vaccine_efficacy <- 0.80
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Total observation period: SFTS has been notifiable in Japan since 2013
# 2026 - 2013 = 13 years of accumulated seroincidence data
years_accumulated <- 13
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Cumulative seroincidence over 13 years (proportion of each group
# who showed serological evidence of past SFTS infection)
# These come from published seroepidemiological studies
baseline_seroprevalence           <- 0.003  # 0.3%  — Healthy individuals living in SFTS-endemic regions of Japan
agriculture_seroprevalence        <- 0.002  # 0.2%  — Agriculture & forestry workers to be changed to 0.042
animal_interaction_seroprevalence <- 0.025  # 2.5%  — People with animal contact
environmental_seroprevalence      <- 0.030  # 3.0%  — Individuals living within 50 m of a forest or agricultural area
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ── OBSERVED SCENARIO RATES (annual) 
# We divide the 13-year cumulative seroincidence by 13 to obtain
# the average annual infection probability for each group.
# This represents the typical infection risk in any given year based on history.
agr_annual_obs <- agriculture_seroprevalence        / years_accumulated  # 0.002 ÷ 13
ani_annual_obs <- animal_interaction_seroprevalence / years_accumulated  # 0.025 ÷ 13
bas_annual_obs <- baseline_seroprevalence           / years_accumulated  # 0.003 ÷ 13 
envi_annual_obs<-environmental_seroprevalence    /  years_accumulated  # 0.030 ÷ 13 
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ── WORST-CASE SCENARIO RATES (no division) ───────────────────────────────────
# Here we assume the entire 13-year cumulative seroincidence could occur
# within a single year. For example, during an outbreak or period of
# intensified tick activity. This is a stress-test of the vaccination programme
# and represents the upper bound of plausible annual infection risk.
agr_annual_wc <- agriculture_seroprevalence         # 0.002 — used as a direct annual rate
ani_annual_wc <- animal_interaction_seroprevalence  # 0.025 — used as a direct annual rate
bas_annual_wc <- baseline_seroprevalence            # 0.003 — used as a direct annual rate
envi_annual_wc<-environmental_seroprevalence        # 0.030 — used as a direct annual rate
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Case fatality rate (CFR): probability of dying given SFTS infection
#5% for adults 18-64 (working-age adults have better outcomes)
#20% for older adults 65+ (older adults face higher mortality due to SFTS)

cfr_young<-0.05
cfr_older<-0.20
# Proportion of infected people who develop symptoms severe enough to be detected/reported
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTANT MODELLING CORRECTION:
# The CFR values above are derived from REPORTED (i.e. clinically confirmed,
# symptomatic) cases only. Most SFTS infections are thought to be mild or
# asymptomatic and are never detected, so applying a hospital-derived CFR to
# every infected person would substantially overstate the number of deaths.
# To correct for this, CFR is only applied to the symptomatic fraction of
# infections (governed by reporting_fraction below); asymptomatic infections
# are assumed to carry a 0% risk of death in this model.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Fraction of infected cases that are detected / symptomatic enough to enter the CFR step
#reporting_fraction_grid <- c(0.04, 0.08, 0.12, 0.20)
reporting_fraction <- 0.08
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# STEP 2: OBSERVED INCIDENCE SCENARIO SIMULATION
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# We loop over each of our 10 synthetic population files (outer loop).
# For each file, we run the simulation 100 times (inner loop).
# 10 files × 100 runs = 1000 total simulation runs.
#
# WHY 1000 RUNS?
# Each run uses rbinom() to randomly assign infections, deaths, and vaccination
# outcomes. Because these are random, each run gives slightly different numbers.
# Running 1000 times gives us a distribution of outcomes, from which we report
# the Median (most typical result) and IQR (how much results vary).
#
# ANNUAL INFECTION PROBABILITY IN THIS SCENARIO:
#   Agriculture & Forestry workers: 0.002 ÷ 13 per year
#   Animal Interaction group:       0.025 ÷ 13 per year
#   General population (baseline):  0.003 ÷ 13 per year
#   Enviromental_50m:               0.030 ÷ 13 per year
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#run_model <- function(reporting_fraction) {  

obs_runs <- data.frame()  # empty container — one row is added per simulation run
  
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #read in the environmental data####
  distance_data <- read.csv("Forest_Ag_Pop_2022.csv")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  distance_data <- distance_data %>%
    mutate(
      
      # Identify individuals living within 50 metres of a forested area
      Forest_50m = case_when(
        Forest_Dis >= 0 & Forest_Dis < 50 ~ 1,
        Forest_Dis == -1 ~ 0,
        Forest_Dis >= 50 ~ 0,
        TRUE ~ NA_real_
      ),
      #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
      # Identify individuals living within 50 metres of agricultural land
      Agri_50m = case_when(
        Agri_Dista >= 0 & Agri_Dista < 50 ~ 1,
        Agri_Dista == -1 ~ 0,
        Agri_Dista >= 50 ~ 0,
        TRUE ~ NA_real_
      ),
      #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
      # Combine both proximity measures into a single environmental exposure
      # variable. Individuals living within 50 metres of either a forested area
      # or agricultural land are classified as environmentally exposed.
      Environmental_50m = if_else(
        Forest_50m == 1 | Agri_50m == 1,
        1,
        0
      )
    ) %>%
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  # Retain only the variables required for merging with each synthetic
  # population dataset
  select(person_id, Environmental_50m) %>%
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  # Ensure each individual contributes only one environmental exposure record
  distinct(person_id, .keep_all = TRUE)
  #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Create a file path to read all synthetic population datasets####
  #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  files <- list.files(
    pattern = "^2015_.*\\.csv$",
    full.names = TRUE
  )
  #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #set random seed for reproducibility####
  set.seed(1802)
  #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  # Loop through each synthetic population dataset, preparing it exactly as in
  # Sections 3-4, before running the observed-incidence simulation on it####
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  for(file in files) {

    # Full pipeline (identical toggles to Section 4): age filter, environmental
    # join, demographic cleaning, behavioural simulation, and risk-flag
    # construction, ready for the observed-incidence simulation below.
    pop_level <- prepare_population(
      file,
      distance_data = distance_data,
      filter_adults = TRUE,
      simulate_behavior = TRUE,
      compute_risk_flags = TRUE
    )
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # RUN SIMULATION 100 TIMES ON THIS POPULATION 
    for (run in 1:100) {
      
      # ASSIGN ANNUAL INFECTION PROBABILITY
      # Every individual initially receives the baseline annual infection probability.
      # Individuals belonging to predefined high-risk subgroups are then assigned the
      # subgroup-specific annual probability derived from published seroprevalence
      # estimates. Where an individual belongs to more than one subgroup, the order of
      # assignment determines which probability is retained.
      pop_level$annual_p <- bas_annual_obs                                    # Baseline annual infection probability
      pop_level$annual_p[pop_level$forest_agric]      <- agr_annual_obs        # Agriculture & forestry workers
      pop_level$annual_p[pop_level$any_interaction]   <- ani_annual_obs        # Individuals with animal interaction
      pop_level$annual_p[pop_level$environmental_50m] <- envi_annual_obs       # Individuals living within 50 m of forest/agricultural area
      
      # COIN FLIP 1 — WHO GETS INFECTED THIS YEAR?
      # rbinom(n, 1, p) flips one coin per person with probability p.
      # Result is 1 (infected) or 0 (not infected).
      pop_level$infected <- rbinom(nrow(pop_level), 1, pop_level$annual_p)
      
      # COIN FLIP 2a — OF THOSE INFECTED, WHO DEVELOPS SYMPTOMS?
      # Most SFTS infections are mild or asymptomatic and never get reported.
      # The CFR from the literature applies only to symptomatic (reported) cases.
      # We first decide who is symptomatic (reporting_fraction), then apply CFR only to those people.
      pop_level$cfr_p          <- ifelse(pop_level$older_adult, cfr_older, cfr_young)
      pop_level$is_symptomatic <- 0L
      idx                      <- which(pop_level$infected == 1)
      pop_level$is_symptomatic[idx] <- rbinom(length(idx), 1, reporting_fraction)
      
      # COIN FLIP 2b — OF THOSE SYMPTOMATIC, WHO WOULD HAVE DIED?
      # would_die starts at 0 for everyone.
      # Only symptomatic people get the death coin flip — asymptomatic cases have 0% chance of dying.
      pop_level$would_die <- 0L
      sym_idx <- which(pop_level$is_symptomatic == 1)
      pop_level$would_die[sym_idx] <- rbinom(length(sym_idx), 1, pop_level$cfr_p[sym_idx])
      
      # COIN FLIP 3 — DOES THE VACCINE WORK FOR THIS PERSON?
      # Each person has an 80% chance that the vaccine is effective.
      # We draw this once and use it across all three coverage scenarios,
      # so vaccine effectiveness is consistent across coverage levels.
      pop_level$ve_works <- rbinom(nrow(pop_level), 1, vaccine_efficacy)
      
      # COIN FLIP 4 — IS EACH PERSON VACCINATED UNDER EACH COVERAGE SCENARIO?
      # Under 25% coverage, each person has a 25% chance of being vaccinated.
      # Under 50% coverage, each person has a 50% chance. And so on.
      # These are independent draws — a person could be vaccinated in one scenario
      # but not in another.
      pop_level$vac_25 <- rbinom(nrow(pop_level), 1, 0.25)
      pop_level$vac_50 <- rbinom(nrow(pop_level), 1, 0.50)
      pop_level$vac_75 <- rbinom(nrow(pop_level), 1, 0.75)
      
      # CALCULATE CASES AVERTED
      # A case is averted if ALL THREE conditions are TRUE:
      #   (1) the person was infected (they were at risk)
      #   (2) the person was vaccinated under this coverage scenario
      #   (3) the vaccine worked for this person (80% chance)
      pop_level$case_s25 <- pop_level$infected == 1 & pop_level$vac_25 == 1 & pop_level$ve_works == 1
      pop_level$case_s50 <- pop_level$infected == 1 & pop_level$vac_50 == 1 & pop_level$ve_works == 1
      pop_level$case_s75 <- pop_level$infected == 1 & pop_level$vac_75 == 1 & pop_level$ve_works == 1
      
      # CALCULATE DEATHS AVERTED
      # A death is averted if the person would have died AND their case was averted.
      # In other words: they were a fatal case that the vaccine prevented.
      pop_level$death_s25 <- pop_level$would_die == 1 & pop_level$vac_25 == 1 & pop_level$ve_works == 1
      pop_level$death_s50 <- pop_level$would_die == 1 & pop_level$vac_50 == 1 & pop_level$ve_works == 1
      pop_level$death_s75 <- pop_level$would_die == 1 & pop_level$vac_75 == 1 & pop_level$ve_works == 1
      
      # FILTER TO EACH SUBGROUP
      # We split the population into 4 subsets so we can count outcomes separately.
      # y_ = younger adults (18-64), o_ = older adults (65+)
      # agr = agriculture & forestry, ani = animal interaction
      # FILTER TO EACH SUBGROUP
      young_agr     <- pop_level[!pop_level$older_adult & pop_level$forest_agric, ]
      young_ani     <- pop_level[!pop_level$older_adult & pop_level$any_interaction, ]
      young_env     <- pop_level[!pop_level$older_adult & pop_level$environmental_50m, ]
      young_int_env <- pop_level[!pop_level$older_adult & pop_level$any_interaction & pop_level$environmental_50m, ]
      
      older_agr     <- pop_level[ pop_level$older_adult & pop_level$forest_agric, ]
      older_ani     <- pop_level[ pop_level$older_adult & pop_level$any_interaction, ]
      older_env     <- pop_level[ pop_level$older_adult & pop_level$environmental_50m, ]
      older_int_env <- pop_level[ pop_level$older_adult & pop_level$any_interaction & pop_level$environmental_50m, ]
      
      # COUNT AND STORE RESULTS FOR THIS RUN
      # We count outcomes in each subgroup and store as one row.
      # Column naming convention:
      #   y_ / o_  = younger / older adults
      #   agr / ani = agriculture / animal interaction subgroup
      #   _cases / _deaths = pre-vaccination burden
      #   _c25 / _d25 etc. = cases/deaths averted at 25%/50%/75% coverage
      one_run <- data.frame(
        # Adults 18-64: Agriculture & Forestry
        y_agr_cases  = sum(young_agr$infected),
        y_agr_deaths = sum(young_agr$would_die),
        y_agr_c25 = sum(young_agr$case_s25),  y_agr_d25 = sum(young_agr$death_s25),
        y_agr_c50 = sum(young_agr$case_s50),  y_agr_d50 = sum(young_agr$death_s50),
        y_agr_c75 = sum(young_agr$case_s75),  y_agr_d75 = sum(young_agr$death_s75),
        
        # Adults 18-64: Animal Interaction
        y_ani_cases  = sum(young_ani$infected),
        y_ani_deaths = sum(young_ani$would_die),
        y_ani_c25 = sum(young_ani$case_s25),  y_ani_d25 = sum(young_ani$death_s25),
        y_ani_c50 = sum(young_ani$case_s50),  y_ani_d50 = sum(young_ani$death_s50),
        y_ani_c75 = sum(young_ani$case_s75),  y_ani_d75 = sum(young_ani$death_s75),
        
        # Adults 18-64: Environmental 50m
        y_env_cases  = sum(young_env$infected),
        y_env_deaths = sum(young_env$would_die),
        y_env_c25 = sum(young_env$case_s25),  y_env_d25 = sum(young_env$death_s25),
        y_env_c50 = sum(young_env$case_s50),  y_env_d50 = sum(young_env$death_s50),
        y_env_c75 = sum(young_env$case_s75),  y_env_d75 = sum(young_env$death_s75),
        
        # Adults 18-64: Animal Interaction + Environmental 50m
        y_int_env_cases  = sum(young_int_env$infected),
        y_int_env_deaths = sum(young_int_env$would_die),
        y_int_env_c25 = sum(young_int_env$case_s25),  y_int_env_d25 = sum(young_int_env$death_s25),
        y_int_env_c50 = sum(young_int_env$case_s50),  y_int_env_d50 = sum(young_int_env$death_s50),
        y_int_env_c75 = sum(young_int_env$case_s75),  y_int_env_d75 = sum(young_int_env$death_s75),
        
        # Older Adults 65+: Agriculture & Forestry
        o_agr_cases  = sum(older_agr$infected),
        o_agr_deaths = sum(older_agr$would_die),
        o_agr_c25 = sum(older_agr$case_s25),  o_agr_d25 = sum(older_agr$death_s25),
        o_agr_c50 = sum(older_agr$case_s50),  o_agr_d50 = sum(older_agr$death_s50),
        o_agr_c75 = sum(older_agr$case_s75),  o_agr_d75 = sum(older_agr$death_s75),
        
        # Older Adults 65+: Animal Interaction
        o_ani_cases  = sum(older_ani$infected),
        o_ani_deaths = sum(older_ani$would_die),
        o_ani_c25 = sum(older_ani$case_s25),  o_ani_d25 = sum(older_ani$death_s25),
        o_ani_c50 = sum(older_ani$case_s50),  o_ani_d50 = sum(older_ani$death_s50),
        o_ani_c75 = sum(older_ani$case_s75),  o_ani_d75 = sum(older_ani$death_s75),
        
        # Older Adults 65+: Environmental 50m
        o_env_cases  = sum(older_env$infected),
        o_env_deaths = sum(older_env$would_die),
        o_env_c25 = sum(older_env$case_s25),  o_env_d25 = sum(older_env$death_s25),
        o_env_c50 = sum(older_env$case_s50),  o_env_d50 = sum(older_env$death_s50),
        o_env_c75 = sum(older_env$case_s75),  o_env_d75 = sum(older_env$death_s75),
        
        # Older Adults 65+: Animal Interaction + Environmental 50m
        o_int_env_cases  = sum(older_int_env$infected),
        o_int_env_deaths = sum(older_int_env$would_die),
        o_int_env_c25 = sum(older_int_env$case_s25),  o_int_env_d25 = sum(older_int_env$death_s25),
        o_int_env_c50 = sum(older_int_env$case_s50),  o_int_env_d50 = sum(older_int_env$death_s50),
        o_int_env_c75 = sum(older_int_env$case_s75),  o_int_env_d75 = sum(older_int_env$death_s75)
      )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
      obs_runs <- rbind(obs_runs, one_run)  # append this run's results to the storage
      
    }  # end inner loop — next of 100 runs for this file
    
  }  # end outer loop — next of 10 synthetic population files
  
  message("Observed scenario done. Total runs stored: ", nrow(obs_runs))
  
  #names(read.csv(files[1]))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # STEP 3: WORST-CASE INCIDENCE SCENARIO SIMULATION
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Identical structure to Step 2, with one key difference:
  # We do NOT divide the seroincidence by 13.
  #
  # RATIONALE: Instead of spreading the 13-year cumulative burden evenly,
  # we assume the full cumulative seroincidence could occur within a single year.
  # This represents a plausible worst-case scenario — for example, an outbreak
  # year, a spike in tick density, or increased human-wildlife contact.
  #
  # ANNUAL INFECTION PROBABILITY IN THIS SCENARIO:
  #   Agriculture & Forestry workers: 0.002 per year (full cumulative rate)
  #   Animal Interaction group:       0.025 per year (full cumulative rate)
  #   General population (baseline):  0.003 per year (full cumulative rate)
  #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  wc_runs <- data.frame()  # fresh empty container for worst-case results
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#read in the environmental data####
distance_data <- read.csv("Forest_Ag_Pop_2022.csv")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
distance_data <- distance_data %>%
  mutate(
    
    # Identify individuals living within 50 metres of a forested area
    Forest_50m = case_when(
      Forest_Dis >= 0 & Forest_Dis < 50 ~ 1,
      Forest_Dis == -1 ~ 0,
      Forest_Dis >= 50 ~ 0,
      TRUE ~ NA_real_
    ),
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    # Identify individuals living within 50 metres of agricultural land
    Agri_50m = case_when(
      Agri_Dista >= 0 & Agri_Dista < 50 ~ 1,
      Agri_Dista == -1 ~ 0,
      Agri_Dista >= 50 ~ 0,
      TRUE ~ NA_real_
    ),
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    # Combine both proximity measures into a single environmental exposure
    # variable. Individuals living within 50 metres of either a forested area
    # or agricultural land are classified as environmentally exposed.
    Environmental_50m = if_else(
      Forest_50m == 1 | Agri_50m == 1,
      1,
      0
    )
  ) %>%
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
# Retain only the variables required for merging with each synthetic
# population dataset
select(person_id, Environmental_50m) %>%
  #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
# Ensure each individual contributes only one environmental exposure record
distinct(person_id, .keep_all = TRUE)
#----------------------------------------------------------------------------------------------------------------------------------------------------------
  for (file in files) {

    # Full pipeline (identical toggles to Section 4 and the observed-scenario
    # loop above): age filter, environmental join, demographic cleaning,
    # behavioural simulation, and risk-flag construction, ready for the
    # worst-case-incidence simulation below.
    pop_level <- prepare_population(
      file,
      distance_data = distance_data,
      filter_adults = TRUE,
      simulate_behavior = TRUE,
      compute_risk_flags = TRUE
    )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
    for (run in 1:100) {
      
      # ASSIGN ANNUAL INFECTION PROBABILITY (worst-case — no division by 13)
      # Each group keeps its own rate but it is applied as a full annual probability
      pop_level$annual_p <- bas_annual_wc                                          # 0.003 for everyone
      pop_level$annual_p[pop_level$forest_agric]      <- agr_annual_wc             # 0.002 for agriculture & forestry
      pop_level$annual_p[pop_level$any_interaction]   <- ani_annual_wc             # 0.025 for animal interaction
      pop_level$annual_p[pop_level$environmental_50m] <- envi_annual_wc            # 0.030 for environmental exposure
      #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
      # COIN FLIP 1 — WHO GETS INFECTED THIS YEAR?
      pop_level$infected <- rbinom(nrow(pop_level), 1, pop_level$annual_p)
      
      # COIN FLIP 2a — OF THOSE INFECTED, WHO DEVELOPS SYMPTOMS?
      # Same correction as observed loop — CFR applies only to symptomatic fraction
      pop_level$cfr_p          <- ifelse(pop_level$older_adult, cfr_older, cfr_young)
      pop_level$is_symptomatic <- 0L
      idx                      <- which(pop_level$infected == 1)
      pop_level$is_symptomatic[idx] <- rbinom(length(idx), 1, reporting_fraction)
      
      # COIN FLIP 2b — OF THOSE SYMPTOMATIC, WHO WOULD HAVE DIED?
      pop_level$would_die <- 0L
      sym_idx <- which(pop_level$is_symptomatic == 1)
      pop_level$would_die[sym_idx] <- rbinom(length(sym_idx), 1, pop_level$cfr_p[sym_idx])
      
      
      # COIN FLIP 3 — DOES THE VACCINE WORK FOR THIS PERSON? (80% chance)
      pop_level$ve_works <- rbinom(nrow(pop_level), 1, vaccine_efficacy)
      
      # COIN FLIP 4 — IS EACH PERSON VACCINATED UNDER EACH COVERAGE SCENARIO?
      pop_level$vac_25 <- rbinom(nrow(pop_level), 1, 0.25)
      pop_level$vac_50 <- rbinom(nrow(pop_level), 1, 0.50)
      pop_level$vac_75 <- rbinom(nrow(pop_level), 1, 0.75)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
      # CALCULATE CASES AVERTED
      # A case is averted if: infected + vaccinated + vaccine worked
      pop_level$case_s25 <- pop_level$infected == 1 & pop_level$vac_25 == 1 & pop_level$ve_works == 1
      pop_level$case_s50 <- pop_level$infected == 1 & pop_level$vac_50 == 1 & pop_level$ve_works == 1
      pop_level$case_s75 <- pop_level$infected == 1 & pop_level$vac_75 == 1 & pop_level$ve_works == 1
      
      # CALCULATE DEATHS AVERTED
      # A death is averted if: would have died + their case was prevented by the vaccine
      pop_level$death_s25 <- pop_level$would_die == 1 & pop_level$vac_25 == 1 & pop_level$ve_works == 1
      pop_level$death_s50 <- pop_level$would_die == 1 & pop_level$vac_50 == 1 & pop_level$ve_works == 1
      pop_level$death_s75 <- pop_level$would_die == 1 & pop_level$vac_75 == 1 & pop_level$ve_works == 1
      
      # FILTER TO EACH SUBGROUP
      young_agr     <- pop_level[!pop_level$older_adult & pop_level$forest_agric, ]
      young_ani     <- pop_level[!pop_level$older_adult & pop_level$any_interaction, ]
      young_env     <- pop_level[!pop_level$older_adult & pop_level$environmental_50m, ]
      young_int_env <- pop_level[!pop_level$older_adult & pop_level$any_interaction & pop_level$environmental_50m, ]
      
      older_agr     <- pop_level[ pop_level$older_adult & pop_level$forest_agric, ]
      older_ani     <- pop_level[ pop_level$older_adult & pop_level$any_interaction, ]
      older_env     <- pop_level[ pop_level$older_adult & pop_level$environmental_50m, ]
      older_int_env <- pop_level[ pop_level$older_adult & pop_level$any_interaction & pop_level$environmental_50m, ]
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
      # COUNT AND STORE RESULTS FOR THIS RUN
      one_run <- data.frame(
        # Adults 18-64: Agriculture & Forestry
        y_agr_cases  = sum(young_agr$infected),
        y_agr_deaths = sum(young_agr$would_die),
        y_agr_c25 = sum(young_agr$case_s25),  y_agr_d25 = sum(young_agr$death_s25),
        y_agr_c50 = sum(young_agr$case_s50),  y_agr_d50 = sum(young_agr$death_s50),
        y_agr_c75 = sum(young_agr$case_s75),  y_agr_d75 = sum(young_agr$death_s75),
        
        # Adults 18-64: Animal Interaction
        y_ani_cases  = sum(young_ani$infected),
        y_ani_deaths = sum(young_ani$would_die),
        y_ani_c25 = sum(young_ani$case_s25),  y_ani_d25 = sum(young_ani$death_s25),
        y_ani_c50 = sum(young_ani$case_s50),  y_ani_d50 = sum(young_ani$death_s50),
        y_ani_c75 = sum(young_ani$case_s75),  y_ani_d75 = sum(young_ani$death_s75),
        
        # Adults 18-64: Environmental 50m
        y_env_cases  = sum(young_env$infected),
        y_env_deaths = sum(young_env$would_die),
        y_env_c25 = sum(young_env$case_s25),  y_env_d25 = sum(young_env$death_s25),
        y_env_c50 = sum(young_env$case_s50),  y_env_d50 = sum(young_env$death_s50),
        y_env_c75 = sum(young_env$case_s75),  y_env_d75 = sum(young_env$death_s75),
        
        # Adults 18-64: Animal Interaction + Environmental 50m
        y_int_env_cases  = sum(young_int_env$infected),
        y_int_env_deaths = sum(young_int_env$would_die),
        y_int_env_c25 = sum(young_int_env$case_s25),  y_int_env_d25 = sum(young_int_env$death_s25),
        y_int_env_c50 = sum(young_int_env$case_s50),  y_int_env_d50 = sum(young_int_env$death_s50),
        y_int_env_c75 = sum(young_int_env$case_s75),  y_int_env_d75 = sum(young_int_env$death_s75),
        
        # Older Adults 65+: Agriculture & Forestry
        o_agr_cases  = sum(older_agr$infected),
        o_agr_deaths = sum(older_agr$would_die),
        o_agr_c25 = sum(older_agr$case_s25),  o_agr_d25 = sum(older_agr$death_s25),
        o_agr_c50 = sum(older_agr$case_s50),  o_agr_d50 = sum(older_agr$death_s50),
        o_agr_c75 = sum(older_agr$case_s75),  o_agr_d75 = sum(older_agr$death_s75),
        
        # Older Adults 65+: Animal Interaction
        o_ani_cases  = sum(older_ani$infected),
        o_ani_deaths = sum(older_ani$would_die),
        o_ani_c25 = sum(older_ani$case_s25),  o_ani_d25 = sum(older_ani$death_s25),
        o_ani_c50 = sum(older_ani$case_s50),  o_ani_d50 = sum(older_ani$death_s50),
        o_ani_c75 = sum(older_ani$case_s75),  o_ani_d75 = sum(older_ani$death_s75),
        
        # Older Adults 65+: Environmental 50m
        o_env_cases  = sum(older_env$infected),
        o_env_deaths = sum(older_env$would_die),
        o_env_c25 = sum(older_env$case_s25),  o_env_d25 = sum(older_env$death_s25),
        o_env_c50 = sum(older_env$case_s50),  o_env_d50 = sum(older_env$death_s50),
        o_env_c75 = sum(older_env$case_s75),  o_env_d75 = sum(older_env$death_s75),
        
        # Older Adults 65+: Animal Interaction + Environmental 50m
        o_int_env_cases  = sum(older_int_env$infected),
        o_int_env_deaths = sum(older_int_env$would_die),
        o_int_env_c25 = sum(older_int_env$case_s25),  o_int_env_d25 = sum(older_int_env$death_s25),
        o_int_env_c50 = sum(older_int_env$case_s50),  o_int_env_d50 = sum(older_int_env$death_s50),
        o_int_env_c75 = sum(older_int_env$case_s75),  o_int_env_d75 = sum(older_int_env$death_s75)
      )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
      wc_runs <- rbind(wc_runs, one_run)  # append this run's results to storage
      
    }  # end inner loop — next of 100 runs for this file
    
  }  # end outer loop — next of 10 synthetic population files
  
  message("Worst-case scenario done. Total runs stored: ", nrow(wc_runs))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
  
  # STEP 4: HELPER FUNCTION — FORMAT AS MEDIAN [IQR]
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # fmt_iqr() takes a column of 1000 numbers (one per simulation run) and
  # summarises it as: "Median [Q25–Q75]"
  #
  # Example: if across 1000 runs, Agriculture workers had between 3 and 6
  # infections with a typical value of 5, the output would be: "5 [3–6]"
  #
  # Q25 = 25th percentile (lower bound of the IQR)
  # Q75 = 75th percentile (upper bound of the IQR)
  # IQR captures the middle 50% of all simulation results
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  fmt_iqr <- function(x) {
    paste0(
      round(median(x)),           " [",   # Median: most typical result
      round(quantile(x, 0.25)),   "\u2013",  # Q25: lower bound (– is an en-dash)
      round(quantile(x, 0.75)),   "]"     # Q75: upper bound
    )
  }
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # STEP 5: BUILD THE FOUR OUTPUT TABLES
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# build_table() assembles one table from simulation results.
# Each table has 12 rows:
#   - Agriculture & Forestry                 at 25%, 50%, 75% coverage (3 rows)
#   - Animal Interaction                     at 25%, 50%, 75% coverage (3 rows)
#   - Environmental 50m                      at 25%, 50%, 75% coverage (3 rows)
#   - Animal Interaction + Environmental 50m at 25%, 50%, 75% coverage (3 rows)
#
# Arguments:
#   results         — the data frame of 1000 simulation runs (obs_runs or wc_runs)
#   age_group       — "young" (18-64) or "older" (65+); determines which columns to use
#   cfr_label       — the CFR shown in the table ("5%" or "20%")
#   agr_prob_label  — the label shown in the Annual Infection Prob. column for Agriculture & Forestry
#   ani_prob_label   — the label shown in the Annual Infection Prob. column for Animal Interaction
#   env_prob_label   — the label shown in the Annual Infection Prob. column for Environmental 50m
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
build_table <- function(results, age_group, cfr_label, agr_prob_label, ani_prob_label, env_prob_label) {
  
  # p is the column name prefix: "y_" for young adults, "o_" for older adults
  # This lets the same function work for both age groups without repeating code
  p <- if (age_group == "young") "y_" else "o_"
  
  # Build 3 rows for Agriculture & Forestry
  agr_rows <- data.frame(
    Risk_Subgroup         = "Agriculture & Forestry",   # risk group label
    Annual_Infection_Prob = agr_prob_label,             # e.g. "0.002/13" or "0.002"
    Vaccine_Coverage      = c("25%", "50%", "75%"),     # one row per coverage level
    CFR                   = cfr_label,                  # e.g. "5%" or "20%"
    Cases                 = fmt_iqr(results[[paste0(p, "agr_cases")]]),    # pre-vaccination cases
    Deaths                = fmt_iqr(results[[paste0(p, "agr_deaths")]]),   # pre-vaccination deaths
    Cases_Averted         = c(
      fmt_iqr(results[[paste0(p, "agr_c25")]]),  # cases averted at 25% coverage
      fmt_iqr(results[[paste0(p, "agr_c50")]]),  # cases averted at 50% coverage
      fmt_iqr(results[[paste0(p, "agr_c75")]])   # cases averted at 75% coverage
    ),
    Deaths_Averted        = c(
      fmt_iqr(results[[paste0(p, "agr_d25")]]),  # deaths averted at 25% coverage
      fmt_iqr(results[[paste0(p, "agr_d50")]]),  # deaths averted at 50% coverage
      fmt_iqr(results[[paste0(p, "agr_d75")]])   # deaths averted at 75% coverage
    )
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
  # Build 3 rows for Animal Interaction
  ani_rows <- data.frame(
    Risk_Subgroup         = "Animal Interaction",
    Annual_Infection_Prob = ani_prob_label,
    Vaccine_Coverage      = c("25%", "50%", "75%"),
    CFR                   = cfr_label,
    Cases                 = fmt_iqr(results[[paste0(p, "ani_cases")]]),
    Deaths                = fmt_iqr(results[[paste0(p, "ani_deaths")]]),
    Cases_Averted         = c(
      fmt_iqr(results[[paste0(p, "ani_c25")]]),
      fmt_iqr(results[[paste0(p, "ani_c50")]]),
      fmt_iqr(results[[paste0(p, "ani_c75")]])
    ),
    Deaths_Averted        = c(
      fmt_iqr(results[[paste0(p, "ani_d25")]]),
      fmt_iqr(results[[paste0(p, "ani_d50")]]),
      fmt_iqr(results[[paste0(p, "ani_d75")]])
    )
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
  # Build 3 rows for Environmental 50m
  env_rows <- data.frame(
    Risk_Subgroup         = "Environmental 50m",
    Annual_Infection_Prob = env_prob_label,
    Vaccine_Coverage      = c("25%", "50%", "75%"),
    CFR                   = cfr_label,
    Cases                 = fmt_iqr(results[[paste0(p, "env_cases")]]),
    Deaths                = fmt_iqr(results[[paste0(p, "env_deaths")]]),
    Cases_Averted         = c(
      fmt_iqr(results[[paste0(p, "env_c25")]]),
      fmt_iqr(results[[paste0(p, "env_c50")]]),
      fmt_iqr(results[[paste0(p, "env_c75")]])
    ),
    Deaths_Averted        = c(
      fmt_iqr(results[[paste0(p, "env_d25")]]),
      fmt_iqr(results[[paste0(p, "env_d50")]]),
      fmt_iqr(results[[paste0(p, "env_d75")]])
    )
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
  # Build 3 rows for Animal Interaction + Environmental 50m
  int_env_rows <- data.frame(
    Risk_Subgroup         = "Animal Interaction + Environmental 50m",
    Annual_Infection_Prob = env_prob_label,
    Vaccine_Coverage      = c("25%", "50%", "75%"),
    CFR                   = cfr_label,
    Cases                 = fmt_iqr(results[[paste0(p, "int_env_cases")]]),
    Deaths                = fmt_iqr(results[[paste0(p, "int_env_deaths")]]),
    Cases_Averted         = c(
      fmt_iqr(results[[paste0(p, "int_env_c25")]]),
      fmt_iqr(results[[paste0(p, "int_env_c50")]]),
      fmt_iqr(results[[paste0(p, "int_env_c75")]])
    ),
    Deaths_Averted        = c(
      fmt_iqr(results[[paste0(p, "int_env_d25")]]),
      fmt_iqr(results[[paste0(p, "int_env_d50")]]),
      fmt_iqr(results[[paste0(p, "int_env_d75")]])
    )
  )
  
  rbind(agr_rows, ani_rows, env_rows, int_env_rows)
}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Generate all 4 tables by calling build_table() with the appropriate inputs:
  #   obs_runs = observed scenario results (÷13)
  #   wc_runs  = worst-case scenario results (no ÷13)
  table_1 <- build_table(obs_runs, "young", "5%",  "0.002/13", "0.025/13", "0.030/13")# Adults 18-64,    Observed
  table_2 <- build_table(wc_runs,  "young", "5%",  "0.002",    "0.025",    "0.030") # Adults 18-64,    Worst-case
  table_3 <- build_table(obs_runs, "older", "20%", "0.002/13", "0.025/13", "0.030/13")# Older Adults 65+, Observed
  table_4 <- build_table(wc_runs,  "older", "20%", "0.002",    "0.025",    "0.030")# Older Adults 65+, Worst-case
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # STEP 6: DISPLAY TABLES IN RSTUDIO VIEWER
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # render_gt() converts a plain data frame into a formatted gt table.
  # - tab_header() adds only a subtitle (no table number — you will add that in Word)
  # - cols_label() renames columns to publication-friendly headings
  # - tab_spanner() groups related columns under a shared header, matching
  #   your supervisor's table template
  # - opt_stylize() applies a clean visual style
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  render_gt <- function(tbl, subtitle) {
    tbl |>
      gt() |>
      tab_header(title = " ", subtitle = subtitle) |>   # blank title, subtitle only (you will number in Word)
      cols_label(
        Risk_Subgroup         = "Risk Subgroup",
        Annual_Infection_Prob = "Annual Infection Prob.",
        Vaccine_Coverage      = "Vaccine Coverage",
        CFR                   = "CFR",
        Cases                 = "Baseline Cases",
        Deaths                = "Baseline Deaths",
        Cases_Averted         = "Cases Averted",
        Deaths_Averted        = "Deaths Averted"
      ) |>
      tab_spanner(
        label   = "Scenario Parameters",           # groups the three input columns
        columns = c(Annual_Infection_Prob, Vaccine_Coverage, CFR)
      ) |>
      tab_spanner(
        label   = "Outcomes \u2014 Median [IQR]",  # groups the four outcome columns
        columns = c(Cases, Deaths, Cases_Averted, Deaths_Averted)
      ) |>
      opt_stylize(style = 1)                       # clean blue-bordered style
  }
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Display each table — appears in the RStudio Viewer pane
  #gtsave(render_gt(table_1, "Adults (18\u201364) | Observed Incidence | CFR: 5%"),"table_1.docx")
  #gtsave(render_gt(table_2, "Adults (18\u201364) | Worst-case Incidence | CFR: 5%"),"table_2.docx")
  #gtsave(render_gt(table_3, "Older Adults (65+) | Observed Incidence | CFR: 20%"),"table_3.docx")
  #gtsave(render_gt(table_4, "Older Adults (65+) | Worst-case Incidence | CFR: 20%"),"table_4.docx")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #Plot box plot to show the impact of the vaccines across ####
  # Prepare data for panelled box plots
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  subgroup_labels <- c(
    agr = "Agriculture & Forestry",
    ani = "Animal Interaction",
    env = "Environmental 50m",
    int_env = "Animal Interaction + Environmental 50m"
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  # make_boxplot_data() reshapes the wide simulation-run columns (e.g.
  # y_agr_c25, y_agr_c50, ...) into a long, tidy data frame with one row per
  # (subgroup x outcome x coverage) combination per simulation run — the
  # format ggplot2 needs for faceted/grouped plotting.
  make_boxplot_data <- function(results, scenario_label, prefix) {
    bind_rows(lapply(names(subgroup_labels), function(sg) {
      bind_rows(lapply(c("c", "d"), function(metric) {
        bind_rows(lapply(c("25", "50", "75"), function(cov) {
          col_name <- paste0(prefix, sg, "_", metric, cov)
          
          tibble(
            scenario = scenario_label,
            subgroup = subgroup_labels[[sg]],
            outcome  = ifelse(metric == "c", "Cases Averted", "Deaths Averted"),
            coverage = factor(paste0(cov, "%"), levels = c("25%", "50%", "75%")),
            value    = results[[col_name]]
          )
        }))
      }))
    }))
  }
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  plot_data_young <- bind_rows(
    make_boxplot_data(obs_runs, "Observed incidence", "y_"),
    make_boxplot_data(wc_runs,  "Worst-case incidence", "y_")
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  plot_data_older <- bind_rows(
    make_boxplot_data(obs_runs, "Observed incidence", "o_"),
    make_boxplot_data(wc_runs,  "Worst-case incidence", "o_")
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  # AGGREGATE SUMMARY DATA FOR BAR PLOTS
  # Rather than plotting all 1,000 raw simulation values as a box plot, we
  # summarise each (scenario x subgroup x outcome x coverage) combination
  # down to its median and IQR up front. This keeps the resulting bar chart
  # clean and avoids near-zero-variance boxes rendering as invisible slivers.
  #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  make_bar_data <- function(df_long) {
    df_long %>%
      group_by(scenario, subgroup, outcome, coverage) %>%
      summarise(
        median_val = median(value, na.rm = TRUE),
        lower_iqr  = quantile(value, 0.25, na.rm = TRUE),
        upper_iqr  = quantile(value, 0.75, na.rm = TRUE),
        .groups    = "drop"
      )
  }
  
  bar_data_young <- make_bar_data(plot_data_young)
  bar_data_older <- make_bar_data(plot_data_older)
  
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  # UNIFIED BAR PANEL PLOTTING FUNCTION
  # Draws grouped bars (median outcome per coverage level, grouped by risk
  # subgroup) with error bars showing the IQR across the 1,000 simulation runs.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  
  make_bar_panel <- function(df, scenario_keep, outcome_keep, panel_title, y_label = "Averted Outcomes (Median [IQR])") {
    ggplot(
      df %>% filter(scenario == scenario_keep, outcome == outcome_keep),
      aes(x = coverage, y = median_val, fill = subgroup)
    ) +
      # Clean grouped bars
      geom_bar(
        stat = "identity",
        position = position_dodge(width = 0.8),
        width = 0.7,
        color = "white",
        linewidth = 0.2
      ) +
      # Error bars representing the interquartile range (IQR) of the 1,000 runs
      geom_errorbar(
        aes(ymin = lower_iqr, ymax = upper_iqr),
        position = position_dodge(width = 0.8),
        width = 0.3,
        color = "black",
        linewidth = 0.5
      ) +
      labs(
        title = panel_title,
        x = "Vaccine Coverage",
        y = y_label,
        fill = "Risk Subgroup"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  }
  
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  # Adults 18–64: Figure 2a–d (Bar Layout)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  p_a <- make_bar_panel(bar_data_young, "Observed incidence", "Cases Averted",
                        "Observed incidence — cases averted")
  
  p_b <- make_bar_panel(bar_data_young, "Worst-case incidence", "Cases Averted",
                        "Worst-case incidence — cases averted")
  
  p_c <- make_bar_panel(bar_data_young, "Observed incidence", "Deaths Averted",
                        "Observed incidence — deaths averted")
  
  p_d <- make_bar_panel(bar_data_young, "Worst-case incidence", "Deaths Averted",
                        "Worst-case incidence — deaths averted")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  fig_young <- ((p_a | p_b) / (p_c | p_d)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      tag_levels = "a"
    ) &
    theme(
      plot.tag = element_text(face = "bold"),
      legend.position = "bottom"
    )
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  #ggsave(
  #"figure_young_vaccine_barplots_panelled.png",
   #plot = fig_young,
   #width = 16,
   #height = 11,
   #dpi = 300
  #)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  # Older adults 65+: Figure 3a–d (Bar Layout)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  p2_a <- make_bar_panel(bar_data_older, "Observed incidence", "Cases Averted",
                         "Observed incidence — cases averted")
  
  p2_b <- make_bar_panel(bar_data_older, "Worst-case incidence", "Cases Averted",
                         "Worst-case incidence — cases averted")
  
  p2_c <- make_bar_panel(bar_data_older, "Observed incidence", "Deaths Averted",
                         "Observed incidence — deaths averted")
  
  p2_d <- make_bar_panel(bar_data_older, "Worst-case incidence", "Deaths Averted",
                         "Worst-case incidence — deaths averted")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  fig_older <- ((p2_a | p2_b) / (p2_c | p2_d)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      tag_levels = "a"
    ) &
    theme(
      plot.tag = element_text(face = "bold"),
      legend.position = "bottom"
    )
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
  #ggsave(
   #"figure_older_vaccine_barplots_panelled.png",
   #plot = fig_older,
   #width = 16,
   #height = 11,
   # dpi = 300
  #)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
#===================================================================================================================
  # SECTION 6: SENSITIVITY ANALYSIS
  #
  # Section 5 fixed the reporting fraction (the proportion of infections that
  # are symptomatic/detected) at 8%. Because this value is a judgement call
  # rather than a directly-measured quantity, this section re-runs the full
  # vaccination simulation across a small grid of plausible reporting
  # fractions (4%, 8%, 12%, 20%) and checks how sensitive the deaths-averted
  # estimates are to that assumption. Everything else (rates, CFR, vaccine
  # efficacy, coverage levels) is held constant.
#===================================================================================================================
  
  set.seed(1802)
  
  # Same Median [IQR] formatter as Section 5's fmt_iqr(), duplicated with an
  # explicit na.rm = TRUE for safety when summing across subgroup columns.
  fmt_iqr_num <- function(x) {
    paste0(
      round(median(x, na.rm = TRUE)), " [",
      round(quantile(x, 0.25, na.rm = TRUE)), "–",
      round(quantile(x, 0.75, na.rm = TRUE)), "]"
    )
  }
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
  # Column names holding "deaths averted" for each risk subgroup x coverage
  # level combination, split by age group. Summing across these columns gives
  # the total deaths averted across all subgroups for a given coverage level.
  death_cols_y <- c(
    "y_agr_d25", "y_agr_d50", "y_agr_d75",
    "y_ani_d25", "y_ani_d50", "y_ani_d75",
    "y_env_d25", "y_env_d50", "y_env_d75",
    "y_int_env_d25", "y_int_env_d50", "y_int_env_d75"
  )
  
  death_cols_o <- c(
    "o_agr_d25", "o_agr_d50", "o_agr_d75",
    "o_ani_d25", "o_ani_d50", "o_ani_d75",
    "o_env_d25", "o_env_d50", "o_env_d75",
    "o_int_env_d25", "o_int_env_d50", "o_int_env_d75"
  )
  
  total_deaths <- function(df, cols) {
    rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
  }
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
  # prepare_population_for_sensitivity() now simply calls the shared
  # prepare_population() function defined in Section 1 with the full set of
  # toggles (age filter + behaviour simulation + risk flags) — the same
  # pipeline used in Sections 4-5. Keeping this thin wrapper (rather than
  # calling prepare_population() directly below) preserves the original
  # function name so nothing else in this section needs to change.
  prepare_population_for_sensitivity <- function(file, distance_data) {
    prepare_population(
      file,
      distance_data = distance_data,
      filter_adults = TRUE,
      simulate_behavior = TRUE,
      compute_risk_flags = TRUE
    )
  }
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
  # simulate_one_scenario() runs the same 100-iterations-per-file stochastic
  # simulation used in Section 5 (infection -> symptom -> death -> vaccine ->
  # coverage coin flips), but takes the annual infection rates as an argument
  # so it can be reused for both the observed and worst-case scenarios, across
  # any reporting fraction, without duplicating the simulation logic again.
  simulate_one_scenario <- function(prepared_populations, annual_rates, reporting_fraction) {
    
    runs <- data.frame()
    
    for (pop_level in prepared_populations) {
      for (run in 1:100) {
        
        pop_level$annual_p <- annual_rates$baseline
        pop_level$annual_p[pop_level$forest_agric]      <- annual_rates$agriculture
        pop_level$annual_p[pop_level$any_interaction]   <- annual_rates$animal
        pop_level$annual_p[pop_level$environmental_50m] <- annual_rates$environmental
        
        pop_level$infected <- rbinom(nrow(pop_level), 1, pop_level$annual_p)
        
        pop_level$cfr_p <- ifelse(pop_level$older_adult, cfr_older, cfr_young)
        pop_level$is_symptomatic <- 0L
        idx <- which(pop_level$infected == 1)
        pop_level$is_symptomatic[idx] <- rbinom(length(idx), 1, reporting_fraction)
        
        pop_level$would_die <- 0L
        sym_idx <- which(pop_level$is_symptomatic == 1)
        pop_level$would_die[sym_idx] <- rbinom(length(sym_idx), 1, pop_level$cfr_p[sym_idx])
        
        pop_level$ve_works <- rbinom(nrow(pop_level), 1, vaccine_efficacy)
        
        pop_level$vac_25 <- rbinom(nrow(pop_level), 1, 0.25)
        pop_level$vac_50 <- rbinom(nrow(pop_level), 1, 0.50)
        pop_level$vac_75 <- rbinom(nrow(pop_level), 1, 0.75)
        
        pop_level$case_s25 <- pop_level$infected == 1 & pop_level$vac_25 == 1 & pop_level$ve_works == 1
        pop_level$case_s50 <- pop_level$infected == 1 & pop_level$vac_50 == 1 & pop_level$ve_works == 1
        pop_level$case_s75 <- pop_level$infected == 1 & pop_level$vac_75 == 1 & pop_level$ve_works == 1
        
        pop_level$death_s25 <- pop_level$would_die == 1 & pop_level$vac_25 == 1 & pop_level$ve_works == 1
        pop_level$death_s50 <- pop_level$would_die == 1 & pop_level$vac_50 == 1 & pop_level$ve_works == 1
        pop_level$death_s75 <- pop_level$would_die == 1 & pop_level$vac_75 == 1 & pop_level$ve_works == 1
        
        young_agr     <- pop_level[!pop_level$older_adult & pop_level$forest_agric, ]
        young_ani     <- pop_level[!pop_level$older_adult & pop_level$any_interaction, ]
        young_env     <- pop_level[!pop_level$older_adult & pop_level$environmental_50m, ]
        young_int_env <- pop_level[!pop_level$older_adult & pop_level$any_interaction & pop_level$environmental_50m, ]
        
        older_agr     <- pop_level[ pop_level$older_adult & pop_level$forest_agric, ]
        older_ani     <- pop_level[ pop_level$older_adult & pop_level$any_interaction, ]
        older_env     <- pop_level[ pop_level$older_adult & pop_level$environmental_50m, ]
        older_int_env <- pop_level[ pop_level$older_adult & pop_level$any_interaction & pop_level$environmental_50m, ]
        
        one_run <- data.frame(
          y_agr_cases  = sum(young_agr$infected),
          y_agr_deaths = sum(young_agr$would_die),
          y_agr_c25 = sum(young_agr$case_s25),  y_agr_d25 = sum(young_agr$death_s25),
          y_agr_c50 = sum(young_agr$case_s50),  y_agr_d50 = sum(young_agr$death_s50),
          y_agr_c75 = sum(young_agr$case_s75),  y_agr_d75 = sum(young_agr$death_s75),
          
          y_ani_cases  = sum(young_ani$infected),
          y_ani_deaths = sum(young_ani$would_die),
          y_ani_c25 = sum(young_ani$case_s25),  y_ani_d25 = sum(young_ani$death_s25),
          y_ani_c50 = sum(young_ani$case_s50),  y_ani_d50 = sum(young_ani$death_s50),
          y_ani_c75 = sum(young_ani$case_s75),  y_ani_d75 = sum(young_ani$death_s75),
          
          y_env_cases  = sum(young_env$infected),
          y_env_deaths = sum(young_env$would_die),
          y_env_c25 = sum(young_env$case_s25),  y_env_d25 = sum(young_env$death_s25),
          y_env_c50 = sum(young_env$case_s50),  y_env_d50 = sum(young_env$death_s50),
          y_env_c75 = sum(young_env$case_s75),  y_env_d75 = sum(young_env$death_s75),
          
          y_int_env_cases  = sum(young_int_env$infected),
          y_int_env_deaths = sum(young_int_env$would_die),
          y_int_env_c25 = sum(young_int_env$case_s25),  y_int_env_d25 = sum(young_int_env$death_s25),
          y_int_env_c50 = sum(young_int_env$case_s50),  y_int_env_d50 = sum(young_int_env$death_s50),
          y_int_env_c75 = sum(young_int_env$case_s75),  y_int_env_d75 = sum(young_int_env$death_s75),
          
          o_agr_cases  = sum(older_agr$infected),
          o_agr_deaths = sum(older_agr$would_die),
          o_agr_c25 = sum(older_agr$case_s25),  o_agr_d25 = sum(older_agr$death_s25),
          o_agr_c50 = sum(older_agr$case_s50),  o_agr_d50 = sum(older_agr$death_s50),
          o_agr_c75 = sum(older_agr$case_s75),  o_agr_d75 = sum(older_agr$death_s75),
          
          o_ani_cases  = sum(older_ani$infected),
          o_ani_deaths = sum(older_ani$would_die),
          o_ani_c25 = sum(older_ani$case_s25),  o_ani_d25 = sum(older_ani$death_s25),
          o_ani_c50 = sum(older_ani$case_s50),  o_ani_d50 = sum(older_ani$death_s50),
          o_ani_c75 = sum(older_ani$case_s75),  o_ani_d75 = sum(older_ani$death_s75),
          
          o_env_cases  = sum(older_env$infected),
          o_env_deaths = sum(older_env$would_die),
          o_env_c25 = sum(older_env$case_s25),  o_env_d25 = sum(older_env$death_s25),
          o_env_c50 = sum(older_env$case_s50),  o_env_d50 = sum(older_env$death_s50),
          o_env_c75 = sum(older_env$case_s75),  o_env_d75 = sum(older_env$death_s75),
          
          o_int_env_cases  = sum(older_int_env$infected),
          o_int_env_deaths = sum(older_int_env$would_die),
          o_int_env_c25 = sum(older_int_env$case_s25),  o_int_env_d25 = sum(older_int_env$death_s25),
          o_int_env_c50 = sum(older_int_env$case_s50),  o_int_env_d50 = sum(older_int_env$death_s50),
          o_int_env_c75 = sum(older_int_env$case_s75),  o_int_env_d75 = sum(older_int_env$death_s75)
        )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------         
        runs <- rbind(runs, one_run)
      }
    }
    
    runs
  }
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  # run_vaccine_sen() runs the full observed + worst-case simulation pair for
  # ONE reporting fraction, using the reusable helpers above. Called once per
  # value in reporting_fraction_grid below.
  run_vaccine_sen <- function(reporting_fraction) {
    
    distance_data <- read.csv("Forest_Ag_Pop_2022.csv") %>%
      mutate(
        Forest_50m = case_when(
          Forest_Dis >= 0 & Forest_Dis < 50 ~ 1,
          Forest_Dis == -1 ~ 0,
          Forest_Dis >= 50 ~ 0,
          TRUE ~ NA_real_
        ),
        Agri_50m = case_when(
          Agri_Dista >= 0 & Agri_Dista < 50 ~ 1,
          Agri_Dista == -1 ~ 0,
          Agri_Dista >= 50 ~ 0,
          TRUE ~ NA_real_
        ),
        Environmental_50m = if_else(Forest_50m == 1 | Agri_50m == 1, 1, 0)
      ) %>%
      select(person_id, Environmental_50m) %>%
      distinct(person_id, .keep_all = TRUE)
    
    files <- list.files(pattern = "^2015_.*\\.csv$", full.names = TRUE)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    # Prepare all 10 synthetic populations once per reporting fraction (the
    # population prep itself does not depend on the reporting fraction, but is
    # kept inside this function so each call is self-contained and reproducible
    # on its own).
    prepared_populations <- lapply(
      files,
      prepare_population_for_sensitivity,
      distance_data = distance_data
    )
    
    obs_rates <- list(
      baseline = bas_annual_obs,
      agriculture = agr_annual_obs,
      animal = ani_annual_obs,
      environmental = envi_annual_obs
    )
    
    wc_rates <- list(
      baseline = bas_annual_wc,
      agriculture = agr_annual_wc,
      animal = ani_annual_wc,
      environmental = envi_annual_wc
    )
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    obs_runs <- simulate_one_scenario(prepared_populations, obs_rates, reporting_fraction)
    message("Observed scenario done. Total runs stored: ", nrow(obs_runs))
    
    wc_runs <- simulate_one_scenario(prepared_populations, wc_rates, reporting_fraction)
    message("Worst-case scenario done. Total runs stored: ", nrow(wc_runs))
    
    list(
      obs_runs = obs_runs,
      wc_runs = wc_runs
    )
  }
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
  # Grid of plausible reporting fractions to stress-test against: 4% (very
  # conservative — most infections silent), 8% (the value used in Section 5),
  # 12%, and 20% (closer to a well-monitored outbreak setting).
  reporting_fraction_grid <- c(0.04, 0.08, 0.12, 0.20)
  
  sensitivity_results <- setNames(
    lapply(reporting_fraction_grid, run_vaccine_sen),
    paste0("rf_", reporting_fraction_grid)
  )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
  # Collapse each reporting fraction's results down to total deaths averted
  # (summed across all four risk subgroups) for young and older adults, under
  # both the observed and worst-case incidence scenarios.
  sensitivity_table <- do.call(rbind, lapply(names(sensitivity_results), function(nm) {
    out <- sensitivity_results[[nm]]
    
    data.frame(
      reporting_fraction = sub("^rf_", "", nm),
      young_observed     = fmt_iqr_num(total_deaths(out$obs_runs, death_cols_y)),
      older_observed     = fmt_iqr_num(total_deaths(out$obs_runs, death_cols_o)),
      young_worst_case   = fmt_iqr_num(total_deaths(out$wc_runs,  death_cols_y)),
      older_worst_case   = fmt_iqr_num(total_deaths(out$wc_runs,  death_cols_o)),
      stringsAsFactors   = FALSE
    )
  }))
  
  print(sensitivity_table)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
  sensitivity_gt <- sensitivity_table %>%
    gt() %>%
    tab_header(
      title = "",
      subtitle = ""
    ) %>%
    cols_label(
      reporting_fraction = "Reporting fraction",
      young_observed     = "Deaths averted, adults 18–64 years (observed)",
      older_observed     = "Deaths averted, adults ≥65 years (observed)",
      young_worst_case   = "Deaths averted, adults 18–64 years (worst-case)",
      older_worst_case   = "Deaths averted, adults ≥65 years (worst-case)"
    )
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  sensitivity_gt
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  #gtsave(sensitivity_gt, "Supplementary_Table_S5.docx")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
