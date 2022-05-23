library(glmmTMB)
library(foreign)
library(tidyverse)
library(TMB)
library(lme4)
source("functions.R")

#####--------------------
##### PIRLS DATA ----
#####--------------------
asg_data <-load.file.type(file.type = "asg.*")
asg_data_long <-data.table::rbindlist(asg_data)
#Student data
s_data <- asg_data_long %>%
  select(IDCNTRY:ITLANG,
         ASBG01, ASDAGE, ASBG03, ASBGHRL, ASDGHRL, ASBGSCR, ASDGSCR,
         ASRREA01:ASRRSI05) #Plausable score values

names(s_data) <- c("Country", "Book", "School", "Class", "Student_ID", "Grade",
                   "Sex", "ADMINI", "Lang", "Gender", "Age", "Lang_home",
                   "Resources", "Resources_cat", "Conf", "Conf_cat",
                   "Overall1", "Overall2", "Overall3", "Overall4", "Overall5",
                   "Lit_purpose1", "Lit_purpose2", "Lit_purpose3", "Lit_purpose4",
                   "Lit_purpose5", "Inf_purpose1", "Inf_purpose2", "Inf_purpose3",
                   "Inf_purpose4", "Inf_purpose5", "Inter1", "Inter2", "Inter3",
                   "Inter4", "Inter5", "Str_for1", "Str_for2", "Str_for3", "Str_for4", "Str_for5")

categorical_vars <- c("Country", "School", "Class",  "Grade",  "Sex", "Lang",
                      "Gender", "Lang_home",  "Resources_cat", "Conf_cat")
score_vars <- names(s_data)[which(names(s_data) == "Overall1"):which(names(s_data) =="Str_for5")]
numeric_vars <- c("Conf")
# # make sure variables are the correct class
s_data <- s_data %>%
  mutate_at(categorical_vars, as.factor) %>%
  mutate_at(score_vars, as.character ) %>%
  mutate_at(score_vars, as.numeric ) %>%
  mutate_at(numeric_vars, as.character ) %>%
  mutate_at(numeric_vars, as.numeric )

# ##### School Data
acg_data <-load.file.type(file.type = "acg.*")
acg_data_long <-data.table::rbindlist(acg_data)
sch_data <- acg_data_long %>%
  select(IDCNTRY, IDSCHOOL,
         ACBG03A, #GEN\STUDENTS BACKGROUND\ECONOMIC DISADVA
         ACDG09 ##SIZE OF SCHOOL LIBRARY
  )
names(sch_data) <- c("Country", "School", "Eco_disad", "Size_lib")
pirls <- left_join(s_data, sch_data, by = intersect(names(s_data), names(sch_data)))

country_id <- unique(pirls$Country)
level_key <- c("36" = "Australia", "40" = "Austria", "31" = "Azerbaijan", "48" = "Bahrain", "956" = "Belgium (Flemish)",
               "957" = "Belgium (French)", "100" = "Bulgaria", "124" = "Canada", "152" = "Chile",
               "203" = "Czech Republic", "208" = "Denmark", "246" = "Finland", "250" = "France",
               "268" = "Georgia", "276" = "Germany", "926" = "England", "724" = "Spain",
               "784" = "United Arab Emirates", "7841" = "Dubai, UAE", "7842" = "Abu Dhabi, UAE",
               "32001" = "Buenos Aires, Argentina", "9132" = "Ontario, Canada", "9133" = "Quebec, Canada",
               "72401" = "Andalusia, Spain", "724005" = "Madrid, Spain")
pirls$Country <-  plyr::revalue(pirls$Country, level_key)

dup <- pirls[duplicated(pirls$Student_ID),]
s129 <- dup[dup$Student_ID == 1290201,]

### Schools are labelled the same in each country
table(pirls[ , c("School", "Country")])

pirls$SchoolsInCountry <- with(pirls, factor(Country:School))
table(pirls[ , c("SchoolsInCountry", "Country")])

save(pirls, file = "Data/pirls_cat.RData")
