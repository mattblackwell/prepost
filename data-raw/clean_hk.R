library(tidyverse)
library(readstata13)

## Load data ----
hk_raw <- read.dta13("data-raw/Horowitz-Klaus-data.dta")

## Relevant variables

# From Stata:

# . codebook treatment support_candidate bloc_before education age close_own_di gender
# kalenjin kikuyu land_security_low, compact
#
# Variable           Obs  Min  Max  Label
# ------------------------------------------------------------------------------------------
# treatment          834  0    2   Treatment group
# support_candidate  818  1    5   How likely would you be to support this candidate?
# bloc_before        834  0    1   Land question battery answered before treatments
# education          834  1    8   Highest level of education
# age                834  18   84  Age
# close_own_di       825  0    1   R feels 'very close' to own ethnic group
# gender             834  0    1   Gender
# kalenjin           834  0    1   Kalenjin
# kikuyu             834  0    1   Kikuyu
# land_security_low  807  0    1   Perceived land security = 'not secure'
# ------------------------------------------------------------------------------------------

# Kalenjin (insiders) and Kikuyu (outsiders) are ethnic groups

# Treatment is a hypothetical candidate speech
# 3 factors
# 0 = control = neutral speech
# 1 = T1 = control + land appeal ("I will make sure you have land for farming")
# 2 = T2 = control + land appeal + ethnic grievance ("migrants have stolen our land")

# Since we need a binary treatment variable, we combine T1 and T2.

hk <-
  hk_raw %>%
  mutate(
    # treatment
    treat = as.factor(treatment),
    # binary treatment
    treat_comb = case_when(
      treatment == 0 ~ 0, 
      treatment == 1 ~ 1, 
      treatment == 2 ~ 1
      ),
    # prepost indicator
    prepost = case_when(bloc_before == 1 ~ 0,
                  bloc_before == 0 ~ 1),
    # outcome
    support_candidate_num =
      case_when(
        support_candidate == "Very unlikely" ~ -2,
        support_candidate == "Somewhat unlikely" ~ -1,
        support_candidate == "Neither unlikely now likely" ~ 0,
        support_candidate == "Somewhat likely" ~ 1,
        support_candidate == "Very likely" ~ 2,
      ),
    support = ifelse(support_candidate_num > 0, 1, 0),
    # binary moderator
    land_insecure = land_security_low,
    # covariates
    educ =
      case_when(
        education == "No formal schooling" ~ 1,
        education == "Some primary school" ~ 2,
        education == "Primary school completed" ~ 3,
        education == "Some secondary school" ~ 4,
        education == "Secondary school completed" ~ 5,
        education == "College" ~ 6,
        education == "Some university" ~ 7,
        education == "University completed" ~ 8,
      ),
    female =
      case_when(female == "Yes" ~ 1, female == "No" ~ 0),
    close_own = close_own_di
  ) %>%
  select(
    treat,
    treat_comb,
    prepost,
    support,
    land_insecure,
    educ,
    age,
    close_own,
    female,
    kalenjin
  ) %>%
  as.data.frame()

# get rid of redundant attributes (e.g. labels for dropped variables)
attributes(hk) <- list(
  names = names(hk),
  class = c("tbl_df", "tbl", "data.frame"),
  row.names = 1:nrow(hk)
)

# subset to Kalenjin only
hk <-
  hk %>%
  filter(hk$kalenjin == 1) %>%
  select(support, land_insecure, treat_comb, prepost, age, female, close_own, educ, treat) %>%
  na.omit() %>%
  as.data.frame()

# save
write_csv(hk,
          "./data-raw/land_experiment_cleaned.csv")

land_experiment <- hk

usethis::use_data(land_experiment,
                  overwrite = TRUE,
                  version = 2)
