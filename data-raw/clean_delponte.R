library(readstata13)
library(tidyverse)


## Data found at supplement of https://doi.org/10.1177/1465116520966653
dat <- read.dta13("./data-raw/EUP_ReplicationDataset.dta")

dat <- dat %>%
  # omit "pure control" group as they
  # were not asked the outcome var
  # omit respondents who failed attention checks
  # (following author)
  filter(vignetteCheck > 0,
         supercheck > 0,
         atcheck != 8,
         atcheck != 4)

# 2x2 treatment represented as 4 dummies
# symcom = Symbolic Commonality
# symdiv = Symbolic Division
# econcom = Economic Commonality
# econdiv = Economic Division
dat <-
  dat %>%
  mutate(treatment = case_when(symcom == 1 ~ "symcom",
                               symdiv == 1 ~ "symdiv",
                               econcom == 1 ~ "econcom",
                               econdiv == 1 ~ "econdiv",
                               control == 1 ~ "control"),
         treatment = as.factor(treatment),
         treatment = relevel(treatment, ref = "control"),
         t_commonality = case_when(symcom == 1 ~ 1,
                                   symdiv == 1 ~ 0,
                                   econcom == 1 ~ 1,
                                   econdiv == 1 ~ 0),
         t_symbolic = case_when(symcom == 1 ~ 1,
                                symdiv == 1 ~ 1,
                                econcom == 1 ~ 0,
                                econdiv == 1 ~ 0)
         )

# Outcome var = anger after reading vignette
# angscale is a scale composed of angry and furiou
# items measured with 4 response categories and scaled to 0-1
dat %>% select(angry, furiou) %>% summary()

# Moderator = Italian identity
# itaid is a scale composed of itaimp, itaeff, itatyp, itaus
# items measured with 4 response categories and scaled to 0-1
dat %>% select(itaimp, itaeff, itatyp, itaus) %>% summary()

# same for European identity
dat %>% select(euimp, eueff, eutyp, euus) %>% summary()
table(dat$euimp)

# Create binary version of outcome and moderator
# using the midpoint of the 4-point scale as the threshold
dat <-
  dat %>%
  mutate(
    angry_bin = ifelse(angry + furiou > 1, 1, 0),
    itaid_bin = ifelse(itaimp + itaeff + itatyp + itaus > 2, 1, 0),
    euid_bin = ifelse(euimp + eueff + eutyp + euus > 2, 1, 0),
    eudic_bin = ifelse(eudic > 0.5, 1, 0)
    )


dat <- dat |>
  select(
    angry_bin, itaid_bin, t_commonality, north, satisf, sopscale, Corriere
  ) |>
  na.omit()


# save
write_csv(dat,
          "./data-raw/dp_replication_data_cleaned.csv")

delponte <- dat

usethis::use_data(delponte,
                  overwrite = TRUE,
                  version = 2)
