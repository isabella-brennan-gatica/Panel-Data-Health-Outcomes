rm(list=ls()) # clear workspace
setwd("C:/Users/.../")

library(plm)
library(lmtest)
library(sandwich)
library(tidyverse)
library(stargazer)
library(readr)
library(haven)
library(dplyr)

# Import Health data
hdta<- read_csv("Assignment/Assignment 2/Health Outcomes/Health Outcomes.csv",show_col_types = FALSE)


# View the structure of the imported data
ls.str(hdta)
view(hdta)
summary(hdta)
##########################################################################
#                        LOG VARIABLES :HEXP & GDPC                      #
##########################################################################

hdta$lhexp <- log(hdta$HEXP)
hdta$lgdpc <- log(hdta$GDPC)

##########################################################################
#             SET UP PANEL DATA  & VERIFY IF DATASET IS BALANCED         #
##########################################################################

# Let R know data is panel
hdta<- pdata.frame(hdta, index=c("COUNTRY","YEAR") )

pdim(hdta)
# RESULT: 
# Unbalanced Panel: n = 191, T = 1-5, N = 754

##########################################################################
#                   VERIFY: HOW UNBALNCED IS THE DATASET?                #
##########################################################################


# Count the number of unique time periods observed for each COUNTRY
country_obs <- hdta %>%
  group_by(COUNTRY) %>%
  summarise(period_reported = n_distinct(YEAR))

# Identify COUNTRIES observed for all time periods
Balanced_data <- country_obs$COUNTRY[country_obs$period_reported == max(country_obs$period_reported)]

# Count the number of COUNTRIES observed for all time periods
Total_countries<-length(country_obs$COUNTRY)
Balanced_countries<-length(Balanced_data)
Unbalanced_countries <- Total_countries-length(Balanced_data)


print(paste("There are",Total_countries, "countries represented in the data.", 
            "Of them,",Balanced_countries,"are observed for all periods and",
            Unbalanced_countries, "are not"))

##########################################################################
#       VERIFY IF VARIABLES ARE TIME-INVARIANT OR TIME VARIENT           #
##########################################################################
 
#CALCULATE STANDARD DEVIATIONS PER COUNTRY PER VARIABLE
summary_sd_stats <- hdta %>%
  group_by(COUNTRY) %>%
  summarise(
    sd_W_DALE = sd(DALE),
    sd_W_HEXP = sd(HEXP),
    sd_W_HC3 = sd(HC3),
    sd_W_GINI = sd(GINI),
    sd_W_TROPICS = sd(TROPICS),
    sd_W_POPDEN = sd(POPDEN),
    sd_W_OECD = sd(OECD),
    sd_W_GEFF = sd(GEFF),
    sd_W_VOICE = sd(VOICE),
    sd_W_PUBTHE = sd(PUBTHE),
    sd_W_GDPC = sd(GDPC),
    sd_W_COMP = sd(COMP),
  )%>%drop_na()

summary(summary_sd_stats)

#CALCULATE VARIANCE PER COUNTRY PER VARIABLE
summary_var_stats <- hdta %>%
  group_by(COUNTRY) %>%
  summarise(
    var_W_DALE = var(DALE),
    var_W_HEXP = var(HEXP),
    var_W_HC3 = var(HC3),
    var_W_GINI = var(GINI),
    var_W_TROPICS = var(TROPICS),
    var_W_POPDEN = var(POPDEN),
    var_W_OECD = var(OECD),
    var_W_GEFF = var(GEFF),
    var_W_VOICE = var(VOICE),
    var_W_PUBTHE = var(PUBTHE),
    var_W_GDPC = var(GDPC),
    var_W_COMP = var(COMP),
  )%>%drop_na()

summary(summary_var_stats)

##########################################################################
#             EXAMINING MOSTLY TIME-INVARIANT VARIABLE: GEFF             #
##########################################################################

Geff_variation <- with(hdta, table(COUNTRY, GEFF)) # Calculate the variation in the specific variable for each country

num_GEFF_variation <- sum(rowSums(Geff_variation > 0) > 1) # Count the number of countries with variation

print(paste("The GEFF variable varies for",num_GEFF_variation, "countries"))

##########################################################################
#                       OLS,POLS  w.o time dummies                       #
##########################################################################


OLS<- lm(DALE ~ HC3 + lhexp + GINI + GEFF + PUBTHE + TROPICS + lgdpc, data=hdta)
POLS<- plm(DALE ~ HC3 + lhexp + GINI + GEFF + PUBTHE + TROPICS + lgdpc, data=hdta, model="pooling")
POLS<- coeftest(POLS, vcov = vcovHC(POLS, type = "HC1", cluster = "group"))

ols_pols_label <- c("OLS", "POOLED OLS")

stargazer(OLS,POLS,type="text", digits = 4, column.labels = ols_pols_label)

##########################################################################
#                        FE & RE  w.o time dummies                       #
##########################################################################

FE<-plm(DALE ~ HC3 + lhexp + GINI + GEFF + PUBTHE + TROPICS + lgdpc, data=hdta, index = c("COUNTRY","YEAR"), model="within")
se_fe <- sqrt(diag(vcovHC(FE, type = "HC1", cluster = "group")))

RE<-plm(DALE ~ HC3 + lhexp + GINI + GEFF + PUBTHE + TROPICS + lgdpc, data=hdta, index = c("COUNTRY","YEAR"), model="random")
se_re <- sqrt(diag(vcovHC(RE, type = "HC1", cluster = "group")))

fe_re_label <- c("FIXED EFFECT", "RANDOM EFFECTS")

stargazer(FE,RE, type="text", digits = 4, column.labels = fe_re_label,se=list(se_fe, se_re))


##########################################################################
#                       OLS & POLS & FE & RE                             #
##########################################################################

all_label <- c("FIXED EFFECT", "RANDOM EFFECTS", "OLS", "POOLED OLS")

stargazer(FE,RE,POLS,OLS, type="text", digits = 4, column.labels = all_label,se=list(se_fe, se_re))

##########################################################################
#                               HAUSMAN TEST                             #
##########################################################################

# H0 is that the preferred model is random effects vs. the alternative the fixed effects 

phtest(FE,RE) 


##########################################################################
#                   FIXED EFFECTS WITH TIME DUMMIES                      #
##########################################################################
FE_timedummy<-plm(DALE ~ HC3 + lhexp + GINI + GEFF + PUBTHE + TROPICS + lgdpc + factor(YEAR), data=hdta, index = c("COUNTRY","YEAR"), model="within")
se_FE_timedummy <- sqrt(diag(vcovHC(FE_timedummy, type = "HC1", cluster = "group")))

fe_re_label <- c("FIXED EFFECT W TIME DUMMY","FIXED EFFECT")
stargazer(FE_timedummy,FE,type="text", digits = 4, column.labels = fe_re_label,se=list(se_FE_timedummy,se_fe))

##########################################################################
#                   F TEST BETWEEN FIXED EFFECTS MODELS                  #
##########################################################################

pFtest(FE_timedummy, FE)

# Bc the p-value is < 0.05 then the fixed effects model with time dummies is a better choice


