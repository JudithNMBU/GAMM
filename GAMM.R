
library(mgcViz)
library(DHARMa)
library(MetBrewer)
library(hms)
library(data.table)
library(tidyverse)
library(lubridate)
library(janitor)
library(kableExtra)
library(papeR)
library(skimr)
library(vtable)
library(RColorBrewer)
library(cowplot)
library(readxl)
library(writexl)
library(suncalc)
library(mgcv)
library(stringr)


## Setup output directory 
output <-"C:/Users/apmc/OneDrive - Norwegian University of Life Sciences/BatLab Norway/Projects/CoastalMonitoring/Analyses/Outputs/Reed" # where you want to save your data

file.name <- "GAMM_JudithF"

todays_date <- Sys.Date()

dir.name <- str_c(output,"/", file.name, "_", todays_date)
dir.name

output_today <- dir.name
output_today

dir.create(output_today)
output_today
# "C:/Users/apmc/OneDrive - Norwegian University of Life Sciences/BatLab Norway/Projects/CoastalMonitoring/Analyses/Outputs/Reed/GAMM_JudithF_2026-03-17"



### PREP FOR GAMM ###


#upload csv files
cm <- read.csv("C:/Users/apmc/OneDrive - Norwegian University of Life Sciences/BatLab Norway/Projects/CoastalMonitoring/Analyses/Judith 6 lakes/JF_inputs/final_cm 1.csv")
# 468400 obs of 39 vars

auto.data <- read.csv("C:/Users/apmc/OneDrive - Norwegian University of Life Sciences/BatLab Norway/Projects/CoastalMonitoring/Analyses/Judith 6 lakes/JF_inputs/JudithSites2024_all.csv")
# 769409 obs of 47 vars

#Logistic regression to determine cut-off, REPEAT WITH PPYG!!!
names(cm)


enil_data <- cm %>%
  dplyr::filter(AUTO.ID == "EPTNIL") %>%
  dplyr::mutate(MATCH.RATIO = as.numeric(MATCH.RATIO),
    correct = ifelse(manual.id == "ENIL", 1, 0)) %>%
  dplyr::filter(!is.na(MATCH.RATIO))
# 7085 obs of 40 vars 

enil_glm <- glm(
  correct ~ MATCH.RATIO,
  family = binomial,
  data = enil_data)

summary(enil_glm)
# Call:
#   glm(formula = correct ~ MATCH.RATIO, family = binomial, data = enil_data)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.1676     0.2057  -0.815    0.415    
# MATCH.RATIO   2.3048     0.2348   9.818   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 5606.5  on 7084  degrees of freedom
# Residual deviance: 5516.0  on 7083  degrees of freedom
# AIC: 5520

### DATA WRESTLE ###

# calculates the match ratio cut-off by inverting the logistic regression
get_cutoff <- function(p, model) {
  beta0 <- coef(model)[1]
  beta1 <- coef(model)[2]
  (log(p / (1 - p)) - beta0) / beta1
}

# target probabilities for which we want to know the match ratio
target_probs <- c(0.5, 0.6, 0.7, 0.8, 0.9)

# create table with probability column and fill with values from get-cutoff + lost and retained calls + percentage
loss_table <- data.frame(Probability = target_probs) %>%
  dplyr::mutate(
    Cutoff_Ratio = pmin(
      sapply(Probability, get_cutoff, model = enil_glm),
      1
    ),
    Calls_retained = sapply(Cutoff_Ratio, function(c) 
      sum(enil_data$MATCH.RATIO >= c)
    ),
    Calls_lost = sapply(Cutoff_Ratio, function(c) 
      sum(enil_data$MATCH.RATIO < c)
    ),
    Percent_loss = round(
      (Calls_lost / nrow(enil_data)) * 100, 
      2
    )
  )

print(loss_table)



# prep auto.data
auto.reduced <- auto.data %>%
  dplyr::select(OUT.FILE.FS, DATE, TIME, HOUR, DATE.12, TIME.12, HOUR.12,
         AUTO.ID., site, MATCH.RATIO) %>%
  dplyr::rename(IN.FILE = OUT.FILE.FS) %>%
  dplyr::mutate(
    MATCH.RATIO = as.numeric(MATCH.RATIO),
    auto.id = dplyr::case_when(
      AUTO.ID. == "PIPPYG" ~ "PPYG",
      AUTO.ID. == "EPTNIL" ~ "ENIL",
      TRUE ~ AUTO.ID.))

summary(auto.reduced)
dim(auto.reduced)


# define cut-offs
cutoff_enil  <- 0.66
cutoff_ppyg  <- 0.76

## Before you filter, you should make sure you have not lost anny manually verified recordings! 

## Note to self - start here and test if any manually verified passes are lost. 

#filter autoIDs for species and cut-off
auto.filtered <- auto.reduced %>%
  dplyr::filter(
    (auto.id == "ENIL" & MATCH.RATIO >= cutoff_enil) |
      (auto.id == "PPYG" & MATCH.RATIO >= cutoff_ppyg))

# remove manual IDs to avoid duplicates, insert remaining taxa from auto.id to final.id
manual.files <- unique(cm$IN.FILE)

## These are the 
auto.new <- auto.filtered %>%
  dplyr::filter(!(IN.FILE %in% manual.files)) %>%
  dplyr::mutate(
    final.id = auto.id)
dim(auto.new)
# 56155    12

# add final.id and insert taxa from manual.id
cm <- cm %>%
  mutate(
    final.id = manual.id)

# combine
combi_cm <- bind_rows(cm, auto.new)

View(combi_cm)
table(combi_cm$auto.id) #newly added auto IDs



### GAMM WITH PNAT ###

# confirm date format and create the column night 
cm$DATE.12 <- as.Date(cm$DATE.12)
cm <- cm %>% 
  mutate(night = as.Date(DATE.12))



# only PNAT, group by site&night, count n of bat passes
pnat_nights <- cm %>%
  filter(manual.id == "PNAT") %>%
  group_by(SITE, night) %>%
  dplyr::summarize(batpass = n(), .groups = "drop")


# create timeline
all_nights <- unique(cm$night)

#create a fully expanded grid with nights x site
full_grid <- expand.grid(
  SITE  = unique(cm$SITE),
  night = all_nights)


# insert pnat data
data_full <- full_grid %>%
  left_join(pnat_nights, by = c("SITE", "night"))
# adds NAs which will be true zeroes


# insert true zeroes into grid
data_full$batpass[
  is.na(data_full$batpass)
] <- 0

View(data_full)



###fitting the model###

# GAMM needs numeric dates
data_full <- data_full %>%
  mutate(
    day = as.numeric(night - min(night)))

data_full$SITE <- as.factor(data_full$SITE)


#scatterplot
ggplot(data_full, aes(x = day, y = batpass)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  labs(
    x = "Night",
    y = "PNAT bat passes per night")


#try Poisson
poisson <- gam(
  batpass ~ 
    s(day, k = 10) +
    s(SITE, bs = "re"),
  family = poisson,
  data = data_full,
  method = "REML")

#check overdipersion
deviance(poisson) / df.residual(poisson)
#that's a lot of overdispersion, move on to NB


#tried NB and adjusted k-value
nb <- gam(
  batpass ~ 
    s(day, k = 28) + s(SITE, bs = "re"), 
  family = nb(), #negative binomial
  data = data_full,
  method = "REML")


#check residuals with Dharma
res_nb <- simulateResiduals(nb)
plot(res_nb)
#left plot looks good but waves on the right side should optimally be straight


#check if k-value is appropriate (edf = effective degrees of freedom)
gam.check(nb)


# test for auto-correlation
acf(residuals(nb))
# no auto-correlation


testZeroInflation(res_nb)
# no zero inflation


# is there a significant seasonal trend?
summary(nb)



# Plot with mgcviz
v_nb <- getViz(nb)

plot(sm(v_nb, 1)) + 
  l_fitLine(colour = "blue") + 
  l_ciPoly(fill = "lightblue", alpha = 0.5) + 
  l_rug() +
  labs(title = "Seasonal Trend of PNAT Activity", 
       x = "Day of Sampling Period", 
       y = "Number of expected Bat Passes")