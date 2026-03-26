library(qgam)
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
# "C:/Users/apmc/OneDrive - Norwegian University of Life Sciences/BatLab Norway/Projects/CoastalMonitoring/Analyses/Outputs/Reed/GAMM_JudithF_2026-03-24"



### PREP FOR GAMM ###
#upload csv files
cm <- read.csv("C:/Users/apmc/OneDrive - Norwegian University of Life Sciences/BatLab Norway/Projects/CoastalMonitoring/Analyses/Judith 6 lakes/JF_inputs/final_cm 1.csv")
# 468400 obs of 39 vars

auto.data <- read.csv("C:/Users/apmc/OneDrive - Norwegian University of Life Sciences/BatLab Norway/Projects/CoastalMonitoring/Analyses/Judith 6 lakes/JF_inputs/JudithSites2024_all.csv")
# 769409 obs of 47 vars

#Logistic regression to determine cut-off, REPEAT WITH PPYG!!!
names(cm)
### PREP FOR GAMM ###



#Logistic regression to determine cut-off, REPEAT WITH PPYG!!!

enil_data <- cm %>%
  filter(AUTO.ID == "EPTNIL") %>%
  mutate(
    MATCH.RATIO = as.numeric(MATCH.RATIO),
    correct = ifelse(manual.id == "ENIL", 1, 0)) %>%
  filter(!is.na(MATCH.RATIO))


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
  mutate(
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
# Probability Cutoff_Ratio Calls_retained Calls_lost Percent_loss
# 1         0.5    0.0727363           7085          0         0.00
# 2         0.6    0.2486550           7085          0         0.00
# 3         0.7    0.4403525           7055         30         0.42
# 4         0.8    0.6742063           6448        637         8.99
# 5         0.9    1.0000000           2809       4276        60.35


# prep auto.data
auto.reduced <- auto.data %>%
  select(OUT.FILE.FS, DATE, TIME, HOUR, DATE.12, TIME.12, HOUR.12,
         AUTO.ID., site, MATCH.RATIO) %>%
  rename(IN.FILE = OUT.FILE.FS) %>%
  mutate(
    MATCH.RATIO = as.numeric(MATCH.RATIO),
    auto.id = case_when(
      AUTO.ID. == "PIPPYG" ~ "PPYG",
      AUTO.ID. == "EPTNIL" ~ "ENIL",
      TRUE ~ AUTO.ID.))

# define cut-offs
cutoff_enil  <- 0.66
cutoff_ppyg  <- 0.76

#filter autoIDs for species and cut-off
auto.filtered <- auto.reduced %>%
  filter(
    (auto.id == "ENIL" & MATCH.RATIO >= cutoff_enil) |
      (auto.id == "PPYG" & MATCH.RATIO >= cutoff_ppyg))

dim(auto.filtered)
# 483228     11

# remove manual IDs to avoid duplicates, insert remaining taxa from auto.id to final.id
manual.files <- unique(cm$IN.FILE)

auto.new <- auto.filtered %>%
  filter(!(IN.FILE %in% manual.files)) %>%
  mutate(
    final.id = auto.id)

# add final.id and insert taxa from manual.id
cm <- cm %>%
  mutate(
    final.id = manual.id)

# combine
combi_cm <- bind_rows(cm, auto.new)
dim(combi_cm)
# 524555     40
# Check that no manual ID rows were lost
#test 
cm$manual.id <- factor(cm$manual.id)
summary(cm)
summary(cm$manual.id)
#     ?   ENIL   MDAU   MYSP   NNOC   NoID  Noise   NYCT   PAUR   PISP   PNAT   PPYG   VMUR   NA's 
#    104  16058      3   6431     26    355   2227   3279     27    407   7287  14340      1 417855 
#417855 NAs - no manual ID
dim(cm)
468400 - 417855
# 50545 manual IDs

str(combi_cm)
# Check that no manual IDs were overwritten
combi_cm$manual.id <- factor(combi_cm$manual.id)
combi_cm$final.id <- factor(combi_cm$final.id)

# compare with the cm object size 
summary(combi_cm$final.id)
#     ?   ENIL   MDAU   MYSP   NNOC   NoID  Noise   NYCT   PAUR   PISP   PNAT   PPYG   VMUR   NA's 
#    104  56803      3   6431     26    355   2227   3279     27    407   7287  29750      1 417855 
summary(combi_cm$manual.id)
#     ?   ENIL   MDAU   MYSP   NNOC   NoID  Noise   NYCT   PAUR   PISP   PNAT   PPYG   VMUR   NA's 
#    104  16058      3   6431     26    355   2227   3279     27    407   7287  14340      1 474010 

474010-417855 # difference between NAs (no manual id)
# 56155
29750-14340 # difference between PPYG obs 
# 15410
56803 - 16058 # difference between ENIL obs 
# 40745

40745 + 15410 # GOOD, this adds up. 
# 56155

## but in the end there should be no NAs at this point. We should drop what is left. 
View(combi_cm)
table(combi_cm$auto.id) #newly added auto IDs



### GAMM WITH PNAT ###

# confirm date format and create the column night 
cm$DATE.12 <- as.Date(cm$DATE.12)
cm <- cm %>% 
  mutate(night = as.Date(DATE.12)) 


# only PNAT, group by site&night, count n of bat passes
pnat_nights <- cm %>%
  dplyr::filter(manual.id == "PNAT") %>%
  dplyr::group_by(SITE, night) %>%
  dplyr::summarize(batpass = dplyr::n(), .groups = "drop")


# create timeline
all_nights <- unique(cm$night)

#create a fully expanded grid with nights x site
full_grid <- expand.grid(
  SITE  = unique(cm$SITE),
  night = all_nights)


# insert pnat data
data_full <- full_grid %>%
  dplyr::left_join(pnat_nights, by = c("SITE", "night"))
# adds NAs which will be true zeroes



# insert true zeroes into grid
data_full$batpass[
  is.na(data_full$batpass)
] <- 0

sum(data_full$batpass)
# 7287


View(data_full)



### Lets do some more data exploration before we go any further

###fitting the model###

# GAMM needs numeric dates
data_full <- data_full %>%
  mutate(
    day = as.numeric(night - min(night)))

data_full$SITE <- as.factor(data_full$SITE)



### DAT EXPLORATION

sum(data_full$batpass )
# 7287
summary(combi_cm$final.id)
# also 7287 PNAT, so that is good. 

## Check the continuous variables
hist(data_full$day)
hist(data_full$batpass)

## Check factor variable(s)
summary(data_full)
# SITE        night               batpass            day       
# CM-04:92   Min.   :2024-07-01   Min.   :  0.00   Min.   : 0.00  
# CM-05:92   1st Qu.:2024-07-23   1st Qu.:  0.00   1st Qu.:22.75  
# CM-06:92   Median :2024-08-15   Median :  1.00   Median :45.50  
# CM-21:92   Mean   :2024-08-15   Mean   : 13.20   Mean   :45.50  
# CM-23:92   3rd Qu.:2024-09-07   3rd Qu.:  9.25   3rd Qu.:68.25  
# CM-42:92   Max.   :2024-09-30   Max.   :718.00   Max.   :91.00  

## There should *not* be 92 days for each detector

## OBS! ## 
## Need to double check that the nights where detectors were NOT working and we know that are not included here. 

ggplot(data_full, aes(x = night, y = batpass)) +
  geom_point() +
  facet_wrap(~SITE) + 
  coord_flip()

# For CM-23 you can see there are dates with data that should not have, for example
data_full %>% dplyr::filter(SITE == "CM-23") %>% droplevels() %>% 
  ggplot() +  geom_point(aes(x=night, y = batpass))

head(data_full)


# Equipment failure dates: 
# CM-04 -> 16 days in July 2024 (02.-17.07); 
# CM-05 -> 9 days in July 2024 (10.-18.07) and 1 in September (30.09)
# CM-06 -> 14 days in July 2024 (05.-18.07) and 2 in September (12.-13.09); 
# CM-21 -> 10 days in July 2024 (9.-18.07)
# CM-23 -> 5 days in August (19.-23.08)
# CM-42 -> 0


## Here you need to remove nights with detector failures
test <- data_full %>% dplyr::filter(
  !(SITE == "CM-04" & night %in% c("2024-07-01", "2024-07-17")),
  !(SITE == "CM-05" & night %in% c("2024-07-10", "2024-07-18")),
  !(SITE == "CM-05" & night == "2024-09-30"),
  !(SITE == "CM-06" & night %in% c("2024-05-07", "2024-05-18")),
  !(SITE == "CM-06" & night %in% c("2024-09-12", "2024-09-13")),
  !(SITE == "CM-21" & night %in% c("2024-09-09", "2024-09-18")),
  !(SITE == "CM-23" & night %in% c("2024-08-09", "2024-08-23")),
)
sum(test$batpass)
#7204 
# 7287- 7204 - 83 bat passes lost 

# find out which dates we need to adjust for the night-1 day effect, 
# then re-build your data object. 

#scatterplot
library(ggplot2)
ggplot(data_full, aes(x = night, y = batpass)) +
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
# 27.08
#that's a lot of overdispersion, move on to NB


#tried NB and adjusted k-value
nb <- gam(
  batpass ~ 
    s(day, k = 28) + s(SITE, bs = "re"), 
  family = nb(), #negative binomial
  data = data_full,
  method = "REML")


#check residuals with Dharma
library(promises)
library(DHARMa)
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

