
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
# "model" is an R object that represents a fitted logistic regression model 
# "p" is a probability ratio (between 0 and 1)
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
# Probability Cutoff_Ratio Calls_retained Calls_lost Percent_loss
# 1         0.5    0.0727363           7085          0         0.00
# 2         0.6    0.2486550           7085          0         0.00
# 3         0.7    0.4403525           7055         30         0.42
# 4         0.8    0.6742063           6448        637         8.99
# 5         0.9    1.0000000           2809       4276        60.35


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
# 769409     11

# define cut-offs
cutoff_enil  <- 0.66
cutoff_ppyg  <- 0.76

## Before you filter, you should make sure you have not lost any manually verified recordings! 
temp <- cm %>% dplyr::filter(manual.id != "NA") 
temp1<- cm %>% dplyr::filter(manual.id != "NA") %>% dplyr::mutate(IN.FILE = factor(IN.FILE)) 
# 50545 obs 
summary(temp1)
# there are only 36920 unique files... where are these duplicates coming from? 

templist <- temp1$IN.FILE
templist <- unique(temp1$IN.FILE)
head(templist)
# 36920 obs unique files... why are there so many duplicates...? 

dup_rows <- duplicated(temp1$IN.FILE) | duplicated(temp1$IN.FILE, fromLast = TRUE)

# Subset the data frame to keep only those rows
duplicates <- temp1[dup_rows, ]
duplicates$SITE <- factor(duplicates$SITE)
# 25004 duplicates? 

summary(duplicates)
# CM-04 CM-05 CM-06 CM-21 CM-23 CM-42 
# 1626   760  6719  4868  1842  9189

dups <- duplicates %>% dplyr::select(c(IN.FILE, SITE, DATE, TIME, AUTO.ID, manual.id, behavior, final.id))
summary(duplicates$SITE)
## double check that these were all observations of multiple bat taxa in single recordings

dups1 <- dups %>% dplyr::select(IN.FILE, AUTO.ID, manual.id) %>% dplyr::distinct()
# There are no duplicates here now, so it is all different bat taxa detected in the manual verification process. Coool beans. 

## Now we want to make sure that these 36920 manually verified files are retained in the following filtering process: 
manual.files <- unique(temp$IN.FILE)
# 36920 unique file names 

## Save these separately with the other relevant data: 

manual.fix <- cm %>% dplyr::filter(IN.FILE %in% manual.files) %>% droplevels()
dim(manual.fix)
# 50545    39

#filter autoIDs for species and cut-off
auto.filtered <- auto.reduced %>%
  dplyr::filter(
    (auto.id == "ENIL" & MATCH.RATIO >= cutoff_enil) |
      (auto.id == "PPYG" & MATCH.RATIO >= cutoff_ppyg))

test <- auto.filtered %>% dplyr::filter(IN.FILE %in% manual.files)
#9218 obs of 11 variables. 
36920-9218 
# with your current pipeline, you lose 27702 manually verified files, which we do not want to have happen. 

# so lets try it another way: 

#first remove all the manually verified files
auto.fix <- auto.reduced %>% dplyr::filter(!IN.FILE %in% manual.files) 

dim(auto.fix)
# 732489     11
dim(auto.reduced)
# 769409     11

dim(auto.filtered)
# 483228     11

## Check that the math makes sense
# N obs in aut.reduced - n files with manual.ids
769409 - 36920
# 732489 - good it makes sense. 

## Now add the cutoffs to auto.fix 
auto.fix1 <- auto.fix %>%  dplyr::filter(
  (auto.id == "ENIL" & MATCH.RATIO >= cutoff_enil) |
    (auto.id == "PPYG" & MATCH.RATIO >= cutoff_ppyg))

dim(auto.fix1)
# 474010     11
dim(auto.filtered)
# 483228     11

## Now combine back with the manual ids 
summary(manual.fix)
summary(auto.fix1)

manual.fix1 <- manual.fix %>% dplyr::select(IN.FILE, DATE, TIME, HOUR, DATE.12, TIME.12, HOUR.12, AUTO.ID., SITE, MATCH.RATIO, auto.id, manual.id) %>% 
  dplyr::rename(site = SITE) %>% droplevels()

summary(manual.fix1)
summary(auto.fix1)
str(manual.fix1)
str(auto.fix1)

dim(manual.fix1)
# 50545

dim(auto.fix1)
# 474010

50545 + 474010
# 524555 
## Now these match 

summary(cm)
#housekeeping
auto.new1 <- auto.new1 %>% dplyr::mutate(AUTO.ID. = factor(AUTO.ID.), 
                                         auto.id = factor(auto.id), 
                                         manual.id = factor(manual.id)) 

summary(auto.new1)
# 474010   NAs in the manual ids, for the 474010 added autoids 
# The math checks out!

# Now final id - remember to not over-write your own manual ids. 
auto.new2 <- auto.new1 %>% dplyr::mutate(
  final.id = factor(dplyr::case_when(
    is.na(manual.id) ~ auto.id, 
    TRUE ~ manual.id)))

summary(auto.new2$manual.id)
summary(auto.new2$final.id)

#     ?   ENIL   MDAU   MYSP   NNOC   NoID  Noise   NYCT   PAUR   PISP   PNAT   PPYG   VMUR   NA's 
#    104  16058      3   6431     26    355   2227   3279     27    407   7287  14340      1 474010 
# > summary(auto.new2$final.id)
#   ENIL   PPYG      ?   MDAU   MYSP   NNOC   NoID  Noise   NYCT   PAUR   PISP   PNAT   VMUR 
# 423131  81277    104      3   6431     26    355   2227   3279     27    407   7287      1 

## Now we need to get back the information from the cm data object that we want for modellling 

# Combine the filtered autoids with the manual ids 
DF <- dplyr::full_join(auto.new2, manual.fix1)
dim(DF)
# 524555     13
# The math checks out!

summary(DF)
summary(DF$final.id)
str(DF)

meta <- cm %>% dplyr::select(c(IN.FILE, mean_temp, mean_wind)) %>% dplyr::distinct()

dim(meta)
# 454775 

summary(meta)

dim(DF)
# 524555     13

summary(DF)

combi_cm <- dplyr::left_join(DF, meta, by = "IN.FILE")
dim(combi_cm)
summary(combi_cm)
# 15435 filenames are not matching, resulting in 56155 NA weather columns...

test <- combi_cm %>% dplyr::filter(is.na(MATCH.RATIO)) 
summary(test) # - these have weather data... 

test2 <- combi_cm %>% dplyr::filter(is.na(mean_temp)) 

testfiles <- test2$IN.FILE
head(testfiles)

test3 <- cm %>% dplyr::filter(IN.FILE %in% testfiles) # how can these not exist in the cm files? 

############# Judith's original work flow: ############# 
# summary(auto.new1)
# ## These are the 
# auto.new <- auto.filtered %>%
#   dplyr::filter(!(IN.FILE %in% manual.files)) %>%
#   dplyr::mutate(
#     final.id = auto.id) 
# dim(auto.new)
# # 56155    12

# add final.id and insert taxa from manual.id # 
# cm <- cm %>%
#   mutate(
#     final.id = manual.id)

# # combine
# combi_cm <- bind_rows(cm, auto.new)
#############  #############  #############  ############# 


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