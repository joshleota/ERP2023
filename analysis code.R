# ERP2023: The Exercise Recovery Project.
# R code for analyses

library(mgcv)
library(emmeans)
library(effsize)
library(lmerTest)
library(lme4)
library(gratia)
library(dplyr)
library(ggplot2)
library(qpdf)

# read data

data.1 <- readRDS("fullsample.rds")

starttime <- proc.time()

# sleep_onset_dec ----------------------------------------------------------------------------------------------

# run GAMM

system.time({
  gam_bio_sleep_onset_dec <- bam(sleep_onset_dec ~ workout_strain 
                                 + s(workout_offset_relative_habitual_new, by = workout_strain) + s(id_factor, bs='re') 
                                 + gender_c + s(age_c),
                                 data=data.1, discrete = TRUE)
})

# GAMM stats

k.check(gam_bio_sleep_onset_dec)
summary(gam_bio_sleep_onset_dec)

# Estimated marginal means comparisons

em_bio_sleep_onset_dec <- emmeans(gam_bio_sleep_onset_dec, consec ~ workout_strain|workout_offset_relative_habitual_new, 
                                  at = list(workout_offset_relative_habitual_new = c(seq(-10, 2, by = 2)), gender_c = 0, age_c = 0), data = data.1, adjust = "mvt")

# estimated marginal means 

em_bio_sleep_onset_dec$emmeans

# contrasts testing the dose-response relationship

em_bio_sleep_onset_dec$contrasts

# contrast effect sizes

eff_size(em_bio_sleep_onset_dec,sigma = sigma(gam_bio_sleep_onset_dec), edf = summary(em_bio_sleep_onset_dec$contrasts)$df[1])
# tst ----------------------------------------------------------------------------------------------

# run GAMM

system.time({
  gam_bio_tst <- bam(tst ~ workout_strain 
                     + s(workout_offset_relative_habitual_new, by = workout_strain) + s(id_factor, bs='re') 
                     + gender_c + s(age_c),
                     data=data.1, discrete = TRUE)
})

# GAMM stats

k.check(gam_bio_tst)
summary(gam_bio_tst)

# Estimated marginal means comparisons

em_bio_tst <- emmeans(gam_bio_tst, consec ~ workout_strain|workout_offset_relative_habitual_new, 
                      at = list(workout_offset_relative_habitual_new = c(seq(-10, 2, by = 2)), gender_c = 0, age_c = 0), data = data.1, adjust = "mvt")

# estimated marginal means 

em_bio_tst$emmeans

# contrasts testing the dose-response relationship

em_bio_tst$contrasts

# contrast effect sizes

eff_size(em_bio_tst,sigma = sigma(gam_bio_tst), edf = summary(em_bio_tst$contrasts)$df[1])
# se ----------------------------------------------------------------------------------------------

# run GAMM

system.time({
  gam_bio_se <- bam(se ~ workout_strain 
                    + s(workout_offset_relative_habitual_new, by = workout_strain) + s(id_factor, bs='re') 
                    + gender_c + s(age_c),
                    data=data.1, discrete = TRUE)
})

# GAMM stats

k.check(gam_bio_se)
summary(gam_bio_se)

# Estimated marginal means comparisons

em_bio_se <- emmeans(gam_bio_se, consec ~ workout_strain|workout_offset_relative_habitual_new, 
                     at = list(workout_offset_relative_habitual_new = c(seq(-10, 2, by = 2)), gender_c = 0, age_c = 0), data = data.1, adjust = "mvt")

# estimated marginal means 

em_bio_se$emmeans

# contrasts testing the dose-response relationship

em_bio_se$contrasts

# contrast effect sizes

eff_size(em_bio_se,sigma = sigma(gam_bio_se), edf = summary(em_bio_se$contrasts)$df[1])

# resting_heart_rate ----------------------------------------------------------------------------------------------

# run GAMM

system.time({
  gam_bio_resting_heart_rate <- bam(resting_heart_rate ~ workout_strain 
                                    + s(workout_offset_relative_habitual_new, by = workout_strain) + s(id_factor, bs='re') 
                                    + gender_c + s(age_c),
                                    data=data.1, discrete = TRUE)
})

# GAMM stats

k.check(gam_bio_resting_heart_rate)
summary(gam_bio_resting_heart_rate)

# Estimated marginal means comparisons

em_bio_resting_heart_rate <- emmeans(gam_bio_resting_heart_rate, consec ~ workout_strain|workout_offset_relative_habitual_new, 
                                     at = list(workout_offset_relative_habitual_new = c(seq(-10, 2, by = 2)), gender_c = 0, age_c = 0), data = data.1, adjust = "mvt")

# estimated marginal means 

em_bio_resting_heart_rate$emmeans

# contrasts testing the dose-response relationship

em_bio_resting_heart_rate$contrasts

# contrast effect sizes

eff_size(em_bio_resting_heart_rate,sigma = sigma(gam_bio_resting_heart_rate), edf = summary(em_bio_resting_heart_rate$contrasts)$df[1])

# HRV ----------------------------------------------------------------------------------------------

# run GAMM

system.time({
  gam_bio_HRV <- bam(HRV ~ workout_strain 
                     + s(workout_offset_relative_habitual_new, by = workout_strain) + s(id_factor, bs='re') 
                     + gender_c + s(age_c),
                     data=data.1, discrete = TRUE)
})

# GAMM stats

k.check(gam_bio_HRV)
summary(gam_bio_HRV)

# Estimated marginal means comparisons

em_bio_HRV <- emmeans(gam_bio_HRV, consec ~ workout_strain|workout_offset_relative_habitual_new, 
                      at = list(workout_offset_relative_habitual_new = c(seq(-10, 2, by = 2)), gender_c = 0, age_c = 0), data = data.1, adjust = "mvt")

# estimated marginal means 

em_bio_HRV$emmeans

# contrasts testing the dose-response relationship

em_bio_HRV$contrasts

# contrast effect sizes

eff_size(em_bio_HRV,sigma = sigma(gam_bio_HRV), edf = summary(em_bio_HRV$contrasts)$df[1])
