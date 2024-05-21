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
                                 + gender_c + s(age_c) + weekday_v_weekend + prior_sleep_onset_dec + strain_mean_c,
                                 data=data.1, discrete = TRUE)
})

# GAMM stats

k.check(gam_bio_sleep_onset_dec)
summary(gam_bio_sleep_onset_dec)

# Estimated marginal means comparisons

em_bio_sleep_onset_dec <- emmeans(gam_bio_sleep_onset_dec, consec ~ workout_strain|workout_offset_relative_habitual_new, 
                                  at = list(workout_offset_relative_habitual_new = c(seq(-10, 2, by = 2)), gender_c = 0, age_c = 0, 
                                            weekday_v_weekend = 0, prior_sleep_onset_dec = 0, strain_mean_c = 0), data = data.1, adjust = "mvt")

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
                     + gender_c + s(age_c) + weekday_v_weekend + prior_tst + strain_mean_c,
                     data=data.1, discrete = TRUE)
})

# GAMM stats

k.check(gam_bio_tst)
summary(gam_bio_tst)

# Estimated marginal means comparisons

em_bio_tst <- emmeans(gam_bio_tst, consec ~ workout_strain|workout_offset_relative_habitual_new, 
                      at = list(workout_offset_relative_habitual_new = c(seq(-10, 2, by = 2)), gender_c = 0, age_c = 0,
                                weekday_v_weekend = 0, prior_tst = 0, strain_mean_c = 0), data = data.1, adjust = "mvt")

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
                    + gender_c + s(age_c) + weekday_v_weekend + prior_se + strain_mean_c,
                    data=data.1, discrete = TRUE)
})

# GAMM stats

k.check(gam_bio_se)
summary(gam_bio_se)

# Estimated marginal means comparisons

em_bio_se <- emmeans(gam_bio_se, consec ~ workout_strain|workout_offset_relative_habitual_new, 
                     at = list(workout_offset_relative_habitual_new = c(seq(-10, 2, by = 2)), gender_c = 0, age_c = 0,
                               weekday_v_weekend = 0, prior_se = 0, strain_mean_c = 0), data = data.1, adjust = "mvt")

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
                                    + gender_c + s(age_c) + weekday_v_weekend + prior_resting_heart_rate + strain_mean_c,
                                    data=data.1, discrete = TRUE)
})

# GAMM stats

k.check(gam_bio_resting_heart_rate)
summary(gam_bio_resting_heart_rate)

# Estimated marginal means comparisons

em_bio_resting_heart_rate <- emmeans(gam_bio_resting_heart_rate, consec ~ workout_strain|workout_offset_relative_habitual_new, 
                                     at = list(workout_offset_relative_habitual_new = c(seq(-10, 2, by = 2)), gender_c = 0, age_c = 0,
                                               weekday_v_weekend = 0, prior_resting_heart_rate = 0, strain_mean_c = 0), data = data.1, adjust = "mvt")

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
                     + gender_c + s(age_c) + weekday_v_weekend + prior_HRV + strain_mean_c,
                     data=data.1, discrete = TRUE)
})

# GAMM stats

k.check(gam_bio_HRV)
summary(gam_bio_HRV)

# Estimated marginal means comparisons

emm_options(opt.digits = FALSE)

em_bio_HRV <- emmeans(gam_bio_HRV, consec ~ workout_strain|workout_offset_relative_habitual_new, 
                      at = list(workout_offset_relative_habitual_new = c(seq(-10, 2, by = 2)), gender_c = 0, age_c = 0,
                                weekday_v_weekend = 0, prior_HRV = 0, strain_mean_c = 0), data = data.1, adjust = "mvt")

# estimated marginal means 

em_bio_HRV$emmeans

# contrasts testing the dose-response relationship

em_bio_HRV$contrasts

# contrast effect sizes

eff_size(em_bio_HRV,sigma = sigma(gam_bio_HRV), edf = summary(em_bio_HRV$contrasts)$df[1])

# Figures ----------------------------------------------------------------------------------------

# example of data creation from GAMMs to visualise associations:

# create new dataset to create predicted values to plot

newdata_sleep_onset_dec <- expand.grid(workout_strain = c("none","light","moderate","high","All out"), 
                                       workout_offset_relative_habitual_new = seq(min(data.1$workout_offset_relative_habitual_new), max(data.1$workout_offset_relative_habitual_new), by = 1/60),
                                       gender_c = 0,
                                       age_c = 0,
                                       id_factor = 0,
                                       weekday_v_weekend = 0, 
                                       prior_sleep_onset_dec = 0,
                                       strain_mean_c = 0)

yhat <- as.data.frame(predict.gam(gam_bio_sleep_onset_dec, newdata = newdata_sleep_onset_dec, se.fit = TRUE)))
newdata_sleep_onset_dec$est <- yhat$fit
newdata_sleep_onset_dec$lower_ci <- yhat$fit - 1.96 * yhat$se.fit
newdata_sleep_onset_dec$upper_ci <- yhat$fit + 1.96 * yhat$se.fit

# creation of figures:

colours <- c("#003C67FF", "#7AA6DCFF", "#FDD262", "#CD534CFF")
colours_r <- c("#003C67FF", "#7AA6DCFF", "#FDD262", "#CD534CFF")

labels <- c("Light", "Moderate", "High", "Maximal")
linetypes <- c(rep("solid", length(unique(newdata_sleep_onset_dec$workout_strain))))

# sleep_onset_dec ----------------------------------------------------------------------------------------------

newdata_sleep_onset_dec <- readRDS("newdata_sleep_onset_dec.rds")

none_sleep_onset_dec <- newdata_sleep_onset_dec$est[1]

newdata_sleep_onset_dec <- newdata_sleep_onset_dec %>%
  filter(workout_strain != "none")

sleep_onset_dec_plot <- ggplot(newdata_sleep_onset_dec, aes(x = workout_offset_relative_habitual_new, y = est, color = workout_strain)) +
  geom_vline(xintercept = 0, color = c("black"), lwd = 1.2, linetype = "dotted") + 
  geom_hline(yintercept = none_sleep_onset_dec, color = c("#868686FF"), lwd = 1.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = workout_strain), alpha = 0.15, color = NA) +  
  geom_line(aes(linetype = workout_strain), lwd = 1.2) +
  labs(x = "", y = "Sleep Onset (Hours)") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linewidth = 2),
        axis.text = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 12),
        # legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  scale_color_manual(values = colours,
                     labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_fill_manual(values = colours_r,
                    labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_linetype_manual(values = linetypes,
                        labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_x_continuous(breaks = seq(min(-10), max(2), by = 2), limits = c(-10,2), labels = c(-10,-8,-6,-4,-2,"0","+2")) +
  scale_y_continuous(breaks = seq(min(23), max(28), by = 1), limits = c(22.8,28.2),
                     labels = c("11 PM","12 AM","1 AM","2 AM","3 AM","4 AM"))
# + coord_fixed(ratio = 1.2)

sleep_onset_dec_plot

pdf("sleep_onset_dec_plot_new.pdf")
print(sleep_onset_dec_plot)
dev.off()

ggsave("sleep_onset_dec_plot.pdf", plot = sleep_onset_dec_plot, width = 8, height = 5)

# tst ----------------------------------------------------------------------------------------------

newdata_tst <- readRDS("newdata_tst.rds")

none_tst <- newdata_tst$est[1]

newdata_tst <- newdata_tst %>%
  filter(workout_strain != "none")

labels <- c("Light", "Moderate", "High", "Maximal")
linetypes <- c(rep("solid", length(unique(newdata_tst$workout_strain))))

tst_plot <- ggplot(newdata_tst, aes(x = workout_offset_relative_habitual_new, y = est, color = workout_strain)) +
  geom_vline(xintercept = 0, color = c("black"), lwd = 1.2, linetype = "dotted") + 
  geom_hline(yintercept = none_tst, color = c("#868686FF"), lwd = 1.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = workout_strain), alpha = 0.15, color = NA) +  
  geom_line(aes(linetype = workout_strain), lwd = 1.2) +
  labs(x = "", y = "Sleep Duration (Hours)") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linewidth = 2),
        axis.text = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 12),
        plot.margin = margin(0, 20, 0, 0, "pt"),
        legend.position = "none") +
  scale_color_manual(values = colours,
                     labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_fill_manual(values = colours_r,
                    labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_linetype_manual(values = linetypes,
                        labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_x_continuous(breaks = seq(min(-10), max(2), by = 2), limits = c(-10,2), labels = c(-10,-8,-6,-4,-2,"0","+2")) +
  scale_y_continuous(breaks = seq(min(round(newdata_tst$lower_ci)), max((round(newdata_tst$upper_ci))), by = 0.5))
# + coord_fixed(ratio = 1.2)

tst_plot

pdf("tst_plot_new.pdf")
print(tst_plot)
dev.off()

ggsave("tst_plot.pdf", plot = tst_plot, width = 8, height = 5)

# se ----------------------------------------------------------------------------------------------

newdata_se <- readRDS("newdata_se.rds")

none_se <- newdata_se$est[1]

newdata_se <- newdata_se %>%
  filter(workout_strain != "none")

labels <- c("Light", "Moderate", "High", "Maximal")
linetypes <- c(rep("solid", length(unique(newdata_se$workout_strain))))

se_plot <- ggplot(newdata_se, aes(x = workout_offset_relative_habitual_new, y = est, color = workout_strain)) +
  geom_vline(xintercept = 0, color = c("black"), lwd = 1.2, linetype = "dotted") + 
  geom_hline(yintercept = none_se, color = c("#868686FF"), lwd = 1.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = workout_strain), alpha = 0.15, color = NA) +  
  geom_line(aes(linetype = workout_strain), lwd = 1.2) +
  labs(x = "Exercise Ending Relative to Habitual Sleep Onset (Hours)",
       y = "Sleep Percentage (%)") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linewidth = 2),
        axis.text = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 12),
        # legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  scale_color_manual(values = colours,
                     labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_fill_manual(values = colours_r,
                    labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_linetype_manual(values = linetypes,
                        labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_x_continuous(breaks = seq(min(-10), max(2), by = 2), limits = c(-10,2), labels = c(-10,-8,-6,-4,-2,"0","+2")) +
  scale_y_continuous(breaks = seq(min(round(newdata_se$lower_ci)), max((round(newdata_se$upper_ci))), by = 1))
# + coord_fixed(ratio = 1.2)

se_plot

pdf("se_plot_new.pdf")
print(se_plot)
dev.off()

ggsave("se_plot.pdf", plot = se_plot, width = 8, height = 5)

# resting_heart_rate ----------------------------------------------------------------------------------------------

newdata_resting_heart_rate <- readRDS("newdata_resting_heart_rate.rds")

none_resting_heart_rate <- newdata_resting_heart_rate$est[1]

newdata_resting_heart_rate <- newdata_resting_heart_rate %>%
  filter(workout_strain != "none")

resting_heart_rate_plot <- ggplot(newdata_resting_heart_rate, aes(x = workout_offset_relative_habitual_new, y = est, color = workout_strain)) +
  geom_vline(xintercept = 0, color = c("black"), lwd = 1.2, linetype = "dotted") + 
  geom_hline(yintercept = none_resting_heart_rate, color = c("#868686FF"), lwd = 1.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = workout_strain), alpha = 0.15, color = NA) +  
  geom_line(aes(linetype = workout_strain), lwd = 1.2) +
  labs(x = "Exercise Ending Relative to Habitual Sleep Onset (Hours)",
       y = "Nocturnal Resting Heart Rate (BPM)") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linewidth = 2),
        axis.text = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 11),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none",
        plot.margin = margin(0, 20, 0, 0, "pt")) +
  scale_color_manual(values = colours,
                     labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_fill_manual(values = colours_r,
                    labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_linetype_manual(values = linetypes,
                        labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_x_continuous(breaks = seq(min(-10), max(2), by = 2), limits = c(-10,2), labels = c(-10,-8,-6,-4,-2,"0","+2")) +
  scale_y_continuous(breaks = seq(min(55), max(75), by = 5),limits = c(55,75))
# + coord_fixed(ratio = 1.2)

resting_heart_rate_plot

pdf("resting_heart_rate_plot_new.pdf")
print(resting_heart_rate_plot)
dev.off()

ggsave("resting_heart_rate_plot.pdf", plot = resting_heart_rate_plot, width = 8, height = 5)

# HRV ----------------------------------------------------------------------------------------------

newdata_HRV <- readRDS("newdata_HRV.rds")

none_HRV <- newdata_HRV$est[1]

newdata_HRV <- newdata_HRV %>%
  filter(workout_strain != "none")

HRV_plot <- ggplot(newdata_HRV, aes(x = workout_offset_relative_habitual_new, y = est, color = workout_strain)) +
  geom_vline(xintercept = 0, color = c("black"), lwd = 1.2, linetype = "dotted") + 
  geom_hline(yintercept = none_HRV, color = c("#868686FF"), lwd = 1.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = workout_strain), alpha = 0.15, color = NA) +  
  geom_line(aes(linetype = workout_strain), lwd = 1.2) +
  labs(x = "Exercise Ending Relative to Habitual Sleep Onset (Hours)",
       y = "Noctural HRV (RMSSD)") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linewidth = 2),
        axis.text = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 11),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 12),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  scale_color_manual(values = colours,
                     labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_fill_manual(values = colours_r,
                    labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_linetype_manual(values = linetypes,
                        labels = labels, guide = guide_legend(title = "Exercise Strain", title.theme = element_text(size = 12, face = "bold"))) +
  scale_x_continuous(breaks = seq(min(-10), max(2), by = 2), limits = c(-10,2), labels = c(-10,-8,-6,-4,-2,"0","+2")) +
  scale_y_continuous(breaks = seq(min(20), max(60), by = 5))
# + coord_fixed(ratio = 1.2)

HRV_plot

pdf("HRV_plot_new.pdf")
print(HRV_plot)
dev.off()

ggsave("HRV_plot.pdf", plot = HRV_plot, width = 8, height = 5)

sleep <- ggarrange(sleep_onset_dec_plot,tst_plot,se_plot, nrow = 3)

ggsave("figure_1.pdf", plot = sleep, width = 8, height = 15)


sleep_onset_dec_plot/tst_plot/se_plot +
  plot_annotation(tag_levels = 'A')

resting_heart_rate_plot/HRV_plot +
  plot_annotation(tag_levels = 'A')

all <- (plot_spacer() | sleep_onset_dec_plot) / 
  (tst_plot | se_plot) /
  (resting_heart_rate_plot | HRV_plot) + plot_annotation(tag_levels = 'A', tag_prefix = '1') & 
  theme(plot.tag = element_text(face = "bold")) 

all

