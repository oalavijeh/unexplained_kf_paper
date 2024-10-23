#reviewer analysis 
#first calculate models and then work out AUC and bootstrp confidence intervals 

library(pscl, lib.loc = "~/re_gecip/renal/Rpackages/")
library(pwr, lib.loc = "~/re_gecip/renal/Rpackages/")
library(ggplot2)
library(dplyr, lib.loc = "~/re_gecip/renal/Rpackages/")
library(logistf)
library(data.table, lib.loc = "~/re_gecip/renal/Rpackages/")
library(FSA)
library(oddsratio, lib.loc = "~/re_gecip/renal/Rpackages/")
library(msm)
library(ggpubr)
#defo need this one for multinomial logistic regression 
library(nnet)
options(bitmapType="cairo")

options(scipen=0)

# Calculate liability threshold heritability
#K=pop prevalence
#P=proportion of cases in study 
#hsq=Heritability estimate (on observed scale)
#bigT = liability threshold
#tau = density of gaussian

h2_liab <- function(x){
  K=0.001
  P=0.04632894
  zv <- dnorm(qnorm(K))
  x * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
}
#read in prs
setwd("/re_gecip/renal/oalavijeh/projects/unexplained_esrf/prs/")
# # Read in the prs file
ckd_prs <- fread("unexplained_esrf_controls_prs_apol1_status_covariates.txt")
colnames(ckd_prs)[14] <- "ckd_prs"
iga_prs <- fread("unexplained_iga_77snp_prs_scores.sscore")
colnames(iga_prs)[1] <- "id"
colnames(iga_prs)[3] <- "iga_prs"
ssns_prs <- fread("ssns/unexplained_ssns_prs_scores.sscore")
colnames(ssns_prs)[1] <- "id"
colnames(ssns_prs)[3] <- "ssns_prs"
membranous_prs <- fread("membranous/unexplained_membranous_prs_scores.sscore")
colnames(membranous_prs)[1] <- "id"
colnames(membranous_prs)[3] <- "membranous_prs"
# 
# #merge all 
prs <- Reduce(merge,list(ckd_prs, iga_prs[,c(1,3)],membranous_prs[,c(1,3)],ssns_prs[,c(1,3)]))
# #get rid of monogenic 
monogenic <- fread("/re_gecip/renal/oalavijeh/projects/unexplained_esrf/apol1/all_solved_cases_with_apol1.txt")
# prs <- subset(prs, !id %in% monogenic$plate_key)
prs$monogenic <- ifelse(prs$id %in% monogenic$plate_key,1,0)
# #standardise for each score
control.prs <- subset(prs, pheno==0)
prs$iga_prs_standardised <- (prs$iga_prs-mean(control.prs$iga_prs))/sd(control.prs$iga_prs)
prs$ssns_prs_standardised <- (prs$ssns_prs-mean(control.prs$ssns_prs))/sd(control.prs$ssns_prs)
prs$membranous_prs_standardised <- (prs$membranous_prs-mean(control.prs$membranous_prs))/sd(control.prs$membranous_prs)
prs$ckd_prs_standardised <- (prs$ckd_prs-mean(control.prs$ckd_prs))/sd(control.prs$ckd_prs)
# apol1 types 
# apol1$status[apol1$gt_final=="G1aHET:G1bHET" & apol1$type == "aa_esrf"] <- "Case - Non Rx APOL1"
# apol1$status[apol1$gt_final=="G1aHET" & apol1$type == "aa_esrf"] <- "Case - Non Rx APOL1"
# apol1$status[apol1$gt_final=="G1aHET:G1bHET:G2HET" & apol1$type == "aa_esrf"] <- "Case - High Rx APOL1"
# apol1$status[apol1$gt_final=="G1aHOM:G1bHOM" & apol1$type == "aa_esrf"] <- "Case - High Rx APOL1"
# apol1$status[apol1$gt_final=="G2HET" & apol1$type == "aa_esrf"] <- "Case - Non Rx APOL1"
# apol1$status[apol1$gt_final=="G2HOM" & apol1$type == "aa_esrf"] <- "Case - High Rx APOL1"

# #make apol1 column 
prs$apol1 <- ifelse(prs$type=="Case - High Rx APOL1"| prs$type=="Control - High Rx APOL1",1,0)
#add in age
age <- fread("/re_gecip/renal/oalavijeh/projects/phenotypes/final_phenotypes/all_agg_ancestry.tsv")
age$age <- 2024-age$year_of_birth
prs_age <- merge(prs, age[,c(2,12)], by.x = "id", by.y= "plate_key", all.x = T, all.y = F)
prs <- prs_age[,c(1,3,2,25,4:13,15,19,20:24)]
# #write_table
#write.table(prs, "unexplained_esrf/prs/unexplained_case_rdcontrol_iga_mem_ssns_ckd_prs_standardised_pcs_age_sex_apol1_massns_monogenic.txt", col.names = T, row.names = F, quote = F, sep = '\t')
####start here with completed files #####
prs <- fread("unexplained_esrf/prs/unexplained_case_rdcontrol_iga_mem_ssns_ckd_prs_standardised_pcs_age_sex_apol1_massns_monogenic.txt")
#switch out apol1s with low risk variants 
prs$apol1 <- ifelse(prs$id=="x"|
                      prs$id=="x"|
                      prs$id=="x"|
                      prs$id=="x", 0,prs$apol1)

#make pheno files for modelling with variants analysis
formodel.prs <- prs[,c(2:14,17:22)]

#modelling - first do cases versus controls only ten multinomial lr
# We can then calculate the null model (model with PRS) using logistic regression 
formodel.null.prs <- formodel.prs[,c(1:13)]
null.model <- glm(pheno~., family=binomial(link='logit'),data=formodel.null.prs, maxit=100)
# And the R2 of the null model is 
null.r2 <- pR2(null.model)
prs.result <- NULL
#fit all patients model
model <- glm(pheno~., family=binomial(link='logit'),data=formodel.prs, maxit=1000)

# model R2 is obtained as 
model.r2 <- pR2(model)[6]
# R2 of PRS is simply calculated as the model R2 minus the null R2
prs.r2 <- model.r2-null.r2
#prs.h2_liab <- h2_liab(prs.r2)
# We can also obtain the coeffcient and p-value of association of PRS as follow
prs.ckd.coef <- summary(model)$coeff["ckd_prs_standardised",]
prs.beta <- as.numeric(prs.ckd.coef[1])
prs.se <- as.numeric(prs.ckd.coef[2])
prs.p <- as.numeric(prs.ckd.coef[4])
# We can then store the results
#prs.result <- rbind(prs.result, data.frame(R2=prs.r2, H2=prs.h2_liab, P=prs.p, BETA=prs.beta,SE=prs.se))
prs.result <- rbind(prs.result, data.frame(R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))

# Best result is:
prs.result[which.max(prs.result$R2),]
#get model data
summary(model)

#forest_plot
#run logistic models
#apol1
model_apol1 <- glm(pheno~apol1 +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+sex+age, data=prs, family = 'binomial')
summary_apol1 <- summary(model_apol1)
OR_apol1<-exp(coef(summary_apol1)["apol1","Estimate"])
CI_apol1<- exp(confint(model_apol1, level=0.95)["apol1",])
CI_lower_apol1 <- CI_apol1[1]
CI_upper_apol1 <- CI_apol1[2]

#iga
model_iga <- glm(pheno~iga_prs_standardised +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+sex+age, data=prs, family = 'binomial')
summary_iga <- summary(model_iga)
OR_iga<-exp(coef(summary_iga)["iga_prs_standardised","Estimate"])
CI_iga<- exp(confint(model_iga, level=0.95)["iga_prs_standardised",])
CI_lower_iga <- CI_iga[1]
CI_upper_iga <- CI_iga[2]

#ssns
model_ssns <- glm(pheno~ssns_prs_standardised +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+sex+age, data=prs, family = 'binomial')
summary_ssns <- summary(model_ssns)
OR_ssns<-exp(coef(summary_ssns)["ssns_prs_standardised","Estimate"])
CI_ssns<- exp(confint(model_ssns, level=0.95)["ssns_prs_standardised",])
CI_lower_ssns <- CI_ssns[1]
CI_upper_ssns <- CI_ssns[2]

#membranous
model_membranous <- glm(pheno~membranous_prs_standardised +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+sex+age, data=prs, family = 'binomial')
summary_membranous <- summary(model_membranous)
OR_membranous<-exp(coef(summary_membranous)["membranous_prs_standardised","Estimate"])
CI_membranous<- exp(confint(model_membranous, level=0.95)["membranous_prs_standardised",])
CI_lower_membranous <- CI_membranous[1]
CI_upper_membranous <- CI_membranous[2]

#CKD
model_ckd <- glm(pheno~ckd_prs_standardised +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+sex+age, data=prs, family = 'binomial')
summary_ckd <- summary(model_ckd)
OR_ckd<-exp(coef(summary_ckd)["ckd_prs_standardised","Estimate"])
CI_ckd<- exp(confint(model_ckd, level=0.95)["ckd_prs_standardised",])
CI_lower_ckd <- CI_ckd[1]
CI_upper_ckd <- CI_ckd[2]

#monogenic
model_monogenic <- glm(pheno~monogenic +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+sex+age, data=prs, family = 'binomial')
summary_monogenic <- summary(model_monogenic)
OR_monogenic<-exp(coef(summary_monogenic)["monogenic","Estimate"])
CI_monogenic<- exp(confint(model_monogenic, level=0.95)["monogenic",])
CI_lower_monogenic <- CI_monogenic[1]
CI_upper_monogenic <- CI_monogenic[2]

#summary_data - with monogenic
summary_stats <- data.frame(
  Variable=c("APOL1","Monogenic","IgA PRS","SSNS PRS","Membranous PRS", "CKD PRS"),
  OR=c(OR_apol1,OR_monogenic, OR_iga, OR_ssns, OR_membranous,OR_ckd),
  CI_Lower=c(CI_lower_apol1,CI_lower_monogenic, CI_lower_iga, CI_lower_ssns, CI_lower_membranous,CI_lower_ckd),
  CI_Upper=c(CI_upper_apol1, CI_upper_monogenic,CI_upper_iga, CI_upper_ssns,CI_upper_membranous,CI_upper_ckd)
  )

#write.table(summary_stats, "with_age_summary_stats_for_prs_analysis_for_forest_monogenic.txt", col.names = T, row.names = F, quote = F, sep ='\t')
write.table(summary_stats,"unexplained_esrf/prs/summary_stats_ofr_prs_analysis_for_forest_protective_apol1_removed.txt",
            col.names = T, row.names = F, sep = '\t', quote = F)

summary_stats <- fread("unexplained_esrf/prs/summary_stats_ofr_prs_analysis_for_forest_protective_apol1_removed.txt")
#dont need monogenic 
summary_stats <- summary_stats[c(1,3,4,5,6),]

#plot in ggplot2 
# Load necessary library
library(ggplot2)

# Add a column for plotting intervals
#summary_stats$Pvalue <- c(5.6e-8,0.11,0.06,0.97,0.73)
summary_stats$Variable <- factor(summary_stats$Variable, levels = summary_stats$Variable[order(summary_stats$OR)])
summary_stats$y <- seq_along(summary_stats$Variable)
#summary_stats$Variable_P <- paste0(summary_stats$Variable, "(", summary_stats$Pvalue,")")
#summary_stats$PvalueFormatted <- sprintf("p = %.3f", summary_stats$Pvalue)
# Create the plot
ggplot(summary_stats, aes(x = OR, y = Variable, color=Variable)) +
  geom_point() +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Forest Plot", x = "Odds Ratio", y = "")
  #geom_text(aes(label = PvalueFormatted, x = Inf), hjust = -0.1, color = "black") # Position text to the right
#theme(plot.margin = margin(1, 100, 1, 1)) # Adjust right margin to ensure P-values are visible

#compare means for SSNS Across african cohorts
african <- fread("unexplained_esrf/apol1/apol1_status_of_esrf_cases_controlsv2.csv")
african_prs <- merge(prs, african[,c(1,3)], by.x = "id", by.y = "sample", all.x = F, all.y = T)
african_prs <- na.omit(african_prs)
ukf_africans_apol1 <- subset(african_prs, african_prs$type.y=="aa_esrf" & african_prs$apol1==1)
aa_control_no_apol1 <- subset(african_prs, african_prs$type.y=="aa_control" & african_prs$apol1==0)
ukf_aa_apol1_aa_control_no_apol1 <- rbind(ukf_africans_apol1, aa_control_no_apol1)

aa_control_apol1 <- subset(african_prs, african_prs$type.y=="aa_control" & african_prs$apol1==1)
ukf_aa_apol1_aa_control_apol1 <- rbind(ukf_africans_apol1, aa_control_apol1)
#run models - appol1 vs controls without apol1
model_ukf_aa_apol1_aa_control_no_apol1 <- ukf_aa_apol1_aa_control_no_apol1[,c(2:14,17:20,22)]
model <- glm(pheno~., family=binomial(link='logit'),data=model_ukf_aa_apol1_aa_control_no_apol1, maxit=1000)
summary(model)

#apol1 vs controls with apol1
model_ukf_aa_apol1_aa_control_apol1 <- ukf_aa_apol1_aa_control_apol1[,c(2:14,17:20,22)]
model <- glm(pheno~., family=binomial(link='logit'),data=model_ukf_aa_apol1_aa_control_apol1, maxit=1000)
summary(model)

#plot african apol1 
african_prs$apol1_type <- NA
african_prs$apol1_type <- ifelse(african_prs$type.y=="aa_esrf" & african_prs$apol1==1, "Case APOL1",african_prs$apol1_type)
african_prs$apol1_type <- ifelse(african_prs$type.y=="aa_esrf" & african_prs$apol1==0, "Case",african_prs$apol1_type)
african_prs$apol1_type <- ifelse(african_prs$type.y=="aa_control" & african_prs$apol1==1, "Control APOL1",african_prs$apol1_type)
african_prs$apol1_type <- ifelse(african_prs$type.y=="aa_control" & african_prs$apol1==0, "Control",african_prs$apol1_type)

for_plotting <- subset(african_prs, african_prs$apol1_type=="Case APOL1" | 
                         african_prs$apol1_type=="Control APOL1" |
                         african_prs$apol1_type=="Control")


my_comparisons <- list(c("Case APOL1","Control APOL1"),c("Case APOL1","Control"),c("Control", "Control APOL1"))
ggviolin(for_plotting, x = "apol1_type", y = "sns", fill = "apol1_type",
         palette = c("#00AFBB", "#FC4E07","#f2241f"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")# Add significance levels
stat_compare_means(label.y = 8) +
  #stat_cor(p.accuracy = 0.01) +
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") 

kruskal.test(ssns ~ apol1_type, data=for_plotting)
pairwise.wilcox.test(for_plotting$ssns, for_plotting$apol1_type, p.adjust.method = "BH")

kruskal.test(ssns_prs_standardised ~ apol1_type, data=for_plotting)
pairwise.wilcox.test(for_plotting$ssns_prs_standardised, for_plotting$apol1_type, p.adjust.method = "BH")

#prs plot with monogenic, APOL1 and unsolved as cohorts 
library(ggplot2)
library(tidyr)
library(dplyr)

# Reshape the data from wide to long format
prs$ckd_apol_adjusted <- ifelse(prs$apol1==1, prs$ckd_prs_standardised-1, prs$ckd_prs_standardised)
df_long <- prs %>%
  pivot_longer(cols = c("iga_prs_standardised", "ssns_prs_standardised",
                        "membranous_prs_standardised", "ckd_apol_adjusted"), 
               names_to = "PRS_score", values_to = "PRS_value") %>%
  mutate(cohort = case_when(
    pheno == 1 & apol1 == 1 ~ "APOL1", # APOL1_Positive takes priority
    pheno == 1 & monogenic == 1 ~ "Monogenic",
    pheno == 1 & monogenic == 0 ~ "Unsolved",    TRUE ~ as.character(NA)
  ))%>%
  filter(cohort %in% c("Unsolved", "Monogenic", "APOL1")) # Filter to include only the specified cohorts

# Calculate means and confidence intervals
df_summary <- df_long %>%
  group_by(PRS_score, cohort) %>%
  summarise(
    mean = mean(PRS_value, na.rm = TRUE),
    lower_ci = mean - qt(0.975, df=n()-1) * sd(PRS_value) / sqrt(n()),
    upper_ci = mean + qt(0.975, df=n()-1) * sd(PRS_value) / sqrt(n()),
    .groups = 'drop'
  )

# Plotting
ggplot(df_summary, aes(x = PRS_score, y = mean, color = cohort)) +
  geom_point(position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, position = position_dodge(width = 0.25)) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Average PRS with 95% CI by Cohort", x = "PRS Type", y = "Standardised Average PRS") +
  theme_minimal() +
  scale_x_discrete(labels = c("iga_prs_standardised" = "IgA PRS", 
                              "ssns_prs_standardised" = "SSNS PRS",
                             "membranous_prs_standardised" = "Membranous PRS",
                             "ckd_apol_adjusted" = "CKD PRS")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"), # Keep axis lines, make them black
        axis.line.x = element_line(color = "black"), # Explicitly set for x-axis if needed
        axis.line.y = element_line(color = "black"))+ # Explicitly set for y-axis if needed
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")


#with ssns
prs2 <- prs
prs2$ssns <- prs2$ssns*-1
df_long <- prs2 %>%
  pivot_longer(cols = c("iga_prs_standardised", "ssns_prs_standardised",
                        "membranous_prs_standardised", "ckd_apol_adjusted"), 
               names_to = "PRS_score", values_to = "PRS_value") %>%
  mutate(cohort = case_when(
    pheno == 1 & apol1 == 1 ~ "APOL1 (n=16)", # APOL1_Positive takes priority
    pheno == 1 & monogenic == 1 ~ "Monogenic (n=38)",
    pheno == 1 & monogenic == 0 ~ "Unsolved (n=180)",    TRUE ~ as.character(NA)
  ))%>%
  filter(cohort %in% c("Unsolved (n=180)", "Monogenic (n=38)", "APOL1 (n=16)")) # Filter to include only the specified cohorts

# Calculate means and confidence intervals
df_summary <- df_long %>%
  group_by(PRS_score, cohort) %>%
  summarise(
    mean = mean(PRS_value, na.rm = TRUE),
    lower_ci = mean - qt(0.975, df=n()-1) * sd(PRS_value) / sqrt(n()),
    upper_ci = mean + qt(0.975, df=n()-1) * sd(PRS_value) / sqrt(n()),
    .groups = 'drop'
  )
  df_summary$PRS_score <- factor(df_summary$PRS_score, levels = c("ssns_prs_standardised",
                                                          "membranous_prs_standardised",
                                                          "ckd_apol_adjusted","iga_prs_standardised"))

# Plotting
  
ggplot(df_summary, aes(x = PRS_score, y = mean, color = cohort)) +
  geom_point(position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, position = position_dodge(width = 0.25)) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Average PRS with 95% CI by Cohort", x = "PRS Type", y = "Standardised Average PRS") +
  theme_minimal() +
  scale_x_discrete(labels = c("iga_prs_standardised" = "IgA PRS", 
                              "ssns_prs_standardised" = "SSNS PRS",
                              "membranous_prs_standardised" = "Membranous PRS",
                              "ckd_apol_adjusted" = "CKD PRS")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"), # Keep axis lines, make them black
        axis.line.x = element_line(color = "black"), # Explicitly set for x-axis if needed
        axis.line.y = element_line(color = "black"))+ # Explicitly set for y-axis if needed
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

#stats between case cohorts
library(dplyr)
library(stats)
library(multcomp)

# Subset data for one PRS type, for example, Iga
df_iga <- df_long %>% filter(PRS_score == "iga_prs_standardised")
df_ssns <- df_long %>% filter(PRS_score == "ssns_prs_standardised")
df_ckd <- df_long %>% filter(PRS_score == "ckd_apol_adjusted")
df_membranous <- df_long %>% filter(PRS_score == "membranous_prs_standardised")

# ANOVA
anova_result <- aov(PRS_value ~ cohort, data = df_iga)
summary(anova_result)
anova_result <- aov(PRS_value ~ cohort, data = df_ssns)#significant for APOL1
summary(anova_result)
anova_result <- aov(PRS_value ~ cohort, data = df_ckd)
summary(anova_result)
anova_result <- aov(PRS_value ~ cohort, data = df_membranous)
summary(anova_result)

# If ANOVA is significant, perform Tukey's HSD test
if (summary(anova_result)[[1]]$'Pr(>F)'[1] < 0.05) {
  tukey_result <- TukeyHSD(anova_result)
  print(tukey_result)
}

anova_result <- aov(PRS_value ~ cohort, data = df_iga)
summary(anova_result)

# If ANOVA is significant, perform Tukey's HSD test
if (summary(anova_result)[[1]]$'Pr(>F)'[1] < 0.05) {
  tukey_result <- TukeyHSD(anova_result)
  print(tukey_result)}
  
# Repeat the process for each PRS type
  
###african mean prs comparison
  # Calculate means and confidence intervals
  df_summary <- df_long %>%
    group_by(PRS_score, cohort) %>%
    summarise(
      mean = mean(PRS_value, na.rm = TRUE),
      lower_ci = mean - qt(0.975, df=n()-1) * sd(PRS_value) / sqrt(n()),
      upper_ci = mean + qt(0.975, df=n()-1) * sd(PRS_value) / sqrt(n()),
      .groups = 'drop'
    )
  df_summary$PRS_score <- factor(df_summary$PRS_score, levels = c("ssns_prs_standardised",
                                                                  "membranous_prs_standardised",
                                                                  "ckd_apol_adjusted","iga_prs_standardised"))
 #plot without non ssns scores
  df_summary <- subset(df_summary, df_summary$PRS_score=="ssns_prs_standardised")
  # Plotting
  
  ggplot(df_summary, aes(x = PRS_score, y = mean, color = cohort)) +
    geom_point(position = position_dodge(width = 0.25)) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, position = position_dodge(width = 0.25)) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Average PRS with 95% CI by Cohort", x = "PRS Type", y = "Standardised Average PRS") +
    theme_minimal() +
    scale_x_discrete(labels = c("iga_prs_standardised" = "IgA PRS", 
                                "ssns_prs_standardised" = "SSNS PRS",
                                "membranous_prs_standardised" = "Membranous PRS",
                                "ckd_apol_adjusted" = "CKD PRS")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"), # Keep axis lines, make them black
          axis.line.x = element_line(color = "black"), # Explicitly set for x-axis if needed
          axis.line.y = element_line(color = "black"))+ # Explicitly set for y-axis if needed
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  
  #stats between case cohorts
  library(dplyr)
  library(stats)
  library(multcomp)
  
  # Subset data for one PRS type, for example, Iga
  df_iga <- df_long %>% filter(PRS_score == "iga_prs_standardised")
  df_ssns <- df_long %>% filter(PRS_score == "ssns_prs_standardised")
  df_ckd <- df_long %>% filter(PRS_score == "ckd_apol_adjusted")
  df_membranous <- df_long %>% filter(PRS_score == "membranous_prs_standardised")
  
  # ANOVA
  anova_result <- aov(PRS_value ~ cohort, data = df_iga)
  summary(anova_result)
  anova_result <- aov(PRS_value ~ cohort, data = df_ssns)#significant for APOL1
  summary(anova_result)
  anova_result <- aov(PRS_value ~ cohort, data = df_ckd)
  summary(anova_result)
  anova_result <- aov(PRS_value ~ cohort, data = df_membranous)
  summary(anova_result)
  
  # If ANOVA is significant, perform Tukey's HSD test
  if (summary(anova_result)[[1]]$'Pr(>F)'[1] < 0.05) {
    tukey_result <- TukeyHSD(anova_result)
    print(tukey_result)
  }
  
  anova_result <- aov(PRS_value ~ cohort, data = df_iga)
  summary(anova_result)
  
  # If ANOVA is significant, perform Tukey's HSD test
  if (summary(anova_result)[[1]]$'Pr(>F)'[1] < 0.05) {
    tukey_result <- TukeyHSD(anova_result)
    print(tukey_result)
  }
  
##analysis per subcohort adjusted and tested against PRS
  prs_analysis <- prs
  prs_analysis <-  prs_analysis %>%
    mutate(cohort = case_when(
    pheno == 1 & apol1 == 1 ~ "APOL1", # APOL1_Positive takes priority
    pheno == 0   ~ "Controls",
    pheno == 1 & apol1 == 0 & monogenic == 1 ~ "Monogenic", 
    pheno == 1 & apol1 == 0 & monogenic == 0 ~ "Unsolved", 
    pheno == 1 & apol1 == 1 & monogenic == 1 ~ "APOL1",    TRUE ~ as.character(NA)
  ))
    prs_analysis$cohort <- as.factor(prs_analysis$cohort)
    prs_analysis$sex <- as.factor(prs_analysis$sex)
    
#model 
model <- lm(ckd_apol_adjusted ~ cohort + sex + age + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+cohort:ckd_apol_adjusted, data = prs_analysis)
    summary(model)
    par(mfrow=c(2,2))
    plot(model)
    
#increase in various things and OR
#calculate SD increase in apol1
sd_increase_in_prs <- exp(model$coefficients["apol1"]*sd(prs$apol1))
sd_increase_in_prs
#95% CI
five_percent <- exp(confint(model)["apol1",1]*sd(prs$apol1))
ninety_five_percent <- exp(confint(model)["apol1",2]*sd(prs$apol1))

#plot difference between controls and cases by apol1 and other PRS
to_plot <- prs

#plot membranous 
ggplot(to_plot, aes(x = membranous_prs_standardised, fill = factor(pheno))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of PRS in cases and controls",
       x = "PRS",
       y = NULL) +
  theme_minimal() + 
  labs(fill = "pheno")

#plot ckd versus apol1 type 
ggplot(to_plot, aes(x = ckd_prs_standardised, fill = factor(type))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of PRS in cases and controls",
       x = "PRS",
       y = NULL) +
  theme_minimal() + 
  labs(fill = "type")

#plot ssns and iga
ggplot(to_plot, aes(x = membranous_prs_standardised, fill = factor(type))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of PRS in cases and controls",
       x = "PRS",
       y = NULL) +
  theme_minimal() + 
  labs(fill = "type")

#plot box and violin plot 
#comparing two distributions - ckd
kruskal.test(ckd_prs_standardised ~ type, data=prs)
pairwise.wilcox.test(prs$ckd_prs_standardised, prs$type, p.adjust.method = "BH")
kruskal.test(apol1 ~ pheno, data=prs)
#iga
kruskal.test(iga_prs_standardised ~ type, data=prs)
pairwise.wilcox.test(prs$iga_prs_standardised, prs$type, p.adjust.method = "BH")
kruskal.test(iga_prs_standardised ~ pheno, data=prs)
#ssns
kruskal.test(ssns_prs_standardised ~ type, data=prs)
pairwise.wilcox.test(prs$ssns_prs_standardised, prs$type, p.adjust.method = "BH")
kruskal.test(ssns_prs_standardised ~ pheno, data=prs)
#membranous
kruskal.test(membranous_prs_standardised ~ type, data=prs)
pairwise.wilcox.test(prs$membranous_prs, prs$type, p.adjust.method = "BH")
kruskal.test(membranous_prs_standardised ~ pheno, data=prs)
#pop <- aov(ldpred_beta_SUM ~ PHENO, data=subgroup.pheno.prs)
#summary(pop)
#t.test(x$ldpred_beta_SUM, y$ldpred_beta_SUM)

#plot with stats 
my_comparisons <- list(c("Case - High Rx APOL1", "Case - Non Rx APOL1 Control",
                         "Control - High Rx APOL1", "Control - Non Rx APOL1"))
ggplot(to_plot, aes(x=type, y=ckd_prs_standardised, fill=type)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  scale_fill_brewer(palette = "Dark2") +
  ylim(-4,5) + 
  theme_classic() +
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") +
  theme(axis.text.x = element_text(angle=30, vjust=0.6, hjust=0.5)) +
  stat_compare_means(comparisons = my_comparisons) + 
  #scale_x_discrete(labels = c("Control", "Solved Case")) + 
  stat_compare_means(label.y =50)


ggviolin(to_plot, x = "type", y = "ckd_prs_standardised", fill = "pheno",
         palette = c("#00AFBB", "#f2241f"),
         add = "boxplot", add.params = list(fill = "white"))+
  #facet_wrap(~type)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 8) + 
  stat_cor(p.accuracy = 0.01) +
  #scale_x_discrete(labels = c("Control", "EEHTN", "Primary HTN")) + 
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") 
#ggsave("slc34a3_subgroup_prs_stones_ggpubr_4groups.violin.png", height = 7, width = 10)
ggviolin(to_plot, x = "pheno", y = "ckd_prs_standardised", fill = "pheno",
         palette = c("#00AFBB", "#f2241f"),
         add = "boxplot", add.params = list(fill = "white"))+
  #facet_wrap(~type)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 8) + 
  stat_cor(p.accuracy = 0.01) +
  #scale_x_discrete(labels = c("Control", "EEHTN", "Primary HTN")) + 
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") 


#compare distibutions 

prs$prs_standardised <- as.numeric(prs$prs_standardised)
comparison_results <- prs %>%
  group_by(type) %>%
  kruskal.test(prs_standardised)
comparison_results 
kruskal.test(prs_standardised ~ pheno, data=pheno.prs)

# Create faceted density plots of PRS by type and pheno_esrf using ggplot2
ggplot(prs, aes(x = prs_standardised, fill = factor(type))) +
  geom_density(alpha = 0.5) +
  labs(title = "Faceted Density Plot of PRS by Type and pheno_esrf",
       x = "PRS",
       y = "Density") +
  scale_fill_discrete(name = "pheno_esrf") +
  facet_grid(type ~ ., scales = "free_y") 

#compare
library(dplyr)
library(stats)

# Assuming your_data_frame is your data frame
# Perform t-test or Wilcoxon rank-sum test for each combination of type and pheno_esrf
comparison_results <- prs %>%
  group_by(factor(type)) %>%
  summarize(p_value = wilcox.test(prs_standardised ~ factor(type))$p.value)

print(comparison_results)

library(dplyr)
library(stats)

# Assuming your_data_frame is your data frame
# Perform Kruskal-Wallis test for each type
comparison_results <- prs %>%
  group_by(type) %>%
  summarise(p_value = kruskal.test(prs_standardised ~ factor(pheno))$p.value)

print(comparison_results)

## Assuming your_data_frame is your data frame
# Perform one-way ANOVA
anova_result <- aov(ckd_prs_standardised ~ type, data = prs)
summary(anova_result)
# Perform post-hoc tests (Tukey's HSD)
posthoc_result <- TukeyHSD(anova_result)

# View ANOVA table
print(summary(anova_result))

# View post-hoc test results
print(posthoc_result)



#AUC 
library(DescTools)
library(boot)
Cstat(model)
#bootstrap confidence intervals 
FUN <- function(x,i) {
  r.glm <- glm(pheno ~ ., data=x[i,], family=binomial)
  Cstat(r.glm)
}
boot.res <- boot(formodel.prs, FUN, R=999) 
boot.ci(boot.res, type="perc") #bca for narrower but need 10,000

# the percentile confidence intervals
boot.ci(boot.res, type="perc")


#calculate SD increase in PRS
sd_increase_in_prs <- exp(model$coefficients["prs_standardised"]*sd(prs$prs_standardised))
sd_increase_in_prs
#95% CI
five_percent <- exp(confint(model)["prs_standardised",1]*sd(prs$prs_standardised))
ninety_five_percent <- exp(confint(model)["prs_standardised",2]*sd(prs$prs_standardised))

#plot centile of model as prediction of getting USD
library(jtools)
library(OddsPlotty)
library(dplyr)
unsolved_plotting <- formodel.prs
unsolved_plotting$rank <- rank(unsolved_plotting$prs_standardised)
unsolved_plotting$rank <- round(unsolved_plotting$rank/length(unsolved_plotting$rank) * 100)
x <- unsolved_plotting %>% group_by(rank) %>% mutate(prevelance = (sum(pheno==1)/length(pheno)*100))
l <- x[,c(15,16)]
l <- unique(l)
l <- l[order(l$rank,decreasing = F),]
ggplot(l, aes(x=rank, y=prevelance)) + geom_smooth()
#fit model with rank 
unsolved_plotting_centile <- unsolved_plotting[,c(1,15,4:13)]
unsolved_centile <- glm(pheno~., family=binomial(link='logit'),data=unsolved_plotting_centile, maxit=10000)
effect_plot(data = unsolved_plotting_centile, unsolved_centile, pred=rank, interval = T, x.label = "PRS Centile",
            y.label = "Prevelance of USD")
#plot_summs(unsolved_model, model, coefs=coef_names)


#plot standardies graphs and distributions 
#create another pheno cohort to subdivide solved, slc34a3 and controls that have slc34a3 
subgroup.pheno.prs <- pheno.prs

subgroup.pheno.prs$PHENO <- as.factor(subgroup.pheno.prs$PHENO)
subgroup.pheno.prs$PHENO <- reorder(subgroup.pheno.prs$PHENO, subgroup.pheno.prs$prs_standardised, mean)

#stats on subgroup 
group_by(prs, type) %>%
  summarise(
    count = n(),
    mean = mean(prs_standardised, na.rm = TRUE),
    median = median(prs_standardised, na.rm = TRUE),
    sd = sd(prs_standardised, na.rm = TRUE)
  )



#plot with stats 
my_comparisons <- list(c("Case - High Rx APOL1", "Control - High Rx APOL1"), c("Case - Non Rx APOL1", "Control - Non Rx APOL1"), c("Case - High Rx APOL1","Case - Non Rx APOL1 Control"))

ggviolin(prs, x = "pheno", y = "membranous_prs_standardised", fill = "pheno",
         palette = c("#00AFBB", "#FC4E07","#f2241f", "black"),
         add = "boxplot", add.params = list(fill = "white"))+
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  #stat_compare_means(label.y = 8) + 
  #stat_cor(p.accuracy = 0.01) +
  labs(x="Cohort", y="Polygenic Risk Score", fill="Cohort") 
ggsave("slc34a3_subgroup_prs_stones_ggpubr_4groups.violin.png", height = 7, width = 10)
