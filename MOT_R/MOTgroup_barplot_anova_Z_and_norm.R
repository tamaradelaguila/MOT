rm(list = ls())
# Restart to clear also loaded packages: Ctrl + Shift + (Fn) +  F10 >>> RESTARTS

# SETTINGS
sourcefile ="group2_Zscore_long_format_forR_group_05filt2_circle_rej1_11roi.xls" # MIND THIS: the first blank sheets have to be deleted (so the first page is that one of interest)

sinkfile =  'Rgroup2_'  # to sink the data in 'Rgroup_R_'  'Rgroup_2roi_'
infotitle=  "05filt2"

analysis ='1way'
# analysis = '2way'

# measure = 'peakZ'
# measure = 'wmeanZ'
measure = 'slopeZ'

# measure = 'peak'
# measure= 'wmean'
# measure= 'slope'

# measure = 'peak_norm' #for old format
# measure = 'wmean_norm' #for old format

# END OF SETTINGS
# get Rejection mode from title
rej <- sub(".*rej", "", sourcefile)  # Extract characters after pattern
rej=substring(rej, 1,1)
tit = paste(infotitle ,measure, 'rej', rej)
sinkfile = paste(sinkfile, infotitle, measure, '_rej', rej)

# GROUP BARPLOT - MEASURES 
library("gdata")
setwd("/home/tamara/Documents/MATLAB/VSDI/MOT")
rawdata <-read.xls(sourcefile)

# subset
rawdata <-  rawdata  %>% filter (rawdata$cond == 2 | rawdata$cond == 3 | rawdata$cond == 4)

# Three-way mixed ANOVA: 2 between- and 1 within-subjects factors
# https://www.datanovia.com/en/lessons/mixed-anova-in-r/#three-way-bww-b

# LOAD PACKAGES
library(tidyverse) #for data manipulation and visualization
library(ggpubr) #for creating easily publication ready plots. contains ggqqplot
library(rstatix) #provides pipe-friendly R functions for easy statistical analyses. contains: shapiro_test and levene_test
library("scales")
library (emmeans)
library(car)

# Change 'tibble' options so all columns and rows are printed:
options(tibble.width = Inf)
options(tibble.print_min = Inf)

# LOAD AND RENAME DATA (wmean) - AND TURN INTO FACTORS
rawdata$id  <-as.factor(rawdata$id)
rawdata$cond  <-as.factor(rawdata$cond)
rawdata$roiname  <- as.factor(rawdata$roi_n)
rawdata$roi  <- as.factor(rawdata$roi)

if (measure == 'peak') {rawdata$value  <-rawdata$peak
} else if (measure == 'wmean'){
  rawdata$value  <-rawdata$wmean
  data  <- rawdata[c(1,2,3,11)] # preserve only cols of interest
  
} else if (measure == 'peakZ') {
    rawdata$value  <-rawdata$peak_z
    data  <- rawdata[c(1,2,3,11)] # preserve only cols of interest
    
} else if (measure == 'wmeanZ'){ rawdata$value  <-rawdata$wmean_z
data  <- rawdata[c(1,2,3,11)] # preserve only cols of interest

} else if (measure == 'slope'){ rawdata$value  <-rawdata$slope
data  <- rawdata[c(1,2,3,11)] # preserve only cols of interest 

} else if (measure == 'slopeZ'){ rawdata$value  <-rawdata$slope_z
data  <- rawdata[c(1,2,3,11)] # preserve only cols of  interest
}

#!!! subset if desired
# rawdata <-  rawdata  %>% filter (rawdata$cond == 1 | data$cond == 2 | data$cond == 3)
# data$roi  <-as.factor(data$roi)
# levels(data$roin) <- 1:11


# LOAD AND RENAME DATA (bradi)

# BARPLOT
library(ggplot2)
library(dplyr)
library(FSA)

# BARPLOT FOR THE MEASURE WMEAN

# First summarize the data 
# sum = Summarize(value ~ roi + cond,
#                 data=data)
sum = Summarize(value ~ roi.n,
                 data=data)
sum$se = sum$sd / sqrt(sum$n)

# Then plot

pd = position_dodge(0.9)    ### How much to jitter the bars on the plot,
###   0.5 for overlapping bars

barplot_mean <-  ggplot(sum,                ### The data frame to use.
                        aes(x= roi,
                            y= mean)) +
  
  geom_bar(stat= "identity",
           color= "black",
           position= pd) +
  
  geom_errorbar(aes(ymin  = mean - sd,
                    ymax  = mean + sd),
                width = 0.1,
                size  = 0.5,
                position = pd,
                color = "black"
  ) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold")) +
  theme(legend.position = "bottom")+
  ylab("peak +- sd") +
ggtitle(tit)
barplot_mean


# BOXPLOT
# bxp <- ggboxplot(
#   datapeak, x = "roi", y = "peak",
#   color = "cond", 
#   ggtheme = ggthemes::theme_fivethirtyeight(), # choosing a different theme
# ) +
#   ggtitle('by roi')
# 
# bxp

# PANEL
library(gridExtra)
grid.arrange(barplot_bradi, barplot_mean, ncol = 1)
barplot_bradi
barplot_mean

#########################################################################
## 2-WAY ANOVA: FORMATED AND SINKING INTO OUTPUT.txt #### 
#########################################################################
library(ggthemes) #https://rafalab.github.io/dsbook/ggplot2.html#add-on-packages
library(wesanderson)


if (analysis == '2way_rm') { file = paste('ANOVA_', sinkfile , '_2waymixed.txt')} else if (analysis == '2way') {file = paste(sinkfile , '_ANOVA_2way.txt') }  else if (analysis == '1way') {file = paste(sinkfile , '_ANOVA_1way.txt') }
  
# 1. SUMMARY STATISTICS ####

# sink(file)

print("SUMMARY STATISTICS")
print("")
print("Summary grouped by factor ' roi'")
data %>%
  group_by(roi.n) %>%
  get_summary_stats(value, type = "mean_sd") #value cannot be a factor


print("")
print("Summary grouped by both factors")
data%>%
  group_by(roi.n) %>%
  get_summary_stats(value, show = c("mean", "sd", "se", "median"))


# sink()

# 2. BOXPLOT
# change and save manually:
cmap = wes_palette(name= "Darjeeling1")
cmap[6] = wes_palette(name= "Royal1")[1]
cmap[7:12] = cmap[1:6]

bxp <- ggboxplot(
  data, x = "roi.n", y = "value",color="roi.n",
  ggtheme = ggthemes::theme_fivethirtyeight() # choosing a different theme
  ,palette = cmap,
) +
  ggtitle(tit)

bxp

jpeg(paste('Rboxplot', sinkfile, 'byroi' ,'.jpg'), width = 500, height = 400)
bxp
dev.off()


#graphics.off()

# 3 . ANOVA: 2-way anova, not repeated measures (when 'Error' term unmuted, then it's repeated measures for that variable)

if (analysis == '2way') { 
# #########################################################################
# # OPTION 1 - NOT REPEATED MEASURES
print("ANOVA: ANOVA SUMMARY")
anova_summary(aov(value ~roi*cond, data=data))
res.aov <- anova_test(data=data, formula = value ~roi*cond,
                      effect.size='pes', type=3)
} else if (analysis == '2way_rm') {
#########################################################################
# OPTION 2 - REPEATED MEASURES
#if 'Error' term is included, it's a mixed effect (r.m.)
#error function: we want to control for that between-participant variation over all of our within-subjects variables.
  
anova_summary(aov(value ~roi*cond + Error(idx/(cond)), data=data))
res.aov <- anova_test(data=data, formula = value ~roi*cond  + Error(id/(cond)),
                      effect.size='pes', type=3)
} else if (analysis == '1way') {
# #########################################################################
# # OPTION 1 - NOT REPEATED MEASURES
print("ANOVA: ANOVA SUMMARY")
anova_summary(aov(value ~roi, data=data))
res.aov <- anova_test(data=data, formula = value ~roi,
                      effect.size='pes', type=3)
#########################################################################
}

sink(file)
print(paste('source file:', sourcefile))
print(paste('analysis= ' ,analysis))
print(paste('measure= ' ,measure))

print(tit)
print("ANOVA: ANOVA TABLE")
#get_anova_table(res.aov)
knitr::kable(get_anova_table(res.aov), caption = 'ANOVA Table (type III tests) - lo mismo, mas bonito')

print("Complete ANOVA object to check Mauchly and Sphericity corrections") --???
  res.aov  

# 4. POST-HOC ANALYSIS

print("POST-HOC with Bonferroni adjustment")

# 4.1: PAIRWISE COMPARISONS each level of the within-s factors

# 4.1.1. Each level of block
# Paired t-test is used because we have repeated measures by time


print("POST-HOC: PAIRWISE COMPARISONS BETWEEN ROI LEVELS  (non paired)")
pwc1 <- data %>%
  pairwise_t_test(value ~roi.n,p.adjust.method = "bonferroni" )
#  %>% select(-df, -statistic, -p) # Remove details
#pwc1

knitr::kable(pwc1, caption = 'POST-HOC: PAIRWISE COMPARISONS BETWEEN ROI LEVELS')


print("POST-HOC: PAIRWISE COMPARISONS BETWEEN mA LEVELS (non paired)")
# non-pair to compare EC? because is not a rm?
pwc2 <- data %>%
  group_by(cond) %>%
  pairwise_t_test(
    formula = value ~ roi, paired = FALSE, 
    p.adjust.method = "bonferroni"
  )
#%>%   select(-df, -statistic, -p) # Remove details
#pwc2

knitr::kable(pwc2, caption = 'POST-HOC: PAIRWISE COMPARISONS BETWEEN cond LEVELS')

sink ()

# ...from here on it has to be pasted manually (otherwise it'll break the format)
#########################################################################
#########################################################################

# 5. TEST ASSUMPTIONS  ####

# 5.1: PRESENCE OF OUTLIERS 

print("ASSUMPTIONS: PRESENCE OF OUTLIERS")
data %>%
  group_by(roi, cond) %>%
  identify_outliers(value)

data %>%
  group_by(cond) %>%
  identify_outliers(value)


# 5.2: NORMALITY ASSUMPTION
# The norcondlity assumption can be checked by computing Shapiro-Wilk test for each time point.
# p > 0.05: normal distribution (H0)

print("ASSUMPTIONS: NORMALITY - SHAPIRO-WILK TEST")
# Shapiro test needs at least three points to run
data %>%
  group_by(roi, cond) %>%
  shapiro_test(value)

# CREATE QQ plot for each cell of design: (for bigger sample sizes)
# ggqqplot(data, "brady", ggtheme = theme_bw()) +
# facet_grid(roi ~ cond, labeller = "label_both")
# If all the points fall approx. along the ref.line, for each cell, we can assume normality of the data.

# 5.3. HOMOCEDASTICITY = homogenity of variance: not needed (not between subjects comparisons)


# 5.4: SPERICITY ASSUMPTION: done and corrected internally by the functions anova_test() and get_anova_table()
# Sphericity is an important assumption of a repeated-measures ANOVA 

print("ASSUMPTIONS: SPERICITY ASSUMPTION is performed and applied internally")

# 5.5 : HOMOGENEITY OF  COVARIANCES  !!! sale diferente
print("ASSUMPTIONS: EQUALITY OF COVARIANCEs - Box's Text")
box_m(data[, "value", drop = FALSE], data$roi)
# The assumption of homogeneity of variance is an assumption of the independent samples t-test and ANOVA stating
#that all comparison groups have the same variance.  The independent samples t-test and ANOVA utilize the t and F statistics respectively,
# which are generally robust to violations of the assumption as long as group sizes are equal

#######################################
library(reshape2)
p.val.test<-pwpm(emmeans(value ~roiname, "roi.n"),means = FALSE, flip = TRUE,reverse = TRUE) # p-values presented compactly in matrix form
p.val.test<-sub("[<>]", "", p.val.test)
p.matx<-matrix(as.numeric((p.val.test)),nrow = length(p.val.test[,1]),ncol = length(p.val.test[,1])) #if your factor has 5 levels ncol and nrow=5
rownames(p.matx) <- colnames(p.matx) <-colnames(p.val.test)
p.matx[upper.tri(p.matx, diag=FALSE)] <- NA
melt(p.matx) %>%
  ggplot(aes(Var1, Var2, fill = value)) + geom_tile() +
  geom_text(aes(label = value))

p.matx
# 
# 
# #########################################################################
# ## ROBUST 2-WAY ANOVA: WRS2
# #########################################################################
# library(WRS2)
# # Median-based (more robust to outliers)
# set.seed(123)
# #check unused factors levels
# b1= nlevels(data$cond)
# b2= length(unique(data$roi)) # if b1 and b2 are different, drop unused levels:
# 
# # we drop unused levels (otherwise med2way will raise an error)
# data$cond <- droplevels(data$cond)
# data$roi <- droplevels(data$roi)
# 
# med2way(value~roi*cond, data = data)
# 
# postHoc1 <- mcp2a(value~roi*cond, data = data,)
# postHoc1
# 
# # see levels to be able to understand the posthoc
# levels(data$roi)
# levels(data$cond)
# 
# # DUNN's TEST (PAIWISE COMPARISONS) for each level of cond
# data  %>%
#   group_by(cond) %>%
#   dunn_test(value ~ roi)
# 
# data  %>%
#   group_by(roi) %>%
#   dunn_test(value ~ cond)