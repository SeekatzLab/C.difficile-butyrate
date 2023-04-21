## Paper: SCFA and C. difficile
## Figure 3: toxin production in presence of scfa
### 3a: toxin production CD630 in different scfa
### Supplementary figure S3: toxin production in VPI10463, R20291; and toxin production over time 630
### 3b: toxin production in presence of butyrate concentration gradient
### 3d: tcdC, tcdR qPCR results


rm(list = ls())
setwd("/Data/figure3_S3/")

## Library

library(readxl)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(nlme)
library(stringi)
library(stringr)
library(plyr)
library(dplyr)
library(growthcurver)
library(DescTools)

# Figure 3a: effect of butyrate on toxin production

tox_630 <- read_tsv(file = "./ToxinData_ForGraph630.txt")

tox_630_2 <- subset(tox_630, !(tox_630$group == "But10" | tox_630$group == "But50"))
tox_630_2$group <- factor(tox_630_2$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))

ggplot(tox_630_2, aes(x=group, y=Toxin, color=group)) + geom_point(size =0.5, position = position_jitter(width = 0.2)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  ylim(0,6) +
  #scale_y_continuous(limits= c(0, 10), breaks = seq(0,10,2))+
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c"))

##STATS

aov_tox_630 <- aov(Toxin~group, data=tox_630_2)
summary(aov_tox_630)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#group        6  43.17   7.195   11.93 7.71e-10 ***
#  Residuals   91  54.88   0.603                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_tox_630)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = Toxin ~ group, data = tox_630_2)

#$group
#diff         lwr       upr     p adj
#Ace5-BHI     1.277778e+00  0.40508924 2.1504663 0.0005380 ***
#Ace25-BHI    1.277778e+00  0.40508924 2.1504663 0.0005380 ***
#Pro5-BHI     1.444444e+00  0.57175591 2.3171330 0.0000578 ***
#Pro25-BHI    1.611111e+00  0.73842257 2.4837996 0.0000053 ***
#But5-BHI     1.444444e+00  0.63986537 2.2490235 0.0000102 ***
#But25-BHI    2.131944e+00  1.32736537 2.9365235 0.0000000 ***
#Ace25-Ace5   4.440892e-15 -0.95598240 0.9559824 1.0000000
#Pro5-Ace5    1.666667e-01 -0.78931573 1.1226491 0.9984065
#Pro25-Ace5   3.333333e-01 -0.62264906 1.2893157 0.9402679
#But5-Ace5    1.666667e-01 -0.72757298 1.0609063 0.9976833
#But25-Ace5   8.541667e-01 -0.04007298 1.7484063 0.0709906
#Pro5-Ace25   1.666667e-01 -0.78931573 1.1226491 0.9984065
#Pro25-Ace25  3.333333e-01 -0.62264906 1.2893157 0.9402679
#But5-Ace25   1.666667e-01 -0.72757298 1.0609063 0.9976833
#But25-Ace25  8.541667e-01 -0.04007298 1.7484063 0.0709906
#Pro25-Pro5   1.666667e-01 -0.78931573 1.1226491 0.9984065
#But5-Pro5    0.000000e+00 -0.89423965 0.8942396 1.0000000
#But25-Pro5   6.875000e-01 -0.20673965 1.5817396 0.2468137
#But5-Pro25  -1.666667e-01 -1.06090632 0.7275730 0.9976833
#But25-Pro25  5.208333e-01 -0.37340632 1.4150730 0.5806217
#But25-But5   6.875000e-01 -0.14040504 1.5154050 0.1701835

DunnettTest(x=tox_630_2$Toxin, g=tox_630_2$group, "BHI")


#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#              diff    lwr.ci   upr.ci    pval    
#Ace5-BHI  1.277778 0.5141749 2.041381 0.00016 ***
#Ace25-BHI 1.277778 0.5141749 2.041381 0.00015 ***
#Pro5-BHI  1.444444 0.6808416 2.208047 1.5e-05 ***
#Pro25-BHI 1.611111 0.8475082 2.374714 1.7e-06 ***
#But5-BHI  1.444444 0.7404374 2.148451 2.6e-06 ***
#But25-BHI 2.131944 1.4279374 2.835951 8.8e-12 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


## Figure 3b: butyrate concentration gradient on C.diff tox
tox_but_630 <- read_xlsx(path = "./toxin_Julian_multibut_time.xlsx")
tox_but_630$group <- factor(tox_but_630$group, levels = c("BHI","But5", "But10", "But25", "But50"))


ggplot(tox_but_630, aes(x=as.factor(time), y=Toxin, color=group)) + geom_point(size =0.5, position = position_jitterdodge(jitter.width = 0.5)) +
  #geom_line() + 
  #facet_grid(~time) +
  geom_boxplot(alpha=0.4, size=0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  ylim(0,6.5)+
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black", hcl.colors(4, palette = "TealGrn")))+
  scale_fill_manual(values=c("black",hcl.colors(4, palette = "TealGrn")))


## STATS

tox_but_630_6 <- subset(tox_but_630, tox_but_630$time == 6)
aov_tox_but_630_6 <- aov(Toxin~group, data = tox_but_630_6)
summary(aov_tox_but_630_6)
TukeyHSD(aov_tox_but_630_6)

## no significance to BHI

DunnettTest(x = tox_but_630_6$Toxin, g = tox_but_630_6$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#             diff     lwr.ci     upr.ci   pval    
#But5-BHI   0.1250 -0.2296574 0.47965742 0.7920    
#But10-BHI -0.1250 -0.4796574 0.22965742 0.7920    
#But25-BHI -0.1875 -0.5421574 0.16715742 0.4890    
#But50-BHI -0.3125 -0.6671574 0.04215742 0.0998 .  

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


tox_but_630_12 <- subset(tox_but_630, tox_but_630$time == 12)
aov_tox_but_630_12 <- aov(Toxin~group, data = tox_but_630_12)
summary(aov_tox_but_630_12)
TukeyHSD(aov_tox_but_630_12)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = Toxin ~ group, data = tox_but_630_12)

#$group
#diff         lwr         upr     p adj
#But5-BHI     0.5625  0.11232153  1.01267847 0.0070035
#But10-BHI    0.0625 -0.38767847  0.51267847 0.9950983
#But25-BHI    0.3125 -0.13767847  0.76267847 0.3054583
#But50-BHI    0.3750 -0.07517847  0.82517847 0.1473963

DunnettTest(x=tox_but_630_12$Toxin, g=tox_but_630_12$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#            diff      lwr.ci    upr.ci   pval    
#But5-BHI  0.5625  0.16070192 0.9642981 0.0030 ** 
#But10-BHI 0.0625 -0.33929808 0.4642981 0.9857    
#But25-BHI 0.3125 -0.08929808 0.7142981 0.1713    
#But50-BHI 0.3750 -0.02679808 0.7767981 0.0741 .  

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

tox_but_630_18 <- subset(tox_but_630, tox_but_630$time == 18)
aov_tox_but_630_18 <- aov(Toxin~group, data = tox_but_630_18)
summary(aov_tox_but_630_18)
TukeyHSD(aov_tox_but_630_18)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = Toxin ~ group, data = tox_but_630_18)

#$group
#diff        lwr       upr     p adj
#But5-BHI    -2.500000e-01 -0.7019828 0.2019828 0.5362737
#But10-BHI   -2.500000e-01 -0.7019828 0.2019828 0.5362737
#But25-BHI   -1.250000e-01 -0.5769828 0.3269828 0.9375805
#But50-BHI    3.750000e-01 -0.0769828 0.8269828 0.1502759
#But10-But5   2.664535e-15 -0.4519828 0.4519828 1.0000000
#But25-But5   1.250000e-01 -0.3269828 0.5769828 0.9375805
#But50-But5   6.250000e-01  0.1730172 1.0769828 0.0021247
#But25-But10  1.250000e-01 -0.3269828 0.5769828 0.9375805
#But50-But10  6.250000e-01  0.1730172 1.0769828 0.0021247
#But50-But25  5.000000e-01  0.0480172 0.9519828 0.0226881

DunnettTest(x=tox_but_630_18$Toxin, g=tox_but_630_18$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#            diff     lwr.ci    upr.ci   pval    
#But5-BHI  -0.250 -0.6534085 0.1534085 0.3483    
#But10-BHI -0.250 -0.6534085 0.1534085 0.3484    
#But25-BHI -0.125 -0.5284085 0.2784085 0.8548    
#But50-BHI  0.375 -0.0284085 0.7784085 0.0759 .  

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

tox_but_630_24 <- subset(tox_but_630, tox_but_630$time == 24)
aov_tox_but_630_24 <- aov(Toxin~group, data = tox_but_630_24)
summary(aov_tox_but_630_24)
TukeyHSD(aov_tox_but_630_24)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = Toxin ~ group, data = tox_but_630_24)

#$group
#diff         lwr       upr     p adj
#But5-BHI     0.500  0.06546072 0.9345393 0.0159635
#But10-BHI    0.250 -0.18453928 0.6845393 0.4969644
#But25-BHI    0.625  0.19046072 1.0595393 0.0012610
#But50-BHI    1.125  0.69046072 1.5595393 0.0000000
#But10-But5  -0.250 -0.68453928 0.1845393 0.4969644
#But25-But5   0.125 -0.30953928 0.5595393 0.9285788
#But50-But5   0.625  0.19046072 1.0595393 0.0012610
#But25-But10  0.375 -0.05953928 0.8095393 0.1233082
#But50-But10  0.875  0.44046072 1.3095393 0.0000029
#But50-But25  0.500  0.06546072 0.9345393 0.0159635

DunnettTest(x=tox_but_630_24$Toxin, g=tox_but_630_24$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#           diff     lwr.ci    upr.ci    pval    
#But5-BHI  0.500  0.1121604 0.8878396 0.00712 ** 
#But10-BHI 0.250 -0.1378396 0.6378396 0.31466    
#But25-BHI 0.625  0.2371604 1.0128396 0.00051 ***
#But50-BHI 1.125  0.7371604 1.5128396 4.1e-10 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


## Figure 3C: qPCR of tcdC and tcdR

qpcr <- read_tsv(file = "./qPCR_full.txt")
qpcr <- subset(qpcr, !is.na(qpcr$time))

qpcr.bhi.log <- qpcr %>% 
  group_by(gene, time) %>% 
  dplyr::summarise(mean = mean(logfc), n = length(logfc), sd = sd(logfc), se = sd / sqrt(n))

qpcr.bhi.log.tox <- subset(qpcr.bhi.log, qpcr.bhi.log$gene == "tcdC" | qpcr.bhi.log$gene == "tcdR")

ggplot(qpcr.bhi.log.tox, aes(x=gene, y=mean, fill=time)) + geom_col(position = position_dodge()) +
  theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("orangered3", "dodgerblue")) +
  ylab("mean logfc") +
  theme(axis.title.x=element_text(size=12, face='bold'),
        axis.title.y=element_text(size=12, face='bold'),
        panel.background = element_rect(fill = 'gray99'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #+
#facet_wrap(~time)

## STATS

#### nothing was found significant
tcdC_qpcr.logfc <- subset(qpcr, qpcr$gene == "tcdC")
aov_tox_qpcr.logfc <- aov(logfc~time, data = tcdC_qpcr.logfc)
summary(aov_tox_qpcr.logfc)
DunnettTest(x=tcdC_qpcr.logfc$logfc, g=tcdC_qpcr.logfc$time)
kruskal.test(logfc~time, data = tcdC_qpcr.logfc)
wilcox.test(logfc~time, data = tcdC_qpcr.logfc)

tcdR_qpcr.logfc <- subset(qpcr, qpcr$gene == "tcdR")
aov_tcdR_qpcr.logfc <- aov(logfc~time, data = tcdR_qpcr.logfc)
summary(aov_tcdR_qpcr.logfc)
kruskal.test(logfc~time, data = tcdR_qpcr.logfc)
wilcox.test(logfc~time, data = tcdR_qpcr.logfc)

qpcr.bhi <- qpcr %>% 
  group_by(gene, time) %>% 
  dplyr::summarise(mean = mean(BHI_delCt), n = length(BHI_delCt), sd = sd(BHI_delCt), se = sd / sqrt(n))

qpcr.bhi$type <- "BHI"

qpcr.but <- qpcr %>% 
  group_by(gene, time) %>% 
  dplyr::summarise(mean = mean(Butyrate_delCt), n = length(Butyrate_delCt), sd = sd(Butyrate_delCt), se = sd / sqrt(n))

qpcr.but$type <- "Butyrate"

qpcr.avg <- rbind(qpcr.bhi, qpcr.but)
qpcr.avg.tcdC <- subset(qpcr.avg, qpcr.avg$gene == "tcdC")
qpcr.avg.tcdC.aov <- aov(mean ~ type, data = qpcr.avg.tcdC)
summary(qpcr.avg.tcdC.aov)
qpcr.early.tcdC <- subset(qpcr, qpcr$time== "early" & qpcr$gene == "tcdC")
kruskal.test(BHI_delCt ~ Butyrate_delCt, data = qpcr.early.tcdC)

qpcr.late.tcdC <- subset(qpcr, qpcr$time== "early" & qpcr$gene == "tcdC")
kruskal.test(BHI_delCt ~ Butyrate_delCt, data = qpcr.early.tcdC)

qpcr.avg.tcdR <- subset(qpcr.avg, qpcr.avg$gene == "tcdR")
qpcr.avg.tcdR.aov <- aov(mean ~ type, data = qpcr.avg.tcdR)
summary(qpcr.avg.tcdR.aov)
qpcr.avg.early.tcdR <- subset(qpcr.avg.tcdR, qpcr.avg.tcdR$time== "early")
qpcr.early.tcdR <- subset(qpcr, qpcr$time== "early" & qpcr$gene == "tcdR")
kruskal.test(BHI_delCt ~ Butyrate_delCt, data = qpcr.early.tcdR)

qpcr.late.tcdR <- subset(qpcr, qpcr$time== "late" & qpcr$gene == "tcdR")
kruskal.test(BHI_delCt ~ Butyrate_delCt, data = qpcr.late.tcdR)


#######################Supplementary Figure S3##################################

## Figure S3A: R20291 toxin

tox_r2 <- read_tsv(file = "./ToxinData_ForGraphR20291.txt")

#tox_r2_2 <- subset(tox_630, !(tox_630$group == "But10" | tox_630$group == "But50"))
tox_r2$group <- factor(tox_r2$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))


ggplot(tox_r2, aes(x=group, y=Toxin, color=group)) + geom_point(size =0.5, position = position_jitter(width = 0.2)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  ylim(0,8) +
  #scale_y_continuous(limits= c(0, 10), breaks = seq(0,10,2))+
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c"))

##STATS
aov_tox_r2 <- aov(Toxin~group, data=tox_r2)
summary(aov_tox_r2)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#group        6  28.90   4.816   7.475 1.89e-05 ***
#  Residuals   41  26.42   0.644                    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=tox_r2$Toxin, g=tox_r2$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#              diff    lwr.ci   upr.ci    pval    
#Ace5-BHI  2.666667 1.3114326 4.021901 3.4e-05 ***
#Ace25-BHI 2.000000 0.7143121 3.285688 0.00093 ***
#Pro5-BHI  3.166667 1.8114326 4.521901 5.6e-07 ***
#Pro25-BHI 2.625000 1.3393121 3.910688 2.3e-05 ***
#But5-BHI  2.625000 1.3393121 3.910688 1.4e-05 ***
#But25-BHI 2.250000 0.9643121 3.535688 0.00021 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


TukeyHSD(aov_tox_r2)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = Toxin ~ group, data = tox_r2)

#$group
#diff        lwr       upr     p adj
#Ace5-BHI     2.66666667  1.0608503 4.2724830 0.0001340
#Ace25-BHI    2.00000000  0.4765888 3.5234112 0.0036585
#Pro5-BHI     3.16666667  1.5608503 4.7724830 0.0000060
#Pro25-BHI    2.62500000  1.1015888 4.1484112 0.0000723
#But5-BHI     2.62500000  1.1015888 4.1484112 0.0000723
#But25-BHI    2.25000000  0.7265888 3.7734112 0.0007957

pairwise.t.test(x=tox_r2$Toxin, g=tox_r2$group, p.adjust.method = "bonferroni")

#Pairwise comparisons using t tests with pooled SD 

#data:  tox_r2$Toxin and tox_r2$group 

#         BHI     Ace5    Ace25   Pro5    Pro25   But5   
#  Ace5  0.00015 -       -       -       -       -      
#  Ace25 0.00440 1.00000 -       -       -       -      
#  Pro5  6.3e-06 1.00000 0.21536 -       -       -      
#  Pro25 7.8e-05 1.00000 1.00000 1.00000 -       -      
#  But5  7.8e-05 1.00000 1.00000 1.00000 1.00000 -      
#  But25 0.00091 1.00000 1.00000 0.85243 1.00000 1.00000

#P value adjustment method: bonferroni 


## Figure S3B: VPI10463 toxin

tox_vpi <- read_tsv(file = "./ToxinData_ForGraphVPI.txt")
tox_vpi$group <- factor(tox_vpi$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))

ggplot(tox_vpi, aes(x=group, y=Toxin, color=group)) + geom_point(size =0.5, position = position_jitter(width = 0.2)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  #ylim(0,15) +
  #scale_y_continuous(limits= c(0, 10), breaks = seq(0,10,2))+
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_y_continuous(breaks = seq(0, 14, by = 2), limits = c(0,14))

## STATS

aov_tox_vpi <- aov(Toxin~group, data=tox_vpi)
summary(aov_tox_vpi)

#Df Sum Sq Mean Sq F value Pr(>F)  
#group        6  48.57   8.096   2.927 0.0126 *
#  Residuals   77 212.99   2.766                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=tox_vpi$Toxin, g=tox_vpi$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#              diff     lwr.ci   upr.ci   pval    
#Ace5-BHI  2.200000  0.1619868 4.238013 0.0297 *  
#Ace25-BHI 1.785714 -0.1185098 3.689938 0.0732 .  
#Pro5-BHI  2.300000  0.2619868 4.338013 0.0210 *  
#Pro25-BHI 2.357143  0.4529188 4.261367 0.0097 ** 
#But5-BHI  2.714286  0.8100616 4.618510 0.0021 ** 
#But25-BHI 2.714286  0.8100616 4.618510 0.0022 ** 

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

TukeyHSD(aov_tox_vpi)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = Toxin ~ group, data = tox_vpi)

#$group
#diff         lwr      upr     p adj
#Ace5-BHI     2.20000000 -0.18851788 4.588518 0.0910759
#Ace25-BHI    1.78571429 -0.44600500 4.017434 0.2034754
#Pro5-BHI     2.30000000 -0.08851788 4.688518 0.0667201
#Pro25-BHI    2.35714286  0.12542357 4.588862 0.0315777
#But5-BHI     2.71428572  0.48256643 4.946005 0.0074811
#But25-BHI    2.71428572  0.48256643 4.946005 0.0074811

pairwise.t.test(x=tox_vpi$Toxin, g=tox_vpi$group, p.adjust.method = "bonferroni")

#Pairwise comparisons using t tests with pooled SD 

#data:  tox_vpi$Toxin and tox_vpi$group 

#BHI   Ace5  Ace25 Pro5  Pro25 But5 
#Ace5  0.140 -     -     -     -     -    
#  Ace25 0.373 1.000 -     -     -     -    
#  Pro5  0.098 1.000 1.000 -     -     -    
#  Pro25 0.042 1.000 1.000 1.000 -     -    
#  But5  0.009 1.000 1.000 1.000 1.000 -    
#  But25 0.009 1.000 1.000 1.000 1.000 1.000

#P value adjustment method: bonferroni

# Figure S3C: CD630 time course with all SCFA

tox_630_time <- read_tsv(file = "./toxin_time_allSCFA_Julian.txt")

tox_630_time$group <- factor(tox_630_time$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))

ggplot(tox_630_time, aes(x=as.factor(time), y=Toxin, color=group)) + geom_point(size =0.5, position = position_jitterdodge(jitter.width = 0.5)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  ylim(0,6) +
  #scale_y_continuous(limits= c(0, 10), breaks = seq(0,10,2))+
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c"))

## STATS

#aov_tox_630_time <- aov(Toxin~group + as.factor(time), data=tox_630_time)
#summary(aov_tox_630_time)

tox_630_time6 <- subset(tox_630_time, tox_630_time$time == 6)
aov_tox_630_time6 <- aov(Toxin~group, data = tox_630_time6)
summary(aov_tox_630_time6)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#group        6  7.842  1.3070   16.21 2.39e-08 ***
#  Residuals   31  2.500  0.0806                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_tox_630_time6)
#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = Toxin ~ group, data = tox_630_time6)

#$group
#diff        lwr       upr     p adj
#Ace5-BHI     1.833333e+00  1.1030174 2.5636492 0.0000001
#Ace25-BHI    2.000000e+00  1.2696841 2.7303159 0.0000000
#Pro5-BHI     2.166667e+00  1.4363508 2.8969826 0.0000000
#Pro25-BHI    2.000000e+00  1.2696841 2.7303159 0.0000000
#But5-BHI     2.000000e+00  1.2696841 2.7303159 0.0000000
#But25-BHI    1.833333e+00  1.1030174 2.5636492 0.0000001

DunnettTest(x=tox_630_time6$Toxin, g=tox_630_time6$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#              diff   lwr.ci   upr.ci    pval    
#Ace5-BHI  1.833333 1.234944 2.431723 2.3e-09 ***
#Ace25-BHI 2.000000 1.401611 2.598389 1.9e-10 ***
#Pro5-BHI  2.166667 1.568277 2.765056 2.9e-11 ***
#Pro25-BHI 2.000000 1.401611 2.598389 1.4e-09 ***
#But5-BHI  2.000000 1.401611 2.598389 2.1e-10 ***
#But25-BHI 1.833333 1.234944 2.431723 2.1e-09 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


tox_630_time12 <- subset(tox_630_time, tox_630_time$time == 12)
aov_tox_630_time12 <- aov(Toxin~group, data = tox_630_time12)
summary(aov_tox_630_time12)
TukeyHSD(aov_tox_630_time12)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = Toxin ~ group, data = tox_630_time12)

#$group
#diff         lwr         upr     p adj
#Ace5-BHI     1.166667e+00  0.08343311  2.24990022 0.0282089
#Ace25-BHI    1.833333e+00  0.75009978  2.91656689 0.0001551
#Pro5-BHI     1.500000e+00  0.41676645  2.58323355 0.0022922
#Pro25-BHI    1.833333e+00  0.75009978  2.91656689 0.0001551
#But5-BHI     1.000000e+00 -0.08323355  2.08323355 0.0856951
#But25-BHI    1.500000e+00  0.41676645  2.58323355 0.0022922
#Ace25-Ace5   6.666667e-01 -0.09929512  1.43262846 0.1214140
#Pro5-Ace5    3.333333e-01 -0.43262846  1.09929512 0.8126916
#Pro25-Ace5   6.666667e-01 -0.09929512  1.43262846 0.1214140
#But5-Ace5   -1.666667e-01 -0.93262846  0.59929512 0.9925116
#But25-Ace5   3.333333e-01 -0.43262846  1.09929512 0.8126916
#Pro5-Ace25  -3.333333e-01 -1.09929512  0.43262846 0.8126916
#Pro25-Ace25 -8.881784e-16 -0.76596179  0.76596179 1.0000000
#But5-Ace25  -8.333333e-01 -1.59929512 -0.06737154 0.0259476
#But25-Ace25 -3.333333e-01 -1.09929512  0.43262846 0.8126916
#Pro25-Pro5   3.333333e-01 -0.43262846  1.09929512 0.8126916
#But5-Pro5   -5.000000e-01 -1.26596179  0.26596179 0.4022293
#But25-Pro5  -4.440892e-16 -0.76596179  0.76596179 1.0000000
#But5-Pro25  -8.333333e-01 -1.59929512 -0.06737154 0.0259476

DunnettTest(x=tox_630_time12$Toxin, g=tox_630_time12$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#              diff    lwr.ci   upr.ci    pval    
#Ace5-BHI  1.166667 0.2791118 2.054222 0.00736 ** 
#Ace25-BHI 1.833333 0.9457784 2.720888 1.4e-05 ***
#Pro5-BHI  1.500000 0.6124451 2.387555 0.00062 ***
#Pro25-BHI 1.833333 0.9457784 2.720888 1.7e-05 ***
#But5-BHI  1.000000 0.1124451 1.887555 0.02393 *  
#But25-BHI 1.500000 0.6124451 2.387555 0.00061 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


tox_630_time18 <- subset(tox_630_time, tox_630_time$time == 18)
aov_tox_630_time18 <- aov(Toxin~group, data = tox_630_time18)
summary(aov_tox_630_time18)

#Df Sum Sq Mean Sq F value Pr(>F)
#group        6  1.208  0.2014   0.657  0.685
#Residuals   25  7.667  0.3067   

#TukeyHSD(aov_tox_630_time18)

DunnettTest(x=tox_630_time18$Toxin, g=tox_630_time18$group, "BHI")

#no sig

tox_630_time24 <- subset(tox_630_time, tox_630_time$time == 24)
aov_tox_630_time24 <- aov(Toxin~group, data = tox_630_time24)
summary(aov_tox_630_time24)
TukeyHSD(aov_tox_630_time24)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = Toxin ~ group, data = tox_630_time24)

#$group
#diff        lwr       upr     p adj
#Ace5-BHI     2.500000e+00  1.1035859 3.8964141 0.0001568
#Ace25-BHI    2.500000e+00  1.1834482 3.8165518 0.0000715
#Pro5-BHI     2.000000e+00  0.5280503 3.4719497 0.0038190
#Pro25-BHI    2.000000e+00  0.6509358 3.3490642 0.0015214
#But5-BHI     2.000000e+00  0.6509358 3.3490642 0.0015214
#But25-BHI    2.666667e+00  1.1947170 4.1386163 0.0001342

DunnettTest(x=tox_630_time24$Toxin, g=tox_630_time24$group, "BHI")

# Significance Ace5-BHI 0.0907




