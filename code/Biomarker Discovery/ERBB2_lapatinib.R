library(PharmacoGx)
library(Biobase)
library(calibrate)
library(forestplot)
library(survcomp)

# load UHN 2019 & GRAY 2017 PSets

load("/data/PSets/GRAY_2017.RData")
load("/data/PSets/UHN_2019.RData")



#### UHNBreast 2019 ####

# get RNAseq data
rnaSeq_UHN <- t(exprs(summarizeMolecularProfiles(UHNBreast2019,mDataType = "rnaseq.exp",fill.missing = F)))

drug <- "lapatinib"

# get AAC sensitivity values for lapatinib
sensitivity_drug <- summarizeSensitivityProfiles(pSet = UHNBreast2019,sensitivity.measure = "auc_recomputed",drugs = drug,fill.missing = F)

commonSamples <- intersect(names(sensitivity_drug),rownames(rnaSeq_UHN))
UHN_common <- commonSamples

UHNdf <- data.frame(name = commonSamples,
                    LAP = sensitivity_drug[commonSamples] ,
                    ERBB2 = rnaSeq_UHN[commonSamples,"ENSG00000141736"])

ci <- survcomp::concordance.index(sensitivity_drug[commonSamples], surv.time = -rnaSeq_UHN[commonSamples,"ENSG00000141736"], surv.event = rep(1,length(sensitivity_drug[commonSamples])),outx = F, method="noether")



UHNci <- ci$c.index
UHNpvalue <- ci$p.value
UHNserr <- ci$se
UHNupper <- ci$upper
UHNlower <- ci$lower
#UHNconpair <- ci$concordant.pairs

gene= "ERBB2"


#### GRAY 2017 ####

# get RNAseq data
rnaSeq_GRAY <- t(exprs(summarizeMolecularProfiles(GRAY2017,mDataType = "rnaseq.exp",fill.missing = F)))

drug <- "lapatinib"

# get AAC sensitivity values for lapatinib
sensitivity_drug <- summarizeSensitivityProfiles(pSet = GRAY2017,sensitivity.measure = "auc_recomputed",drugs = drug,fill.missing = F)

commonSamples <- intersect(names(sensitivity_drug),rownames(rnaSeq_GRAY))
GRAY_common <- commonSamples

GRAYdf <- data.frame(name = commonSamples,
                     LAP = sensitivity_drug[commonSamples] ,
                     ERBB2 = rnaSeq_GRAY[commonSamples,"ENSG00000141736"])

ci <- survcomp::concordance.index(sensitivity_drug[commonSamples], surv.time = -rnaSeq_GRAY[commonSamples,"ENSG00000141736"], surv.event = rep(1,length(sensitivity_drug[commonSamples])),outx = F, method="noether")


GRAYci <- ci$c.index
GRAYpvalue <- ci$p.value
GRAYserr <- ci$se
GRAYupper <- ci$upper
GRAYlower <- ci$lower
#GRAYconpair <- ci$concordant.pairs


combined_ci <- combine.est(
  c(
    GRAYci, UHNci
  ),
  c(
    GRAYserr, UHNserr
  ),na.rm = TRUE,hetero = TRUE)


combined_ci_lower <- combined_ci$estimate + qnorm(0.025, lower.tail=TRUE) *  combined_ci$se
combined_ci_upper <- combined_ci$estimate + qnorm(0.025, lower.tail=FALSE) *  combined_ci$se
combined_ci_p <- pnorm((combined_ci$estimate - 0.5)/combined_ci$se, lower.tail = combined_ci$estimate < 0.5) * 2



######CREATE FOREST PLOT######

c_indices <- structure(
  list(
    mean  = c(NA, GRAYci, UHNci, combined_ci$estimate),
    lower = c(NA, GRAYlower, UHNlower, combined_ci_lower),
    upper = c(NA, GRAYupper, UHNupper, combined_ci_upper)
  ),
  .Names = c("C-index    ", "lower", "upper"),
  row.names = c(NA, -2L), 
  class = "data.frame"
)

c_tabletext <- cbind(
  c("PSet", "GRAY", "UHNBreast", "Meta analysis"),
  c("N", 39, 50, 89), #common samples for each dataset
  c("C-index", formatC(GRAYci, format = "e", digits = 2), formatC(UHNci, format = "e", digits = 2), formatC(combined_ci$estimate, format = "e", digits = 2)),
  c("P-value", formatC(GRAYpvalue, format = "e", digits = 2), formatC(UHNpvalue, format = "e", digits = 2), formatC(combined_ci_p, format = "e", digits = 2))
)



fileName = "../results/Figure_3.pdf"
pdf(fileName, width=7, height=3, onefile=FALSE)

forestplot(c_tabletext, c_indices, new_page = TRUE, boxsize = 0.3, is.summary=c(T,F,F), xlab="Concordance Index", 
           title="", zero=c(.49, .51),hrzl_lines=list("2"=gpar(lty=2, columns=1:4, col = "#000044")),
           txt_gp=fpTxtGp(label=gpar(fontfamily = "", cex = 0.8, fontface=2),
                          ticks=gpar(fontfamily = "", cex=.5, fontface=1),
                          xlab=gpar(fontfamily = "", cex=0.8, fontface=2),
                          legend=gpar(fontfamily = "", cex = 1, fontface=1)),
           col=fpColors(box=RColorBrewer::brewer.pal(n=4, name="Set2"),
                        line=RColorBrewer::brewer.pal(n=4, name="Set2"),
                        summary="blue"),
           xticks= c(.4, .5, .6, .7, .8)
)

dev.off()