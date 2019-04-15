library(PharmacoGx)
library(Biobase)
library(calibrate)
library(mCI)
library(forestplot)

# load UHN & GRAY PSets

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

ci <- paired.concordance.index(sensitivity_drug[commonSamples], rnaSeq_UHN[commonSamples,"ENSG00000141736"], delta.pred=0, delta.obs=0)

UHNci <- ci$cindex
UHNpvalue <- ci$p.value
UHNserr <- ci$sterr
UHNupper <- ci$upper
UHNlower <- ci$lower
UHNconpair <- ci$concordant.pairs

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

ci <- paired.concordance.index(sensitivity_drug[commonSamples], rnaSeq_GRAY[commonSamples,"ENSG00000141736"], delta.pred=0, delta.obs=0)

GRAYci <- ci$cindex
GRAYpvalue <- ci$p.value
GRAYserr <- ci$sterr
GRAYupper <- ci$upper
GRAYlower <- ci$lower
GRAYconpair <- ci$concordant.pairs

gene= "ERBB2"



#COMBINED FOREST PLOT

c_indices <- structure(
  list(
    mean  = c(NA, GRAYci, UHNci),
    lower = c(NA, GRAYlower, UHNlower),
    upper = c(NA, GRAYupper, UHNupper)
  ),
  .Names = c("C-index    ", "lower", "upper"),
  row.names = c(NA, -2L), 
  class = "data.frame"
)

c_tabletext <- cbind(
  c("PSet", "GRAY", "UHNBreast"),
  c("N", 39, 50), #common samples for each dataset
  c("C-index", format(round(GRAYci, 4), nsmall = 4), format(round(UHNci, 4), nsmall = 4)),
  c("P-value", formatC(GRAYpvalue, format = "e", digits = 2), formatC(UHNpvalue, format = "e", digits = 2))
)


c_fn <- local({
  i = 0
  b_clrs =  c("red","red")
  l_clrs =  c("red","red")
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
  }
})

c_fn1 <- local({
  i = 0
  s_clrs =c()
  function(..., col){
    i <<- i + 1
    fpDrawSummaryCI(...,col=s_clrs[i])
  }
})

fileName = "../results/Figure_3.pdf"
pdf(fileName, width = 10, height = 10, onefile=FALSE)

forestplot(c_tabletext, c_indices, new_page = TRUE, is.summary=c(T,F,F), xlab="C-index", 
           title="", clip=c(0,1), xlog=FALSE, col=fpColors(text="black"),
           zero=0.5, align='l', boxsize=0.2, txt_gp = fpTxtGp(ticks = gpar(cex=0.8), xlab  = gpar(cex = 1)), 
           graphwidth=unit(4, "inches"), fn.ci_norm=c_fn,  fn.ci_sum=c_fn1
)

dev.off()
