## PGxCWL: Creating reproducible pharmacogenomic analysis pipelines

Anthony Mammoliti, Petr Smirnov, Zhaleh Safikhani, Wail Ba-Alawi, Benjamin Haibe-Kains


This capsule contains our reproducible CWL workflows for generating PSets for the following breast cancer pharmacogenomic datasets:

**1.** GRAY (2013, 2017)

**2.** UHNBreast (2017, 2019)

PSets available on Harvard Dataverse at the following DOI: https://doi.org/10.7910/DVN/BXIY5W

## CWL Workflow Execution Instructions:

**1.** Install cwltool (version 1.0.20190228155703)

**2.** Unzip /code/Compiled CWL Pipelines/**CWL_PSets.zip**

**3.** Proceed to a workflow directory for a PSet you would like to generate 
                        (e.g. `getUHN2017_Workflow`)

**4.** Run cwltool on the .cwl and .yml files within the workflow directory: 

        cwltool getUHN2017_Workflow.cwl getUHN2017_Workflow.yml

*To utilize the data provenance option to create a Research Object, use the --provenance flag:*

 `cwltool --provenance /outputdir getUHN2017_Workflow.cwl getUHN2017_Workflow.yml`
 

## Workflow Raw Data:

The raw data used in our CWL workflows are available under /data/CWL Raw Data. 

**NOTE:** CWL_PSets.zip already contains this data. 



## Biomarker Discovery

To investigate the gene-drug association between lapatinib and ERBB2 amplification in GRAY 2017 and UHN 2019, please execute `run.sh`. 

To run this PSet biomarker discovery locally:

**1.** Copy the GRAY 2017 and UHN 2019 PSets into one directory

**2.** Set this directory to the working directory in a R session

**3.** `Execute ERBB2_lapatinib.R`
