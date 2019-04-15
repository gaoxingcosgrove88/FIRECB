FROM registry.codeocean.com/codeocean/r-base:3.4.4-ubuntu16.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'options(warn=2); devtools::install_version("calibrate", version = "1.7.2")'
RUN Rscript -e 'options(warn=2); devtools::install_version("forestplot", version = "1.7.2")'

RUN Rscript -e ' \
      source("http://bioconductor.org/biocLite.R"); \
      biocLite(c( \
        "Biobase", \
        "PharmacoGx" \
      ), suppressUpdates = TRUE)'

RUN Rscript -e 'devtools::install_github("bhklab/mCI", \
      upgrade_dependencies = FALSE, \
      ref = "v0.1.0")'

COPY postInstall /
RUN /postInstall
