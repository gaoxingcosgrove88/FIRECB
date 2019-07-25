FROM registry.codeocean.com/codeocean/r-base:3.6.0-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && apt-key adv --keyserver "hkp://keyserver.ubuntu.com:80" --recv-keys \
      "0xAD2323F17326AE31401037733E05EBFF05441C52" \
    && add-apt-repository -y "deb http://deb.codeocean.com/rstudio-server-trusty/ ubuntu main" \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
      "pandoc=1.19.2.4~dfsg-1build4" \
      "rstudio-server=1.2.1335" \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(warn=2); devtools::install_version("calibrate", version = "1.7.2")'
RUN Rscript -e 'options(warn=2); devtools::install_version("forestplot", version = "1.7.2")'

RUN Rscript -e 'devtools::install_github("bhklab/wCI", \
      upgrade_dependencies = FALSE, \
      ref = "v0.1.0")'

COPY postInstall /
RUN /postInstall
