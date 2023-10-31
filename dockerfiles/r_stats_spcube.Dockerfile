FROM r-base:4.2.1

RUN apt-get update -y && \
    apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev

RUN R -e "install.packages(c('BiocManager', 'curl', 'XML'))"

RUN R -e "BiocManager::install(c('DESeq2', 'EnrichmentBrowser', 'clusterProfiler'))"
