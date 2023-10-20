FROM r-base:4.2.1
WORKDIR /usr/local/src
RUN apt-get update \
&& apt-get install libxml2-dev -y \
&& apt-get install libssl-dev -y \
&& apt-get install libcurl4-openssl-dev -y 


RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('curl')"
RUN R -e "install.packages('XML')"

RUN R -e "library(BiocManager)"
RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('EnrichmentBrowser')"
RUN R -e "BiocManager::install('clusterProfiler')"

RUN R -e "rownames(installed.packages())"

