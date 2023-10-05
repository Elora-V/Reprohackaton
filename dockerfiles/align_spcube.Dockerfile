FROM ubuntu:latest as builder

RUN apt-get update  \
    && apt-get install -y wget unzip

RUN wget -O bowtie-0.12.7-linux-x86_64.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download \
    && unzip bowtie-0.12.7-linux-x86_64.zip \
    && rm bowtie-0.12.7-linux-x86_64.zip

FROM ubuntu:latest as base

RUN apt-get update && apt-get install -y samtools

COPY --from=builder bowtie-0.12.7/bowtie /bin
COPY --from=builder bowtie-0.12.7/bowtie-build /bin
