FROM ubuntu:latest as builder

RUN apt-get update  \
    && apt-get install -y wget libxml-libxml-perl

RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.7/sratoolkit.3.0.7-ubuntu64.tar.gz \
    && tar zxvf sratoolkit.3.0.7-ubuntu64.tar.gz  \
    && rm sratoolkit.3.0.7-ubuntu64.tar.gz

FROM ubuntu:latest as base

RUN apt-get update && apt-get install -y \
    wget pigz

COPY --from=builder sratoolkit.3.0.7-ubuntu64/bin /bin


