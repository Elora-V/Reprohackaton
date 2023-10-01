FROM ubuntu:latest as builder

RUN apt-get update && apt-get -y install git

RUN git clone https://github.com/FelixKrueger/TrimGalore.git

FROM ubuntu:latest as base

RUN apt-get update && apt-get -y install python3-pip fastqc
RUN pip3 install cutadapt==4.4 #difference avec l'article : v1.11

COPY --from=builder /TrimGalore/trim_galore /bin/