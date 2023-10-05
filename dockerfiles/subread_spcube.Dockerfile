FROM ubuntu:latest as builder

RUN apt-get update && apt-get install -y wget

RUN wget -O subread-1.4.6-Linux-x86_64.tar.gz https://sourceforge.net/projects/subread/files/subread-1.4.6/subread-1.4.6-Linux-x86_64.tar.gz/download \
    && tar zxvf subread-1.4.6-Linux-x86_64.tar.gz  \
    && rm subread-1.4.6-Linux-x86_64.tar.gz

FROM ubuntu:latest as base

COPY --from=builder subread-1.4.6-Linux-x86_64/bin /bin
