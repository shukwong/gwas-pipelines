FROM ubuntu:20.04

MAINTAINER Wendy Wong <wendy.wong@gmail.com>

RUN apt-get update && apt-get install -y wget libopenblas-base

RUN wget https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.6-Ubuntu16.04-x86_64.tgz && \
    tar -zxvf qctool_v2.0.6-Ubuntu16.04-x86_64.tgz
    
RUN cp qctool_v2.0.6-Ubuntu16.04-x86_64/qctool /usr/local/bin/ &&  export PATH=/usr/local/bin/:$PATH 

