FROM python:3.6-slim-stretch

MAINTAINER Wendy Wong "wendy.wong@gmail.com"
#modified from https://hub.docker.com/r/crukcibioinformatics/crossmap/dockerfile
#and https://github.com/shaze/py3plink/blob/master/Dockerfile

RUN apt-get update && apt-get install wget unzip gcc liblzo2-dev libz-dev -y

RUN pip3 install CrossMap

RUN wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip && \ 
    unzip plink2_linux_x86_64.zip -d /usr/local/bin && \
    rm plink2_linux_x86_64.zip && chmod a+x /usr/local/bin/plink2
