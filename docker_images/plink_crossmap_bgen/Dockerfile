FROM python:3.6-slim-stretch

MAINTAINER Wendy Wong "wendy.wong@gmail.com"
#modified from https://hub.docker.com/r/crukcibioinformatics/crossmap/dockerfile
#and https://github.com/shaze/py3plink/blob/master/Dockerfile

RUN apt-get update && apt-get install wget unzip gcc liblzo2-dev libz-dev -y

RUN pip3 install CrossMap

RUN wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip && \
    unzip plink2_linux_x86_64.zip -d /usr/local/bin && \
    rm plink2_linux_x86_64.zip && chmod a+x /usr/local/bin/plink2

RUN apt-get update && apt-get install -y wget build-essential zlib1g-dev unzip

RUN wget  https://enkre.net/cgi-bin/code/bgen/tarball/6ac2d582f9/BGEN-6ac2d582f9.tar.gz && \
    tar -xzvf BGEN-6ac2d582f9.tar.gz

RUN cd BGEN-6ac2d582f9/ && ./waf configure && ./waf && \
    cp ./build/apps/bgenix /usr/local/bin/ &&  export PATH=/usr/local/bin/:$PATH    