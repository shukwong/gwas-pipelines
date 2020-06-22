FROM python:3.8-slim-buster

MAINTAINER Wendy Wong <wendy.wong@gmail.com>

RUN apt-get update && apt-get install -y wget build-essential zlib1g-dev unzip 

RUN wget  https://bitbucket.org/gavinband/bgen/get/tip.zip && \
    unzip tip.zip 

RUN cd gavinband-bgen-fe6d17aa6933/ && ./waf configure && ./waf && \
    cp build/apps/bgenix /usr/local/bin/ &&  export PATH=/usr/local/bin/:$PATH 
 