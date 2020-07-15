FROM centos:centos8

MAINTAINER Wendy Wong <wendy.wong@gmail.com>

RUN yum -y install wget gcc-c++ blas lapack

RUN wget https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz && \
    tar -zxvf qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64.tgz 

RUN cp qctool_v2.0.8-CentOS\ Linux7.6.1810-x86_64/qctool /usr/local/bin/ &&  export PATH=/usr/local/bin/:$PATH 

