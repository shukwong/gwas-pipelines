FROM ubuntu:20.04

#modified from https://github.com/manning-lab/metal/blob/master/Dockerfile

RUN apt-get update && apt-get install -y build-essential zlib1g-dev wget

RUN cd /bin
RUN wget -O generic-metal-2011-03-25.tar.gz http://csg.sph.umich.edu/abecasis/Metal/download/generic-metal-2011-03-25.tar.gz
RUN tar -xzf ./generic-metal-2011-03-25.tar.gz && rm generic-metal-2011-03-25.tar.gz && cd generic-metal && make all && mv executables/metal /bin && cd /
