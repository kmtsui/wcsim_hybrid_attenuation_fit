# Container image that runs your code
FROM rootproject/root:6.22.08-centos7

# Install neccessary software
#RUN yum -y install wget git make cmake3 gcc-c++ gcc binutils libX11-devel libXpm-devel libXft-devel libXext-devel python openssl-devel

# Copies your code file from your action repository to the filesystem path `/` of the container
#COPY ../../../../ /

RUN git clone https://github.com/kmtsui/wcsim_hybrid_attenuation_fit.git optical_fit

RUN pwd

RUN ls -lh

# Compile the code
RUN cmake3 -B build optical_fit
RUN cmake3 --build build
