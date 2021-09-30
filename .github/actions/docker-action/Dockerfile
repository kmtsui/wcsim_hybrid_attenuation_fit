# Container image that runs your code
FROM rootproject/root:6.22.08-centos7

# Install neccessary software
RUN yum -y install wget git make cmake3 gcc-c++ gcc binutils libX11-devel libXpm-devel libXft-devel libXext-devel python openssl-devel

# Copies your code file from your action repository to the filesystem path `/` of the container
COPY . .

# Compile the code
RUN cmake -B build
RUN cmake --build build
