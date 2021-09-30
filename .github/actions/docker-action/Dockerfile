# Container image that runs your code
FROM rootproject/root:6.22.08-centos7

# Copies your code file from your action repository to the filesystem path `/` of the container
COPY . .

# Compile the code
RUN cmake -B build
RUN cmake --build build
