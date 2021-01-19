FROM ubuntu:20.04

# Set timezone in tzdata
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Berlin"

# Install required packages
RUN apt-get update
RUN apt-get install -y build-essential wget git autoconf

# Install dependencies for AUGUSTUS comparative gene prediction mode (CGP)
RUN apt-get install -y libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev
RUN apt-get install -y libsqlite3-dev libmysql++-dev

# Install dependencies for the optional support of gzip compressed input files
RUN apt-get install -y libboost-iostreams-dev zlib1g-dev

# Install dependencies for bam2hints and filterBam 
RUN apt-get install -y libbamtools-dev

# Install additional dependencies for bam2wig
RUN apt-get install -y samtools libhts-dev

# Install additional dependencies for homGeneMapping and utrrnaseq
RUN apt-get install -y libboost-all-dev

# Install additional dependencies for scripts
RUN apt-get install -y cdbfasta diamond-aligner libfile-which-perl libparallel-forkmanager-perl libyaml-perl libdbd-mysql-perl
RUN apt-get install -y --no-install-recommends python3-biopython

# Install hal - required by homGeneMapping 
# execute the commented out code if you want to use this program - see auxprogs/homGeneMapping/Dockerfile
#RUN apt-get install -y libhdf5-dev
#RUN git clone https://github.com/benedictpaten/sonLib.git
#WORKDIR /root/sonLib
#RUN make
#WORKDIR /root
#RUN git clone https://github.com/ComparativeGenomicsToolkit/hal.git
#WORKDIR /root/hal
#ENV RANLIB=ranlib
#RUN make
#ENV PATH="${PATH}:/root/hal/bin"

# Clone AUGUSTUS repository
ADD / /root/augustus

# Build AUGUSTUS
WORKDIR "/root/augustus"
RUN make clean
RUN make
RUN make install
ENV PATH="/root/augustus/bin:/root/augustus/scripts:${PATH}"

# Test AUGUSTUS
RUN make unit_test
