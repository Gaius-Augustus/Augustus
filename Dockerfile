FROM ubuntu:20.04


# Install required packages
RUN apt-get update
RUN apt-get install -y build-essential wget git autoconf

# Install dependencies for AUGUSTUS comparative gene prediction mode (CGP)
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev
RUN apt-get install -y libsqlite3-dev libmysql++-dev

# Install dependencies for  the optional support of gzip compressed input files
RUN apt-get install -y libboost-iostreams-dev zlib1g-dev

# Install dependencies for bam2hints and filterBam 
RUN apt-get install -y libbamtools-dev

# Install samtools
RUN apt-get install samtools

# Install additional dependencies for bam2wig
RUN apt-get install -y libhts-dev

# Install additional dependencies for homgenemapping and utrrnaseq
RUN apt-get install -y libboost-all-dev

# Install hal
WORKDIR /root
RUN wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz
RUN tar xzf hdf5-1.10.1.tar.gz
WORKDIR /root/hdf5-1.10.1
RUN ./configure --enable-cxx
RUN make && make install
ENV PATH="${PATH}:/root/hdf5-1.10.1/hdf5/bin"
WORKDIR /root
RUN git clone https://github.com/benedictpaten/sonLib.git
WORKDIR /root/sonLib
RUN make
WORKDIR /root
RUN git clone https://github.com/ComparativeGenomicsToolkit/hal.git
WORKDIR /root/hal
ENV RANLIB=ranlib
RUN make
ENV PATH="${PATH}:/root/hal/bin"

# Clone AUGUSTUS repository
ADD / /root/augustus

# Build AUGUSTUS
WORKDIR "/root/augustus"
RUN make HTSLIBS='-lhts'
RUN make install
ENV PATH="/root/augustus/bin:${PATH}"

# Test AUGUSTUS
RUN make unit_test
