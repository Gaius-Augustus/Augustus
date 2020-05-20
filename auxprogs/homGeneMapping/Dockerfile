FROM ubuntu:18.04


# Install required packages
RUN apt-get update
RUN apt-get install -y build-essential wget git autoconf
RUN apt-get install -y zlib1g-dev

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

# Install dependencies for homGeneMapping
RUN apt-get update
RUN apt-get install -y libgsl-dev liblpsolve55-dev
RUN apt-get install -y libsqlite3-dev  libmysql++-dev
RUN apt-get install -y libboost-all-dev

# Clone AUGUSTUS repository
WORKDIR /root
RUN git clone https://github.com/Gaius-Augustus/Augustus.git
WORKDIR /root/Augustus

# Build load2sqlitedb
WORKDIR /root/Augustus/src/
RUN mkdir -p /root/Augustus/bin
RUN make load2sqlitedb
ENV PATH="${PATH}:/root/Augustus/bin"

# Build homGeneMapping
ADD / /root/Augustus/auxprogs/homGeneMapping
WORKDIR /root/Augustus/auxprogs/homGeneMapping/
RUN make clean
RUN make

# Test homGeneMapping
RUN make test
