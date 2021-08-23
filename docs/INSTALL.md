
# Building AUGUSTUS


## Download sources

To obtain the most recent complete version, clone the GitHub repository

    git clone https://github.com/Gaius-Augustus/Augustus.git

or, alternatively, [download](http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz) and unpack the AUGUSTUS source package with

    tar -xzf augustus.current.tar.gz


## Install dependencies

- Dependencies

  The following dependencies are required for AUGUSTUS:
  - for gzip compressed input:
   (set ZIPINPUT = false in [common.mk](../common.mk) if this feature is not required or the required libraries are not available)
    - libboost-iostreams-dev
    - zlib1g-dev
  - for [comparative AUGUSTUS](README-cgp.md) (multi-species, CGP):
    (set COMPGENEPRED = false in [common.mk](../common.mk) if the libraries required by the CGP version are not available. Augustus can then only be run in single-genome mode, which is what most users need.)
    - libgsl-dev
    - libboost-all-dev
    - libsuitesparse-dev
    - liblpsolve55-dev
    - libsqlite3-dev (add SQLITE = false to [common.mk](../common.mk) if this feature is not required or the required library is not available)
    - libmysql++-dev (add MYSQL = false to [common.mk](../common.mk) if this feature is not required or the required library is not available)
  - for compiling utilities bam2hints and filterBam:
    - libbamtools-dev
  - for compiling utility utrrnaseq:
    - libboost-all-dev (version must be >Boost_1_49_0)
  - for compiling utility bam2wig:
    - Follow [these instructions](../auxprogs/bam2wig/README.md). Note that it shouldn't be a problem to compile AUGUSTUS without bam2wig. In practice, you can simply use `bamToWig.py` to accomplish the same task.
  - For compiling homgenemapping
      (set BOOST = FALSE in [auxprogs/homgenemapping/src/Makefile](../auxprogs/homgenemapping/src/Makefile) if the option --printHomologs is not required or the required libraries are not available)
    - libboost-all-dev
  - for scripts:
    - Perl
    - Python3

- Install dependencies
  - as root use your package manager to install the desired dependencies

        apt-get update
        apt-get install build-essential wget git autoconf

        # Install dependencies for AUGUSTUS comparative gene prediction mode (CGP)
        apt-get install libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev
        apt-get install libsqlite3-dev libmysql++-dev

        # Install dependencies for the optional support of gzip compressed input files
        apt-get install libboost-iostreams-dev zlib1g-dev

        # Install dependencies for bam2hints and filterBam
        apt-get install libbamtools-dev

        # Install additional dependencies for bam2wig
        apt-get install samtools libhts-dev

        # Install additional dependencies for homGeneMapping and utrrnaseq
        apt-get install libboost-all-dev

        # Install additional dependencies for scripts
        apt-get install cdbfasta diamond-aligner libfile-which-perl libparallel-forkmanager-perl libyaml-perl libdbd-mysql-perl
        apt-get install --no-install-recommends python3-biopython

  - as a non-root user you can either
    - download package archives for the dependencies and install them into a specified directory (for an example see here: [Troubleshooting/MySQL](#MySQL) )
    - or install the dependencies from source
      
      in both cases, you need to set the appropriate include and library paths, e.g. in [common.mk](../common.mk) (see here: [Troubleshooting](#Troubleshooting) ).

## Compile from sources

Once all dependencies are available, you can compile AUGUSTUS using make.

    make augustus

After compilation has finished, the command bin/augustus should be executable and print a usage message.

For AUGUSTUS and utilities use

    make



## Troubleshooting

The proposed solutions have been tested with Ubuntu 20.04. For other systems/architectures, the paths and commands given may have to be adapted.  
With root rights, you can install the libraries in the default folders. To do this, please remove the prefix argument in the configure, cmake and bootstrap commands and INSTALL from the make call.

### [zlib](https://zlib.net) - library for compression methods

- possible error messages

        /usr/bin/ld: cannot find -lz
        /usr/bin/ld: genbank.o: in function boost::detail::sp_counted_impl_p<boost::iostreams::symmetric_filter ...
        genbank.cc: ... undefined reference to boost::iostreams::detail::zlib_base::reset(bool, bool)'

- solutions
  - switch off zlib usage by setting `ZIPINPUT = false` in [common.mk](../common.mk) and do not make utility programs **bam2wig**, **bam2hints** and **filterBam**
  - or install package `zlib1g-dev`
  - or install from source
    - execute

            mkdir -p /your/path/to/zlib/zlib_build
            cd /your/path/to/zlib/zlib_build
            wget -O zlib-1.2.11.tar.gz https://sourceforge.net/projects/libpng/files/zlib/1.2.11/zlib-1.2.11.tar.gz/download
            tar xzf zlib-1.2.11.tar.gz
            cd /your/path/to/zlib/zlib_build/zlib-1.2.11
            ./configure --prefix=/your/path/to/zlib/zlib_install
            make
            make install

    - add to common.mk

            INCLUDE_PATH_ZLIB := -I/your/path/to/zlib/zlib_install/include
            LIBRARY_PATH_ZLIB := -L/your/path/to/zlib/zlib_install/lib -Wl,-rpath,/your/path/to/zlib/zlib_install/lib

    - check: the missing file should be here

           /your/path/to/zlib/zlib_install/lib/libz.so


### [Boost](https://www.boost.org/) - C++ libraries

- possible error messages

        ../include/types.hh:16:10: fatal error: boost/archive/text_oarchive.hpp: No such file or directory
        ../include/genome.hh:24:10: fatal error: boost/graph/adjacency_list.hpp: No such file or directory
        /usr/bin/ld: cannot find -lboost_iostreams
        /usr/bin/ld: cannot find -lboost_serialization

- solutions
  - switch off boost usage by setting `BOOST = false`, `COMPGENEPRED = false` and `ZIPINPUT = false` in [common.mk](../common.mk) do not make utility program **utrrnaseq**
  - or install package `libboost-all-dev`
  - or install from source
    - execute (don't set ZLIB variables if zlib is installed on standard locations)

            mkdir -p /your/path/to/boost/boost_build
            cd /your/path/to/boost/boost_build
            wget -O boost_1_76_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.76.0/boost_1_76_0.tar.gz/download
            tar xzf boost_1_76_0.tar.gz
            cd /your/path/to/boost/boost_build/boost_1_76_0
            export ZLIB_INCLUDE=/your/path/to/zlib/zlib_install/include
            export ZLIB_LIBRARY_PATH=/your/path/to/zlib/zlib_install/lib
            ./bootstrap.sh --prefix=/your/path/to/boost/boost_install --with-libraries=all
            ./b2 install --prefix=/your/path/to/boost/boost_install

    - add to common.mk

            INCLUDE_PATH_BOOST := -I/your/path/to/boost/boost_install/include
            LIBRARY_PATH_BOOST := -L/your/path/to/boost/boost_install/lib -Wl,-rpath,/your/path/to/boost/boost_install/lib

    - check: the missing files should be here

           /your/path/to/boost/boost_install/include/boost/archive/text_oarchive.hpp
           /your/path/to/boost/boost_install/lib/libboost_iostreams.so


### [MySQL](https://tangentsoft.com/mysqlpp)

- possible error messages

        ../include/randseqaccess.hh:21:10: fatal error: mysql++/mysql++.h: No such file or directory
        /usr/bin/ld: cannot find -lmysqlclient
        /usr/bin/ld: cannot find -lmysqlpp

- solutions
  - switch off MySQL usage by setting `MYSQL = false` in [common.mk](../common.mk)
  - or install package `libmysql++-dev`
  - or as non-root user install package archives into a specified directory (beware of further recursive dependencies at link and runtime)
    - execute (replace the archives mentioned here with the ones necessary for your system/architecture)

            mkdir /your/path/to/mysql
            cd /your/path/to/mysql
            # libmysql++3v5:
            wget        libmysql++3v5_3.2.5-1build1_amd64.deb http://de.archive.ubuntu.com/ubuntu/pool/universe/m/mysql++/libmysql++3v5_3.2.5-1build1_amd64.deb
            dpkg-deb -x libmysql++3v5_3.2.5-1build1_amd64.deb /your/path/to/mysql/mysql_install
            # libmysql++-dev:
            wget        libmysql++-dev_3.2.5-1build1_amd64.deb http://de.archive.ubuntu.com/ubuntu/pool/universe/m/mysql++/libmysql++-dev_3.2.5-1build1_amd64.deb
            dpkg-deb -x libmysql++-dev_3.2.5-1build1_amd64.deb /your/path/to/mysql/mysql_install
            # libmysqlclient21:
            wget        libmysqlclient21_8.0.23-0ubuntu0.20.04.1_amd64.deb http://security.ubuntu.com/ubuntu/pool/main/m/mysql-8.0/libmysqlclient21_8.0.23-0ubuntu0.20.04.1_amd64.deb
            dpkg-deb -x libmysqlclient21_8.0.23-0ubuntu0.20.04.1_amd64.deb /your/path/to/mysql/mysql_install
            # libmysqlclient-dev:
            wget        libmysqlclient-dev_8.0.23-0ubuntu0.20.04.1_amd64.deb http://security.ubuntu.com/ubuntu/pool/main/m/mysql-8.0/libmysqlclient-dev_8.0.23-0ubuntu0.20.04.1_amd64.deb
            dpkg-deb -x libmysqlclient-dev_8.0.23-0ubuntu0.20.04.1_amd64.deb /your/path/to/mysql/mysql_install


    - add to common.mk

            INCLUDE_PATH_MYSQL := -I/your/path/to/mysql/mysql_install/usr/include -I/your/path/to/mysql/mysql_install/usr/include/mysql
            LIBRARY_PATH_MYSQL := -L/your/path/to/mysql/mysql_install/usr/lib/x86_64-linux-gnu -Wl,-rpath,/your/path/to/mysql/mysql_install/usr/lib/x86_64-linux-gnu

    - check: the missing files should be here

           /your/path/to/mysql/mysql_install/usr/include/mysql++/mysql++.h
           /your/path/to/mysql/mysql_install/usr/lib/x86_64-linux-gnu/libmysqlclient.so

  - installation from source code is quite time-consuming and requires many architecture-specific dependencies and is therefore not listed here


### [SQLite](https://www.sqlite.org)

- possible error messages

        ../include/sqliteDB.hh:13:10: fatal error: sqlite3.h: No such file or directory
        /usr/bin/ld: cannot find -lsqlite3

- solutions
  - switch off SQLite usage by setting `SQLITE = false` in [common.mk](../common.mk)
  - or install package `libsqlite3-dev`
  - or install from source
    - execute

            mkdir -p /your/path/to/sqlite3
            cd /your/path/to/sqlite3
            wget sqlite-autoconf-3350500.tar.gz https://www.sqlite.org/2021/sqlite-autoconf-3350500.tar.gz
            tar zxf sqlite-autoconf-3350500.tar.gz
            cd /your/path/to/sqlite3/sqlite-autoconf-3350500
            ./configure --prefix=/your/path/to/sqlite3/sqlite3_install
            make
            make install


    - add to common.mk

            INCLUDE_PATH_SQLITE := -I/your/path/to/sqlite3/sqlite3_install/include
            LIBRARY_PATH_SQLITE := -L/your/path/to/sqlite3/sqlite3_install/lib -Wl,-rpath,/your/path/to/sqlite3/sqlite3_install/lib

    - check: the missing files should be here

           /your/path/to/qlite3/sqlite3_install/include/sqlite3.h
           /your/path/to/sqlite3/sqlite3_install/lib/libsqlite3.so


### [GSL](https://www.gnu.org/software/gsl/) - GNU Scientific Library

- possible error messages

        parser/../../include/contTimeMC.hh:18:10: fatal error: gsl/gsl_matrix.h: No such file or directory
        /usr/bin/ld: cannot find -lgsl

- solutions
  - switch off gsl usage by setting `COMPGENEPRED = false` in [common.mk](../common.mk)
  - or install package `libgsl-dev`
  - or install from source
    - execute

            mkdir /your/path/to/gsl
            cd /your/path/to/gsl
            wget gsl-latest.tar.gz https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz
            tar zxf gsl-latest.tar.gz
            cd /your/path/to/gsl/gsl-2.6/
            ./configure --prefix=/your/path/to/gsl/gsl_install
            make
            make install

    - add to common.mk

            INCLUDE_PATH_GSL := -I/your/path/to/gsl/gsl_install/include
            LIBRARY_PATH_GSL := -L/your/path/to/gsl/gsl_install/lib -Wl,-rpath,/your/path/to/gsl/gsl_install/lib

    - check: the missing files should be here

           /your/path/to/gsl/gsl_install/include/gsl/gsl_matrix.h
           /your/path/to/gsl/gsl_install/lib/libgsl.so


### [lp_solve](http://lpsolve.sourceforge.net/) - Mixed Integer Linear Programming (MILP) solver

- possible error messages

        alignment.cc:16:10: fatal error: lp_lib.h: No such file or directory
        /usr/bin/ld: cannot find -llpsolve55

- solutions
  - switch off lp_solve usage by setting `COMPGENEPRED = false` in [common.mk](../common.mk)
  - or install package `liblpsolve55-dev`
  - or install from source
    - execute

            mkdir /your/path/to/lpsolve
            cd /your/path/to/lpsolve
            wget -O lp_solve_5.5.2.11_source.tar.gz https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz/download
            tar zxf lp_solve_5.5.2.11_source.tar.gz
            cd /your/path/to/lpsolve/lp_solve_5.5
            chmod +x ./lpsolve55/ccc ./lp_solve/ccc
            cd /your/path/to/lpsolve/lp_solve_5.5/lpsolve55
            ./ccc
            cd /your/path/to/lpsolve/lp_solve_5.5/lp_solve
            ./ccc
            mkdir -p /your/path/to/lpsolve/lpsolve_install/include
            cp /your/path/to/lpsolve/lp_solve_5.5/*.h /your/path/to/lpsolve/lpsolve_install/include/
            mkdir -p /your/path/to/lpsolve/lpsolve_install/lib
            cp /your/path/to/lpsolve/lp_solve_5.5/lpsolve55/bin/ux64/liblpsolve55.* /your/path/to/lpsolve/lpsolve_install/lib/

    - add to common.mk

            INCLUDE_PATH_LPSOLVE := -I/your/path/to/lpsolve/lpsolve_install/include
            LIBRARY_PATH_LPSOLVE := -L/your/path/to/lpsolve/lpsolve_install/lib -Wl,-rpath,/your/path/to/lpsolve/lpsolve_install/lib

    - check: the missing files should be here

           /your/path/to/lpsolve/lpsolve_install/include/lp_lib.h
           /your/path/to/lpsolve/lpsolve_install/lib/liblpsolve55.so


### [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html) - A Suite of Sparse matrix software

- possible error messages

        /usr/bin/ld: cannot find -lcolamd

- solutions
  - switch off SuiteSparse usage by setting `COMPGENEPRED = false` in [common.mk](../common.mk)
  - or install package `libsuitesparse-dev`
  - or install from source
    - execute

            git clone https://github.com/DrTimothyAldenDavis/SuiteSparse.git /your/path/to/suitesparse/suitesparse
            cd /your/path/to/suitesparse/suitesparse/SuiteSparse_config
            make install INSTALL=/your/path/to/suitesparse/suitesparse_install
            cd /your/path/to/suitesparse/suitesparse/COLAMD
            make install INSTALL=/your/path/to/suitesparse/suitesparse_install

    - add to common.mk

            INCLUDE_PATH_SUITESPARSE := -I/your/path/to/suitesparse/suitesparse_install/include
            LIBRARY_PATH_SUITESPARSE := -L/your/path/to/suitesparse/suitesparse_install/lib -Wl,-rpath,/your/path/to/suitesparse/suitesparse_install/lib

    - check: the missing file should be here

           /your/path/to/suitesparse/suitesparse_install/lib/libcolamd.so


### [Bamtools](https://github.com/pezmaster31/bamtools) - handling BAM files

- possible error messages

        bam2hints.cc:16:10: fatal error: api/BamReader.h: No such file or directory
        /usr/bin/ld: cannot find -lbamtools

- solutions
  - do not make utility programs **bam2hints** and **filterBam**
  - or install package `libbamtools-dev`
  - or install from source
    - cmake or cmake3 must be installed (use the installed one in next steps)

               apt-get install cmake

      (or see https://cmake.org)

    - download sources

            git clone https://github.com/pezmaster31/bamtools.git /your/path/to/bamtools/bamtools

      or alternatively download the zipped version from https://github.com/pezmaster31/bamtools/zipball/master and unpack:

            mkdir -p /your/path/to/bamtools/bamtools
            unzip pezmaster31-bamtools-v2.5.1-23-g2391b1a.zip -d /your/path/to/bamtools/bamtools
            mv /your/path/to/bamtools/bamtools/pezmaster31-bamtools-2391b1a/* /your/path/to/bamtools/bamtools/

    - execute (don't set ZLIB variables if zlib is installed on standard locations)

            mkdir -p /your/path/to/bamtools/bamtools_build
            cd /your/path/to/bamtools/bamtools_build
            cmake \
              -DCMAKE_INSTALL_PREFIX=/your/path/to/bamtools/bamtools_install \
              -DZLIB_LIBRARY=/your/path/to/zlib/zlib_install/lib/libz.so \
              -DZLIB_INCLUDE_DIR=/your/path/to/zlib/zlib_install/include \
              /your/path/to/bamtools/bamtools
            make
            make install

    - add to common.mk

            INCLUDE_PATH_BAMTOOLS := -I/your/path/to/bamtools/bamtools_install/include/bamtools
            LIBRARY_PATH_BAMTOOLS := -L/your/path/to/bamtools/bamtools_install/lib -Wl,-rpath,/your/path/to/bamtools/bamtools_install_dir/lib

    - check: the missing files should be here

           /your/path/to/bamtools/bamtools_install/include/bamtools/api/BamReader.h
           /your/path/to/bamtools/bamtools_install/lib/libbamtools.a


### [HTSlib](https://github.com/samtools/htslib) - accessing common file formats, such as SAM, CRAM and VCF

- possible error messages

        bam2wig.c:12:10: fatal error: bgzf.h: No such file or directory
        /usr/bin/ld: cannot find -lhts

- solutions
  - do not make utility program **bam2wig**
  - or install packages `samtools libhts-dev`
  - or install from source
    - execute (don't set ZLIB variables if zlib is installed on standard locations, set variables for bz2 and lzma if available and used
      or don't set the disable options if bz2 and lzma are installed on standard locations)

            git clone https://github.com/samtools/htslib.git /your/path/to/htslib/htslib_build
            cd /your/path/to/htslib/htslib_build
            git submodule update --init --recursive
            autoreconf -i
            ./configure --prefix=/your/path/to/htslib/htslib_install \
              CPPFLAGS=-I/your/path/to/zlib/zlib_install/include \
              LDFLAGS=-L/your/path/to/zlib/zlib_install/lib \
              --disable-bz2 \
              --disable-lzma
            make
            make install

    - add to common.mk

            INCLUDE_PATH_HTSLIB   := -I/your/path/to/htslib/htslib_install/include/htslib
            LIBRARY_PATH_HTSLIB   := -L/your/path/to/htslib/htslib_install/lib -Wl,-rpath,/your/path/to/htslib/htslib_install/lib

    - check: the missing files should be here

           /your/path/to/htslib/htslib_install/include/htslib/bgzf.h
           /your/path/to/htslib/htslib_install/lib/libhts.a

### Test failures

#### [SAMtools](https://github.com/samtools/samtools)

- possible error messages

        samtools: not found

- solutions
  - install package `apt-get install samtools`
  - or install from source
    - execute (don't set ZLIB and htslib variables and flags if both are installed on standard locations, remove the --without-curses if this library is installed on standard locations)


            git clone https://github.com/samtools/samtools.git /your/path/to/samtools/samtools_build
            cd /your/path/to/samtools/samtools_build
            autoheader
            autoconf -Wno-syntax
            ./configure \
              --prefix=/your/path/to/samtools/samtools_install \
              --with-htslib=/your/path/to/htslib/htslib_install \
              --without-curses \
              CPPFLAGS=-I/your/path/to/zlib/zlib_install/include \
              LDFLAGS="-L/your/path/to/zlib/zlib_install/lib -Wl,-rpath,/your/path/to/zlib/zlib_install/lib -L/your/path/to/htslib/htslib_install/lib -Wl,-rpath,/your/path/to/htslib/htslib_install/lib"
            make
            make install

    - change PATH in your shell

            export PATH="$PATH:/your/path/to/samtools/samtools_install/bin/"


#### [HAL](https://github.com/ComparativeGenomicsToolkit/hal) - Hierarchical Alignment (HAL) Format API

- possible error messages

        halLiftover is not executable

- solution
  - install from source
    - install package `apt-get install libhdf5-dev`
    - execute (don't set ZLIB variables if zlib is installed on standard locations)

            # install hdf5 from source if libhdf5-dev is not installed:
            mkdir -p /your/path/to/hdf5/hdf5_build
            cd /your/path/to/hdf5/hdf5_build
            wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz
            tar xzf hdf5-1.10.1.tar.gz
            cd /your/path/to/hdf5/hdf5_build/hdf5-1.10.1
            ./configure --enable-cxx --prefix=/your/path/to/hdf5/hdf5_install \
               --with-zlib="/your/path/to/zlib/zlib_install/include,/your/path/to/zlib/zlib_install/lib"
            make
            make install
            export PATH="$PATH:/your/path/to/hdf5/hdf5_install/bin"
            
            # install sonLib
            git clone https://github.com/benedictpaten/sonLib.git /your/path/to/sonLib
            cd /your/path/to/sonLib
            export CPPFLAGS="-I/your/path/to/zlib/zlib_install/include -L/your/path/to/zlib/zlib_install/lib"
            make
            unset CPPFLAGS
            
            # install hal
            git clone https://github.com/ComparativeGenomicsToolkit/hal.git /your/path/to/hal
            cd /your/path/to/hal
            export RANLIB=ranlib
            make

    - change PATH in your shell

            export PATH="$PATH:/your/path/to/hdf5/hdf5_install/bin"
            export PATH="$PATH:/your/path/to/hal/bin/"
