# bam2wig

Convert bam files to wiggle files for usage with AUGUSTUS (exonpart hints)

## INSTALLATION

### Install dependencies
    
- HTSlib version 1.10 or higher
  
  - For Ubuntu 20.04 and Debian 11 or later:
    ```
      sudo apt install samtools libhts-dev
    ```

  - **OR** install HTSlib from https://github.com/samtools
    
    ```
      apt install -y samtools libbz2-dev liblzma-dev libncurses5-dev libssl-dev libcurl3-dev

      git clone https://github.com/samtools/htslib.git /path/where/HTSlib/reside/htslib
      cd /path/where/HTSlib/reside/htslib
      autoheader
      autoconf
      ./configure
      make
    
      export TOOLDIR=/path/where/HTSlib/reside

    ```

### Compile
Once all dependencies are available, you can compile bam2wig using make.

```
  make
```

### Test

Examples:
A file tests/auxprogs/bam2wig/test_files/test.s.bam has been included.  
You can try bam2wig by trying out the following examples:

1) ``` ./bam2wig test.s.bam ```
2) ``` ./bam2wig -t "my_specified_track" -r chr3L test.s.bam ```

Example (2) can only be done if an index file for test.s.bam exists.  
Do "samtools index test.s.bam" and a file "test.s.bam.bai" will be generated.
