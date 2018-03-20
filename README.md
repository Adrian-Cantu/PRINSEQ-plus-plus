![Prinseq++](prinseq_logo.png)

PRINSEQ++ is a C++ implementation of the prinseq-lite.pl program. It can be used to filter, reformat or trim genomic and metagenomic sequence data. It is 5X faster than prinseq-lite.pl and uses less RAM thanks to the use of multi-threading and the cboost libraries. It can read and write compressed (gzip) files, drastically reducing the use of hard drive.

## Requierments
1. automake
2. g++
3. make
4. boost-devel
5. pthread

## Download
If you are just interested in compiling and using PRINSEQ++, download the latest [version](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus/releases/download/v1.0/prinseq++-1.0.tar.gz).
If you want to edit the source code, clone this repository

## To install
1. tar -xvf pprinseqc-0.9.1.tar.gz
2. cd pprinseqc-0.9.1
3. ./configure
4. make
4. sudo make install

## To use the repository
1. ./autogen.sh
2. ./configure
3. make
4. make test 



