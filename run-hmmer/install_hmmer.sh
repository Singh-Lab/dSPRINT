#!/usr/bin/env bash

# HMMER 2.3.2 (released October 2003) and HMMER 3.1b2 (released February 2015)
# can both be downloaded from http://hmmer.org/download.html.
# I recommend installing HMMER 2.3.2 first (note that the scripts in this repository require the HMMER 2.3.2 binaries
# to be renamed with the "232" suffix as follows):

# NOTE: I got a "Can't locate getopts.pl in @INC" failure when I tried to install HMMER 2.3.2 and needed to run:
# sudo apt-get install libperl4-corelibs-perl

wget http://eddylab.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz
tar -xvzf hmmer-2.3.2.tar.gz
cd hmmer-2.3.2
./configure
make
make check
for file in hmmalign hmmbuild hmmcalibrate hmmconvert hmmemit hmmfetch hmmindex hmmpfam hmmsearch ; do
  sudo cp src/$file /usr/local/bin/${file}232;
done
cd ../
rm hmmer-2.3.2.tar.gz


# And download and install HMMER 3.1b2 (as the default!)

wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2.tar.gz
tar -xvzf hmmer-3.1b2.tar.gz
cd hmmer-3.1b2
./configure
make
make check
sudo make install
cd ../
rm hmmer-3.1b2.tar.gz