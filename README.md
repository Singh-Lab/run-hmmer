# Running HMMER
Run both HMMER 2.3.2 and HMMER 3 on a set of Pfam HMMs, parse and combine results.

## Downloading and installing HMMER
HMMER 2.3.2 (released October 2003) and HMMER 3.0 (released March 2010) can both be downloaded from http://hmmer.org/download.html. I recommend installing HMMER 2.3.2 first (note that the scripts in this repository require the HMMER 2.3.2 binaries to be renamed with the "232" suffix as follows):

```bash
wget http://eddylab.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz
tar -xvzf hmmer-2.3.2.tar.gz
cd hmmer-2.3.2
./configure
make
make check
for file in hmmalign hmmbuild hmmcalibrate hmmconvert hmmemit hmmfetch hmmindex hmmpfam hmmsearch ; do
  cp src/$file /usr/local/bin/${file}232;
done
cd ../
rm -rf hmmer-2.3.2*
```

To download and install HMMER 3.0:

```bash
wget http://eddylab.org/software/hmmer3/3.0/hmmer-3.0.tar.gz
tar -xvzf hmmer-3.0.tar.gz
cd hmmer-3.0
./configure
make
make check
make install
cd ../
rm -rf hmmer-3.0*
```

## Downloading and parsing Pfam HMMs
You may only want to download a subset of Pfam HMMs. The most recent set of Pfam HMMs can always be found at ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz; and the release number can be found in the second line of ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt. 






