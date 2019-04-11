# Running HMMER
Run both HMMER 2.3.2 and HMMER 3 on a set of Pfam HMMs, parse and combine results. We run both versions of HMMER because in practice we have found that the sets of domains returned are slightly different, particularly for *short* domains; running both versions guarantees the most exhaustive recovery of all potential domain matches.

## Downloading and installing HMMER
HMMER 2.3.2 (released October 2003) and HMMER 3.1b2 (released February 2015) can both be downloaded from http://hmmer.org/download.html. I recommend installing HMMER 2.3.2 first (note that the scripts in this repository require the HMMER 2.3.2 binaries to be renamed with a "232" suffix):

**NOTE:** *If you do not have sudo access* wherever you are installing (e.g., your institution's cluster), you *must* edit the script below to install locally. Specifically, (1) change where HMMER 2.3.2 files are copied, and (2) run `./configure --prefix=/somewhere/else/than/usr/local` before `make` and `make install` for HMMER 3.1b2.

```bash
sh install_hmmer.sh
```


## Downloading and formatting Pfam HMMs
The most current set of Pfam-A HMMs can be found at: 
ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz. To download this file then parse 
it into individual files to use with HMMER 2.3.2 and HMMER 3.1b2 (as the scripts in this repository 
require), run the following. Note that the default **pfam_path** is "respository_directory/pfam/".

```bash
python download_pfam.py --pfam_path pfam/
```

## Running HMMER and parsing results
The general format to run is as follows, where the ending value is 1 more than the actual ending 0-index. 
For instance, to run on HMMs 0, 1, and 2, we would call process_hmmer.py with --start 0 and --end **3**.

```bash
python process_hmmer.py --pfam_path pfam/ 
                        --pfam_version 32 
                        --fasta_infile <full path to NON-GZIPPED fasta-formatted sequence file> 
                        --results_path domains/ 
                        --start 0 
                        --end 10
```

**NOTE:** This step should be run in parallel on the cluster!

## Combining HMMER output

We run both HMMER 2.3.2 and HMMER 3.1b2 and therefore expect lots of duplicate domains. We combine all output across 
all Pfam HMMs with the following call, which will produce a file called `all-hmmer-results-by-prot-v32.txt.gz`. You can
change this name using the `--hmmer_results` argument, and you can change the subset of HMMER results to be included 
(if desired) by editing the script directly. 

```bash
python create_domain_output.py --concatenate_hmmer_results 
                               --pfam_path pfam/
                               --pfam_version 32
                               --fasta_infile <full path to NON-COMPRESSED fasta-formatted sequence file>
                               --results_path domains/processed-v32/
                               --hmmer_results domains/all-hmmer-results-by-prot-v32.txt.gz
```

We **must** run the previous function before calling the following, as the following depends on the intermediate 
results (i.e., found in `all-hmmer-results-by-prot-v32.txt.gz`). We run this step to restrict to domains that:

* are complete (i.e., matched from the very start to the very end of the HMM)
* passed the gathering threshold (taking into account both domain- and sequence-based cutoffs)
* have the appropriate residue at high information content positions (to remove "deprecated" domains)

```bash
python create_domain_output.py --filter_domains
                               --pfam_path pfam/
                               --pfam_version 32
                               --fasta_infile <full path to NON-COMPRESSED fasta-formatted sequence file>
                               --results_path domains/processed-v32/
                               --hmmer_results domains/all-hmmer-results-by-prot-v32.txt.gz
                               --processed_results domains/all-domains-by-prot-v32.txt.gz
```

This produces a single, beautiful file (called `all-domains-by-prot-v32.txt.gz`) with all the domains from your 
original FASTA file.
