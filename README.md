# Running HMMER
Run both HMMER 2.3.2 and HMMER 3 on a set of Pfam HMMs, parse and combine results.

## Downloading and installing HMMER
HMMER 2.3.2 (released October 2003) and HMMER 3.0 (released March 2010) can both be downloaded from http://hmmer.org/download.html. I recommend installing HMMER 2.3.2 first (note that the scripts in this repository require the HMMER 2.3.2 binaries to be renamed with the "232" suffix as follows):

```bash
sh install_hmmer.sh
```


## Downloading and formatting Pfam HMMs
The most current set of Pfam-A HMMs can be found at: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz. To download this file then parse it into individual files to use with HMMER 2.3.2 and HMMER 3.1b2 (as the scripts in this repository require), run:

```bash
python download_pfam.py
```

## Running HMMER and parsing results
The general format to run is as follows, where the ending value is 1 more than the actual ending 0-index. For instance, to run on HMMs 0, 1, and 2, we would call processhmmer.py with --start 0 and --end **3**.

```bash
python processhmmer.py --start 0 --end 10
```

**NOTE:** You *will* need to edit this file with your hardcoded directories, protein input file (in FASTA format), and the version of Pfam you are running on.

## Combining HMMER output

We run both HMMER 2.3.2 and HMMER 3.0 and therefore expect lots of duplicate domains. We combine all output from individual Pfam HMMs with the following call, which will produce a file called allhmmresbyprot-v31.txt.gz. (You can change this name in the script.)

```bash
python createdomainoutput.py --function concatenate_hmmer_results
```

We **must** run this function first before calling the following to restrict to domains that:

* are complete (i.e., matched from the very start to the very end of the HMM)
* passed the gathering threshold (taking into account both domain- and sequence-based cutoffs)
* have the appropriate residue at high information content positions (to remove "deprecated" domains)

```bash
python createdomainoutput.py --function filter_domains
```

This produces a single, beautiful file (called domsbyprot.txt.gz) with all the domains from your original FASTA file.