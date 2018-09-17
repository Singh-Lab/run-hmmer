# Running HMMER
Run both HMMER 2.3.2 and HMMER 3 on a set of Pfam HMMs, parse and combine results.

## Downloading and installing HMMER
HMMER 2.3.2 (released October 2003) and HMMER 3.1b2 (released February 2015) can both be downloaded from http://hmmer.org/download.html. I recommend installing HMMER 2.3.2 first (note that the scripts in this repository require the HMMER 2.3.2 binaries to be renamed with a "232" suffix):

```bash
sh install_hmmer.sh
```


## Downloading and formatting Pfam HMMs
The most current set of Pfam-A HMMs can be found at: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz. To download this file then parse it into individual files to use with HMMER 2.3.2 and HMMER 3.1b2 (as the scripts in this repository require), run the following. Note that the default **path_to_pfam** is "current_working_directory/pfam/".

```bash
python download_pfam.py --path_to_pfam pfam/
```

## Running HMMER and parsing results
The general format to run is as follows, where the ending value is 1 more than the actual ending 0-index. For instance, to run on HMMs 0, 1, and 2, we would call processhmmer.py with --start 0 and --end **3**.

```bash
python process_hmmer.py --path_to_pfam pfam/ 
                        --pfam_version 31 
                        --fasta_infile <full path to NON-GZIPPED fasta-formatted sequence file> 
                        --path_to_results domains/ 
                        --start 0 
                        --end 10
```

**NOTE:** This step should be run in parallel on the cluster!

## Combining HMMER output

We run both HMMER 2.3.2 and HMMER 3.1b2 and therefore expect lots of duplicate domains. We combine all output across all Pfam HMMs with the following call, which will produce a file called allhmmresbyprot-v31.txt.gz. (You can change this name, and perhaps also what subset of HMM results you want to include, in the script.)

```bash
python create_domain_output.py --concatenate_hmmer_results 
                               --path_to_pfam pfam/
                               --pfam_version 31
                               --fasta_infile <full path to NON-GZIPPED fasta-formatted sequence file>
                               --path_to_results domains/
```

We **must** run the previous function before calling the following, as the following depends on the intermediate results. We run this to restrict to domains that:

* are complete (i.e., matched from the very start to the very end of the HMM)
* passed the gathering threshold (taking into account both domain- and sequence-based cutoffs)
* have the appropriate residue at high information content positions (to remove "deprecated" domains)

```bash
python create_domain_output.py --filter_domains
                               --path_to_pfam pfam/
                               --pfam_version 31
                               --fasta_infile <full path to NON-GZIPPED fasta-formatted sequence file>
                               --path_to_results domains/
```

This produces a single, beautiful file (called domsbyprot.txt.gz) with all the domains from your original FASTA file.
