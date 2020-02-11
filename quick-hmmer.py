#!/usr/bin/python
__author__ =  "Javed M. Aman"
#%%
import sys
sys.path.append('../ppi-dome/util')
import process_hmmer
import create_domain_output
import fasta_util
import argparse
import os
#import pygtrie as trie

#%%
def chain_domain_dict(dom_dict_file) :
    domain_dict = {}
    with open(dom_dict_file, 'rt') as csvfile :
        for line in csvfile :
            fields = line.strip().split(',')
            pdb_id = fields[0]
            pfam_ids = (fields[1].split('#'))[1:]
            domain_dict[pdb_id] = pfam_ids
    return domain_dict
#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quickly execute run-hmmer")
    script_path = os.path.dirname(os.path.abspath(__file__)) + '/'
    parser.add_argument('--fasta_infile', type=str, help='Full path to fasta-formatted sequence file to run HMMER on.')
    parser.add_argument('--dom_dict_infile', type=str, help='Full path to dictionary csv with ')
    parser.add_argument('--pfam_path', type=str, default=script_path + 'pfam/',
                        help='Full path to a directory where Pfam HMMs are stored')
    parser.add_argument('--results_path', type=str, help='directory to store domain results')
    #parser.add_argument('--pfam_version', type=str, default='32', choices=[str(n) for n in range(28, 33)],
    #                    help='Pfam version we are running on.')

    args = parser.parse_args()
    domain_dict = chain_domain_dict(args.dom_dict_infile)

    if not args.pfam_path.endswith('/'):
        args.pfam_path = args.pfam_path + '/'
    if not args.results_path.endswith('/'):
        args.results_path = args.results_path + '/'
    hmm_dir = args.pfam_path

    (_, _, hmm_files) = next(os.walk(hmm_dir))

    hmm_file_dict = {}
    for hmm_file in hmm_files : 
        hmm_id = hmm_file.split('_')[0]
        hmm_file_dict[hmm_id] = hmm_file
    
    print(hmm_file_dict['PF00595'])
    fasta = fasta_util.read_fasta(args.fasta_infile)
    for (comment, sequence) in fasta :
        (chain, _, _, _, _) = fasta_util.split_comment_chain(comment)[0]
        if not chain in domain_dict:
            print(chain + " is not in pdb dictionary of domains")
            continue
        for hmm in domain_dict[chain] :
            outfile = args.results_path+hmm+'_'+chain+'.hmmres.gz'
            process_hmmer.find_domain_matches(hmm_dir+'/'+hmm_file_dict[hmm], outfile, args.fasta_infile, 'error.log', [chain])


    


# %%
