#!/usr/bin/python

"""
Given formatted by-HMM HMMER hits, concatenate all results to remove duplications and order results, 
then also create a file listing ALL complete domains that passed their gathering thresholds.

Contact snadimpa@princeton.edu with questions.
"""

import os
import sys
import gzip
import math
import argparse

####################################################################################################
# CONSTANTS -- UPDATE THESE!!!
####################################################################################################

# Path to where data is stored
DATAPATH = os.getcwd() + '/'

PFAMVERSION = '31'  # Version of Pfam we are using

# Full path to where HMMs have been saved. We assume files are named as
#  PfamID_PfamName.hmm (e.g., PF00096_zf-C2H2.hmm)
HMMPATH = DATAPATH + 'pfam/hmms-v' + PFAMVERSION + '/'

# Full paths to the FASTA file of protein sequences that we supposedly ran on:
#  (for the sake of filling in the header properly)
PROTFILE = DATAPATH + 'human_test_sequences.fa'

# Output processed from the processhmmer.py script is here, for each domain (e.g., v31/PF00096_zf-C2H2-v31.hmmres.gz)
PROCESSED_HMMER_RESULTS = DATAPATH + 'domains/processed-v' + PFAMVERSION + '/'

# Release dates and number of entries can always be found for all Pfam releases at
# ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt
HMMINFO = {'27': {'date': 'March 2013', 'entries': '14,831'},
           '28': {'date': 'May 2015', 'entries': '16,230'},
           '29': {'date': 'December 2015', 'entries': '16,295'},
           '30': {'date': 'June 2016', 'entries': '16,306'},
           '31': {'date': 'March 2017', 'entries': '16,712'}}


####################################################################################################

def create_allhmmresbyprot(pfamversion=PFAMVERSION,
                           inputdir=PROCESSED_HMMER_RESULTS,
                           outfile=DATAPATH + 'domains/allhmmresbyprot-v31.tsv.gz'):
  """
  :param pfamversion: Which Pfam version to use? 31 is our standard
  :param inputdir: location of all .hmmres.gz files generated from the processhmmer.py script
  :param outfile: final file of all non-redundant domain hits
  :return: concatenate all HMMER results into a single file, removing duplicates as best we can
           between the two versions of HMMER
  """
  output_handle = gzip.open(outfile, 'w') if outfile.endswith('gz') else open(outfile, 'w')
  output_handle.write('# HMMER 2.3.2 and HMMER 3.1b2 results on all protein sequences found in '+PROTFILE+'\n')
  output_handle.write('# Pfam version '+pfamversion+'.0, released '+HMMINFO[pfamversion]['date']+
                      ', containing '+HMMINFO[pfamversion]['entries']+' entries\n')

  output_handle.write('\t'.join(['#TargetID', 'HMM_Name', 'E-value', 'BitScore', 'TargetStart', 'TargetEnd', 'HMM_Seq',
                                 'Target_Seq', 'HMM_Pos', 'Description'])+'\n')

  totaldomainhits = 0  # total number of domains processed

  for fname in sorted([a for a in os.listdir(inputdir)
                       if a.startswith('PF') and a.endswith('-v' + pfamversion + '.hmmres.gz')]):

    # (targetID, start, end, pfamHMM, desc) -> (bitscore, evalue, HMM sequence, target sequence, HMM match states)
    results = {}

    for l in gzip.open(inputdir+fname):
      if not l.startswith('#') and len(l.strip()) > 0:
        targ_id, hmm_id, evalue, bitscore, tstart, tend, hmmseq, tseq, hmmpos, desc = l[:-1].split('\t')

        # This is a new result; add it:
        if (targ_id, int(tstart), int(tend), hmm_id, desc) not in results:
          results[(targ_id, int(tstart), int(tend), hmm_id, desc)] = [float(bitscore), evalue, hmmseq, tseq, hmmpos]

        # This is a DUPLICATE result; update if it is better:
        if float(bitscore) > results[(targ_id, int(tstart), int(tend), hmm_id, desc)][0]:
          results[(targ_id, int(tstart), int(tend), hmm_id, desc)] = [float(bitscore), evalue, hmmseq, tseq, hmmpos]

    # Write out all sorted results:
    for ((targ_id, tstart, tend, hmm_id, desc), (bitscore, evalue, hmmseq, tseq, hmmpos)) in sorted(results.items()):
      output_handle.write('\t'.join([targ_id, hmm_id, evalue, str(bitscore), str(tstart), str(tend),
                                     hmmseq, tseq, hmmpos, desc])+'\n')

    totaldomainhits += len(results.keys())  # keep track of grand total across all domain types
  output_handle.close()

  sys.stderr.write('Concatenated all files from '+inputdir+' to '+outfile+'\n')
  sys.stderr.write(str(totaldomainhits)+' total domains!\n')


####################################################################################################

def find_domains_from_file(hmmresfile=PROCESSED_HMMER_RESULTS+'allhmmresbyprot-v31.tsv.gz'):
  """
  :param hmmresfile: full path to a formatted HMMER result file
  :return: set of all domains found in that HMMER result file
  """

  x = gzip.open(hmmresfile) if hmmresfile.endswith('.gz') else open(hmmresfile)
  allhmms = set([a.strip().split('\t')[1] for a in x if not a.startswith('#')])
  x.close()

  return allhmms


####################################################################################################

def get_high_information_content(hmmfile):
  """
  :param hmmfile: full path to a Pfam HMM file
  :return: a dictionary of 1-indexed match state -> required amino acid assignment if the 
     information content at that match state is >=4 (corresponding to ~95% the same amino acid at
     that position)
  """

  requiredstates = {}

  reachhmm = False
  with open(hmmfile) as x:
    for l in x:
      if l.startswith('HMM') and not l.startswith('HMMER'):
        aas = l.strip().split()[1:]  # amino acids in order
        x.next()  # description of transition probabilities
        x.next()  # begin state match state emission probabilities (unnecessary)
        x.next()  # begin state insertion state emission probabilities (also unnecessary)
        x.next()  # begin state transition probabilities
        reachhmm = True  # Have we reached the actual HMM model yet?
        continue
      elif reachhmm and len(l.strip().split()) > len(aas):
        matchstate = l.strip().split()[0]
        probabilities = map(lambda var: math.exp(-1*float(var)), l.strip().split()[1:len(aas)+1])
        if math.log(20, 2) + sum(j*math.log(j, 2) for j in probabilities) >= 4:  # required amino acid!
          requiredstates[matchstate] = sorted(zip(probabilities, aas), reverse=True)[0][1]
  return requiredstates


####################################################################################################

def return_passing_hits(currenthits, gatheringthresholds):
  """
  In the case of repetitive domains, multiple domain hits is evidence of biological importance,
  so domains are included even if individually they did not pass the gathering threshold:
  :param currenthits: hmmID -> start_position -> (tuple(hmm_match_states), bitscore)
  :param gatheringthresholds: hmmID -> (sequence gathering threshold, domain gathering threshold)
  :return: the set of HMM matches that passed the gathering threshold
  """

  passing_hits = []

  for hmm in currenthits.keys():
    allpass = False  # all the domains pass if the sum of bit scores meets the threshold

    if sum([bitscore for _, bitscore in currenthits[hmm].values()]) >= gatheringthresholds[hmm][0]:
      allpass = True

    for startpos in sorted(currenthits[hmm].keys()):

      matchstates, bitscore = currenthits[hmm][startpos]

      if allpass or bitscore >= gatheringthresholds[hmm][1]:
        passing_hits.append((hmm, matchstates))

  return passing_hits


####################################################################################################

def create_domsbyprot(pfamversion=PFAMVERSION,
                      hmmresfile=DATAPATH + 'domains/allhmmresbyprot-v31.tsv.gz',
                      outfile=DATAPATH + 'domains/domsbyprot-v31.txt.gz'):
  """
  :param pfamversion: version of the Pfam database we are using (default is 31)
  :param hmmresfile: full path to a file generated by create_allhmmresbyprot()
  :param outfile: map the HMM match states to the sequence indices and residues at those indices
                  across all proteins for which we have results; note that these results are
                  guaranteed to be complete, pass default gathering thresholds, and have the 
                  required amino acid at match states with high information content
  :return: 
  """

  path_to_hmms = HMMPATH

  # Find all HMMs to consider:
  allhmms = find_domains_from_file(hmmresfile)

  # Determine the set of "required" states for all HMMs, where relevant
  requiredstates = {}  # hmm_id -> matchstate -> amino acid required
  for hmm in allhmms:
    currentrstates = get_high_information_content(path_to_hmms + hmm + '.hmm')
    if len(currentrstates.keys()) > 0:
      requiredstates[hmm] = currentrstates

  # Determine the lengths and instance/sequence gathering threshold cutoffs:
  hmmlens = {}
  gacutoff = {}
  for hmm in allhmms:
    for l in open(path_to_hmms + hmm + '.hmm'):
      if l.startswith('LENG'):
        hmmlens[hmm] = l[:-1].split()[1]
      if l.startswith('GA'):
        gacutoff[hmm] = (float(l[:-1].split()[1].replace(';', '')),
                         float(l[:-1].split()[2].replace(';', '')))
        break

  domsbyprot = {}  # prot_id -> list of domains (hmm_id, set([(matchstate,index),...])) in that protein
  current_protid = ''  # current protein being processed
  current_sum = {}  # pfamID -> sum

  # Process all domains!
  x = gzip.open(hmmresfile) if hmmresfile.endswith('.gz') else open(hmmresfile)
  for l in x:
    if l.startswith('#') or len(l.strip().split('\t')) < 10:
      continue

    prot_id, hmm_id, evalue, bit_score, targ_start, targ_end, hmm_seq, targ_seq, hmm_pos, desc = l[:-1].split('\t')[:10]

    # Save all relevant domain hits for the previous protein:
    if prot_id != current_protid:
      passing_hits = return_passing_hits(current_sum, gacutoff)

      if len(passing_hits) > 0:
        if current_protid not in domsbyprot:
          domsbyprot[current_protid] = []
        for hit in passing_hits:
          domsbyprot[current_protid].append(hit)

      current_protid = prot_id
      current_sum = {}

    # (1) Include only COMPLETE domains:
    if not hmm_pos.startswith('1,') or not hmm_pos.endswith(',' + hmmlens[hmm_id]):
      continue

    # (2) Make sure that the first and last positions are ungapped:
    mstatetoseq = zip(hmm_pos.split(','), list(targ_seq))
    if mstatetoseq[-1][1] == '-' or mstatetoseq[0][1] == '-':
      continue

    # (3) High information-content sites must have the appropriate residue assignment:
    if hmm_id in requiredstates:
      badmatch = False
      for hmmState, seqAA in mstatetoseq:
        if hmmState in requiredstates[hmm_id] and seqAA != requiredstates[hmm_id][hmmState]:
          #sys.stderr.write('In HMM ' + hmm_id + ', we expected ' + requiredstates[hmm_id][hmmState] +
          #                 ' in state ' + hmmState + ' but saw ' + seqAA + '\n')
          badmatch = True
          break
      if badmatch:
        continue

    # Map HMM match state -> the HMM sequence expected -> the actual sequence hit
    seq = range(int(targ_start) - 1, int(targ_end))  # HMM sequence
    i = 0  # keep track of ungapped matches
    currhmmstates = []
    for hmmState, seqAA in mstatetoseq:
      if seqAA != '-':
        currhmmstates.append((hmmState, seq[i], seqAA.upper()))  # THIS HAS BEEN CHECKED!! REAL (0) INDICES!
        i += 1

    # Store this domain hit (along with match state mapping) to potentially include later on.
    if hmm_id not in current_sum:
      current_sum[hmm_id] = {}
    startpos = min([seqindex for _, seqindex, _ in currhmmstates])
    if startpos not in current_sum[hmm_id] or float(bit_score) > current_sum[hmm_id][startpos][1]:
      current_sum[hmm_id][startpos] = (tuple(currhmmstates), float(bit_score))

  x.close()

  # Store information for the very last protein:
  if len(current_sum.keys()) > 0:
    passing_hits = return_passing_hits(current_sum, gacutoff)

    if len(passing_hits) > 0:
      if current_protid not in domsbyprot:
        domsbyprot[current_protid] = []
      for hit in passing_hits:
        domsbyprot[current_protid].append(hit)

  # Write out all sorted results:
  output_handle = gzip.open(outfile, 'w') if outfile.endswith('gz') else open(outfile, 'w')
  output_handle.write('# All COMPLETE Pfam domains (version '+pfamversion+'.0, ' +
                      HMMINFO[pfamversion]['date']+', '+HMMINFO[pfamversion]['entries']+' entries) that ' +
                      'passed the gathering threshold found in all amino acid sequences from\n')
  output_handle.write('# '+PROTFILE+'\n')
  output_handle.write('# Original, unfiltered by-domain HMMER results found in ' +
                      PROCESSED_HMMER_RESULTS+' directory\n')
  output_handle.write('\t'.join(['#Protein_Sequence_ID', 'Pfam_HMM_ID', 'matchstate:AA-0-index:AA-value']) + '\n')

  final_domain_count = 0
  for pID in sorted(domsbyprot.keys()):
    for m in sorted(domsbyprot[pID]):
      output_handle.write(pID + '\t' + m[0] + '\t' +
                          ','.join([str(a[0]) + ':' + str(a[1]) + ':' + str(a[2]) for a in m[1]]) + '\n')
      final_domain_count += 1

  output_handle.close()
  sys.stderr.write('Condensed ' + hmmresfile + ' into ' + outfile + '\n')
  sys.stderr.write(str(final_domain_count)+' total (1) complete, (2) non-deprecated domains that ' +
                   '(3) passed the gathering threshold!\n')


####################################################################################################

if __name__ == "__main__":

  """
  PROCESS COMMAND-LINE ARGUMENTS
  """

  # Option of setting the Pfam version from the command line if you want to change defaults:
  parser = argparse.ArgumentParser(description='Concatenate and filter all processed HMMER domain hit results.')

  parser.add_argument('--pfam_version', type=int,
                      help='Pfam database version (e.g., 31)',
                      default=PFAMVERSION,
                      choices=range(27, 32))
  parser.add_argument('--function', type=str,
                      help='Processing step to run (e.g., concatenate_hmmer_results, filter_domains)',
                      default='concatenate_hmmer_results',
                      choices=['concatenate_hmmer_results', 'filter_domains'])

  args = parser.parse_args()

  PFAMVERSION = str(args.pfam_version)

  HMMER_RESULTS = PROCESSED_HMMER_RESULTS
  HMMER_RESULT_FILE = DATAPATH + 'domains/allhmmresbyprot-v' + PFAMVERSION + '.tsv.gz'
  DOMAIN_RESULT_FILE = DATAPATH + 'domains/domsbyprot-v' + PFAMVERSION + '.txt.gz'

  if args.function == 'concatenate_hmmer_results':
    # Concatenate all results from multiple files, removing duplicates as well as we can
    create_allhmmresbyprot(PFAMVERSION, HMMER_RESULTS, HMMER_RESULT_FILE)
    sys.stderr.write('For final, filtered output, remember to run:\n' + 
                     'python '+sys.argv[1]+' --function filter_domains\n')

  if args.function == 'filter_domains':
    # Restrict to domains that:
    # (1) are complete (i.e., matched from the very start to the very end of the HMM)
    # (2) passed the gathering threshold (taking into account both domain- and sequence-based cutoffs)
    # (3) have the appropriate residue at high information content positions (to remove "deprecated" domains)
    if not os.path.isfile(HMMER_RESULT_FILE):
      sys.stderr.write('No such file: '+HMMER_RESULT_FILE+'!\n' + 
                       'Please run: python '+sys.argv[1]+' --function concatenate_hmmer_results\n')
           
    create_domsbyprot(PFAMVERSION, HMMER_RESULT_FILE, DOMAIN_RESULT_FILE)
