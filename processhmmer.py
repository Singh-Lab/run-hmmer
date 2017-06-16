#!/usr/bin/python

"""
Call relevant functions from the hmmer.py script in order to find matches to 
particular Pfam domains (and process the output from HMMER).

Contact snadimpa@princeton.edu with questions.
"""

import os
import sys
import gzip
import argparse
import time
from string import ascii_letters
from string import digits
from random import choice
from Bio import pairwise2
from hmmer import finddom

####################################################################################################
# CONSTANTS -- UPDATE THESE!!!
####################################################################################################

# Path to where data is stored
DATAPATH = os.getcwd() + '/'

PFAMVERSION = '31'  # Version of Pfam we are using

# Full path to where HMMs have been saved. We assume files are named as
#  PfamID_PfamName.hmm (e.g., PF00096_zf-C2H2.hmm)
HMMPATH = DATAPATH + 'pfam/hmms-v' + PFAMVERSION + '/'

# Full paths to the FASTA file of protein sequences that we want to search
#  for domains in. Note that hmmsearch causes a fuss if this files are zipped
#  at all.. (so make sure it is not)
PROTFILE = DATAPATH + 'human_test_sequences.fa'

# Unprocessed HMMER results will go here (good idea not to delete, in case of a crash or debugging)
TEMPORARY_HMMER_RESULTS = DATAPATH + 'domains/hmmres-v' + PFAMVERSION + '/'

# Final, processed output will go here, for each domain (e.g., PF00096_zf-C2H2-v31.hmmres.gz)
PROCESSED_HMMER_RESULTS = DATAPATH + 'domains/processed-v' + PFAMVERSION + '/'

# In case we had a match that we couldn't understand nor rescue, keep track of it:
LOGFILE = PROCESSED_HMMER_RESULTS + 'problems.log'


####################################################################################################

def idgen(size=20, chars=ascii_letters + digits):
  """
  :param size: length of string to generate
  :param chars: characters to randomly draw from in resulting string
  :return: string of length size containing random characters drawn from chars with replacement 
           (intended use as a temporary file name); lowercase & uppercase letters and digits 0-9
  """

  return ''.join(choice(chars) for _ in xrange(size))


####################################################################################################

def process_hmmer_output(hmmresfile, outputfile, pfam):
  """
  :param hmmresfile: full path to the results file from running finddom from hmmer.py
  :param outputfile: full path to the output file containing a list of complete matches
  :param pfam: dictionary of Pfam HMM ID -> the representative sequence for that HMM (obtained from the 
               raw .hmm files)
  :return: None. Output formatted results to the specified outfile. This entire process is necessary
           to determine WHICH match states aligned to which sequence positions (older version of HMMER
           makes this non-trivial)
  """

  output_handle = gzip.open(outputfile, 'w') if outputfile.endswith('gz') else open(outputfile, 'w')

  input_handle = gzip.open(hmmresfile) if hmmresfile.endswith('gz') else open(hmmresfile)

  # Skip the infile header and write out a new header:
  input_handle.next()
  output_handle.write('\t'.join(['#Target', 'HMM', 'E-value', 'BitScore', 'TargetStart', 'TargetEnd', 'HMM_Seq',
                                 'Ali_Seq', 'Pos_Seq', 'Description'])+'\n')

  for current_line in input_handle:
    # Tab-delineated file generated from the "finddom" function in hmmer.py:
    prot_id, hmm_id, evalue, bitscore, tstart, tend, hmmmatch, alimatch, desc = current_line[:-1].split('\t')[:9]
    
    if hmm_id not in pfam:
      sys.stderr.write('Error: Could not find '+hmm_id+' in pfam dictionary!\n')
      continue

    # Reformat the hmmmatch sequence off the bat; usually the "+" means the match was positive at
    #  that position in the HMM alignment (says nothing about the complete match). Full description
    #  of HMMER output can be found in the User's Guide:
    # http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf
    hmmmatch = hmmmatch.replace('+', '')

    # Sometimes only part of the domain matched. Figure out the starting index of the match:
    startingindex = pfam[hmm_id][0].find(hmmmatch.replace('.', '').replace('+', ''))

    # If the resulting HMM match doesn't match what we expect for the overall Pfam sequence match, then
    #  we will not be able to find the starting index. In this case, we need to "translate" our HMM match
    #  into what we expected to see via a quick alignment:
    if startingindex < 0 or len(hmmmatch) != len(alimatch):
        
      # Let's align our hmmmatch to what we SHOULD have, the Pfam sequence.
      # A mismatch is only -0.5 whereas gaps are -2, so mismatch is preferred.
      newh, alim = pairwise2.align.globalms(hmmmatch.replace('.', ''), pfam[hmm_id][0], 2, -0.5, -1.5, -1)[0][:2]
        
      newhmmmatch = []
      i, j, k = 0, 0, 0  # indices into original HMM match, aligned HMM match, and aligned expected Pfam sequence

      # ....Simple quick alignment adjustment:
      while i < len(hmmmatch) or j < len(newh) or k < len(alim):  # alim must ALWAYS be >= in length.
        if newh[j] == '-' and alim[k] != '-' and i >= len(hmmmatch):  # gap at the END, don't add anything, move on.
          j += 1
          k += 1
        elif hmmmatch[i] == '.':
          newhmmmatch.append('.')
          i += 1
        elif newh[j] == '-' and alim[k] != '-' and i == 0:  # gap at the BEGINNING, don't add anything, just move on.
          j += 1
          k += 1
        elif newh[j] == '-' and alim[k] != '-' and i < len(hmmmatch):  # ADD THIS CHARACTER
          newhmmmatch.append(alim[k])
          j += 1
          k += 1
        elif newh[j] != '-' and newh[j] == hmmmatch[i] and alim[k] == '-':  # DELETE THIS CHARACTER
          i += 1
          j += 1
          k += 1
        elif newh[j] != '-' and newh[j] == hmmmatch[i] and alim[k] != '-':  # REPLACE this character
          newhmmmatch.append(alim[k])
          i += 1
          j += 1
          k += 1
        elif newh[j] == alim[k] == hmmmatch[i]:
          newhmmmatch.append(hmmmatch[i])
          i += 1
          j += 1
          k += 1
        else:  # No idea what this case is! Let's check it out...
          sys.stderr.write(str(i) + ':' + hmmmatch[i] + '\t' + hmmmatch + '\n')  # Original HMM match
          sys.stderr.write(str(j) + ':' + newh[j] + '\t' + newh + '\n')  # Aligned HMM match
          sys.stderr.write(str(k) + ':' + alim[k] + '\t' + alim + '\n')  # Aligned Pfam sequence match

      # Did this fix the problem!?
      newhmmmatch = ''.join(newhmmmatch)
      startingindex = pfam[hmm_id][0].find(newhmmmatch.replace('.', '').replace('+', ''))

      if startingindex < 0 or len(newhmmmatch) != len(alimatch):
        logfile_handle = open(LOGFILE, 'a')
        logfile_handle.write(hmm_id+'\n')
        logfile_handle.write('Orig HMM Match:\t'+hmmmatch+'\n')
        logfile_handle.write('Align HMM Matc:\t'+newh+'\n')
        logfile_handle.write('Align Pfam Mat:\t'+alim+'\n')
        logfile_handle.write('New HMM Match: \t'+newhmmmatch+'\n')
        logfile_handle.write('Orig Ali Match:\t'+alimatch+'\n\n')
        logfile_handle.close()
        sys.stderr.write('Problem in '+hmm_id+'\n')
        continue
        
      else:
        hmmmatch = newhmmmatch

    # Actual mapping of Pfam indices:
    mapindices = {i + 1: int(k) for i, k in enumerate(pfam[hmm_id][1].split(','))}
            
    startingindex += 1  # Results should be 1-indexed, not 0-indexed, for consistency with HMMER.

    findex = []
    currgap = True
    currgapi = 0

    # Determine the actual match state -> sequence mapping. Insertion states are prefixed with an "a"
    for i in xrange(len(hmmmatch)):
      if hmmmatch[i] != '.':
        currgap = False
        findex.append(str(mapindices[startingindex]))
        startingindex += 1
      elif hmmmatch[i] == '.':
        if currgap:
          currgapi += 1
        else:
          currgap = True
          currgapi = 0
        findex.append('a'+str(mapindices[startingindex]-1)+'-'+str(currgapi))

    output_handle.write('\t'.join([prot_id, hmm_id, evalue, bitscore, tstart, tend,
                                   hmmmatch, alimatch, ','.join(findex), desc])+'\n')

  input_handle.close()
  output_handle.close()


####################################################################################################

def get_pfamseq(hmm_file):
  """
  :param hmm_file: full path to a Pfam HMM file
  :return: consensus sequence representing the HMM, needed to parse hmmsearch results (v2.3.2) to
           get the match state -> sequence index -> amino acid mapping
  """

  pfamseq = []
  for current_line in open(hmm_file):
    i = current_line.strip().split()
    if len(i) == 26 and i[21].isdigit():
      pfamseq.append(i[22])
  pfamseq = ''.join(pfamseq)

  return pfamseq


####################################################################################################

def find_domain_matches(hmm_file, output_file, infile, ids=set()):
  """
  :param hmm_file: full path to an HMM. Remember our naming convention MUST be PATH + PfamID_PfamName.hmm
  :param output_file: full path to an output file to store the output
  :param infile: full path to a FASTA file with protein sequences that we want to find domains in
  :param ids: a subset of sequence identifiers that we care about from infile
  :return: find all sequence matches to the HMM profile in hmmfile for the specified sequence IDs found in 
           the input fasta file. Write the *complete* (not necessarily all confident), easily-parsed results
           to the output file. 
  """

  # Pfam "consensus" sequence is needed to parse hmmsearch results; obtained from the raw .hmm file:
  pfamseq = get_pfamseq(hmm_file)

  chmm = hmm_file.split('/')[-1].replace('.hmm', '')  # PfamID_PfamName

  output_handle = gzip.open(output_file, 'w') if output_file.endswith('gz') else open(output_file, 'w')
  output_handle.write('# All matching hits from the HMM found in ' + hmm_file + '\n')
  output_handle.write('# on the protein sequences found in ' + PROTFILE + '\n')
  output_handle.write('\t'.join(['#TargetID', 'HMM_Name', 'E-value', 'BitScore', 'TargetStart', 'TargetEnd', 'HMM_Seq',
                                 'Target_Seq', 'HMM_Pos', 'Description'])+'\n')

  # NOW, let's search, seq by seq...  
  input_handle = gzip.open(infile) if infile.endswith('gz') else open(infile)
  for current_line in input_handle:
    # Check if this is a sequence we actually want to run on:
    if current_line.startswith('>'):
      sequence_id = current_line[1:-1].split()[0]
      if len(ids) > 0 and sequence_id not in ids:
        continue

      tmpfiles = ['/tmp/'+idgen()+'.seq',      # sequence input file (single sequence)
                  '/tmp/'+idgen()+'.hmmres',   # original unparsed output
                  '/tmp/'+idgen()+'.hmmres1']  # complete parsed output
        
      # Write out current header and following sequence as input for hmmsearch
      # NOTE: this assumes that each sequence in the FASTA file is on one continuous line (without breaks)
      sequence_handle = open(tmpfiles[0], 'w')
      sequence_handle.write(current_line + input_handle.next())
      sequence_handle.close()
  
      # Now, RUN hmmsearch for EACH hmmID in our list.
      finddom([hmm_file], tmpfiles[0], tmpfiles[1])

      os.system('rm '+tmpfiles[0])
    
      # CHECK IF THIS FAILED!!!
      if sum(1 for _ in open(tmpfiles[1])) < 2:
        sys.stderr.write('Failed to find any HMM matches for ' + hmm_file + ' in ' + sequence_id + '.\n')
        os.system('rm '+tmpfiles[1])
        continue

      # we need the "pfam sequence" for our particular HMM so that we can match positions to match states.
      process_hmmer_output(tmpfiles[1], tmpfiles[2],
                           {chmm[chmm.find('_')+1:]: (pfamseq, ','.join(map(str, range(1, len(pfamseq) + 1))))})

      os.system('rm '+tmpfiles[1])
        
      with open(tmpfiles[2]) as y:
        y.next()
        for l2 in y: 
          # Write out the exact same line, but update the HMM name to include the Pfam ID (e.g., PF00096) also!
          output_handle.write(l2[:-1].split('\t')[0] + '\t' + chmm + '\t' + '\t'.join(l2[:-1].split('\t')[2:]) + '\n')
      os.system('rm '+tmpfiles[2])

  input_handle.close()
  output_handle.close()


####################################################################################################

if __name__ == "__main__":

  # All possible HMMs that we can run on:
  hmms = sorted([a.replace('.hmm', '') for a in os.listdir(HMMPATH) if a.startswith('PF') and a.endswith('.hmm')])

  # In order to parallelize, we can specify a *subset* of these domains to run HMMER on.
  # We have the option of specifying via the command line:
  parser = argparse.ArgumentParser(description='Run, parse, and process results from HMMER 3.0 and HMMER 2.3.2.')

  parser.add_argument('--start', '-s', type=int,
                      help='Starting 0-index of subset of domains from ' + HMMPATH + ' to run on.',
                      default=0,
                      choices=range(len(hmms)))

  parser.add_argument('--end', '-e', type=int,
                      help='Ending 0-index of subset of domains from ' + HMMPATH + ' to run on.',
                      default=len(hmms),
                      choices=range(len(hmms)))

  parser.add_argument('--pfam_version', '-p', type=int,
                      help='Pfam version we are running on.',
                      default=31,
                      choices=range(28,32))

  args = parser.parse_args()

  # Customize log file name if need be:
  if not (args.start == 0 and args.end == len(hmms)):
    LOGFILE = LOGFILE.replace('.log', '-'+str(args.start)+'-'+str(args.end)+'.log')
    hmms = hmms[args.start:args.end]  # Subset to a range of HMMs if it was specified

  # Get set of all sequence IDs in input file:
  allseqs = set()
  x = gzip.open(PROTFILE) if PROTFILE.endswith('gz') else open(PROTFILE)
  for l in x:
    if l.startswith('>'):
      allseqs.add(l[1:-1].split()[0])
  x.close()

  # Clear the logfile to keep track of only the most recent results:
  x = gzip.open(LOGFILE, 'w') if LOGFILE.endswith('gz') else open(LOGFILE, 'w')
  x.close()

  for hmm in hmms:
    hmmfile = HMMPATH + hmm + '.hmm'

    # (1) check if HMM exists
    if not os.path.isfile(hmmfile):
      sys.stderr.write('\n'.join(['Could not find HMM file: ' + hmmfile,
                                  'Download latest version from Pfam:',
                                  'wget -O ' + hmmfile + ' http://pfam.xfam.org/family/' +
                                  hmm.split('_')[0] + '/hmm']) + '\n')
      continue

    # (2) run hmmsearch to subset the genes we want to actually find domain hits in
    #     NOTE: this is just an efficiency step, as this process is FAST but gives us less information
    os.system(' '.join(['hmmsearch',
                        '--domtblout ' + TEMPORARY_HMMER_RESULTS + hmm + '.hmmres-orig',
                        '-o /dev/null',
                        '-T 0 --domT 0 --incT 0 --incdomT 0',  # No cutoffs guarantees more thorough hits
                        hmmfile,
                        PROTFILE]))
    
    # (3) Set the subset of sequence IDs to look through (to speed up process of finding matches):
    whichseqs = set()
    if os.path.isfile(TEMPORARY_HMMER_RESULTS + hmm + '.hmmres-orig'):
      for l in open(TEMPORARY_HMMER_RESULTS + hmm + '.hmmres-orig'):
        if l.startswith('#'):
          continue
        whichseqs.add(l.strip().split()[0])

    if False and len(whichseqs) < 1:
      sys.stderr.write('Failed to get preliminary results for '+hmm+'! ' +
                       'Running on all ' + str(len(allseqs)) + ' sequences...\n')
      whichseqs = allseqs

    # (4) Find the full matches (match state -> sequence index -> sequence residue information):
    if len(whichseqs) > 0:
      sys.stderr.write('Processing '+hmm+'-v'+PFAMVERSION+'...')

      # start the clock to measure performance
      start = time.time()

      outfile = PROCESSED_HMMER_RESULTS + hmm + '-v' + PFAMVERSION + '.hmmres.gz'
      find_domain_matches(hmmfile, outfile, PROTFILE, whichseqs)

      # end the clock and print total elapsed time:
      end = time.time()
      total_seconds = end - start
      m, s = divmod(total_seconds, 60)
      h, m = divmod(m, 60)
      d, h = divmod(h, 24)
      sys.stderr.write('Finished '+str(len(whichseqs))+' sequences in '+
                       ':'.join(map(lambda x: str(int(x)).zfill(2), [d, h, m, s]))+'!\n')
