#!/usr/bin/python

"""
Download latest version of Pfam and set up standard directories to run processhmmer.py
Contact snadimpa@princeton.edu with questions.
"""

import os
import gzip
import sys
from subprocess import call

####################################################################################################

def pfam_release_info(current_release_loc='ftp://ftp.ebi.ac.uk/pub/databases/Pfam/' +
                                          'current_release/relnotes.txt'):
  """
  :param current_release_loc: complete URL to Pfam's current release notes (.txt format)
  :return: current Pfam release (e.g., 31.0), month and year of release (e.g., March 2017), and the total
           number of entries (e.g., 16,712)
  """

  # Download the current release notes and save to a temporary directory
  os.system('wget ' + current_release_loc + ' -O pfam-current-release.txt')

  with open('pfam-current-release.txt') as x:
    x.next()
    release = x.next().strip().split()[1]  # RELEASE information always on the second line

    for l in x:
      if l.strip().startswith(release):
        date = l.strip().split()[1]
        num_entries = l.strip().split()[2]
        break

  # Remove the temporary file:
  os.system('rm pfam-current-release.txt')

  month_abbreviations = {'01': 'January', '02': 'February', '03': 'March', '04': 'April', '05': 'May',
                         '06': 'June', '07': 'July', '08': 'August', '09': 'September', '10': 'October',
                         '11': 'November', '12': 'December'}

  month = month_abbreviations[date.split('/')[0]]
  year = ('20' if int(date.split('/')[-1]) < 40 else '19') + date.split('/')[-1]

  return release, month + ' ' + year, "{:,}".format(int(num_entries))


####################################################################################################

def pfam_hmm_download(output_directory='hmms-v31/', pfam_version='31',
                      hmms_link='ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'):
  """
  :param output_directory: full path to a directory where all separate HMM files will be written to
  :param hmms_link: complete URL to Pfam's current set of "A" HMMs
  :return: None. Download the Pfam-A HMM file, parse each individual HMM into a separate file.
  """

  # Create output directory if need be:
  if not os.path.exists(output_directory):
    os.makedirs(output_directory)

  # Download the file containing all Pfam-A HMMs of interest:
  temporary_hmm_file = 'Pfam-A.current-v'+pfam_version+'-release.hmm' + ('.gz' if hmms_link.endswith('gz') else '.txt')
  os.system('wget ' + hmms_link + ' -O ' + temporary_hmm_file)

  # Keep track of current HMM accession and name (to print out at the appropriate time)
  current_name, current_accession = '', ''

  # Keep track of all lines belonging to this HMM:
  current_lines = []

  hmm_infile = gzip.open(temporary_hmm_file) if temporary_hmm_file.endswith('gz') else open(temporary_hmm_file)

  for l in hmm_infile:
    current_lines.append(l)

    if l.strip() == '//':
      if len(current_lines) > 0:
        hmm_outfile = open(output_directory + current_accession.split('.')[0] + '_' + current_name + '.hmm', 'w')
        map(hmm_outfile.write, current_lines)
        hmm_outfile.close()
      current_name, current_accession = '', ''
      current_lines = []

    else:
      if l.startswith('NAME '):
        current_name = l.strip().split()[-1]
      if l.startswith('ACC '):
        current_accession = l.strip().split()[-1]
  hmm_infile.close()

  os.system('rm ' + temporary_hmm_file)

  sys.stderr.write('Finished processing all HMMs from ' + hmms_link + ' into ' + output_directory + '\n')


####################################################################################################

if __name__ == "__main__":

  # Get information about the most recent version of Pfam:
  version, date, entries = pfam_release_info()
  sys.stderr.write('Pfam version ' + version + ', released ' + date + ', contains ' + entries + ' total HMMs.\n')

  # Set up all directories here as needed:
  DATAPATH = os.getcwd() + '/'
  pfam_version = version[:version.rfind('.')]

  for directory in ['pfam',
                    'pfam/hmms-v'+pfam_version,
                    'domains',
                    'domains/hmmres-v'+pfam_version,
                    'domains/processed-v'+pfam_version]:
    call(['mkdir', DATAPATH+directory])

  # Download the most recent version of Pfam, storing HMMs in the proper directory:
  pfam_hmm_download(DATAPATH+'pfam/hmms-v'+pfam_version+'/', pfam_version)
