#!/usr/bin/env python3

import argparse
import io
import logging
import os
import subprocess
import sys
import tempfile
import vcf

from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path


def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True)
    parser.add_argument('--consensus', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--fail_vcf', required=False)
    parser.add_argument('--max', required=False, default=20, type=int, help='Maximum number of Ns before erroring out')
    parser.add_argument('--force', required=False, default=False, action='store_true')

    return parser


def parse_fasta(input_fasta):
    '''
    Parse the input fasta sequence that is input

    INPUTS:
        input_fasta --> `path` to fasta file from argparse --consensus and --reference input

    RETURNS:
        Uppercase string of the fasta sequence
        Record string of the fasta sequence name
    '''
    record = SeqIO.read(input_fasta, 'fasta')

    return str(record.seq).upper(), record.name


def correct_for_indels(ref_in, con_in, del_dict={}):
    '''
    If there is a deletion, determine the location(s) (through mafft) and return a dictionary containing the starting location
    and the length of the deletion

    INPUTS:
        ref_in      --> `string` of uppercase reference sequence
        con_in      --> `string` of uppercase consensus sequence
        del_dict    --> `dictionary` to be filled

    RETURNS:
        del_dict    --> `dictionary` populated with the starting genomic location of the deletion and the length of the deletion
            {6031: 6, 12033: 18}
    '''

    # List of sequences as Biopython Seq objects (needed!) to be used for mafft input
    sequence_list = [SeqIO.read(ref_in, 'fasta'), SeqIO.read(con_in, 'fasta')]

    # Use tempfile as we don't want to keep the combined fasta file
    fd, concat_temp = tempfile.mkstemp()

    try:
        with open(concat_temp, 'w') as handle:
            SeqIO.write(sequence_list, handle, 'fasta')

        # Run Mafft, adjustdirection sets it to map to reference sequence specifically
        mafft_cline=MafftCommandline(input=concat_temp)
        mafft_cline.adjustdirection = True

        # MAFFT output goes to stdout
        stdout, stderr = mafft_cline()
    
    # Remove temp file as MAFFT output is in stdout and we no longer need the file
    finally:
        os.remove(concat_temp)
    

    alignment = AlignIO.read(io.StringIO(stdout), 'fasta')

    # Find the deletions in the alignment
    # Start enumerate with one to get Genomic positions and not the python index
    deletions = [index for index, character in enumerate(str(alignment[1].seq), start=1) if character == '-']
    
    # Find the non-consecutive numbers in the deletions list which are the gaps between the deletions
    gaps = [[start, end] for start, end in zip(deletions, deletions[1:]) if start+1 < end]

    # Using the gaps, create list iterator object of the start and end of consecutive ranges
    edges = iter(deletions[:1] + sum(gaps, []) + deletions[-1:])

    # Generate list of the start and end (use +1 to get the full range)
    indel_groups = [[start, end+1] for start, end in zip(edges, edges)]
    logging.info('Deletion ranges include: {}'.format(indel_groups))

    for group in indel_groups:
        # Generate a dictionary containing the deletion start sites as keys and the deletion length as values 
        # {5738: 9}
        del_dict[group[0]] = group[1] - group[0]

    return del_dict


def _find_N(sequence, start=True):
    '''
    Find the First (start=True) non-N position or the last non-N (start=False) position

    INPUTS:
        sequence    --> `string` uppercase consensus nucleotide sequence
        start       --> `boolean` switch to control if looking at the start of the sequence
                        (true) or the end of the sequence (false)
    
    RETURNS:
        `integer` position of the first non-N
    '''

    if start:
        # Nucleotide location of first non-N is equal to number of N's at the start
        for position, char in enumerate(sequence):
            if char != 'N':
                return position
    else:
        # Iteration from the end of the sequence, position of first N is equal to the length
        # of the sequence minus the first location that is not an N when iterating from the end
        for position, char in enumerate(sequence[::-1]):
            if char != 'N':
                return len(sequence) - position


def find_all_n(consensus_seq, del_dict):
    '''
    Find the genomic location of all NON-starting N's if no fail.vcf is given

    INPUTS:
        consensus_seq   --> `string` uppercase consensus nucleotide sequence
        del_dict        --> `dictionary` containing starting genomic location of the deletion and the length of the deletion
                            {6031: 6, 12033: 18}

    RETURNS:
        `list` of reference location of the N's found (reference = genomic if no indel)
        `dictionary` of {reference location: genomic location} if there is a deletion, if no deletion, then None
            {9026: 9023, 9027: 9024}
    '''

    first_non_n = _find_N(consensus_seq)
    last_non_n = _find_N(consensus_seq, start=False)

    if last_non_n < first_non_n:
        logging.error('ERROR: Maximum reference location is less than the first non-N character')
        quit()

    if del_dict != {}:
        adjusted_pos_list = []
        tracking_dict = {}
        # Add one to go from index to genomic
        initial_positions = [position+1 for position, char in enumerate(consensus_seq[first_non_n:last_non_n], start=first_non_n) if char == 'N']

        for pos in initial_positions:
            # Track initial positions to modify them if they are greater than the deletion
            initial_position = pos
            for del_start in del_dict.keys():

                if pos > del_start:
                    # Here we add the deletion as we have the genomic location already (from previous enumerate +1 to get genomic)
                    pos = pos + del_dict[del_start]

            # List is the reference positions to check with bcftools
            adjusted_pos_list.append(pos)

            # Tracking dict is a dict that is structured reference_pos: genomic position
            tracking_dict[pos] = initial_position

        return adjusted_pos_list, tracking_dict

    else:
        # We add 1 to the position to go from the index number to its genomic location and None as we have no deletion changes to keep track of
        return [position+1 for position, char in enumerate(consensus_seq[first_non_n:last_non_n], start=first_non_n) if char == 'N'], None


def find_fail_locations(input_vcf, del_dict, output_list=[]):
    '''
    Use fail.vcf file to get the N locations to correct

    INPUTS:
        input_vcf   --> `path` from argparse to the input fail vcf file
        del_dict    --> `dictionary` containing starting genomic location of the deletion and the length of the deletion
                        {6031: 6, 12033: 18}

    RETURNS:
        `list` of reference location of the N's found (reference = genomic if no indel)
        `dictionary` of {reference location: genomic location} if there is a deletion, if no deletion, then None
            {9026: 9023, 9027: 9024}
    '''

    vcf_reader = vcf.Reader(open(input_vcf, 'rb'))

    for rec in vcf_reader:
        pos = rec.POS

        # Add length of the reference to make sure that we check all of the N locations, otherwise it may miss one
        pos_range = list(range(pos, pos + len(rec.REF)))

        output_list.append(pos_range)

    if del_dict != {}:
        # Tracking dict is a dict that is structured reference_pos: genomic position
        tracking_dict = {}
        flat_location_list = [location for sublist in output_list for location in sublist]

        for pos in flat_location_list:
            position_modification = 0

            for del_start in del_dict.keys():
                if pos > del_start:
                    position_modification = position_modification + del_dict[del_start]

            # Here, different from all N's we have the reference position but need the genomic so we subtract
            # Dictionary is formatted the same though
            tracking_dict[pos] = pos - position_modification
        return flat_location_list, tracking_dict

    else:
        # Flattens list to match if no fail vcf given
        return [location for sublist in output_list for location in sublist], None


def generate_location_input(n_locations, ref_name, cmd_out=[]):
    '''
    Generates input location for bcftools

    INPUTS:
        n_locations --> `list` of reference N locations 
        ref_name    --> `string` name of the reference from the reference file
        cmd_out     --> `list` of the command to put into bcftools

    RETURNS:
        `list` of the ref_name:N_locations to be used in bcftools

    '''
    for n_loc_int in n_locations:
        cmd_out.append('{}:{}'.format(ref_name, n_loc_int))

    return cmd_out


def get_ref_from_vcf(filtered_vcf, tracking_dict, out=[]):
    '''
    Takes the bcftools vcf output to get the passing reference called bases and their reference location
    Returns the reference base and the genomic location in a tuple

    INPUTS:
        filtered_vcf    --> `path` to the output vcf file (doesn't change)
        tracking_dict   --> `dictionary` of {reference location: genomic location} if there is a deletion, if no deletion, then None
                                {9026: 9023, 9027: 9024}
        out             --> `list` to be populated with tuples of the (genomic location, reference base)

    RETURNS:
        `list` populated with tuples of (genomic location, reference base)
    '''

    vcf_reader = vcf.Reader(open(filtered_vcf, 'rb'))

    if tracking_dict:
        for rec in vcf_reader:
            # Create tuple of the locations corrected for differences in genome length from the deletion
            # and the reference alleles that pass our selections.
            # If there is more than one base in REF, ignore as its an indel and that isn't supported at the moment
            # It is an indel if it is longer than one as each base input is checked individually
            pos = tracking_dict[rec.POS]
            if len(rec.REF) != 1:
                logging.warning('{} reference position was corrected to an indel of "{}". This correction is not supported. Correct manually or skip.\n'.format(pos, rec.REF))
                continue
            out.append((pos, rec.REF))

    else:
        for rec in vcf_reader:
            # Create tuple of the location and the reference alleles that pass our selections criteria
            # If there is more than one base in REF, ignore as its an indel and that isn't supported at the moment
            # It is an indel if it is longer than one as each base input is checked individually
            if len(rec.REF) != 1:
                logging.warning('{} reference position was corrected to an indel of "{}". This correction is not supported. Correct manually or skip.\n'.format(rec.POS, rec.REF))
                continue
            out.append((rec.POS, rec.REF))

    return out


def generate_fasta(changes, consensus_seq, sample_name, con_name, ref_name):
    '''
    Make changes to the fasta file and outputs a corrected version if the genomic location matches the N in the sequence only
    Also needs the changes to pass the bcftools filtering steps

    INPUTS:
        changes         --> `list` of changes as (genomic location, reference base)
        consensus_seq   --> `string` uppercase of the consensus sequence
        sample_name     --> `string` file sample name to keep the names the same
        con_name        --> `string` header name of the consensus sequence fasta file
        ref_name        --> `string` reference header

    RETURNS:
        None --> outputs the new consensus file
    '''

    list_seq = list(consensus_seq)
    for change in changes:
        # change tuple structured as (Position, Reference Allele)
        # We -1 to the position to go back to index values from genomic ones to change the sequence
        if list_seq[change[0]-1] == 'N':
            logging.info('{} at position {} changed to {}'.format(list_seq[change[0]-1], change[0]-1, change[1]))
            list_seq[change[0]-1] = change[1]
    
    new_seq = SeqRecord(Seq(''.join(list_seq)), id='{}-updated'.format(con_name), description=ref_name, name=sample_name)

    logging.info('Generating fasta file called {}.corrected.consensus.fasta'.format(sample_name))
    with open('{}.corrected.consensus.fasta'.format(sample_name), 'w') as output_handle:
        SeqIO.write(new_seq, output_handle, 'fasta')


def main():
    parser = init_parser()
    args = parser.parse_args()

    ref_seq, ref_name = parse_fasta(args.reference)
    con_seq, con_name = parse_fasta(args.consensus)
    sample_name = os.path.splitext(Path(args.consensus).stem)[0]

    # Logging stuff
    log_name = '{}.corrected.log'.format(sample_name)
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        handlers=[logging.FileHandler(log_name), logging.StreamHandler()]
        )
    logging.info('Starting correction on {}'.format(sample_name))


    # Check for indels that will change locations
    # Blank deletions dict to pass if no deletions
    del_dict = {}

    if len(ref_seq) < len(con_seq):
        logging.error('Insertion detected. Insertions are unsupported at the moment sorry')
        quit()

    elif len(ref_seq) != len(con_seq):
        logging.warning('''
WARNING: Consensus sequence length ({}) does not match the reference sequence ({}).

This means that there is an indel detected and this script is not fully
ready to handle indels. Please note that the correction will not continue unless you pass the "--force" command.

If you do, please double check the BAM file in a viewer such as IGV to make sure everything is accurate!
            '''.format(len(con_seq), len(ref_seq)))

        if args.force:
            logging.info('Forced')
            del_dict = correct_for_indels(args.reference, args.consensus)
        else:
            quit()

    # Find N's based on either the consensus sequence or the failed variants
    if not args.fail_vcf:
        n_locations, tracking_dict = find_all_n(con_seq, del_dict)

        if len(n_locations) > args.max:
            logging.error('WARNING: Greater than {} Ns detected. Change the "--max" command to at least "--max {}" to allow this analysis to run'.format(args.max, len(n_locations)))
            quit()
        logging.info('N locations include {}'.format(n_locations))

    else:
        n_locations, tracking_dict = find_fail_locations(args.fail_vcf, del_dict)
        logging.info('N locations include {}'.format(n_locations))

    # Generate command for bcftools mpileup targeting those sites
    cmd_input = ','.join(generate_location_input(n_locations, ref_name))
    logging.info('Command locations formatting: {}'.format(cmd_input))
  
    if cmd_input == '':
        if con_seq[300:29700].count('N') != 0: # temporary solution for this as we need to see if no internal N's or all N's
            logging.error('Sequence is all Ns!')
            quit()
        logging.info('No internal Ns in the input, no changes needed! :)')
        quit()

    try:
        cmd = '''
    bcftools mpileup {} -Ou -f {} --max-depth 10000 -r {} \
    | bcftools call -c --ploidy 1 \
    | bcftools filter -i "%QUAL>100 && AF1==0 && DP>=100" -o ./{}.filtered.vcf
            '''.format(args.bam, args.reference, cmd_input, sample_name)

        logging.info('Running the following commands: {}'.format(cmd))
        subprocess.run(cmd, shell=True)
        subprocess.run('bgzip -f ./{}.filtered.vcf'.format(sample_name), shell=True)

    except OSError:
        logging.error("OSError: Argument list too long, too many Ns in sample")
        quit()

    changes = get_ref_from_vcf('{}.filtered.vcf.gz'.format(sample_name), tracking_dict)

    if changes == []:
        logging.info('No N changes to be made based! Exiting, have a nice day! :)')
        quit()
    
    else:
        generate_fasta(changes, con_seq, sample_name, con_name, ref_name)

if __name__ == "__main__":
    main()
