#!/usr/bin/env python3

from Bio import SeqIO
from functools import reduce
from pybedtools import BedTool
import subprocess
import vcf
import re
import pandas as pd
import matplotlib.pyplot as plt
import shlex

"""
This script can incorporate as many QC checks as required
as long as it outputs a csv file containing a final column
headed with 'qc_pass' and rows for each sample indcating
'TRUE' if the overall QC check has passed or 'FALSE' if not.
"""

def make_qc_plot(depth_pos, n_density, samplename, window=200):
    depth_df = pd.DataFrame( { 'position' : [pos[1] for pos in depth_pos], 'depth' : [dep[2] for dep in depth_pos] } )
    depth_df['depth_moving_average'] = depth_df.iloc[:,1].rolling(window=window).mean()

    n_df = pd.DataFrame( { 'position' : [pos[0] for pos in n_density], 'n_density' : [dens[1] for dens in n_density] } )

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_xlabel('Position')

    ax1.set_ylabel('Depth', color = 'g')
    ax1.set_ylim(top=10**5, bottom=1)
    ax1.set_yscale('log')
    ax1.plot(depth_df['depth_moving_average'], color = 'g')

    ax2.set_ylabel('N density', color = 'r')  
    ax2.plot(n_df['n_density'], color = 'r')
    ax2.set_ylim(top=1)

    plt.title(samplename)
    plt.savefig(samplename + '.depth.png')

def read_depth_file(bamfile):
    p = subprocess.Popen(['samtools', 'depth', '-a', '-d', '0', bamfile],
                       stdout=subprocess.PIPE)
    out, err = p.communicate()
    counter = 0

    pos_depth = []
    for ln in out.decode('utf-8').split("\n"):
       if ln:
          pos_depth.append(ln.split("\t"))

    return pos_depth


def get_covered_pos(pos_depth, min_depth):
    counter = 0
    for contig, pos,depth in pos_depth:
        if int(depth) >= min_depth:
            counter = counter + 1
    
    return counter

def get_depth_coverage(pos_depth, ref_length):
    depth_total = 0
    for contig, pos, depth in pos_depth:
        depth_total = depth_total + int(depth)

    return round(depth_total/ref_length)

def get_N_positions(fasta):
    n_pos =  [i for i, letter in enumerate(fasta.seq.lower()) if letter == 'n']

    return n_pos

def get_pct_N_bases(fasta):
    
    count_N = len(get_N_positions(fasta))

    pct_N_bases = count_N / len(fasta.seq) * 100

    return pct_N_bases, count_N

def get_largest_N_gap(fasta):
    n_pos = get_N_positions(fasta)

    n_pos = [0] + n_pos + [len(fasta.seq)]

    n_gaps = [j-i for i, j in zip(n_pos[:-1], n_pos[1:])]

    return sorted(n_gaps)[-1]

def get_ref_length(ref):
    record = SeqIO.read(ref, "fasta")
    return len(record.seq)

def sliding_window_N_density(sequence, window=10):

    sliding_window_n_density = []
    for i in range(0, len(sequence.seq), 1):
        window_mid = i + ( window / 2)
        window_seq = sequence.seq[i:i+window]
        n_count = window_seq.lower().count('n')
        n_density = n_count / window

        sliding_window_n_density.append( [ window_mid, n_density ] )

    return sliding_window_n_density

def get_num_reads(bamfile):

    st_filter = '0x900'
    command = 'samtools view -c -F{} {}'.format(st_filter, bamfile)
    what = shlex.split(command)

    return subprocess.check_output(what).decode().strip()

def get_vcf_variants(variants_vcf, variants_list=[], locations=[]):
    '''
    Use pipeline pass.vcf.gz file to get passing variants and their locations
    INPUTS:
        variants_vcf   --> `path` from argparse to the input vcf.gz file
        variants_list  --> `list` to append the variants formated as RefPosAlt
                        ['G231T', 'A3423C']
        locations      --> `list` to append only the variant location ints to
    RETURNS:
        'None' or str variants separated by a ;
        location list or None for the variant int location
    '''
    vcf_reader = vcf.Reader(open(variants_vcf, 'rb'))
    for rec in vcf_reader:
        variant = '{}{}{}'.format(rec.REF, rec.POS, rec.ALT[0])

        # Checking for duplicate variants that have been an issue
        if variant in variants_list:
            pass
        # Removal of N variants
        elif rec.ALT[0] == 'N':
            pass
        else:
            variants_list.append(variant)
            locations.append(rec.POS)


    variants = (';'.join(variants_list))

    if variants == '':
        return 'None', None

    return variants, locations

def get_tsv_variants(variants_tsv, variants_list=[], locations=[]):
    '''
    Use pipeline variants.tsv file to get passing variants and their locations
    INPUTS:
        variants_tsv   --> `path` from argparse to the input tsv file generated by signal
        variants_list  --> `list` to append the variants formated as RefPosAlt
                        ['G231T', 'A3423C']
        locations      --> `list` to append only the variant location ints to
    RETURNS:
        `str` of 'None' or variants separated by ;
        `list` of locations or None if no variants found
    '''
    with open(variants_tsv) as input_handle:

        for index, line in enumerate(input_handle):

            # Pass the header line so as to not include it!
            if index == 0:
                continue

            row = line.strip('\n').split('\t') # Order is [REGION, POS, REF, ALT, ...]

            variant = '{}{}{}'.format(row[2], row[1], row[3])

            # Checking for duplicate variants that have been an issue
            if variant in variants_list:
                pass
            else:
                variants_list.append(variant)
                locations.append(row[1])
    
    variants = (';'.join(variants_list))

    if variants == '':
        return 'None', None
        
    return variants, locations

def find_primer_mutations(pcr_bed, genomic_locations, primer_mutations=[]):
    '''
    Use variant info to check for mutations in currently used primers
    INPUTS:
        pcr_bed            --> `path` from argparse to input bed file of PCR primers
        genomic_locations  --> `list` of integer locations to check or `None` if none found
        primer_mutations   --> `list` to append any primer mutations found
    RETURNS:
        `str` primer mutations statement
    '''
    if genomic_locations is None:
        return 'None'
    
    input_bed = BedTool(pcr_bed)

    for gene in input_bed:
        location = range(gene.start, gene.stop + 1) # Plus one to make sure that we get mutations in the final location of the range

        for variant_pos in genomic_locations:
            if variant_pos in location:
                primer_mutations.append('Variant position {} overlaps PCR primer {}'.format(variant_pos, gene.name))
    
    if primer_mutations != []:
        statement = '; '.join(primer_mutations)
        return 'Warning: {}'.format(statement)

    return 'None'

def get_lineage(pangolin_csv, sample_name):
    '''
    Check Pangolin for the pangoLEARN version as that isn't captured by ncov-tools
    INPUTS:
        pangolin_csv --> `path` from argparse to input pangolin csv file
        sample_name  --> `str` sample name from argparse
    RETURNS:
        `str` pangoLEARN version
    '''
    df = pd.read_csv(pangolin_csv)
    df_slice = df.loc[df['taxon'] == sample_name]

    if not df_slice.empty:
        pangoV = df_slice['pangoLEARN_version'].any()
        return pangoV
    else:
        return 'Unknown'

def get_protein_variants(aa_table):
    '''
    Parse ncov-tools output to report its amino acid mutations along with finding any problem consequences
    INPUTS:
        aa_table  --> `path` from argparse to input aa_table.csv
    RETURNS:
        `str` amino acid mutations separated by a ;
        `str` amino acid gisaid problem sites separated by a ;
    '''
    df = pd.read_csv(aa_table, sep='\t')
    df = df[~df.alt.str.contains("N")] # remove alt allele Ns the N nucleotide means we didn't seq it
    protein_mutations = ';'.join(df['aa'].dropna().tolist())

    # Check for consequences that may cause frame issues. http://pcingola.github.io/SnpEff/se_inputoutput/#eff-field-vcf-output-files
    consequences_to_check = ['stop_gained', 'gene_fusion', 'frameshift_variant', 'start_list', 'feature_ablation']
    df = df[df['Consequence'].isin(consequences_to_check)]
    if df.empty:
        found_consequences = 'NA'
    else:
        df['cons_pro'] = df['Consequence'] + '-' + df['gene']
        found_consequences = ';'.join(set(df['cons_pro']))

    return protein_mutations, found_consequences

def parse_ncov_tsv(file_in, sample, negative=False):
    '''
    Parse ncov-tools output tsv files (summary and negative) to grab data for the sample
    INPUTS:
        file_in   --> `path` from argparse to input ncov-tools tsv file to parse
        sample    --> `str` sample name from argparse
        negative  --> `boolean` for if the tsv table is for the negative control or not
    RETURNS:
        Populated `df`
    '''

    # Try to read file (as negative control may not have data in it)
    try:
        df = pd.read_csv(file_in, sep='\t')

    # If no data, we set up how it should be and then pass it through
    # Could also make is such that runs without negative ctrls just don't have the columns
    except pd.errors.EmptyDataError:
        negative_df = pd.DataFrame(columns=['sample', 'qc', 'genome_covered_bases', 'genome_total_bases', 'genome_covered_fraction', 'amplicons_detected'])
        negative_df.loc[1, 'sample'] = sample
        negative_df.fillna('NA', inplace=True)
    
        return negative_df

    # If the column headers get changed in the input just replace the new sample column below
    # First column is the file name
    if negative:
        new_columns = df.columns.values
        new_columns[0] = 'sample'
        df.columns = new_columns
    # Input is summary_df, drop run_name from this output as its always `nml` and we have our own name that is put through
    else:
        df.drop(columns=['run_name'], inplace=True)

    # Set which column contains the sample
    sample_column = 'sample'

    # Finding the data, all samples will be in ncov-tools summary output (as they had data generated)
    for index, name in enumerate(df[sample_column].tolist()):
        if re.search(sample, name):
            df.loc[index, sample_column] = sample
            df.fillna('NA', inplace=True)
            return df.iloc[[index]]
    
    # If nothing is found, input is not a negative control and we need to keep negative columns
    negative_df = pd.DataFrame(columns=new_columns)
    negative_df.loc[1, sample_column] = sample
    negative_df.fillna('NA', inplace=True)
    
    return negative_df

def get_samplesheet_info(sample_tsv, sample_name, project_id):
    '''
    Parse samplesheet info to allow for IRIDA uploads and adding whatever data wanted to output qc file
    INPUTS:
        sample_tsv   --> `path` from argparse to input samplesheet.tsv file
        sample_name  --> `str` sample name from argparse
        project_id   --> `str` project ID from CL to add if there isn't one
    RETURNS:
        `df` populated with data from samplesheet
    '''
    df = pd.read_csv(sample_tsv, sep='\t', dtype=object)
    samplesheet_columns = df.columns.values
    # Rename run to run_identifier as that is what is already in IRIDA
    if 'run' in samplesheet_columns:
        df.rename(columns={'run': 'run_identifier'}, inplace=True)

    # ncov-tools captures ct and date from this file so remove, scheme is sometimes here and we capture it ourselves so remove
    columns_to_remove = set(['ct', 'date', 'scheme', 'primer_scheme']).intersection(samplesheet_columns)
    if columns_to_remove:
        df.drop(columns=columns_to_remove, inplace=True)

    if 'project_id' not in samplesheet_columns:
        df['project_id'] = project_id
    
    # Get only the sample row and if empty, fill it in to match other rows
    df = df.loc[df['sample'] == sample_name]
    if df.empty:
        df.loc[1, 'sample']  = sample_name
    
    df.fillna('NA', inplace=True)

    return df

def go(args):
    if args.illumina:
        depth = 10
    elif args.nanopore:
        depth = 20

    ## Depth calcs
    ref_length = get_ref_length(args.ref)
    depth_pos = read_depth_file(args.bam)

    depth_covered_bases = get_covered_pos(depth_pos, depth)

    depth_coverage = get_depth_coverage(depth_pos, ref_length)

    pct_covered_bases = depth_covered_bases / ref_length * 100

    ## Number of aligned reads calculaton
    num_reads = get_num_reads(args.bam)

    # Unknown base calcs
    fasta = SeqIO.read(args.fasta, "fasta")

    pct_N_bases   = 0
    largest_N_gap = 0
    qc_pass       = "FALSE"

    if len(fasta.seq) != 0:

        pct_N_bases, count_N = get_pct_N_bases(fasta)
        largest_N_gap = get_largest_N_gap(fasta)

    	# QC PASS / FAIL
        if largest_N_gap >= 10000 or pct_N_bases < 50.0:
                qc_pass = "TRUE"

    ### Added checks ###
    ####################
    if args.nanopore:
        # Vcf passing variants from nanopore pipeline
        variants, variant_locations = get_vcf_variants(args.vcf)
    
    elif args.illumina:
        # Tsv variants from Illumina pipeline
        if args.tsv_variants:
            variants, variant_locations = get_tsv_variants(args.tsv_variants)
        # VCF variants from FREEBAYES
        else:
            variants, variant_locations = get_vcf_variants(args.vcf)

    # Find any overlap of variants in the pcr primer regions
    primer_statement = find_primer_mutations(args.pcr_bed, variant_locations)

    # PangoLEARN version
    pangoLearn = get_lineage(args.pangolin, args.sample)

    # snpEFF output
    protein_variants, found_consequences = get_protein_variants(args.snpeff_tsv)

    # NCOV-Tools Results
    summary_df = parse_ncov_tsv(args.ncov_summary, args.sample)
    negative_df = parse_ncov_tsv(args.ncov_negative, args.sample, negative=True)

    # If we have a samplesheet, use its values to create final output
    if args.sample_sheet:
        sample_sheet_df = get_samplesheet_info(args.sample_sheet, args.sample, args.project_id)

        qc_line = {  'sample' : [args.sample],
           'num_aligned_reads': [num_reads],
                   'variants' : [variants],
            'protein_variants': [protein_variants],
'snpeff_frameshift_consequence' : [found_consequences],
'diagnostic_primer_mutations' : [primer_statement],
                     'scheme' : [args.scheme],
      'sequencing_technology' : [args.sequencing_technology],
         'pangoLEARN_version' : [pangoLearn],
                'script_name' : ['nml-ncov2019-artic-nf'],
                   'revision' : [args.revision]}
        
        qc_df = pd.DataFrame.from_dict(qc_line)
        data_frames = [sample_sheet_df, qc_df, summary_df, negative_df]

    else:
        run_name = 'NA'
        if args.nanopore:
            barcode = re.search(r'\d+', args.sample).group(0)
        else:
            barcode = 'NA'

        qc_line = {      'sample' : [args.sample],
                     'project_id' : [args.project_id],
                        'barcode' : [barcode],
               'num_aligned_reads': [num_reads],
                       'variants' : [variants],
                'protein_variants': [protein_variants],
  'snpeff_frameshift_consequence' : [found_consequences],
    'diagnostic_primer_mutations' : [primer_statement],
                         'scheme' : [args.scheme],
          'sequencing_technology' : [args.sequencing_technology],
             'pangoLEARN_version' : [pangoLearn],
                 'run_identifier' : [run_name],
                    'script_name' : ['nml-ncov2019-artic-nf'],
                       'revision' : [args.revision]}

        qc_df = pd.DataFrame.from_dict(qc_line)

        data_frames = [qc_df, summary_df, negative_df]

    # Merge all dataframes together
    out_df = reduce(lambda left,right: pd.merge(left,right,on='sample', how='left'), data_frames)

    # Remove comma's as some of the ncov-tools fields have commas :(
    out_df.replace(',',';', regex=True, inplace=True)

    # Add final column as the nextflow_qc_pass designation to not error rest of pipeline till I figure it out to change it
    out_df['nextflow_qc_pass'] = qc_pass

    # Output
    out_df.to_csv(args.outfile, sep=',', index=False)
    N_density = sliding_window_N_density(fasta)
    make_qc_plot(depth_pos, N_density, args.sample)

def main():
    import argparse

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--nanopore', action='store_true')
    group.add_argument('--illumina', action='store_true')
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--ref', required=True)
    parser.add_argument('--bam', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--pangolin', required=True)
    parser.add_argument('--ncov_summary', required=True)
    parser.add_argument('--ncov_negative', required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--vcf', required=False)
    parser.add_argument('--tsv_variants', required=False)
    parser.add_argument('--sequencing_technology', required=True)
    parser.add_argument('--scheme', required=True)
    parser.add_argument('--snpeff_tsv', required=True)
    parser.add_argument('--pcr_bed', required=True)
    parser.add_argument('--project_id', required=False, default="NA")
    parser.add_argument('--sample_sheet', required=False)

    args = parser.parse_args()
    go(args)

if __name__ == "__main__":
    main()
