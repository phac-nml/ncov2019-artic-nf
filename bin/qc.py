#!/usr/bin/env python3

from Bio import SeqIO
import csv
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

def get_variants(variants_vcf, variants_list=[]):
    vcf_reader = vcf.Reader(open(variants_vcf, 'rb'))
    for rec in vcf_reader:
        variants_list.append('{}{}{}'.format(rec.REF, rec.POS, rec.ALT[0]))
    
    variants = (';'.join(variants_list))

    if variants == '':
        return 'None'
        
    return variants

def get_lineage(pangolin_csv, sample_name):
    with open(pangolin_csv, 'r') as input_handle:
        reader = csv.reader(input_handle)

        for row in reader: # Row format is ['taxon', 'lineage', 'SH-alrt', 'UFbootstrap', 'lineages_version', 'status', 'note']

            if re.search(sample_name, row[0]):
                return str(row[1])
    
    return 'Unknown'

def get_ncovtools_qc(ncovtools_tsv, sample_name):
    with open(ncovtools_tsv) as input_handle:

        for line in input_handle:
            row = line.strip('\n').split('\t') # Order is [sample, ..., qc_pass/qc_fail]

            if re.search(sample_name, row[0]):
                return str(row[-1])


def get_samplesheet_info(sample_tsv, sample_name):
    with open(sample_tsv) as input_handle:

        for line in input_handle:
            row = line.strip('\n').split('\t') # Order is [sample, run, barcode, project_id, ct]

            if re.search(sample_name, row[0]):
                return str(row[1]), str(row[2]), str(row[3]), str(row[4])

    return 'Unknown', 'Unknown', 'Unknown', 'Unknown'
    
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

    # Vcf passing variants
    variants = get_variants(args.vcf)

    # Pangolin Lineages
    lineage = get_lineage(args.pangolin, args.sample)

    # ncov-tools
    ncov_tools_status = get_ncovtools_qc(args.ncovtools, args.sample).replace(',', '|')

    if args.sample_sheet:
        run_name, barcode, project_id, ct = get_samplesheet_info(args.sample_sheet, args.sample)
    
    else:
        run_name = 'N/A'
        barcode = re.search(r'\d+', args.sample).group(0)
        project_id = 'N/A' 
        ct = 'N/A'


    qc_line = { 'sample_name' : args.sample,
                 'project_id' : project_id,
                    'barcode' : barcode,
                    'count_N' : count_N,
                'pct_N_bases' : "{:.2f}".format(pct_N_bases),
          'pct_covered_bases' : "{:.2f}".format(pct_covered_bases), 
           'longest_no_N_run' : largest_N_gap,
             'depth_coverage' : depth_coverage,
          'num_aligned_reads' : num_reads,
                    'lineage' : lineage,
                   'variants' : variants,
                         'ct' : ct,
                   'run_name' : run_name,
                'script_name' : 'nml-ncov2019-artic-nf',
                   'revision' : args.revision,
              'ncov-tools-qc' : ncov_tools_status,
                    'qc_pass' : qc_pass}


    with open(args.outfile, 'w') as csvfile:
        header = qc_line.keys()
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        writer.writerow(qc_line)

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
    parser.add_argument('--ncovtools', required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--sample_sheet', required=False)
    parser.add_argument('--vcf', required=True)

    args = parser.parse_args()
    go(args)

if __name__ == "__main__":
    main()
