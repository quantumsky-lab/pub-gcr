# Copyright (c) Quantum Sky, Inc. and Jorna, Inc. and their affiliates
#
# This script intends to act like a wrapper for running aligner (BWA) and then SAMTools
# with common parameters that Jorna uses for internal purposes set as default
# for machine image rnax-e2 (aspen)

# Dependencies:
# 1. Trimmomatic
# 2. Blast+ 

import sys, re, os
import argparse
import shutil

from scipy.stats import ttest_ind
import numpy as np
import pandas as pd

from Bio import SeqIO, SearchIO
from Bio.SearchIO import BlastIO
from subprocess import run, PIPE, Popen

def mark_homopolymer_sites(seq, min_length):
    # Return a list of tuples (start, end) of homopolymer sites
    # min_length: minimum length of homopolymer to be considered
    # seq: a string of nucleotides
    # Return: a list of tuples (start, end) of homopolymer sites
    homopolymer_sites = []
    start = 0
    end = 0
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            end = i
        else:
            if end - start + 1 >= min_length:
                homopolymer_sites.append((start, end))
            start = i
            end = i
    if end - start + 1 >= min_length:
        homopolymer_sites.append((start, end))
    return homopolymer_sites

def cal_dist_to_homopolymer(p, homopolymer_sites):
    # Return the distance from position p to the nearest homopolymer site
    # p: a position
    # homopolymer_sites: a list of tuples (start, end) of homopolymer sites
    # Return: the distance from position p to the nearest homopolymer site
    dist = 1000000000
    for site in homopolymer_sites:
        if p >= site[0] and p <= site[1]:
            return 0
        else:
            dist = min(dist, abs(p - site[0]), abs(p - site[1]))
    return dist

def run_makeblastdb(infile, reference, makeblastdb='makeblastdb'):
    makeblastdb_cmd = [makeblastdb, '-in', infile, '-dbtype', 'nucl', '-out', reference]
    p1 = Popen(makeblastdb_cmd)
    p1.communicate()
    return 0

def run_blastn(infile, outfile, reference, blastn='blastn'):
    blastn_cmd = [blastn, '-query', infile, '-db', reference, '-out', outfile, 
                  '-outfmt', '5', '-max_target_seqs', '1', '-evalue', '1e-6']
    p1 = Popen(blastn_cmd, stdout=PIPE, stderr=PIPE)
    p1.communicate()
    return 0

def fastq2fasta(infile, outfile):
    ifh = open(infile, 'r')
    ofh = open(outfile, 'w')
    while 1:
        line = ifh.readline()
        if not line:
            break
        if line.startswith('@'):
            ofh.write('>'+line[1:])
            seq = ifh.readline().rstrip()
            ofh.write(seq+'\n')
            ifh.readline()
            ifh.readline()
    ofh.close()

def test_jar(jar_file):
    try:
        run(["java", "-jar", jar_file], stdout=PIPE, stderr=PIPE)
        return 1
    except FileNotFoundError:
        return 0

def run_trimmomatic(infiles, prefix, temp_dir, trimmomatic='trimmomatic.jar'):
    if len(infiles) == 1:
        mode = 'SE'
    elif len(infiles) == 2:
        mode = 'PE'
    else:
        print('Error: more than 2 infiles provided. Abort')
        exit()
    if mode == 'SE':
        clip_file = os.path.join(os.path.dirname(trimmomatic), 'adapters/TruSeq3-SE.fa')
        if not os.path.isfile(clip_file):
            print('Error: cannot find TruSeq3-SE.fa in the same folder as trimmomatic.jar. Abort')
            exit()
        outfile = os.path.join(temp_dir, prefix + '.trimmed.fastq')
        param = f'ILLUMINACLIP:{clip_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'.split(' ')
        #cmd = ['java', '-jar', trimmomatic, 'SE', '-phred33', infiles[0], outfile]+param
        cmd = ['cp', infiles[0], outfile]
        print(f'  ## Sample {prefix}:\n'+'\t'+' '.join(cmd))
        run(cmd, stdout=sys.stdout, stderr=sys.stderr)
    else:
        print('PE mode not implemented yet. Abort')
        exit()

def reverse_complement(seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', '-':'-'}
    return ''.join([complement[base] for base in seq[::-1]])

# define a few globals
#trimmomatic = '/Applications/bioinformatics/Trimmomatic-0.39/trimmomatic-0.39.jar'

def main():
    parser = argparse.ArgumentParser(description = "Run reads mapping with Jorna default settings for genome editing data.")
    parser.add_argument('--infile-list', '-i', dest='infile_list', help="A list of input files in fastq format in text file;\
                        if paired-end, they should be in the same line, separated by comma. NOTE: no directory should be supplied.", required=True)
    parser.add_argument('--infile-vol', '-d', dest='infile_vol', help="directory where the infiles are stored", required=True)

    parser.add_argument('--reference', '-r', dest='reference', help="Reference sequence in fasta format", required=True)
    parser.add_argument('--temp-dir', '-t', dest='temp_dir', help="Where the intermediate files should live.", default=os.path.join(os.getcwd(), 'temp'))
    parser.add_argument('--keep-temp', '-k', dest='keep_temp', help="Turns on temp dir keeping when specified.", action='store_true')
    parser.add_argument('--prefix', '-o', dest='prefix', help="Output file prefix", required=True)
    parser.add_argument('--blastn', dest='blastn', help="Path to blastn", default='blastn')
    parser.add_argument('--makeblastdb', dest='makeblastdb', help="Path to makeblastdb", default='makeblastdb')
    parser.add_argument('--trimmomatic', dest='trimmomatic', help="Path to trimmomatic jar", default='/root/trimmomatic-0.39.jar')
    parser.add_argument('--homopolymer', dest='homopolymer', help="Homopolymer length threshold", default=3, type=int)
    parser.add_argument('--min-map', dest='min_map', help="Mininum match length threshold", default=200, type=int)
    
    
    args = parser.parse_args()

    # check if all the dependencies are available
    print(f'Checking dependencies...')
    if shutil.which(args.makeblastdb) is None:
        print(f'[FATAL] Blast+:makeblastdb not found.')
        exit()
    if shutil.which(args.blastn) is None:
        print(f'[FATAL] Blast+:blastn not found.')
        exit()
    if test_jar(args.trimmomatic) == 0:
        print(f'[FATAL] trimmomatic not found.')
        exit()
    print(f'All dependencies are available.')
    
    # make sure the scratch directory exists
    if not os.path.exists(args.temp_dir):
        try:
            os.mkdir(args.temp_dir)
        except:
            print(f'[FATAL] Cannot create temporary directory.')
            exit()

    # step 0, check if the input files are in the correct format
    print(f'[Step 0] Checking input files...')
    infiles = {}
    infile_ind = []
    for line in open(args.infile_list, 'r'):
        if len(line.rstrip()) == 0: continue
        if ',' in line:
            fx = [os.path.join(args.infile_vol, a) for a in line.rstrip().split(',')]
        else:
            fx = [os.path.join(args.infile_vol, a) for a in [line.rstrip()]]
        for f in fx:
            if not os.path.exists(f):
                print(f'[FATAL] {f} does not exist.')
                exit()
            if not f.endswith('.fastq'):
                print(f'[FATAL] {f} is not in fastq format.')
                exit()
        fx_infile = os.path.basename(fx[0]).replace('.fastq', '')
        infiles[fx_infile] = fx
        infile_ind.append(fx_infile)

    print(f'Found {len(infiles)} input files in {args.infile_list}.')
    print(f'And they are:')
    for fx_infile in infiles:
        print(f'{fx_infile}: {infiles[fx_infile]}')
    

    # step 1, mark all the homopolymer sites on the reference
    print(f'[Step 1] Marking homopolymer sites on the reference ...')
    refseq = ''
    for rec in SeqIO.parse(args.reference, 'fasta'):
        refseq = str(rec.seq).upper()
    homopolymer_sites = mark_homopolymer_sites(refseq, args.homopolymer)
    for s, e in homopolymer_sites:
        print(f'  {s+1}-{e+1}: {refseq[s:e+1]}')
    print(f'  ## Done.\n') 

    # step 2, trim the fastq files using Trimmomatic
    print(f'[Step 2] Running Trimmomatic...')
    for ind in infile_ind:
        print(f'  [Sample: {ind}]')
        fastqs = infiles[ind]
        if not os.path.exists(fastqs[0]):
            run_trimmomatic(fastqs, ind, args.temp_dir, trimmomatic = args.trimmomatic)
    print(f'  ## Done.\n')
    
    # step 3, run Blast+ to map the reads to the reference
    print(f'[Step 3] Running Blast+ to generate alignments ...')
    print(f'  Using reference: {args.reference}')
    blast_index_prefix = os.path.join(args.temp_dir, '.'.join(os.path.basename(args.reference).split('.')[:-1]))
    blast_index_file = blast_index_prefix+'.nsq'
    blastfiles = {}
    if not os.path.exists(blast_index_file):
        print(f'  ## Building blast index ...')
        run_makeblastdb(args.reference, blast_index_prefix, makeblastdb=args.makeblastdb)
        print(f'  ## Done.\n')

    for ind in infile_ind:
        print(f'  [Sample: {ind}]')
        fastqs = infiles[ind]
        fastas = [os.path.join(args.temp_dir, os.path.basename(f).replace('.fastq', '.fasta')) for f in fastqs]
        blastfiles[ind] = []
        for fastq, fasta in zip(fastqs, fastas):
            fastq2fasta(fastq, fasta)
            blastfile = fasta.replace('.fasta', '.blast')
            if not os.path.exists(blastfile):
                run_blastn(fasta, blastfile, blast_index_prefix, blastn=args.blastn) 
            blastfiles[ind].append(blastfile)
    print(f'  ## Done.\n')
    
    # step 4, parse the blast results, count indels and mutations
    print(f'[Step 4] Parsing blast results ...')
    mpileups = []
    # initialize the mpileup file
    for position, ref_base in enumerate(list(refseq.upper()), 1):
        mpileups.append([position, ref_base, {}]) # position, ref_base
        # for each position, we have a dict keyed by sample name, and valued by list as  [A, C, G, T, insertion, deletion]
        for sample in infile_ind:
            mpileups[-1][2][sample] = [0, 0, 0, 0, 0, 0]
    # parse the blast results
    for sample in infile_ind:
        blast_outputs = blastfiles[sample]
        i = 0
        for blast_output in blast_outputs:
            if os.stat(blast_output).st_size == 0:
                print(f'[WARNING] Skipping {blast_output}, size zero.')
                continue
            for rec in SearchIO.parse(blast_output, 'blast-xml'):
                i+=1
                try:
                    hsp = rec.hsps[0]
                except:
                    continue
                query_start, query_end = hsp.query_start, hsp.query_end  # all 0-based
                hit_start, hit_end = hsp.hit_start, hsp.hit_end # all 0-based
                query_match_length = query_end - query_start + 1
                if query_match_length < args.min_map:
                    continue
                if hsp.hit_frame == -1:
                    hit_seq = reverse_complement(hsp.hit.seq)
                    query_seq = reverse_complement(hsp.query.seq)
                else:
                    hit_seq = hsp.hit.seq
                    query_seq = hsp.query.seq
                query_start, query_end = hsp.query_start, hsp.query_end
                hit_start, hit_end = hsp.hit_start, hsp.hit_end
                xind = hit_start
                for r, q in zip(hit_seq, query_seq):
                    if r == '-':
                        mpileups[xind][2][sample][4] += 1
                    else:
                        if q == '-':
                            mpileups[xind][2][sample][5] += 1
                        if q == 'A':
                            mpileups[xind][2][sample][0] += 1
                        elif q == 'C':
                            mpileups[xind][2][sample][1] += 1
                        elif q == 'G':
                            mpileups[xind][2][sample][2] += 1
                        elif q == 'T':
                            mpileups[xind][2][sample][3] += 1
                        xind += 1
    print(f'  ## Done.\n')

    # step 5, output the mpileup file
    print(f'[Step 5] Outputting the mpileup file ...')
    outfile = args.prefix+'.mpileup'
    print(f'  Output file: {outfile}')
    ofh = open(outfile, 'w')
    sample_str = '\t'.join(infile_ind)
    ofh.write(f'#position\tref_base\tdist_to_homopolymer\t{sample_str}\n')
    ofh.write(f'#per position: (A, C, G, T, insertion, deletion)\n')
    ofh.write(f'#homopolymer length limit: {args.homopolymer}nt\n')
    for position_ind in range(len(mpileups)):
        position, ref_base, xref = mpileups[position_ind]
        dist_to_homopolymer = cal_dist_to_homopolymer(position, homopolymer_sites)
        ofh.write(f'{position}\t{ref_base}\t{dist_to_homopolymer}\t')
        xref_str = []
        for sample in infile_ind:
            xref_str.append(f'{xref[sample]}')
        ofh.write('\t'.join(xref_str)+'\n')
    ofh.close()
    # remove the temp dir
    if not args.keep_temp:
        shutil.rmtree(args.temp_dir)
    print(f'  ## Done.\n')

    # step 6, generate the coverage and indel rate matrices
    print(f'[Step 6] Generating the coverage and indel rate matrices ...')
    cov_outfile = args.prefix + '.cov'
    indel_outfile = args.prefix + '.indel'

    pos_list = []
    indel_rate = []
    cov = []
    for idx, line in enumerate(open(outfile, 'r')):
        if idx == 0:
            cols = line.rstrip().split('\t')
            samples = [a.split('_')[0] for a in cols[3:]]

        if line.startswith('#'): continue
        cols = line.rstrip().split('\t')
        pos_st, ref_base, dist_st = cols[:3]
        dx = []
        dc = []
        for s in cols[3:]:
            sg = re.search('\[(.+)\]', s).group(1)
            a, c, g, t, i, d = [int(n) for n in sg.split(',')]
            if a+c+g+t == 0:
                r = -0.0001
            else:
                r = float(i+d)/(a+c+g+t)
            dx.append(r)
            dc.append(a+c+g+t)
        indel_rate.append(dx)
        cov.append(dc)
        pos_list.append(pos_st)

    indel_df = pd.DataFrame(indel_rate, columns=samples, index=pos_list)
    cov_df = pd.DataFrame(cov, columns=samples, index=pos_list)

    # write to file
    indel_df.to_csv(indel_outfile, sep='\t')
    cov_df.to_csv(cov_outfile, sep='\t')

    print(f'  ## Done.\n')

if __name__ == '__main__':
    main()
