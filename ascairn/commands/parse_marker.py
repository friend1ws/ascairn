import click
import os, statistics

from ascairn.logger import get_logger
from ascairn.my_seq import reverse_complement

logger = get_logger(__name__)

@click.command()
@click.argument("cen_fasta_list", type=click.Path(exists=True))
@click.argument("output_prefix", type=click.Path())
@click.argument("kmer_size", type = int)
@click.option("--alpha_satellite_margin_size", type = int, default = 0)
def parse_marker_command(cen_fasta_list, output_prefix, kmer_size, alpha_satellite_margin_size): 

    # make directory for the output prefix
    output_dir = os.path.dirname(output_prefix)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    logger.info("Parsing raw rare kmers")
    rare_kmer_parse(cen_fasta_list, output_prefix + ".rare_kmer.txt", kmer_size)

    logger.info("Filtering redundant rare kmers via minimizer")
    filt_kmer(output_prefix + ".rare_kmer.txt", cen_fasta_list, output_prefix + ".rare_kmer.pruned.txt")

    logger.info("Adding positional information to rare kmers")
    add_rare_kmer_info(output_prefix + ".rare_kmer.pruned.txt", cen_fasta_list, output_prefix + ".rare_kmer.pruned.annot.txt", alpha_satellite_margin_size)

    proc_rare_kmer_table(output_prefix + ".rare_kmer.pruned.annot.txt", output_prefix + ".rare_kmer.pruned.annot.long.txt")

    logger.info("Completed.")


def gather_kmer(fasta_file, kmer_size): 

    kmer2count = {}

    def rare_kmer_check(qid, qseq):

        for i in range(len(qseq) - kmer_size + 1):
            kmer = qseq[i:(i + kmer_size)]
            rkmer = reverse_complement(kmer)

            if kmer < rkmer:
                if kmer not in kmer2count: kmer2count[kmer] = 0
                kmer2count[kmer] = kmer2count[kmer] + 1
            else:
                if rkmer not in kmer2count: kmer2count[rkmer] = 0
                kmer2count[rkmer] = kmer2count[rkmer] + 1


    temp_rid = None
    temp_rseq = None        
    with open(fasta_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if temp_rid is not None:
                    rare_kmer_check(temp_rid, temp_rseq)

                temp_rid = line.lstrip('>')
                temp_rseq = ''
            else:
                temp_rseq = temp_rseq + line

        if temp_rid is not None:
            rare_kmer_check(temp_rid, temp_rseq)

    return(kmer2count)


def rare_kmer_parse(cen_fasta_list, output_file, kmer_size):

    black_list_kmer = {}
    unique_kmer_list_raw = {}
    with open(cen_fasta_list, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            seq_name = F[0]
            fasta_file = F[1]

            temp_kmer2count = gather_kmer(fasta_file, int(kmer_size))
            # print(f'Read sequence {seq_name}, kmer count: {len(temp_kmer2count)}')
    
            for kmer in temp_kmer2count:
                if temp_kmer2count[kmer] > 2: black_list_kmer[kmer] = 1


            for kmer in temp_kmer2count:    
                if kmer not in black_list_kmer and temp_kmer2count[kmer] == 1: unique_kmer_list_raw[kmer] = 1

            # print(f'Current unique kmer count: {len(unique_kmer_list_raw)}')


    unique_kmer_list = {}
    for kmer in unique_kmer_list_raw:
        if kmer not in black_list_kmer: unique_kmer_list[kmer] = 1

    logger.info(f'Rare kmer count: {len(unique_kmer_list)}')
    with open(output_file, 'w') as hout:
        for kmer in unique_kmer_list:
            print(kmer, file = hout)



def filt_kmer_from_fasta_file(fasta_file, kmer_list, filt_kmer_list, kmer_size = 27):

    def get_filt_kmer_from_seq(qid, qseq):
        cur_kmer2hash = {} 
        for i in range(len(qseq) - kmer_size + 1):
            kmer = qseq[i:(i + kmer_size)]
            rkmer = reverse_complement(kmer)

            if kmer in kmer_list:
                cur_kmer2hash[kmer] = hash(kmer)
            elif rkmer in kmer_list:
                cur_kmer2hash[rkmer] = hash(rkmer)
            elif len(cur_kmer2hash) > 0:
                min_kmer = min(cur_kmer2hash, key = cur_kmer2hash.get)
                filt_kmer_list[min_kmer] = 1
                cur_kmer2hash = {}

        if len(cur_kmer2hash) > 0:
            min_kmer = min(cur_kmer2hash, key = cur_kmer2hash.get)
            filt_kmer_list[min_kmer] = 1

                    
    temp_rid = None
    temp_rseq = None
    with open(fasta_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if temp_rid is not None:
                    get_filt_kmer_from_seq(temp_rid, temp_rseq)

                temp_rid = line.lstrip('>')
                temp_rseq = ''
            else:
                temp_rseq = temp_rseq + line

        if temp_rid is not None:
            get_filt_kmer_from_seq(temp_rid, temp_rseq)



def filt_kmer(kmer_file, cen_fasta_list, output_file):

    kmer_list = {}
    kmer_size = None
    filt_kmer_list = {}
    with open(kmer_file, 'r') as hin:
        for line in hin:
            kmer = line.rstrip('\n')
            kmer_list[kmer] = 1
            kmer_size = len(kmer)


    with open(cen_fasta_list, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            seq_name = F[0]
            fasta_file = F[1]

            filt_kmer_from_fasta_file(fasta_file, kmer_list, filt_kmer_list) 

    logger.info(f'Filtered rare kmer count: {len(filt_kmer_list)}')
    hout = open(output_file, 'w') 
    for kmer in filt_kmer_list:
        print(kmer, file = hout)

    hout.close()


def add_rare_kmer_info(kmer_file, cen_fasta_list, output_file, alpha_satellite_margin_size):

    def kmer_pos_check(qid, qseq, contig_len):

        for i in range(len(qseq) - kmer_size + 1):

            kmer = qseq[i:(i + kmer_size)]
            rkmer = reverse_complement(kmer)

            if kmer in kmer2info:
                kmer2info[kmer].append(f'{qid},{contig_len},{i},+')

            if rkmer in kmer2info:
                kmer2info[rkmer].append(f'{qid},{contig_len},{i},-')


    kmer_size = None
    kmer2info = {}
    with open(kmer_file, 'r') as hin:
        for line in hin:
            kmer2info[line.rstrip('\n')] = []
            kmer_size = len(line.rstrip('\n'))

    with open(cen_fasta_list, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            seq_name = F[0]
            fasta_file = F[1]

            # is_inv_contig, inv_start, inv_end = False, None, None
            # contig_len, strand = None, None

            temp_rid = None
            temp_rseq = None
            with open(fasta_file, 'r') as hin:

                for line in hin:
                    line = line.rstrip('\n')
                    if line.startswith('>'):
                        if temp_rid is not None:
                            kmer_check(temp_rid, temp_rseq)

                        temp_rid = line.lstrip('>')
                        temp_rseq = ''
                    else:
                        temp_rseq = temp_rseq + line

                if temp_rid is not None:
                    contig_len = len(temp_rseq)
                    kmer_pos_check(seq_name, temp_rseq, contig_len)



    with open(output_file, 'w') as hout:
        for kmer in kmer2info:

            rel_pos, num, num_kmer_is_inv, num_is_inv_contig = [], 0, 0, 0
            for FFF in kmer2info[kmer]:
                seq_name, contig_len, pos, strand = FFF.split(',')
                # rel_pos.append( (float(pos) - 100000) / (float(contig_len) - 200000) )
                rel_pos.append( (float(pos) - alpha_satellite_margin_size) / float(int(contig_len) - 2 * alpha_satellite_margin_size) )
                num = num + 1

            rel_pos_median = statistics.median(rel_pos)
            rel_pos_var = statistics.variance(rel_pos) if len(rel_pos) > 1 else None

            tinfo = ';'.join(kmer2info[kmer])
            print(f'{kmer}\t{rel_pos_median}\t{rel_pos_var}\t{num}\t{tinfo}', file = hout)



def proc_rare_kmer_table(kmer_file, output_file):

    with open(kmer_file, 'r') as hin, open(output_file, 'w') as hout:
        print("Marker\tRelative_pos_mean\tRelative_pos_std\tMarker_num\tHaplotype\t\tContig_len\tMarker_pos\tMarker_strand", file = hout)
        for line in hin:
            F = line.rstrip('\n').split('\t')
            FF = F[4].split(';')
            for elm in FF:
                FFF = elm.split(',')
                pos_mean = round(float(F[1]), 8)
                pos_sqrt = "None" if F[3] == "1" else round(float(F[2]), 8)
                print(f'{F[0]}\t{pos_mean}\t{pos_sqrt}\t{F[3]}\t{FFF[0]}\t{FFF[1]}\t{FFF[2]}\t{FFF[3]}', file = hout)

