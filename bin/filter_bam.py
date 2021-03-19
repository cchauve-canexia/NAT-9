#!/usr/bin/env python3

"""
Extracting from a BAM file data (mappings and reads) overlaping the amplicons of
a given manifest.

usage: filter_bam.py input_bam_file input_tsv_manifest_file output_directory \
  -q/--ext_qual <int, default=41>

Assume input_bam_file = .../BAM_FILE.bam

The script reads all entries of the input BAM file and output in a BAM file
in output_directory named "BAM_FILE-filtered.bam" all entries that overlap with
at least one amplicon of the manifest.

The script also outputs two FASTQ files (one file of forward reads and one file
of reverse reads, both in output_directory and named respectively
"BAM_FILE_L001_R1_001_annotated.fq" and "BAM_FILE_L001_R2_001_annotated.fq").
From a given primary  mapping overlapping an amplicon, any part of the read that
maps outside of the amplicon is trimmed out.
There are then two options; if the argument ext_qual > 0 then if the mapping
partially overlaps the amplicon, it is completed by a prefix and a suffix in
order to cover the full amplicon and the Phred qyuality of the extensions is
set to ext_qual. If ext_qual==0, the read is not extended.

In the FASTQ files, the reads header follows the format of annotated FASTQ files
with an added information field OVERLAP that contains the coordinates (in 0-base
amplicon coordinates) of the overlap between the read and the amplicon as
defined by the primary mapping of the read.

A log file in output_directory named "BAM_FILE.log" records statistics about
the extracted reads and mappings.
"""

import argparse
from collections import defaultdict
import os
import subprocess

import pandas as pd
import pysam

# Phred quality for extension of reads not covering fully their amplicon
EXT_QUAL = 41

# Auxiliary functions

def phred_to_str(quality):
    """
    Return the string version of Phred quality
    :param: quality (list(int)): numerical Phred quality
    :return: str: string version of Phred quality
    """
    return ''.join([chr(x + 33) for x in quality])

REVERSE_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

def rev_comp(seq):
    """
    :param seq (str): DNA sequence
    :return: (str): the reverse complement of seq
    taken from cg.pipelineutils.qanexus.codeword_from_fq
    """
    return ''.join(map(lambda x: REVERSE_COMPLEMENT[x], seq))[::-1]

def trim_extend_read(
    ref_start, ref_end,
    read_seq, read_qual,
    amp_seq, amp_start, amp_end,
    ext_qual=EXT_QUAL
):
    """
    Given a mapping of a read over an amplicon, trim or extend the mapping to
    cover the full amplicon.
    If the mapping extends to the left of the amplicon, the extension is
    trimmed.
    If the mapping extends to the right, the extension is trimmed.
    If the mapping starts after the first amplicon base and ext_qual>0, the
    read  is prefixed by the missing prefix of the amplicon.
    If the mapping ends before the end of the amplicon and ext_qual>0, the read
    is suffixed by the missing suffix of the amplicon.
    In case of prefix or suffix extension, the quality of added
    bases is set to ext_qual
    :param: ref_start (int): position in the amplicon (in chromosomal
    coordinates of the first aligned base of the read)
    :param: ref_end (int): same for the last aligned base
    :param: read_seq (str): sequence of the read that is aligned
    :param: amp_seq (str): amplicon sequence
    :param: amp_start, amp_end (int): chromosomal coordinates of the amplicon
    :param: ext_qual (int): if >0 extend in case the read does not cover the
    whole amplicon and assigns ext_qual as Phred quality to the extensions
    :return: (str, list(int), int, int):
    - trimmed/extended read sequence
    - trimmed/extended Phred quality (list of numerical values)
    - start, end of the aligned part of the read in 0-base amplicon coordinates
    """
    if ref_start < amp_start:
        prefix_start = amp_start - ref_start
        seq1, qual1 = read_seq[prefix_start:], read_qual[prefix_start:]
        start = 0
    elif ext_qual > 0:
        prefix_seq = amp_seq[0:ref_start - amp_start]
        prefix_qual = [ext_qual for i in range(len(prefix_seq))]
        seq1, qual1 = f"{prefix_seq}{read_seq}", prefix_qual + read_qual
        start = ref_start - amp_start
    else:
        seq1, qual1 = read_seq, read_qual
        start = ref_start - amp_start
    # Trimming/extending the suffix of the alignment
    if ref_end > amp_end:
        seq, qual = seq1[:-(ref_end - amp_end)], qual1[:-(ref_end - amp_end)]
        end = amp_end - amp_start
    elif ext_qual > 0:
        suffix_seq = amp_seq[ref_end - amp_start + 1:]
        suffix_qual = [ext_qual for i in range(len(suffix_seq))]
        seq, qual = f"{seq1}{suffix_seq}", qual1 + suffix_qual
        end = ref_end - amp_start
    else:
        seq, qual = seq1, qual1
        end = ref_end - amp_start
    return (seq, qual, start, end)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Natera project: preprocess BAM file'
    )
    parser.add_argument('bam_file', type=str, help='Input BAM file')
    parser.add_argument('manifest_tsv_file', type=str, help='Amplicon manifest')
    parser.add_argument('out_dir', type=str, help='Output directory')
    parser.add_argument(
        '-q', '--ext_qual', type=int, default=EXT_QUAL,
        help='Phred quality extension (0=no extension)'
    )
    args = parser.parse_args()

    # Reading the amplicon manifest
    manifest_df = pd.read_csv(args.manifest_tsv_file, sep='\t')

    # Name of the input BAM file
    _, bam_file_name =  os.path.split(args.bam_file)
    # Path to sorted input BAM file
    sorted_bam_file = args.bam_file.replace('.bam', '-sorted.bam')
    # Path to index of sorted BAM file
    index_sorted_bam_file = f"{sorted_bam_file}.bai"
    # Sorting and indexing the input BAM file if not already done
    if not os.path.isfile(sorted_bam_file):
        print('Sorting input BAM file')
        pysam.sort('-o', sorted_bam_file, args.bam_file)
    if not os.path.isfile(index_sorted_bam_file):
        print('Indexing sorted input BAM file')
        subprocess.call(['samtools', 'index', sorted_bam_file])
    # Reading the sorted input BAM file
    pysam_in_bam_file = pysam.AlignmentFile(sorted_bam_file, 'rb')

    # Output FASTQ files
    out_fq_file_prefix = bam_file_name.replace('.bam', '_L001_')
    out_fq_fwd_file_name = os.path.join(
        args.out_dir, f"{out_fq_file_prefix}R1_001_annotated.fq"
    )
    out_fq_rev_file_name = os.path.join(
        args.out_dir, f"{out_fq_file_prefix}R2_001_annotated.fq"
    )
    out_fq_fwd_file = open(out_fq_fwd_file_name, 'w')
    out_fq_rev_file = open(out_fq_rev_file_name, 'w')
    # Strings to be written in FASTQ files
    out_fwd_str, out_rev_str = '', ''
    # Output BAM file
    out_bam_file_name = os.path.join(
        args.out_dir, bam_file_name.replace('.bam', '_filtered.bam')
    )
    pysam_out_bam_file = pysam.AlignmentFile(
        out_bam_file_name, mode='wb', template=pysam_in_bam_file
    )
    # Log file
    log_file = open(out_bam_file_name.replace('_filtered.bam', '.log'), 'w')
    log_file.write(f"INPUT:BAM_FILE\t{args.bam_file}\n")
    log_file.write(f"INPUT:MANIFEST_FILE\t{args.manifest_tsv_file}\n")
    log_file.write(f"INPUT:EXT_QUAL\t{args.ext_qual}\n")
    log_file.write(f"OUTPUT:BAM_FILE\t{out_bam_file_name}\n")
    log_file.write(
        f"OUTPUT:FASTQ_FILES\t{out_fq_fwd_file_name}\t{out_fq_rev_file_name}\n"
    )
    # Dictionary recording if a mapping has already been written in the output
    # BAM file
    out_bam_reads_written = defaultdict(bool)
    # Loop on amplicons
    for _, amplicon in manifest_df.iterrows():
        # Log statistics
        log_nb_bam_entries = 0 # Number of exported BAM entries
        log_nb_fwd_reads = 0 # Number of forward reads
        log_nb_fwd_bases = 0 # Number of bases from forward reads overlaping
        log_nb_rev_reads = 0 # Number of reverse reads
        log_nb_rev_bases = 0 # Number of bases from reverse reads overlaping
        # Amplicon features
        amp_id = amplicon['Amplicon_ID']
        amp_chr = amplicon['Chr']
        amp_start, amp_end = int(amplicon['Start']), int(amplicon['End'])
        amp_seq = amplicon['Amplicon']
        amp_read_count = pysam_in_bam_file.count(
            contig=amp_chr, start=amp_start, stop=amp_end
        )
        if amp_read_count > 0:
            log_file.write(f"STAT:NB_READS\t{amp_id}:{amp_read_count}\n")
        # Loop over all mappings overlaping the current amplicon
        for read in pysam_in_bam_file.fetch(
            contig=amp_chr, start=amp_start, stop=amp_end
        ):
            # Mapping features
            read_name = read.query_name
            ref_start, ref_end = read.reference_start, read.reference_end
            query_start = read.query_alignment_start
            query_end = read.query_alignment_end
            cigar = read.cigarstring
            # Exporting the mapping into a BAM file, ensuring no mapping is
            # recorded twice; idx = index of the mapping
            # Two mappings are identical if they concern the same read, have the
            # same coordinates and the same CIGAR string
            idx = (read_name, ref_start, ref_end, query_start, query_end, cigar)
            if not out_bam_reads_written[idx]:
                log_nb_bam_entries += 1
                out_bam_reads_written[idx] = True
                pysam_out_bam_file.write(read)
            # Exporting primary mappings into FASTQ files
            if not read.is_secondary:
                # Mapping features
                query_seq = read.query_alignment_sequence
                query_qual = list(read.query_alignment_qualities)
                # Computing the overlap of the mapping with the amplicon,
                # trimming it if needed and extensing it if args.ext_qual > 0
                (seq, qual, start, end) = trim_extend_read(
                    ref_start, ref_end, query_seq, query_qual,
                    amp_seq, amp_start, amp_end,
                    ext_qual=args.ext_qual
                )
                # Exporting the read as a forward or reverse read
                header_suffix = (
                    f"N:0:2:AMPLICON:{amp_id}:CONT_BARCODE::CODEWORD::"
                    f"OVERLAP:{start}:{end}"
                )
                if read.is_read1:
                    # Forward read
                    log_nb_fwd_reads += 1
                    log_nb_fwd_bases += end - start + 1
                    header = f"@{read_name} 1:{header_suffix}"
                    qual_str, seq_str = phred_to_str(qual), seq
                    out_fwd_str +=  f"{header}\n{seq_str}\n+\n{qual_str}\n"
                elif read.is_read2:
                    # reverse read
                    log_nb_rev_reads += 1
                    log_nb_rev_bases += end - start + 1
                    header = f"@{read_name} 2:{header_suffix}"
                    qual.reverse()
                    qual_str, seq_str = phred_to_str(qual), rev_comp(seq)
                    out_rev_str += f"{header}\n{seq_str}\n+\n{qual_str}\n"
        log_file.write(f"STAT:NB_BAM_ENTRIES\t{amp_id}:{log_nb_bam_entries}\n")
        log_file.write(f"STAT:NB_FWD_READS\t{amp_id}:{log_nb_fwd_reads}\n")
        log_file.write(f"STAT:NB_FWD_BASES\t{amp_id}:{log_nb_fwd_bases}\n")
        log_file.write(f"STAT:NB_REV_READS\t{amp_id}:{log_nb_rev_reads}\n")
        log_file.write(f"STAT:NB_REV_BASES\t{amp_id}:{log_nb_rev_bases}\n")
    # Writing two annotated FASTQ files, for forward and reverse reads
    out_fq_fwd_file.write(out_fwd_str[:-1]) # [:-1] to exclude the last \n
    out_fq_rev_file.write(out_rev_str[:-1])
    # Closing opened files
    pysam_out_bam_file.close()
    out_fq_fwd_file.close()
    out_fq_rev_file.close()
    log_file.close()
