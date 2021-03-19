
# NAT-9: extract from bam files the reads and mappings that overlap our panel

This task is done by the script `bin/filter_bam.py`.

Usage:
```
filter_bam.py input_bam_file input_tsv_manifest_file output_directory \
  -q/--ext_qual <int, default=41>
```

Assume `input_bam_file` is `XXX/BAM_FILE.bam`

The script reads all entries of the input BAM file and output in a BAM file
in `output_directory` named `BAM_FILE-filtered.bam` all entries that overlap
with at least one amplicon of the manifest.

The script also outputs two FASTQ files (one file of forward reads and one file
of reverse reads, both in output_directory and named respectively
`BAM_FILE_L001_R1_001_annotated.fq` and `BAM_FILE_L001_R2_001_annotated.fq`).
From a given primary  mapping overlapping an amplicon, any part of the read that
maps outside of the amplicon is trimmed out.
There are then two options; if the argument `ext_qual` is > 0 then if the
mapping only partially overlaps the amplicon, it is completed by a prefix and a
suffix in order to cover the full amplicon and the Phred quality of the
extensions is set to `ext_qual`. If `ext_qual` takes value 0, the read is not
extended.

In the FASTQ files, the reads header follows the format of annotated FASTQ files
with an added information field OVERLAP that contains the coordinates (in 0-base
amplicon coordinates) of the overlap between the read and the amplicon as
defined by the primary mapping of the read.

A log file in output_directory named   `BAM_FILE.log` records statistics about
the extracted reads and mappings.

### Example

```
mkdir data results
./tools/sratoolkit.2.11.0-ubuntu64/bin/sam-dump SRR11850665 | samtools view -bS - > ./data/SRR11850665.bam
./bin/filter_bam.py data/SRR11850665.bam assets/CG001v5.1_Amplicon_Manifest_Panel5.1.14_20201119.tsv results
```
